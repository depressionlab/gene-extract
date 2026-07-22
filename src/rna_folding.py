#!/usr/bin/env python3
"""
3'UTR structural divergence pipeline

Hypothesis: differential mRNA localization across species correlates with
differences in 3'UTR secondary structure, not just sequence divergence.

Pipeline:
  1. Align orthologous 3'UTRs (MAFFT)
  2. Fold each sequence individually (RNAfold/ViennaRNA)
  3. Simultaneous align + fold comparison (LocARNA)
  4. Quantitative structural distance (RNAforester)
  5. Covariation support check (R-scape)
"""
import subprocess
import tempfile
from pathlib import Path

import RNA

def align_with_mafft(input_fasta: Path, output_fasta: Path) -> None:
    """Stage 1: MAFFT multiple sequence alignment of orthologous 3'UTRs."""
    with open(output_fasta, "w") as out:
        subprocess.run(["mafft", "--auto", str(input_fasta)], stdout=out, check=True)

def fold_sequence(seq: str) -> tuple[str, float]:
    """Stage 2: RNAfold MFE structure + energy for a single sequence."""
    structure, mfe = RNA.fold(seq)
    return structure, mfe

def fold_with_partition_function(seq: str) -> dict:
    """Stage 2b: MFE + partition function/ensemble diversity."""
    fc = RNA.fold_compound(seq)
    mfe_structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    pf_structure, pf_energy = fc.pf()
    ensemble_diversity = fc.mean_bp_distance()
    return {
        "mfe_structure": mfe_structure,
        "mfe": mfe,
        "pf_structure": pf_structure,
        "pf_energy": pf_energy,
        "ensemble_diversity": ensemble_diversity,
    }

def run_locarna(
    seq_a_fasta: Path,
    seq_b_fasta: Path,
    clustal_out: Path,
    stockholm_out: Path | None = None,
    extra_args: list[str] | None = None,
) -> dict:
    """Stage 3: LocARNA simultaneous alignment + folding of two sequences.

    LocARNA jointly aligns two RNA sequences using both sequence and
    predicted-structure information. This is the point of the pipeline
    extension: sequence alignment (MAFFT, stage 1) ignores structure
    entirely, so two 3'UTRs that are structurally conserved but
    sequence-diverged would look unrelated to MAFFT but LocARNA should
    recover the shared fold.

    `stockholm_out`, if given, also writes a Stockholm-format alignment
    (with the consensus structure as `#=GC SS_cons`). This is the
    input format `run_rscape` (stage 5) needs, so passing this through
    for a whole ortholog group chains stages 3 and 5 directly.

    Returns a dict with the raw stdout, the alignment score parsed from
    it, and the paths to the alignment file(s) LocARNA wrote.
    """
    cmd = ["locarna", str(seq_a_fasta), str(seq_b_fasta), "--clustal", str(clustal_out)]
    if stockholm_out:
        cmd.extend(["--stockholm", str(stockholm_out)])
    if extra_args:
        cmd.extend(extra_args)
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    score = None
    for line in result.stdout.splitlines():
        if line.strip().lower().startswith("score:"):
            try:
                score = float(line.split(":", 1)[1].strip())
            except ValueError:
                pass

    return {
        "score": score,
        "stdout": result.stdout,
        "clustal_file": clustal_out,
        "stockholm_file": stockholm_out,
    }


def run_rnaforester(
    name_a: str,
    structure_a: str,
    name_b: str,
    structure_b: str,
) -> dict:
    """Stage 4: RNAforester quantitative structural (tree-edit) distance.

    Takes two dot-bracket secondary structures (e.g. the `mfe_structure`
    values from `fold_with_partition_function`) and computes their
    pairwise tree-alignment similarity score. This complements LocARNA:
    LocARNA jointly aligns sequence+structure, while RNAforester scores
    the structures alone, so a large sequence-alignment score with a low
    RNAforester score (or vice versa) is itself a signal worth flagging.
 
    NOTE: `--score` only ever compares the *first two* structures in the
    input, even if more are supplied. It is not an all-pairs matrix tool
    the way the name might suggest. This function is deliberately
    pairwise-only (one ortholog pair per call) to avoid silently
    mis-scoring if more sequences are ever passed in.
 
    Returns both the raw tree-alignment score and the length-normalized
    relative score (0-1), the latter being what you want when comparing
    across orthologs of different lengths.
    """
    payload = f">{name_a}\n{structure_a}\n>{name_b}\n{structure_b}\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
        tmp.write(payload)
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            ["RNAforester", f"-f={tmp_path}", "--score", "-r"],
            capture_output=True,
            text=True,
            check=True,
        )
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    raw_score = float(lines[0]) if len(lines) >= 1 else None
    relative_score = float(lines[1]) if len(lines) >= 2 else None

    return {
        "raw_score": raw_score,
        "relative_score": relative_score,
        "stdout": result.stdout,
    }


def run_rscape(alignment_stockholm: Path, outdir: Path) -> dict:
    """Stage 5: R-scape covariation support check.

    R-scape tests whether base pairs in a proposed consensus secondary
    structure show statistically significant covariation across an
    alignment, beyond what phylogeny alone would predict. This is the
    strongest evidence a shared fold is functionally conserved, not just
    coincidentally similar. It needs a real multi-sequence alignment
    with an annotated consensus structure (a Stockholm file with an
    `#=GC SS_cons` line), not a single pairwise comparison. LocARNA's
    `--stockholm` output (from `run_locarna`, or a version of it re-run
    across all species at once) is exactly this format, so the natural
    place to call this is on the LocARNA alignment/consensus of a whole
    ortholog group, once stages 3-4 have flagged a structurally
    interesting region worth testing.

    Returns the list of covarying base pairs found (each with position,
    E-value, and substitution count), plus paths to R-scape's full
    output directory for anyone who wants the plots/full report.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["R-scape", "--outdir", str(outdir), str(alignment_stockholm)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    cov_files = sorted(outdir.glob("*.cov"))
    covarying_pairs = []
    if cov_files:
        with open(cov_files[0]) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split()
                covarying_pairs.append({
                    "in_given_structure": fields[0] == "*",
                    "left_pos": int(fields[1]),
                    "right_pos": int(fields[2]),
                    "score": float(fields[3]),
                    "e_value": float(fields[4]),
                    "p_value": float(fields[5]),
                    "substitutions": int(fields[6]),
                    "power": float(fields[7]),
                })

    return {
        "covarying_pairs": covarying_pairs,
        "cov_file": cov_files[0] if cov_files else None,
        "stdout": result.stdout,
        "outdir": outdir,
    }


def demo():
    """Smoke-test stages 2-4 with placeholder sequences, just to confirm
    the environment is set up correctly. Replace with real orthologous
    3'UTR sequences once we know what region to pull."""
    example_seqs = {
        "dmel_example": "AUGGCUAGCUAGCUACGUACGUAGCUAGCUACGUACGUAGCUAGCUAGCUACGUA",
        "dsim_example": "AUGGCUAGCUAGCUACGUACGUAGCUAGCUACGUACGUAGCUAGCUAGCUACGCA",
    }

    print("=== Stage 2 smoke test: individual folding ===")
    folded = {}
    for name, seq in example_seqs.items():
        result = fold_with_partition_function(seq)
        folded[name] = result
        print(f"{name}: MFE={result['mfe']:.2f}  structure={result['mfe_structure']}")
        print(f"  ensemble diversity: {result['ensemble_diversity']:.3f}")

    print("\n=== Stage 3 smoke test: LocARNA simultaneous align+fold ===")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        fasta_paths = {}
        for name, seq in example_seqs.items():
            p = tmpdir / f"{name}.fa"
            p.write_text(f">{name}\n{seq}\n")
            fasta_paths[name] = p

        names = list(example_seqs)
        stockholm_path = tmpdir / "locarna_out.sto"
        locarna_result = run_locarna(
            fasta_paths[names[0]],
            fasta_paths[names[1]],
            tmpdir / "locarna_out.aln",
            stockholm_out=stockholm_path,
        )
        print(f"LocARNA score: {locarna_result['score']}")
        print(f"Alignment written to: {locarna_result['clustal_file']}")

        print("\n=== Stage 4 smoke test: RNAforester structural distance ===")
        rnaforester_result = run_rnaforester(
            names[0], folded[names[0]]["mfe_structure"],
            names[1], folded[names[1]]["mfe_structure"],
        )
        print(f"RNAforester raw score: {rnaforester_result['raw_score']}")
        print(f"RNAforester relative score: {rnaforester_result['relative_score']}")

        print("\n=== Stage 5 smoke test: R-scape covariation check ===")
        print("(only 2 near-identical placeholder sequences here. Zero")
        print(" statistical power expected; this just confirms the tool")
        print(" runs and the .cov output parses, not a real result)")
        rscape_result = run_rscape(stockholm_path, tmpdir / "rscape_out")
        print(f"Covarying pairs found: {len(rscape_result['covarying_pairs'])}")
        for pair in rscape_result["covarying_pairs"]:
            print(f"  {pair}")


if __name__ == "__main__":
    demo()
