# gene-extract

Pull ortholog gene sequences out of FlyBase genome FASTA files!

## What does this do?

You give it:

1. A FlyBase ortholog TSV (maps *D. melanogaster* genes -> the equivalent gene ID in another *Drosophila* species)
2. A species-specific FASTA genome file (all the gene sequences for that species)
3. A list of the Dmel `FBgn...` gene IDs you want to isolate

...and it hands you back a clean FASTA file containing only the matching ortholog sequences for that species.

---

## Quickstart

### 1. Install Rust

Grab it from [rustup.rs](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### 2. Clone this repo

```bash
git clone https://github.com/depressionlab/gene-extract
cd gene-extract
```

### 3. Get your input data (or just let the tool fetch it automatically)

You need two files per species:

| File | What it is | Auto-downloaded? |
| --- | --- | --- |
| **Ortholog TSV** | Maps Dmel gene IDs -> other species' gene IDs | ✅ Yes, by default (FlyBase FB2019_03 release) |
| **Species FASTA** | The actual gene sequences for the species you want | ✅ Yes, by default (D. pseudoobscura r3.04 which matches the built-in default `--fasta` path) |

So out of the box, with zero flags, `cargo run --release` will download both files itself the first time and use them from then on (they're cached at `./data/` and it won't redownload on later runs).

Want a **different species**? Point `--fasta` at a new local path, and pass a matching `--fasta-url` so the tool knows where to fetch it:

```bash
cargo run --release -- \
  --fasta ./data/dsim-all-gene-r2.02.fasta \
  --fasta-url https://flybase.org/genomes/Drosophila_simulans/dsim_r2.02_FB2017_04/fasta/dsim-all-gene-r2.02.fasta.gz \
  --output ./filtered/dsim_filtered.fasta
```

Browse [flybase.org/genomes/](https://flybase.org/genomes/) to find the right species folder, release, and `-all-gene-*.fasta(.gz)` filename. `.gz` URLs are decompressed automatically.

Already have the files downloaded another way? Just drop them into [`data/`](./data) with the expected filenames (or point `--tsv`/`--fasta` at wherever they are) and the tool will use them as-is, skipping the download step entirely.

Working somewhere without internet access? Pass `--offline` and the tool will fail fast with a clear message instead of trying to reach the network.

### 4. Run it

For the normal, default run:

```bash
cargo run
```

If you want to use custom, non-default data:

```bash
cargo run --release -- \
  --tsv ./data/dmel_orthologs_in_drosophila_species_fb_2019_03.tsv \
  --fasta ./data/dpse-all-gene-r3.04.fasta \
  --output ./filtered/dpse_filtered.fasta
```

That's it! The `filtered/` output folder is created automatically if it doesn't exist yet.

---

## All the options

| Flag | What it does | Default |
| --- | --- | --- |
| `--tsv <path>` | Path to the FlyBase ortholog TSV | `./data/dmel_orthologs_in_drosophila_species_fb_2019_03.tsv` |
| `--tsv-url <url>` | Where to download the TSV from if missing | FlyBase FB2019_03 release |
| `--fasta <path>` | Path to the species genome FASTA | `./data/dpse-all-gene-r3.04.fasta` |
| `--fasta-url <url>` | Where to download the FASTA from if missing | D. pseudoobscura r3.04 release |
| `--output <path>` | Where to write the filtered FASTA | `./filtered/my_filtered.fasta` |
| `--gene <FBgn...>` | A Dmel gene ID to look up. Repeatable: pass it multiple times for multiple genes | built-in list of 11 genes (see below) |
| `--species <code>` | FlyBase species code to filter orthologs for (e.g. `Dpse`, `Dsim`, `Dyak`) — must match your `--fasta` species | `Dpse` |
| `--print` | Also print matched sequences to your terminal | off |
| `--offline` | Never auto-download; fail immediately if an input file is missing | off |

### Filtering multiple species

Just run the command again with different `--fasta`, `--fasta-url`, `--species`, and `--output`. The same `--tsv` and `--gene` list can be reused:

```bash
cargo run --release -- --fasta ./data/dsim-all-gene-r2.02.fasta --species Dsim --output ./filtered/dsim_filtered.fasta
cargo run --release -- --fasta ./data/dyak-all-gene-r1.05.fasta --species Dyak --output ./filtered/dyak_filtered.fasta
```

⚠️ **Always keep `--species` in sync with `--fasta`.** The ortholog TSV maps each Dmel gene to orthologs across *many* species at once. `--species` tells the tool which of those to actually look for in your FASTA. Get this wrong (or forget it) and you'll either miss real matches or chase down a pile of "not found" warnings for the wrong species' IDs.

### Reading the output

After scanning, you'll get a per-gene summary so you can see at a glance whether each gene you asked for was actually resolved for your chosen species:

```bash
▶ Gene summary
  ✓ FBgn0065109
  ✓ FBgn0022981
  ✗ FBgn0261722
  ...
  ⚠ 10/11 requested gene(s) resolved
```

A `✗` means that specific Dmel gene has no ortholog for your `--species` in either the TSV or the FASTA. It's worth double-checking the gene ID and species code, but not necessarily an issue.

### Picking your own genes

By default the tool looks up 11 built-in Dmel genes (`ppk11`, `rpk`, `fwe`, `hts`, `NO66`, `CG32069`, `yki`, `Xpc`, `asp`, `mu2`, `grk`). To use your own instead:

```bash
cargo run --release -- --gene FBgn0000140 --gene FBgn0261722
```

You can find unofficial-name -> Flybase's `FBgn` ID mappings in FlyBase's XML export.

## Troubleshooting

- **"no ortholog mapping found for: ..."**: the Dmel gene ID you asked for isn't in your TSV file. Double-check the `FBgn` ID, or make sure you downloaded the right ortholog TSV.
- **"not found in fasta file: ..."**: the ortholog ID exists in the TSV mapping, but the matching sequence wasn't found in the FASTA file you supplied. You may be pointing at the wrong species' FASTA, or an older/newer genome release with different IDs.
- **"no ortholog mapping at all for: ..."**: the Dmel gene ID you asked for isn't in your TSV file at all. Double-check the `FBgn` ID, or make sure you downloaded the right ortholog TSV.
- **"no `<species>` ortholog found for: ..."**: FlyBase has *some* ortholog data for that gene, just not for the species you specified with `--species`. Double-check the species code matches your `--fasta`.
- **"not found in fasta file: ..."**: the ortholog ID exists in the TSV mapping for your species, but wasn't found in the FASTA file you supplied. You may be pointing at the wrong species' FASTA, an older/newer genome release with different IDs, or `--species` doesn't actually match `--fasta`.
- All three warnings above are non-fatal: the tool still writes out everything it did successfully match, and the **Gene summary** step tells you exactly which requested genes made it through.
- **"failed to open fasta file" / "failed to open tsv file"**: check the path you passed to `--tsv` / `--fasta` is correct and the file actually exists there.
- **Using a non-default species and forgot `--fasta-url`?** The tool will try to download the *default* D. pseudoobscura FASTA to whatever path you gave `--fasta`, which is almost certainly not what you want. Always pass a matching `--fasta-url` when you change `--fasta` to a different species. Find yours at [flybase.org/genomes/](https://flybase.org/genomes/).
- **"server returned an error response for ..."**: the download URL is wrong, moved, or temporarily unreachable. FlyBase occasionally reorganizes its file layout; re-check the path in a browser and update `--tsv-url`/`--fasta-url` accordingly.
- **Download seems stuck / no internet access**: pass `--offline` to skip network calls entirely and get an immediate, clear error instead of a hang.

---

## Data assumptions

- The ortholog TSV follows FlyBase's standard `dmel_orthologs_in_drosophila_species` format (tab-separated, Dmel ID in column 1, ortholog ID in column 6, species name in column 7).
- Gene/sequence IDs are unique within a single genome/species — true for FlyBase gene IDs.
- The FASTA file has a parseable name on each sequence header — true for essentially all FASTA-spec files.

See the [`data/`](./data) directory for example file structures.

---

## 📚 Citing this tool

If you use `gene-extract` in published work, please cite it. A machine-readable [`CITATION.cff`](./CITATION.cff) is included at the repo root. GitHub surfaces this automatically via the **"Cite this repository"** button on the repo page.

BibTeX:

```bibtex
@article{lastname2026geneextract,
  author  = {LASTNAME, FIRSTNAME},
  title   = {TITLE OF YOUR PAPER HERE},
  journal = {JOURNAL NAME},
  year    = {2026},
  note    = {Manuscript in preparation}
}

@software{depressionlab2026geneextractsoftware,
  author  = {{depressionlab}},
  title   = {gene-extract: extracting Drosophila ortholog sequences from FlyBase genome FASTA files},
  year    = {2026},
  url     = {https://github.com/depressionlab/gene-extract},
  version = {0.1.0}
}
```

## License

This project is released under **[The Unlicense](https://unlicense.org)**. The full text is available [`LICENSE`](./LICENSE):

<details>
<summary>Click to expand full license text!</summary>

```text
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org>
```

</details>
