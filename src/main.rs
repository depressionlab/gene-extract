use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read};
use std::path::{Path, PathBuf};

use clap::Parser;
use colored::Colorize;
use eyre::{Context, Result, bail};
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use reqwest::blocking::Client;

/// Default Dmel `FlyBase` gene IDs used when no `--gene` flags are given.
/// You can find the unofficial names -> Dmel `FlyBase` ID mappings in
/// `FlyBase`'s XML export.
const DEFAULT_GENES: [&str; 11] = [
	"FBgn0065109", // 'ppk11'
	"FBgn0022981", // 'rpk'
	"FBgn0261722", // 'fwe'
	"FBgn0263391", // 'hts'
	"FBgn0266570", // 'NO66'
	"FBgn0052069", // 'CG32069'
	"FBgn0034970", // 'yki'
	"FBgn0004698", // 'Xpc'
	"FBgn0000140", // 'asp'
	"FBgn0002872", // 'mu2'
	"FBgn0001137", // 'grk'
];

/// Default download URL for the ortholog TSV, matching the default `--tsv`
/// path. `FlyBase` distributes this file gzipped; it's decompressed
/// automatically on download.
const DEFAULT_TSV_URL: &str = "https://s3ftp.flybase.org/releases/FB2019_03/precomputed_files/orthologs/dmel_orthologs_in_drosophila_species_fb_2019_03.tsv.gz";

/// Default download URL for the species FASTA, matching the default
/// `--fasta` path (D. pseudoobscura, release r3.04). If you point `--fasta`
/// at a different species/version, pass a matching `--fasta-url` too.
const DEFAULT_FASTA_URL: &str = "https://flybase.org/genomes/Drosophila_pseudoobscura/dpse_r3.04_FB2016_05/fasta/dpse-all-gene-r3.04.fasta.gz";

/// Extract Drosophila ortholog sequences from a `FlyBase` genome FASTA file,
/// using a `FlyBase` Dmel-ortholog TSV mapping to resolve gene IDs.
///
/// Missing input files are downloaded automatically when a URL is known
/// (see `--tsv-url` / `--fasta-url`).
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
	/// Path to the `FlyBase` Dmel ortholog TSV mapping file.
	#[arg(
		long,
		default_value = "./data/dmel_orthologs_in_drosophila_species_fb_2019_03.tsv"
	)]
	tsv: PathBuf,

	/// URL to auto-download the ortholog TSV from if it's not already at
	/// `--tsv`. Defaults to the `FlyBase` `FB2019_03` release, matching the
	/// default `--tsv` path. `.gz` URLs are decompressed automatically.
	#[arg(long, default_value = DEFAULT_TSV_URL)]
	tsv_url: String,

	/// Path to the species-specific `FlyBase` genome FASTA file.
	/// Download this from `FlyBase` and pass its path here.
	#[arg(long, default_value = "./data/dpse-all-gene-r3.04.fasta")]
	fasta: PathBuf,

	/// URL to auto-download the species FASTA from if it's not already at
	/// `--fasta`. Defaults to the D. pseudoobscura r3.04 release, matching
	/// the default `--fasta` path. If you change `--fasta` to a different
	/// species/version, pass a matching `--fasta-url` too. Find yours
	/// under <https://flybase.org/genomes/>. `.gz` URLs are decompressed
	/// automatically.
	#[arg(long, default_value = DEFAULT_FASTA_URL)]
	fasta_url: String,

	/// Where to write the filtered FASTA output. Parent directories are created
	/// if needed.
	#[arg(long, default_value = "./filtered/my_filtered.fasta")]
	output: PathBuf,

	/// Dmel `FlyBase` gene ID (`FBgn`...) to look up. Repeat for multiple
	/// genes. If omitted, a built-in default gene list is used.
	#[arg(long = "gene", value_name = "FBGN_ID")]
	genes: Vec<String>,

	/// Also print the matched sequences to stdout.
	#[arg(long)]
	print: bool,

	/// `FlyBase` species code to filter orthologs for (e.g. `Dpse`, `Dsim`,
	/// `Dyak`). The ortholog TSV contains mappings to *many* species at
	/// once. This must match whichever species your `--fasta` file is
	/// for, or you'll get a flood of spurious "not found" warnings for
	/// every other species' ortholog IDs.
	#[arg(long, default_value = "Dpse")]
	species: String,

	/// Never auto-download missing input files; fail immediately instead.
	#[arg(long)]
	offline: bool,
}

fn main() -> Result<()> {
	let cli = Cli::parse();

	print_banner();

	// Dmel is the reference genome the ortholog TSV is keyed *from*, not
	// keyed *to* -- there's no "Dmel ortholog of a Dmel gene" row to look
	// up. Requested Dmel FBgn IDs already ARE the target IDs, so we skip
	// the ortholog-resolution step entirely and go straight to the FASTA.
	// (Previously this silently produced zero mappings for `--species
	// Dmel`, since the TSV's species column never contains "Dmel".)
	let is_dmel = cli.species.eq_ignore_ascii_case("dmel");

	step("Checking input data");
	if !is_dmel {
		ensure_file(
			&cli.tsv,
			Some(cli.tsv_url.as_str()),
			"ortholog TSV",
			!cli.offline,
		)?;
	}
	ensure_file(
		&cli.fasta,
		Some(cli.fasta_url.as_str()),
		"genome FASTA",
		!cli.offline,
	)?;

	let default_genes: Vec<String>;
	let requested_genes: &[String] = if cli.genes.is_empty() {
		default_genes = DEFAULT_GENES.iter().map(ToString::to_string).collect();
		&default_genes
	} else {
		&cli.genes
	};

	// Dmel gene -> ortholog IDs for `cli.species` (usually just one, occasionally
	// a couple if FlyBase records more than one candidate ortholog).
	let mut gene_to_ids: HashMap<String, Vec<String>> = HashMap::new();
	let mut unmapped_genes = Vec::new();
	let mut no_species_match = Vec::new();

	if is_dmel {
		step("Resolving genes (Dmel direct lookup, no ortholog step)");
		info(&format!(
			"species is Dmel -- treating {} requested gene ID(s) as direct FlyBase IDs",
			requested_genes.len().to_string().bold()
		));
		for gene in requested_genes {
			gene_to_ids.insert(gene.clone(), vec![gene.clone()]);
		}
	} else {
		step("Loading ortholog mappings");
		info(&format!(
			"reading {}",
			cli.tsv.display().to_string().dimmed()
		));
		let (mut mappings, all_dmel_genes) = read_tsv(&cli.tsv, &cli.species)?;
		success(&format!(
			"loaded {} ortholog mapping(s) for {} ({} Dmel gene(s) total in file)",
			mappings.len().to_string().bold(),
			cli.species.cyan().bold(),
			all_dmel_genes.len().to_string().bold()
		));

		step("Resolving genes to ortholog IDs");
		info(&format!(
			"resolving {} requested gene(s) for species {}",
			requested_genes.len().to_string().bold(),
			cli.species.cyan().bold()
		));

		for gene in requested_genes {
			match mappings.remove(gene) {
				Some(ortholog_ids) => {
					gene_to_ids.insert(gene.clone(), ortholog_ids);
				}
				None if all_dmel_genes.contains(gene) => no_species_match.push(gene.as_str()),
				None => unmapped_genes.push(gene.as_str()),
			}
		}
	}
	if !unmapped_genes.is_empty() {
		warn(&format!(
			"no ortholog mapping at all for: {}",
			unmapped_genes.join(", ")
		));
	}
	if !no_species_match.is_empty() {
		warn(&format!(
			"no {} ortholog found for: {}",
			cli.species,
			no_species_match.join(", ")
		));
	}
	// Built once, from data we already own. There's no need to clone every ortholog
	// ID a second time on top of the clone already in `gene_to_ids`.
	let flybase_ids: HashSet<String> = gene_to_ids.values().flatten().cloned().collect();
	success(&format!(
		"{} ortholog ID(s) to search for",
		flybase_ids.len().to_string().bold()
	));

	step("Scanning genome FASTA");
	info(&format!(
		"scanning {}",
		cli.fasta.display().to_string().dimmed()
	));
	let sequences = extract_sequences(&cli.fasta, &flybase_ids)?;

	let missing_ids: Vec<&String> = flybase_ids
		.iter()
		.filter(|id| !sequences.contains_key(*id))
		.collect();
	if !missing_ids.is_empty() {
		let missing = missing_ids
			.iter()
			.map(|s| s.as_str())
			.collect::<Vec<_>>()
			.join(", ");
		warn(&format!("not found in fasta file: {missing}"));
	}
	success(&format!(
		"matched {} of {} ortholog ID(s)",
		sequences.len().to_string().bold(),
		flybase_ids.len().to_string().bold()
	));

	step("Gene summary");
	let mut hits = 0usize;
	for gene in requested_genes {
		let found = gene_to_ids
			.get(gene)
			.is_some_and(|ids| ids.iter().any(|id| sequences.contains_key(id)));
		if found {
			hits += 1;
			println!("  {} {}", "✓".green().bold(), gene);
		} else {
			println!("  {} {}", "✗".red().bold(), gene.red());
		}
	}
	println!(
		"  {} {}/{} requested gene(s) resolved",
		if hits == requested_genes.len() {
			"✅"
		} else {
			"⚠"
		},
		hits.to_string().bold(),
		requested_genes.len().to_string().bold()
	);

	step("Writing output");
	if let Some(parent) = cli.output.parent().filter(|p| !p.as_os_str().is_empty()) {
		fs::create_dir_all(parent)
			.with_context(|| format!("failed to create output directory: {}", parent.display()))?;
	}

	let mut file_writer = noodles::fasta::io::Writer::new(BufWriter::with_capacity(
		1 << 20,
		File::create(&cli.output)
			.with_context(|| format!("failed to create output file: {}", cli.output.display()))?,
	));

	let stdout = std::io::stdout();
	let mut stdout_writer = cli
		.print
		.then(|| noodles::fasta::io::Writer::new(stdout.lock()));

	for record in sequences.values() {
		if let Some(writer) = stdout_writer.as_mut() {
			writer.write_record(record)?;
		}
		file_writer.write_record(record)?;
	}

	println!();
	success(&format!(
		"wrote {} sequence(s) to {}",
		sequences.len().to_string().bold(),
		cli.output.display().to_string().green().bold()
	));
	println!("{}", "🧬 all done!".magenta().bold());

	Ok(())
}

/// Reads records from `fasta_file`, keeping only those whose name is in `ids`.
fn extract_sequences(
	fasta_file: &Path,
	ids: &HashSet<String>,
) -> Result<HashMap<String, noodles::fasta::record::Record>> {
	let file = File::open(fasta_file)
		.with_context(|| format!("failed to open fasta file: {}", fasta_file.display()))?;
	// 1 MiB read buffer: genome FASTA files can be large, and the default
	// 8 KiB buffer means far more syscalls than necessary.
	let mut reader = noodles::fasta::io::Reader::new(BufReader::with_capacity(1 << 20, file));

	let mut sequences = HashMap::with_capacity(ids.len());

	for result in reader.records() {
		let record = result.context("failed to read fasta record")?;
		// Only allocate an owned String for records that actually match.
		// Most records in a whole-genome FASTA won't, so there's no point
		// paying for `.into_owned()` on every single one.
		let name = String::from_utf8_lossy(record.name());
		if ids.contains(name.as_ref()) {
			sequences.insert(name.into_owned(), record);
		}
	}

	Ok(sequences)
}

/// Reads the ortholog TSV, keeping only orthologs for `species`, grouped by
/// Dmel flybase ID. Also returns the full set of Dmel gene IDs seen in the
/// file (regardless of species) so callers can distinguish "no ortholog
/// data at all for this gene" from "no ortholog for this species
/// specifically" without having to keep every other species' rows around.
fn read_tsv(
	tsv_path: &Path,
	species: &str,
) -> Result<(HashMap<String, Vec<String>>, HashSet<String>)> {
	let mut reader = csv::ReaderBuilder::new()
		.delimiter(b'\t')
		.from_path(tsv_path)
		.with_context(|| format!("failed to open tsv file: {}", tsv_path.display()))?;

	let mut mappings: HashMap<String, Vec<String>> = HashMap::new();
	let mut all_dmel_genes: HashSet<String> = HashSet::new();

	for result in reader.records() {
		let record = result.context("failed to read tsv record")?;
		let dmel_id = &record[0];
		all_dmel_genes.insert(dmel_id.to_string());

		// Skip rows for any species we don't care about. There's no point paying for
		// the allocation or the HashMap entry if we'll never look at it.
		let Some(record_species) = record[6].get(0..4) else {
			continue;
		};
		if !record_species.eq_ignore_ascii_case(species) {
			continue;
		}

		mappings
			.entry(dmel_id.to_string())
			.or_default()
			.push(record[5].to_string());
	}

	Ok((mappings, all_dmel_genes))
}

fn print_banner() {
	println!(
		"\n{}",
		"🧬  gene-extract: FlyBase ortholog FASTA filter"
			.magenta()
			.bold()
	);
}

fn step(msg: &str) {
	println!("\n{} {}", "▶".cyan().bold(), msg.bold());
}

fn info(msg: &str) {
	println!("  {} {msg}", "ℹ".blue().bold());
}

fn success(msg: &str) {
	println!("  ✅ {msg}");
}

fn warn(msg: &str) {
	eprintln!("  {} {}", "⚠".yellow().bold(), msg.yellow());
}

/// Makes sure `path` exists, downloading it from `url` if it's missing and
/// `allow_download` is true. `what` is a short human-readable label used in
/// log messages (e.g. `"ortholog TSV"`).
fn ensure_file(path: &Path, url: Option<&str>, what: &str, allow_download: bool) -> Result<()> {
	if path.exists() {
		success(&format!(
			"found {what} at {}",
			path.display().to_string().dimmed()
		));
		return Ok(());
	}

	warn(&format!("{what} not found at {}", path.display()));

	if !allow_download {
		bail!(
			"{what} is missing and --offline was set. Download it manually and place it at {}.",
			path.display()
		);
	}

	match url {
		Some(url) => download_file(url, path, what),
		None => {
			bail!(
				"no download URL known for {what}. Pass --fasta-url <URL> (find yours at \
				 https://flybase.org/genomes/), or download it manually and place it at {}. \
				 See README.md for details.",
				path.display()
			);
		}
	}
}

/// Downloads `url` to `dest`, showing a colorful progress bar. Transparently
/// gunzips the content if `url` ends in `.gz`.
fn download_file(url: &str, dest: &Path, what: &str) -> Result<()> {
	println!(
		"  📥 downloading {} from:\n      {}",
		what.cyan().bold(),
		url.dimmed()
	);

	if let Some(parent) = dest.parent().filter(|p| !p.as_os_str().is_empty()) {
		fs::create_dir_all(parent)
			.with_context(|| format!("failed to create directory: {}", parent.display()))?;
	}

	let client = Client::builder()
		.user_agent(concat!("gene-extract/", env!("CARGO_PKG_VERSION")))
		.build()
		.context("failed to build HTTP client")?;

	let response = client
		.get(url)
		.send()
		.with_context(|| format!("failed to reach {url}"))?
		.error_for_status()
		.with_context(|| format!("server returned an error response for {url}"))?;

	let total_size = response.content_length().unwrap_or(0);

	let pb = ProgressBar::new(total_size);
	pb.set_style(
		#[expect(clippy::literal_string_with_formatting_args)]
		ProgressStyle::with_template(
			"      {spinner:.magenta} [{bar:32.cyan/blue}] {bytes}/{total_bytes} ({bytes_per_sec}, eta {eta})",
		)
		.unwrap_or_else(|_| ProgressStyle::default_bar())
		.progress_chars("█▉▊▋▌▍▎▏  "),
	);

	let tracked = ProgressRead {
		inner: response,
		pb: &pb,
	};

	let tmp_path = dest.with_extension("part");
	let mut out = File::create(&tmp_path)
		.with_context(|| format!("failed to create temp file: {}", tmp_path.display()))?;

	if Path::new(url)
		.extension()
		.is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
	{
		let mut decoder = GzDecoder::new(tracked);
		std::io::copy(&mut decoder, &mut out)
			.with_context(|| format!("failed to decompress/write {}", dest.display()))?;
	} else {
		let mut tracked = tracked;
		std::io::copy(&mut tracked, &mut out)
			.with_context(|| format!("failed to write {}", dest.display()))?;
	}

	pb.finish_and_clear();

	fs::rename(&tmp_path, dest)
		.with_context(|| format!("failed to finalize {}", dest.display()))?;

	success(&format!("saved to {}", dest.display().to_string().green()));
	Ok(())
}

/// A `Read` wrapper that ticks a progress bar forward on every chunk read.
struct ProgressRead<'a, R> {
	inner: R,
	pb: &'a ProgressBar,
}

impl<R: Read> Read for ProgressRead<'_, R> {
	fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
		let n = self.inner.read(buf)?;
		self.pb.inc(n as u64);
		Ok(n)
	}
}
