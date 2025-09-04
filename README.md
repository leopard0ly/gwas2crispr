# gwas2crispr

> **GWAS‑to‑CRISPR**: streamlined extraction of significant GWAS SNPs, metadata aggregation and optional FASTA/BED/CSV export for downstream CRISPR design (GRCh38/hg38).

## Overview

Genome‑wide association studies (**GWAS**) link traits to genetic variants, but raw summary statistics are not directly usable for guide design.  `gwas2crispr` bridges this gap.  It retrieves significant single‑nucleotide polymorphisms (**SNPs**) for a given **Experimental Factor Ontology** (**EFO**) trait, annotates them with gene and study metadata, and returns in‑memory summaries.  When requested, it also writes ready‑to‑use **CSV**, **BED** and **FASTA** files for high‑throughput CRISPR target design.  All genomic coordinates are mapped to GRCh38/hg38.

### Core functions

* `fetch_gwas(efo_id, p_cut = 5e-8)`: fetches significant associations for an EFO trait via `gwasrapidd` with a REST API fallback.
* `run_gwas2crispr(efo_id, p_cut = 5e-8, flank_bp = 200, out_prefix = NULL)`: end‑to‑end pipeline that calls `fetch_gwas()`, aggregates variant/gene/study metadata, and returns an object with summaries.  If you provide `out_prefix`, it will also write `CSV`, `BED` and optional `FASTA` files.

> **CRAN‑safe examples:** the package does **not** write files by default.  Examples that perform network operations or file writing are wrapped in `\donttest{}`.  When you supply `out_prefix`, outputs are written only to paths you specify — in documentation we use `tempdir()`.

---

## Installation

### Requirements (read first)

* **R ≥ 4.1**
* **Required runtime client:** **gwasrapidd** — used to query the GWAS Catalog for associations. *Install it before running `gwas2crispr`.*
* Core CRAN stack: `httr`, `dplyr`, `purrr`, `readr`, `tibble`, `tidyr`, `methods`, `utils` (pulled automatically)
* **FASTA output requires (Bioconductor):** `Biostrings`, `BSgenome.Hsapiens.UCSC.hg38`
  If these are missing, CSV/BED are still produced; FASTA is skipped gracefully.
* (Optional for CLI) `optparse`

#### Install the core prerequisite (GWAS Catalog client)

```r
install.packages("gwasrapidd")
```

#### Install Bioconductor dependencies (for FASTA)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "BSgenome.Hsapiens.UCSC.hg38"))
```

### Install from GitHub

Until the package is on CRAN, install the development version directly:

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("leopard0ly/gwas2crispr")
```

After CRAN release you will be able to run:

```r
install.packages("gwas2crispr")
```

---

## Quick start (primary workflow)

Use a clear prefix and write outputs (CSV/BED/FASTA) to your **current working directory**:

```r
library(gwas2crispr)

run_gwas2crispr(
  efo_id    = "EFO_0000707",  # lung disease (example)
  p_cut     = 1e-6,
  flank_bp  = 300,
  out_prefix = "lung"         # produces: lung_snps_full.csv / lung_snps_hg38.bed / lung_snps_flank300.fa
)
```

**Outputs**

* `lung_snps_full.csv` — harmonised SNP metadata from the GWAS Catalog (GRCh38).
* `lung_snps_hg38.bed` — intervals suitable for genomic intersection.
* `lung_snps_flank300.fa` — sequences for CRISPR guide design (**requires** `Biostrings` + `BSgenome.Hsapiens.UCSC.hg38`).

### A) Object‑only (no files written)

```r
library(gwas2crispr)

res <- run_gwas2crispr(
  efo_id     = "EFO_0001663",  # Prostate cancer
  p_cut      = 5e-8,
  flank_bp   = 200,
  out_prefix = NULL            # <- no writing; returns objects only
)

res$summary   # one‑row tibble: n_SNPs, SNPs_w_gene, unique_genes, n_studies
res$chr_freq  # table of chromosomes by SNP count
```

### B) Write files to a safe temporary directory (secondary)

```r
out <- file.path(tempdir(), "prostate")  # CRAN‑friendly
res <- run_gwas2crispr(
  efo_id     = "EFO_0001663",
  p_cut      = 5e-8,
  flank_bp   = 200,
  out_prefix = out
)

res$csv    # path to <prefix>_snps_full.csv
res$bed    # path to <prefix>_snps_hg38.bed
res$fasta  # path to <prefix>_snps_flank<bp>.fa (only if BSgenome installed)
```

**Output file names**

* `<prefix>_snps_full.csv` — unified metadata table
* `<prefix>_snps_hg38.bed` — BED intervals
* `<prefix>_snps_flank<bp>.fa` — FASTA sequences (requires `BSgenome.Hsapiens.UCSC.hg38`)

---

## Command‑line interface (CLI)

A portable Rscript is installed in the package under `inst/scripts/gwas2crispr.R`.  Use it to run the pipeline from the shell.  The script relies on the `optparse` package; install it if missing.

### Windows (Command Prompt)

**Version-agnostic (recommended):**

```bat
Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))" ^
  | Rscript -- -e EFO_0001663 -p 5e-8 -f 200 -o "%CD%\prostate"
```

**Fixed output (current folder):**

```bat
"C:\Program Files\R\R-4.4.1\bin\Rscript.exe" ^
  "C:\Users\ZAD ECT\AppData\Local\R\win-library\4.4\gwas2crispr\scripts\gwas2crispr.R" ^
  -e EFO_0001663 -p 5e-8 -f 200 -o "%CD%\prostate"
```

**Temporary output (system temp):**

```bat
Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))" ^
  | Rscript -- -e EFO_0001663 -p 5e-8 -f 200 -o "%TEMP%\prostate"
```

### Linux/macOS (Bash)

**Fixed output (current folder):**

```bash
Rscript "$(Rscript -e 'cat(system.file("scripts","gwas2crispr.R", package="gwas2crispr"))')" -e EFO_0001663 -p 5e-8 -f 200 -o "$PWD/prostate"
```

**Temporary output (system temp):**

```bash
Rscript "$(Rscript -e 'cat(system.file("scripts","gwas2crispr.R", package="gwas2crispr"))')" -e EFO_0001663 -p 5e-8 -f 200 -o "$(mktemp -d)/prostate"
```

### Windows note (PATH issues):

On some Windows systems, `Rscript` may not be recognized directly in the shell.
In this case, either use the full path to the executable:

```cmd
"C:\Program Files\R\R-4.4.3\bin\Rscript.exe" "C:\Gwas2Crispr\inst\scripts\gwas2crispr.R" -e EFO_0001663 -p 5e-8 -f 200 -o prostate
```

Or add R’s `bin` directory (e.g., `C:\Program Files\R\R-4.4.3\bin`) to your system `PATH` variable,
so that `Rscript` can be invoked directly:

```cmd
Rscript "C:\Gwas2Crispr\inst\scripts\gwas2crispr.R" -e EFO_0001663 -p 5e-8 -f 200 -o prostate
```

This ensures the CLI runs smoothly across Linux, macOS, and Windows.

### Options

* `-e, --efo` (required) — EFO trait ID, e.g. `EFO_0001663`
* `-p, --pthresh` — P‑value cut‑off (default `5e-8`)
* `-f, --flank` — number of flanking bases for FASTA (default `200`)
* `-o, --out` — output file prefix *(optional; omit to run object‑only without writing files)*
* `-v, --verbose` — print progress messages and, when `--out` is omitted, a concise summary

If you omit the `-o/--out` option, no files are written.  Use `-v/--verbose` to emit a concise summary of the run.

---

## Function reference

### `fetch_gwas(efo_id, p_cut = 5e-8)`

Fetch significant associations for an EFO trait.  Tries `gwasrapidd::get_associations()` first; if no rows or an error is returned, falls back to the EBI GWAS REST API.

* **Returns:** an S4 object of class `"associations"` with slots `associations` and `risk_alleles` (compatible with `gwasrapidd`).
* **Notes:** performs network requests and may be rate‑limited.

### `run_gwas2crispr(efo_id, p_cut = 5e-8, flank_bp = 200, out_prefix = NULL)`

Runs the full pipeline: fetches GWAS data, merges gene and study annotations, and returns a list with `summary` and `chr_freq`.  When `out_prefix` is provided, the list also contains file paths to the written `csv`, `bed` and optional `fasta` files.

* **Genome build:** GRCh38/hg38 (requires `BSgenome.Hsapiens.UCSC.hg38` for FASTA export)
* **Return value:** list with components `summary`, `chr_freq` and, if writing, `csv`, `bed`, `fasta` paths.

---

## Reproducibility & file layout

* Large outputs are **not** bundled in the package tarball; they are excluded via `.Rbuildignore`.
* Small example files (if needed) should live under `inst/extdata/` and can be accessed with:

  ```r
  system.file("extdata", "your_example.csv", package = "gwas2crispr")
  ```

---

## Testing

Automated tests live in `tests/testthat/` and avoid network calls on CRAN via `skip_on_cran()`.  To run the test suite locally:

```r
devtools::test()
```

---

## Notes on resources

* **FASTA export** is optional.  If `BSgenome.Hsapiens.UCSC.hg38` is not installed, the FASTA step is skipped gracefully; CSV and BED files will still be produced.
* **Required runtime client:** `gwasrapidd` must be installed to ensure smooth data retrieval from the GWAS Catalog.
* Informational output uses `message()` so that you can silence it with `suppressMessages()` when running scripts or examples.

---

## Citation

Please cite `gwas2crispr` and the resources it builds upon.  To see the formatted citation:

```r
citation("gwas2crispr")
```

Additional background: Sudlow et al. (2015) *UK Biobank: An open access resource for identifying the causes of a wide range of complex diseases of middle and old age* [doi:10.1093/nar/gkv1256](https://doi.org/10.1093/nar/gkv1256).

---

## Getting help

* Report issues or request features at [https://github.com/leopard0ly/gwas2crispr/issues](https://github.com/leopard0ly/gwas2crispr/issues).

## License

MIT © Othman S. I. Mohammed — see the [`LICENSE`](LICENSE) file for details.

## Acknowledgments

This package builds upon **gwasrapidd** and the **EBI GWAS** REST API.  Sequence handling and genome data are powered by **Biostrings** and **BSgenome**.
