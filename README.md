# gwas2crispr

> GWAS-to-CRISPR data pipeline for high-throughput SNP target extraction (GRCh38/hg38)

<!-- Optional badges (safe for GitHub; harmless in CRAN tarball) -->

<!--
[![R-CMD-check](https://github.com/leopard0ly/gwas2crispr/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/leopard0ly/gwas2crispr/actions/workflows/R-CMD-check.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](#license)
-->

## Overview

`gwas2crispr` retrieves significant **GWAS** SNPs for an **EFO** trait, aggregates variant / gene / study
metadata, and exports **CSV**, **BED**, and **FASTA** files suitable for downstream functional genomics and
CRISPR guide design. The package targets **GRCh38/hg38**.

**Core functions**

* `fetch_gwas(efo_id, p_cut)` — fetch significant associations via `gwasrapidd` with a REST fallback.
* `run_gwas2crispr(efo_id, p_cut, flank_bp, out_prefix = NULL)` — end-to-end pipeline; returns in-memory summaries and (optionally) writes CSV/BED/FASTA if you supply an output prefix.

> **CRAN-safe**: Networked / file-writing examples in the package docs use `\donttest{}`. The function does **not** write files by default; you opt-in by providing `out_prefix`.

---

## Installation

### Requirements

* **R ≥ 4.1**
* CRAN: `dplyr`, `httr`, `purrr`, `readr`, `tibble`, `tidyr`, `methods`, `utils`
* Bioconductor: `Biostrings`, `BSgenome.Hsapiens.UCSC.hg38` *(FASTA export is optional; if not installed, the step is skipped gracefully)*
* (Optional, CLI) `optparse`

### Install Bioconductor dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "BSgenome.Hsapiens.UCSC.hg38"))
```

### Install the development version from GitHub

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("leopard0ly/gwas2crispr")
```

---

## Quick start

### A) Object-only (no files written)

```r
library(gwas2crispr)

res <- run_gwas2crispr(
  efo_id     = "EFO_0001663",  # Prostate cancer
  p_cut      = 5e-8,
  flank_bp   = 200,
  out_prefix = NULL            # <- no writing; returns objects
)

res$summary   # one-row tibble: n_SNPs, SNPs_w_gene, unique_genes, n_studies
res$chr_freq  # top chromosomes by SNP count
```

### B) Write files to a safe temp directory

```r
out <- file.path(tempdir(), "prostate")  # CRAN-friendly
res <- run_gwas2crispr(
  efo_id     = "EFO_0001663",
  p_cut      = 5e-8,
  flank_bp   = 200,
  out_prefix = out
)

res$csv
res$bed
res$fasta   # exists only if BSgenome.Hsapiens.UCSC.hg38 is installed
```

**Outputs**

* `<prefix>_snps_full.csv` — unified metadata table
* `<prefix>_snps_hg38.bed` — BED intervals
* `<prefix>_snps_flank<bp>.fa` — FASTA sequences (optional; requires `BSgenome.Hsapiens.UCSC.hg38`)

---

## Command-line interface (CLI)

A portable Rscript is installed under the package `scripts/` directory.

**Windows (Command Prompt):**

```bat
Rscript ^
  "%USERPROFILE%\AppData\Local\R\win-library\%R_MAJOR%.%R_MINOR%\gwas2crispr\scripts\gwas2crispr.R" ^
  -e EFO_0001663 -p 5e-8 -f 200 -o "%TEMP%\prostate"
```

**Linux/macOS (Bash):**

```bash
Rscript "$(Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))")" \
  -e EFO_0001663 -p 5e-8 -f 200 -o "$(mktemp -d)/prostate"
```

**Options**

* `-e, --efo` (required) — EFO trait ID, e.g. `EFO_0001663`
* `-p, --pthresh` — P-value cutoff (default `5e-8`)
* `-f, --flank` — Flanking bases for FASTA (default `200`)
* `-o, --out` — Output file prefix *(required to write files; omit to run object-only)*

---

## Function reference

### `fetch_gwas(efo_id, p_cut = 5e-8)`

Fetch significant associations for an EFO trait. Tries `gwasrapidd::get_associations()` first; if it returns no rows or errors, falls back to the EBI GWAS REST API.

* **Returns**: S4 object of class `"associations"` with slots `associations` and `risk_alleles` (compatible with `gwasrapidd`).
* **Notes**: Performs network I/O and may be rate-limited.

### `run_gwas2crispr(efo_id, p_cut = 5e-8, flank_bp = 200, out_prefix = NULL)`

End-to-end helper that aggregates variant / gene / study metadata, returns an object with summaries, and **optionally** writes CSV/BED/FASTA if `out_prefix` is provided.

* **Genome**: GRCh38/hg38 (`BSgenome.Hsapiens.UCSC.hg38`)
* **Return value**: list with `summary`, `chr_freq`, and, when writing is enabled, `csv`, `bed`, `fasta` paths.

---

## Reproducibility & file layout

* Large outputs are **not** shipped; they’re excluded via `.Rbuildignore`.
* Ship small examples (if needed) under `inst/extdata/` and access them with:

```r
system.file("extdata", "your_example.csv", package = "gwas2crispr")
```

---

## Testing

Minimal tests live in `tests/testthat/` and avoid network on CRAN via `skip_on_cran()`. Locally:

```r
devtools::test()
```

---

## Notes on resources

* **FASTA export** is optional and uses `BSgenome.Hsapiens.UCSC.hg38` (Bioconductor). If it’s not installed, the FASTA step is **skipped gracefully**; CSV and BED still produce.
* Informational output uses `message()` so you can silence it with `suppressMessages(...)`.

---

## Citation

Please cite `gwas2crispr` and the underlying resources. See:

```r
citation("gwas2crispr")
```

Background reference: Sudlow et al. (2015) [doi:10.1093/nar/gkv1256](doi:10.1093/nar/gkv1256).

---

## Getting help

* Issues & feature requests: [https://github.com/leopard0ly/gwas2crispr/issues](https://github.com/leopard0ly/gwas2crispr/issues)

---

## License

MIT © Othman S. I. Mohammed — see [`LICENSE`](LICENSE).

---

## Acknowledgments

* Built atop **gwasrapidd** and the **EBI GWAS** REST API.
* Sequence handling powered by **Biostrings** and **BSgenome** (hg38).
