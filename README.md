# gwas2crispr

> GWAS‑to‑CRISPR data pipeline for high‑throughput SNP target extraction (GRCh38/hg38)

<!-- Optional badges (safe for GitHub; harmless in CRAN tarball) -->

<!--
[![R-CMD-check](https://github.com/leopard0ly/gwas2crispr/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/leopard0ly/gwas2crispr/actions/workflows/R-CMD-check.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](#license)
-->

## Overview

`gwas2crispr` retrieves significant GWAS SNPs for an EFO trait, aggregates variant/gene/study metadata, and exports ready‑to‑use **CSV**, **BED**, and **FASTA** files suitable for downstream functional genomics and CRISPR guide design tools. The package targets **GRCh38/hg38**.

Core functions:

* **`fetch_gwas(efo_id, p_cut)`** — fetch significant associations via `gwasrapidd` with a REST fallback.
* **`run_gwas2crispr(efo_id, p_cut, flank_bp, out_prefix)`** — end‑to‑end pipeline that writes CSV/BED/FASTA.

> **Networking & CRAN**: Examples in the docs are wrapped in `\dontrun{}` because the functions perform network requests and write files. This keeps CRAN checks clean.

---

## Installation

### Requirements

* **R >= 4.1**
* CRAN: `dplyr`, `httr`, `purrr`, `readr`, `tibble`, `tidyr`, `methods`, `utils`
* Bioconductor: `Biostrings`, `BSgenome.Hsapiens.UCSC.hg38`
* (Optional, CLI) `optparse`

### Install dependencies

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

```r
library(gwas2crispr)

# Prostate cancer, genome build GRCh38/hg38
run_gwas2crispr(
  efo_id     = "EFO_0001663",
  p_cut      = 5e-8,
  flank_bp   = 200,
  out_prefix = "prostate"
)
```

This will create the following files in your working directory:

* `prostate_snps_full.csv` — unified metadata table.
* `prostate_snps_hg38.bed` — BED intervals for genomic intersection.
* `prostate_snps_flank200.fa` — FASTA sequences with ±200 bp around each SNP (requires `BSgenome.Hsapiens.UCSC.hg38`).

---

## Command‑line interface (CLI)

A portable Rscript is installed under the package `scripts/` directory.

**Windows (Command Prompt):**

```bat
"C:\\Program Files\\R\\R-4.4.1\\bin\\Rscript.exe" ^
  "%USERPROFILE%\\AppData\\Local\\R\\win-library\\4.4\\gwas2crispr\\scripts\\gwas2crispr.R" ^
  -e EFO_0001663 -p 5e-8 -f 200 -o prostate
```

**Linux/macOS (Bash):**

```bash
Rscript "$(Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))")" \
  -e EFO_0001663 -p 5e-8 -f 200 -o prostate
```

Options:

* `-e, --efo` (required) — EFO trait ID, e.g. `EFO_0001663`.
* `-p, --pthresh` — P‑value cutoff (default `5e-8`).
* `-f, --flank` — Flanking bases for FASTA (default `200`).
* `-o, --out` — Output file prefix (default `output`).

---

## Function reference

### `fetch_gwas(efo_id, p_cut = 5e-8)`

Fetch significant associations for an EFO trait. Tries `gwasrapidd::get_associations()` first; if it returns no rows or errors, falls back to the EBI GWAS REST API.

* **Returns**: S4 object of class `"associations"` with slots `associations` and `risk_alleles` (compatible with `gwasrapidd`).
* **Notes**: Performs network I/O; may be rate‑limited.

### `run_gwas2crispr(efo_id, p_cut = 5e-8, flank_bp = 200, out_prefix = "output")`

End‑to‑end helper that writes CSV/BED/FASTA and returns a list of output paths and summary tables.

* **Genome**: GRCh38/hg38 only (`BSgenome.Hsapiens.UCSC.hg38`).
* **Output**: `*_snps_full.csv`, `*_snps_hg38.bed`, `*_snps_flank<bp>.fa`.

---

## Reproducibility & file layout

* Large example outputs produced by the functions are **not** shipped in the package build; they are ignored via `.Rbuildignore`.
* If you wish to ship small example files, place them under `inst/extdata/` and access them with:

```r
system.file("extdata", "your_example.csv", package = "gwas2crispr")
```

---

## Testing

Minimal tests are provided under `tests/testthat/` and avoid network on CRAN via `skip_on_cran()`. Locally you can run:

```r
devtools::test()
```

---

## Roadmap

* Add GRCh37 support via an optional BSgenome package.
* Expand unit tests and add mocks for network calls.
* Provide precomputed example datasets under `inst/extdata/`.

---

## Getting help

* **Bug reports & feature requests**: [https://github.com/leopard0ly/gwas2crispr/issues](https://github.com/leopard0ly/gwas2crispr/issues)
* **How to cite**: see [`inst/CITATION`](inst/CITATION) (also available via `citation("gwas2crispr")`).

---

## License

MIT © Othman S. I. Mohammed — see [`LICENSE`](LICENSE).

---

## Acknowledgments

* Built on top of the excellent GWAS infrastructure provided by **gwasrapidd** and the EBI GWAS REST API.
* Sequence handling powered by **Biostrings** and **BSgenome** (hg38).
