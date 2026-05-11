# gwas2crispr

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16878244.svg)](https://doi.org/10.5281/zenodo.16878244)

> **GWAS-to-CRISPR**: direct GWAS Catalog REST API v2 retrieval and GRCh38/hg38 CSV, BED, and optional FASTA preparation for downstream CRISPR guide-design workflows.

## Overview

`gwas2crispr` retrieves significant genome-wide association study (**GWAS**) associations for an **Experimental Factor Ontology** (**EFO**) trait directly from the EMBL-EBI GWAS Catalog REST API v2.

The package prepares:

- harmonised SNP metadata in CSV format
- genomic intervals in BED format
- optional FASTA sequences with user-defined flanking regions

All genomic outputs are prepared for GRCh38/hg38.

The package is a computational preparation workflow. It does not perform wet-lab validation, therapeutic interpretation, or biological efficacy testing.

## Core functions

- `fetch_gwas(efo_id, p_cut = 5e-8, verbose = interactive())`
- `run_gwas2crispr(efo_id, p_cut = 5e-8, flank_bp = 200, out_prefix = NULL, verbose = interactive())`

## Installation

### Requirements

- R >= 4.1
- Direct internet access to the EMBL-EBI GWAS Catalog REST API v2
- Core R packages listed in `DESCRIPTION`

### Optional FASTA requirements

FASTA extraction requires:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Biostrings",
  "GenomeInfoDb",
  "BSgenome.Hsapiens.UCSC.hg38"
))
```

CSV and BED outputs can still be produced without the optional FASTA packages.

### Install from GitHub

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("leopard0ly/gwas2crispr")
```

## Quick start

```r
library(gwas2crispr)

res <- run_gwas2crispr(
  efo_id     = "EFO_0000707",
  p_cut      = 1e-6,
  flank_bp   = 300,
  out_prefix = "lung",
  verbose    = TRUE
)

res$summary
res$written
```

Expected output files:

- `lung_snps_full.csv`
- `lung_snps_hg38.bed`
- `lung_snps_flank300.fa`

The FASTA file is written only when the hg38 BSgenome and Biostrings packages are installed.

## Thesis case-study example

```r
library(gwas2crispr)

res <- run_gwas2crispr(
  efo_id     = "EFO_0001663",
  p_cut      = 5e-8,
  flank_bp   = 200,
  out_prefix = "prostate",
  verbose    = TRUE
)

res$summary
res$written
```

## Object-only mode

No files are written when `out_prefix = NULL`.

```r
res <- run_gwas2crispr(
  efo_id     = "EFO_0001663",
  p_cut      = 5e-8,
  flank_bp   = 200,
  out_prefix = NULL,
  verbose    = FALSE
)

res$summary
res$bed
```

## Output files

When `out_prefix` is supplied, the package writes:

- `<prefix>_snps_full.csv`
- `<prefix>_snps_hg38.bed`
- `<prefix>_snps_flank<bp>.fa`

## Notes

- Genome build is fixed to GRCh38/hg38.
- GWAS retrieval uses the EMBL-EBI GWAS Catalog REST API v2 directly.
- Results may change when GWAS Catalog content is updated.
- Network availability and rate limits may affect retrieval.
- FASTA export is optional.
- CSV and BED preparation remain available without optional genome packages.
- The package prepares computational outputs for downstream CRISPR guide-design workflows only.

## Command-line interface

A portable script is available under:

```text
inst/scripts/gwas2crispr.R
```

Example:

```bash
Rscript inst/scripts/gwas2crispr.R -e EFO_0001663 -p 5e-8 -f 200 -o prostate -v
```

Options:

- `-e, --efo` — EFO trait ID, for example `EFO_0001663`
- `-p, --pthresh` — p-value threshold
- `-f, --flank` — number of flanking bases for FASTA extraction
- `-o, --out` — output file prefix
- `-v, --verbose` — print progress messages

## Testing

```r
devtools::test()
```

Network-dependent tests are skipped on CRAN.

## Citation

If you use `gwas2crispr`, cite the Zenodo release:

[https://doi.org/10.5281/zenodo.16878244](https://doi.org/10.5281/zenodo.16878244)

```r
citation("gwas2crispr")
```

## Getting help

Report issues at:

[https://github.com/leopard0ly/gwas2crispr/issues](https://github.com/leopard0ly/gwas2crispr/issues)

## License

MIT © Othman S. I. Mohammed — see the `LICENSE` file.
