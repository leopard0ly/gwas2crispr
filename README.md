# gwas2crispr

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20129602.svg)](https://doi.org/10.5281/zenodo.20129602)
[![CRAN status](https://www.r-pkg.org/badges/version/gwas2crispr)](https://CRAN.R-project.org/package=gwas2crispr)

> **GWAS-to-CRISPR**: direct EMBL-EBI GWAS Catalog REST API v2 retrieval and GRCh38/hg38 CSV, BED, and optional FASTA preparation for downstream CRISPR guide-design workflows.

## Overview

`gwas2crispr` retrieves significant genome-wide association study (**GWAS**) associations for an **Experimental Factor Ontology** (**EFO**) trait directly from the EMBL-EBI GWAS Catalog REST API v2.

The package prepares:

- harmonised SNP metadata in CSV format
- genomic intervals in BED format
- optional FASTA sequences with user-defined flanking regions

All genomic outputs are prepared for **GRCh38/hg38**.

The package is a computational preparation workflow. It does not perform wet-lab validation, therapeutic interpretation, biological causality testing, or biological efficacy testing.

## Core functions

```r
fetch_gwas(
  efo_id,
  p_cut = 5e-8,
  verbose = interactive()
)

run_gwas2crispr(
  efo_id,
  p_cut = 5e-8,
  flank_bp = 200,
  out_prefix = NULL,
  verbose = interactive()
)
```

## Installation

### Requirements

- R >= 4.1
- Direct internet access to the EMBL-EBI GWAS Catalog REST API v2
- Core R packages listed in `DESCRIPTION`

### Install from CRAN

```r
install.packages("gwas2crispr")
```

### Install from GitHub

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("leopard0ly/gwas2crispr")
```

## Optional FASTA requirements

FASTA extraction requires Bioconductor sequence packages:

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

```text
lung_snps_full.csv
lung_snps_hg38.bed
lung_snps_flank300.fa
```

The FASTA file is written only when the hg38 BSgenome and Biostrings packages are installed.

## Thesis case-study example

The prostate cancer case study can be reproduced with:

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

Expected output files:

```text
prostate_snps_full.csv
prostate_snps_hg38.bed
prostate_snps_flank200.fa
```

## Object-only mode

No files are written when `out_prefix = NULL`.

```r
library(gwas2crispr)

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

```text
<prefix>_snps_full.csv
<prefix>_snps_hg38.bed
<prefix>_snps_flank<bp>.fa
```

Example with `out_prefix = "prostate"` and `flank_bp = 200`:

```text
prostate_snps_full.csv
prostate_snps_hg38.bed
prostate_snps_flank200.fa
```

## Command-line interface

`gwas2crispr` includes a portable command-line interface (**CLI**) for running the same GWAS-to-CRISPR preparation workflow outside an interactive R session.

The CLI script is stored in the package source tree at:

```text
inst/scripts/gwas2crispr.R
```

After package installation, the script can be located from R using:

```r
system.file("scripts", "gwas2crispr.R", package = "gwas2crispr")
```

### CLI options

- `-e, --efo` — EFO trait ID, for example `EFO_0001663`
- `-p, --pthresh` — p-value threshold, for example `5e-8`
- `-f, --flank` — number of flanking bases for FASTA extraction
- `-o, --out` — output file prefix
- `-v, --verbose` — print progress messages

### Linux and macOS

If running from a cloned GitHub source folder:

```bash
Rscript inst/scripts/gwas2crispr.R \
  -e EFO_0001663 \
  -p 5e-8 \
  -f 200 \
  -o prostate \
  -v
```

If the package is already installed:

```bash
SCRIPT=$(Rscript -e "cat(system.file('scripts', 'gwas2crispr.R', package = 'gwas2crispr'))")

Rscript "$SCRIPT" \
  -e EFO_0001663 \
  -p 5e-8 \
  -f 200 \
  -o prostate \
  -v
```

### Windows CMD

If running from a cloned GitHub source folder:

```bat
Rscript inst\scripts\gwas2crispr.R -e EFO_0001663 -p 5e-8 -f 200 -o prostate -v
```

If the package is already installed, use `system.file()` to find the installed CLI script automatically:

```bat
for /f "delims=" %i in ('Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))"') do Rscript "%i" -e EFO_0001663 -p 5e-8 -f 200 -o prostate -v
```

If `Rscript` is not available in the Windows PATH, use the full path to `Rscript.exe`.

Example:

```bat
"C:\Program Files\R\R-4.4.3\bin\x64\Rscript.exe" "C:/Users/hp/AppData/Local/R/win-library/4.4/gwas2crispr/scripts/gwas2crispr.R" -e EFO_0001663 -p 5e-8 -f 200 -o prostate -v
```

### Windows PowerShell

If the package is already installed and `Rscript` is available in the PATH:

```powershell
$script = Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))"

Rscript $script -e EFO_0001663 -p 5e-8 -f 200 -o prostate -v
```

If `Rscript` is not available in the PATH, replace `Rscript` with the full path to `Rscript.exe`.

Example:

```powershell
$Rscript = "C:\Program Files\R\R-4.4.3\bin\x64\Rscript.exe"
$script = & $Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))"

& $Rscript $script -e EFO_0001663 -p 5e-8 -f 200 -o prostate -v
```

### Expected CLI output files

For the prostate cancer example above, the CLI writes:

```text
prostate_snps_full.csv
prostate_snps_hg38.bed
prostate_snps_flank200.fa
```

The FASTA file is written only when the optional hg38 sequence packages are installed.

## Testing

```r
devtools::test()
```

Network-dependent tests are skipped on CRAN.

## Notes

- Genome build is fixed to GRCh38/hg38.
- GWAS retrieval uses the EMBL-EBI GWAS Catalog REST API v2 directly.
- Results may change when GWAS Catalog content is updated.
- Network availability and rate limits may affect retrieval.
- FASTA export is optional.
- CSV and BED preparation remain available without optional genome packages.
- The package prepares computational outputs for downstream CRISPR guide-design workflows only.
- The package does not perform therapeutic interpretation, wet-lab validation, or biological efficacy testing.

## Citation

If you use `gwas2crispr`, cite the Zenodo release:

[https://doi.org/10.5281/zenodo.20129602](https://doi.org/10.5281/zenodo.20129602)

You can also run:

```r
citation("gwas2crispr")
```

## Getting help

Report issues at:

[https://github.com/leopard0ly/gwas2crispr/issues](https://github.com/leopard0ly/gwas2crispr/issues)

## License

MIT © Othman S. I. Mohammed — see the `LICENSE` file.
