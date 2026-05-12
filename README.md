# gwas2crispr

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20129602.svg)](https://doi.org/10.5281/zenodo.20129602)
[![CRAN status](https://www.r-pkg.org/badges/version/gwas2crispr)](https://CRAN.R-project.org/package=gwas2crispr)

> GWAS-to-CRISPR: streamlined retrieval of significant GWAS SNPs, metadata aggregation, and optional FASTA/BED/CSV export for downstream CRISPR guide-design workflows using GRCh38/hg38.

## Overview

Genome-wide association studies (GWAS) link traits to genetic variants, but raw GWAS Catalog association records are not directly usable for downstream CRISPR guide-design preparation.

`gwas2crispr` bridges this gap. It retrieves significant single-nucleotide polymorphisms (SNPs) for a given Experimental Factor Ontology (EFO) trait directly from the EMBL-EBI GWAS Catalog REST API v2, aggregates variant, gene, and study metadata, and returns in-memory summaries. When requested, it also writes ready-to-use CSV, BED, and optional FASTA files for high-throughput downstream CRISPR target-design preparation.

All genomic coordinates are prepared for GRCh38/hg38.

The package is a computational preparation workflow. It does not perform wet-lab validation, therapeutic interpretation, biological causality testing, or biological efficacy testing.

## Core functions

- `fetch_gwas(efo_id, p_cut = 5e-8, verbose = interactive())`: fetches significant associations for an EFO trait directly from the EMBL-EBI GWAS Catalog REST API v2.
- `run_gwas2crispr(efo_id, p_cut = 5e-8, flank_bp = 200, out_prefix = NULL, verbose = interactive())`: end-to-end pipeline that calls `fetch_gwas()`, aggregates variant/gene/study metadata, and returns an object with summaries. If you provide `out_prefix`, it also writes CSV, BED, and optional FASTA files.

> CRAN-safe examples: the package does not write files by default. Examples that perform network operations or file writing should use `tempdir()` or user-defined output paths.

---

## Installation

### Requirements (read first)

- R >= 4.1
- Direct internet access to the EMBL-EBI GWAS Catalog REST API v2
- Core CRAN stack listed in `DESCRIPTION`
- FASTA output requires optional Bioconductor sequence packages
- Optional for CLI: `optparse`

### Install from CRAN

```r
install.packages("gwas2crispr")
```

### Install Bioconductor dependencies (for FASTA)

FASTA output requires `Biostrings`, `GenomeInfoDb`, and `BSgenome.Hsapiens.UCSC.hg38`.

If these packages are missing, CSV and BED outputs are still produced, while FASTA is skipped.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Biostrings",
  "GenomeInfoDb",
  "BSgenome.Hsapiens.UCSC.hg38"
))
```

### Install from GitHub

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("leopard0ly/gwas2crispr")
```

---

## Quick start (primary workflow)

Use a clear prefix and write outputs to your current working directory:

```r
library(gwas2crispr)

run_gwas2crispr(
  efo_id     = "EFO_0000707",  # lung disease example
  p_cut      = 1e-6,
  flank_bp   = 300,
  out_prefix = "lung",
  verbose    = TRUE
)
```

Outputs:

- `lung_snps_full.csv` — harmonised SNP metadata from GWAS Catalog associations
- `lung_snps_hg38.bed` — GRCh38/hg38 intervals suitable for genomic intersection
- `lung_snps_flank300.fa` — sequence windows for downstream CRISPR guide-design preparation

The FASTA file is written only when the optional hg38 sequence packages are installed.

### A) Object-only (no files written)

```r
library(gwas2crispr)

res <- run_gwas2crispr(
  efo_id     = "EFO_0001663",  # prostate cancer
  p_cut      = 5e-8,
  flank_bp   = 200,
  out_prefix = NULL,           # no writing; returns objects only
  verbose    = FALSE
)

res$summary
res$bed
```

### B) Write files to a safe temporary directory (secondary)

```r
library(gwas2crispr)

out <- file.path(tempdir(), "prostate")

res <- run_gwas2crispr(
  efo_id     = "EFO_0001663",
  p_cut      = 5e-8,
  flank_bp   = 200,
  out_prefix = out,
  verbose    = TRUE
)

res$summary
res$written
```

Expected output files:

- `<tempdir>/prostate_snps_full.csv`
- `<tempdir>/prostate_snps_hg38.bed`
- `<tempdir>/prostate_snps_flank200.fa`

The FASTA file is written only when the optional hg38 sequence packages are installed.

---

## Output files

When `out_prefix` is supplied, the package writes:

- `<prefix>_snps_full.csv`
- `<prefix>_snps_hg38.bed`
- `<prefix>_snps_flank<bp>.fa`

Example with `out_prefix = "prostate"` and `flank_bp = 200`:

- `prostate_snps_full.csv`
- `prostate_snps_hg38.bed`
- `prostate_snps_flank200.fa`

---

## Command-line interface

A portable command-line interface (CLI) script is available under:

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

If the package is already installed and `Rscript` is available in the Windows PATH:

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

- `prostate_snps_full.csv`
- `prostate_snps_hg38.bed`
- `prostate_snps_flank200.fa`

The FASTA file is written only when the optional hg38 sequence packages are installed.

---

## Testing

```r
devtools::test()
```

Network-dependent tests are skipped on CRAN.

---

## Notes

- Genome build is fixed to GRCh38/hg38.
- GWAS retrieval uses the EMBL-EBI GWAS Catalog REST API v2 directly.
- Results may change when GWAS Catalog content is updated.
- Network availability and rate limits may affect retrieval.
- FASTA export is optional.
- CSV and BED preparation remain available without optional genome packages.
- The package prepares computational outputs for downstream CRISPR guide-design workflows only.
- The package does not perform therapeutic interpretation, wet-lab validation, or biological efficacy testing.

---

## Citation

If you use `gwas2crispr`, cite the Zenodo release:

[https://doi.org/10.5281/zenodo.20129602](https://doi.org/10.5281/zenodo.20129602)

You can also run:

```r
citation("gwas2crispr")
```

---

## Getting help

Report issues at:

[https://github.com/leopard0ly/gwas2crispr/issues](https://github.com/leopard0ly/gwas2crispr/issues)

---

## License

MIT © Othman S. I. Mohammed — see the `LICENSE` file.
