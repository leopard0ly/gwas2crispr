# gwas2crispr

GWAS-to-CRISPR: streamlined retrieval of significant GWAS Catalog associations, metadata aggregation, and optional CSV/BED/FASTA export for downstream CRISPR guide-design workflows using GRCh38/hg38.

## Overview

Genome-wide association studies (GWAS) link traits, diseases, and phenotypes to genetic variants, but raw GWAS Catalog association records are not directly usable for downstream CRISPR guide-design preparation.

`gwas2crispr` bridges this gap. It retrieves significant GWAS Catalog associations for a supported trait identifier, aggregates variant, gene, and study metadata, and returns in-memory summaries. When requested, it also writes ready-to-use CSV, BED, and optional FASTA files for downstream CRISPR target-design preparation.

All genomic coordinates are prepared for GRCh38/hg38.

The package is a computational preparation workflow. It does not perform wet-lab validation, therapeutic interpretation, biological causality testing, biological efficacy testing, guide scoring, or off-target prediction.

## Core functions

- `fetch_gwas(efo_id, p_cut = 5e-8, verbose = interactive())`: fetches significant GWAS Catalog associations for a supported trait identifier.
- `run_gwas2crispr(efo_id, p_cut = 5e-8, flank_bp = 200, out_prefix = NULL, verbose = interactive())`: end-to-end workflow that calls `fetch_gwas()`, aggregates variant/gene/study metadata, creates GRCh38/hg38 intervals, and returns an object with summaries. If `out_prefix` is supplied, it also writes CSV, BED, and optional FASTA files.

The argument name `efo_id` is retained for backward compatibility. Starting from `gwas2crispr` 0.1.5, selected non-EFO GWAS Catalog trait identifiers are accepted through the same argument when supported by the GWAS Catalog API.

CRAN-safe examples: the package does not write files by default in examples. Examples that perform network operations or file writing should use `tempdir()` or user-defined output paths.

## Supported trait identifiers

### Officially supported in 0.1.5

- `EFO`
- `MONDO`
- `NCIT`

### Compatibility accepted in 0.1.5

- `HP`
- `Orphanet`
- `ORPHA`

### Not supported as primary GWAS Catalog trait identifiers in 0.1.5

- `GO`

Accepted input forms include both underscore and colon syntax:

```text
EFO_0000000
EFO:0000000

MONDO_0000000
MONDO:0000000

NCIT_C0000
NCIT:C0000

HP_0000000
HP:0000000

Orphanet_0000
Orphanet:0000

ORPHA_0000
ORPHA:0000
```

Colon syntax is normalized internally to underscore syntax.

Data availability depends on the GWAS Catalog API. Different identifier systems may return different association sets, even when disease or phenotype concepts are related. Users should verify returned studies and trait context before interpretation.

## Installation

### Requirements

- R >= 4.1
- Direct internet access to the EMBL-EBI GWAS Catalog REST API v2
- Optional internet access to the Ensembl REST API for non-fatal rsID coordinate recovery
- Core CRAN stack listed in `DESCRIPTION`
- FASTA output requires optional Bioconductor sequence packages
- Optional for CLI: `optparse`

### Install from CRAN

```r
install.packages("gwas2crispr")
```

### Install Bioconductor dependencies for FASTA output

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

## Quick start

Use a supported GWAS Catalog cancer trait identifier, choose a p-value threshold, choose a flank size, and set an output prefix.

```r
library(gwas2crispr)

trait_id <- "MONDO_0007254"  # breast cancer
```

### Option 1: run directly and write files

Use this form when you only want the output files.

```r
run_gwas2crispr(
  efo_id     = trait_id,
  p_cut      = 1e-6,
  flank_bp   = 300,
  out_prefix = "breast_cancer_run",
  verbose    = TRUE
)
```

Expected output files:

```text
breast_cancer_run_snps_full.csv
breast_cancer_run_snps_hg38.bed
breast_cancer_run_snps_flank300.fa
```

The FASTA file is written only when the optional hg38 sequence packages are installed.

### Option 2: save the returned object and write files

Use this form when you want the output files and also want to inspect the returned R object.

```r
res <- run_gwas2crispr(
  efo_id     = trait_id,
  p_cut      = 1e-6,
  flank_bp   = 300,
  out_prefix = "breast_cancer_run",
  verbose    = TRUE
)

res$summary
res$written
```

`res <-` is not required for file writing. It is used only when the user wants to inspect the returned summary, tables, and written file paths inside R.

## Cancer trait identifier examples

Use the same workflow with any supported GWAS Catalog cancer-related identifier.

```r
trait_id <- "MONDO_0007254"  # breast cancer
trait_id <- "NCIT_C4872"     # breast carcinoma / breast cancer terminology
trait_id <- "EFO_0001663"    # prostate cancer
```

Colon syntax can also be used:

```r
trait_id <- "MONDO:0007254"  # breast cancer
trait_id <- "NCIT:C4872"     # breast carcinoma / breast cancer terminology
trait_id <- "EFO:0001663"    # prostate cancer
```

Then run either Option 1 or Option 2 above.

The identifier determines which GWAS Catalog records are retrieved. Results may differ across identifier systems because GWAS Catalog annotations and ontology mappings are not always equivalent.

## A) Object-only mode: no files written

Set `out_prefix = NULL` to return R objects only.

```r
trait_id <- "MONDO_0007254"  # breast cancer

res <- run_gwas2crispr(
  efo_id     = trait_id,
  p_cut      = 1e-6,
  flank_bp   = 300,
  out_prefix = NULL,
  verbose    = FALSE
)

res$summary
res$snps_full
res$bed
```

## B) Write files to a safe temporary directory

Use `tempdir()` when you want written files without cluttering the working directory.

```r
trait_id <- "NCIT_C4872"  # breast carcinoma / breast cancer terminology
out <- file.path(tempdir(), "breast_carcinoma_run")

res <- run_gwas2crispr(
  efo_id     = trait_id,
  p_cut      = 1e-6,
  flank_bp   = 300,
  out_prefix = out,
  verbose    = TRUE
)

res$summary
res$written
```

Expected output files:

```text
<tempdir>/breast_carcinoma_run_snps_full.csv
<tempdir>/breast_carcinoma_run_snps_hg38.bed
<tempdir>/breast_carcinoma_run_snps_flank300.fa
```

The FASTA file is written only when the optional hg38 sequence packages are installed.

## C) Fetch GWAS records only

Use `fetch_gwas()` when only GWAS Catalog retrieval is needed.

```r
trait_id <- "EFO_0001663"  # prostate cancer

gwas <- fetch_gwas(
  efo_id  = trait_id,
  p_cut   = 5e-8,
  verbose = TRUE
)

names(gwas)
gwas$associations
```

`fetch_gwas()` returns:

```text
associations
risk_alleles
cache
```

## Returned object

`run_gwas2crispr()` returns a list containing:

```text
summary
chr_freq
snps_full
bed
fasta
written
```

Example:

```r
names(res)

res$summary
res$chr_freq
head(res$snps_full)
head(res$bed)
res$fasta
res$written
```

## Output files

When `out_prefix` is supplied, the package writes:

```text
<prefix>_snps_full.csv
<prefix>_snps_hg38.bed
<prefix>_snps_flank<bp>.fa
```

Example with `out_prefix = "breast_cancer_run"` and `flank_bp = 300`:

```text
breast_cancer_run_snps_full.csv
breast_cancer_run_snps_hg38.bed
breast_cancer_run_snps_flank300.fa
```

### CSV output

The CSV file contains harmonised SNP and association metadata.

Typical columns include:

```text
variant_id
chromosome_name
chromosome_position
genes
association_id
pvalue
study_accession
```

### BED output

The BED file contains genomic intervals around each variant under GRCh38/hg38.

BED output is intended for genomic interval operations and downstream CRISPR guide-design preparation.

### FASTA output

The FASTA file contains sequence windows around each variant using the selected flank size.

FASTA output is optional and depends on the availability of the required hg38 sequence packages.

## Inspecting outputs

After a file-writing run:

```r
csv <- read.csv("breast_cancer_run_snps_full.csv")
bed <- read.delim("breast_cancer_run_snps_hg38.bed", header = FALSE)

dim(csv)
head(csv)

dim(bed)
head(bed)
```

If FASTA was produced:

```r
fa <- readLines("breast_cancer_run_snps_flank300.fa")
sum(grepl("^>", fa))
head(fa)
```

## Coordinate recovery

`gwas2crispr` uses coordinates provided by GWAS Catalog records when available.

When some rsID records lack complete coordinates, `gwas2crispr` may attempt non-fatal coordinate recovery through additional metadata routes, including an optional Ensembl REST fallback.

If coordinate recovery is unavailable or incomplete, unresolved variants may be skipped from coordinate-based BED and FASTA outputs. CSV metadata preparation can still proceed when association records are available.

## Command-line interface

A portable command-line interface script is available under:

```text
inst/scripts/gwas2crispr.R
```

After package installation, the script can be located from R using:

```r
system.file("scripts", "gwas2crispr.R", package = "gwas2crispr")
```

### CLI options

```text
-e, --efo      GWAS Catalog trait identifier. The option name is retained for backward compatibility.
-p, --pthresh  p-value threshold, for example 5e-8.
-f, --flank    number of flanking bases for FASTA extraction.
-o, --out      output file prefix.
-v, --verbose  print progress messages.
```

The `--efo` option accepts selected supported identifiers, including EFO, MONDO, and NCIT identifiers, when supported by the GWAS Catalog API.

### Linux and macOS

If running from a cloned GitHub source folder:

```bash
# breast cancer
Rscript inst/scripts/gwas2crispr.R \
  -e MONDO_0007254 \
  -p 1e-6 \
  -f 300 \
  -o breast_cancer_run \
  -v
```

If the package is already installed:

```bash
SCRIPT=$(Rscript -e "cat(system.file('scripts', 'gwas2crispr.R', package = 'gwas2crispr'))")

# breast cancer
Rscript "$SCRIPT" \
  -e MONDO_0007254 \
  -p 1e-6 \
  -f 300 \
  -o breast_cancer_run \
  -v
```

### Windows CMD

If running from a cloned GitHub source folder:

```cmd
REM breast cancer
Rscript inst\scripts\gwas2crispr.R -e MONDO_0007254 -p 1e-6 -f 300 -o breast_cancer_run -v
```

If the package is already installed and `Rscript` is available in the Windows PATH:

```cmd
REM breast cancer
for /f "delims=" %i in ('Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))"') do Rscript "%i" -e MONDO_0007254 -p 1e-6 -f 300 -o breast_cancer_run -v
```

If `Rscript` is not available in the Windows PATH, use the full path to `Rscript.exe`.

Example:

```cmd
REM breast cancer
"C:\Program Files\R\R-4.4.3\bin\x64\Rscript.exe" "C:/Users/hp/AppData/Local/R/win-library/4.4/gwas2crispr/scripts/gwas2crispr.R" -e MONDO_0007254 -p 1e-6 -f 300 -o breast_cancer_run -v
```

### Windows PowerShell

If the package is already installed and `Rscript` is available in the PATH:

```powershell
$script = Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))"

# breast cancer
Rscript $script -e MONDO_0007254 -p 1e-6 -f 300 -o breast_cancer_run -v
```

If `Rscript` is not available in the PATH, replace `Rscript` with the full path to `Rscript.exe`.

Example:

```powershell
$Rscript = "C:\Program Files\R\R-4.4.3\bin\x64\Rscript.exe"
$script = & $Rscript -e "cat(system.file('scripts','gwas2crispr.R', package='gwas2crispr'))"

# breast cancer
& $Rscript $script -e MONDO_0007254 -p 1e-6 -f 300 -o breast_cancer_run -v
```

### Expected CLI output files

For the example above, the CLI writes:

```text
breast_cancer_run_snps_full.csv
breast_cancer_run_snps_hg38.bed
breast_cancer_run_snps_flank300.fa
```

The FASTA file is written only when the optional hg38 sequence packages are installed.

## Testing

```r
devtools::test()
```

Some local tests may take several minutes because they exercise live retrieval and fallback behavior.

Network-dependent tests are skipped on CRAN.

## Notes

- Genome build is fixed to GRCh38/hg38.
- GWAS retrieval uses the EMBL-EBI GWAS Catalog REST API v2 directly.
- Optional coordinate recovery may use Ensembl REST when GWAS Catalog records lack complete coordinates.
- Results may change when GWAS Catalog content is updated.
- Network availability and rate limits may affect retrieval.
- Different supported identifiers may return different association sets.
- FASTA export is optional.
- CSV and BED preparation remain available without optional genome packages.
- The package prepares computational outputs for downstream CRISPR guide-design workflows only.
- The package does not perform therapeutic interpretation, wet-lab validation, biological causality testing, guide scoring, off-target prediction, or biological efficacy testing.

## Citation

If you use `gwas2crispr`, cite the CRAN package:

```r
citation("gwas2crispr")
```

CRAN package DOI:

```text
https://doi.org/10.32614/CRAN.package.gwas2crispr
```

## Getting help

Report issues at:

```text
https://github.com/leopard0ly/gwas2crispr/issues
```

## License

MIT © Othman S. I. Mohammed — see the `LICENSE` file.
