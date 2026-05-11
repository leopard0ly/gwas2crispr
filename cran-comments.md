## Submission (0.1.4)

This submission updates the GWAS Catalog retrieval backend and removes the former third-party GWAS retrieval dependency.

### Key changes since 0.1.2

- **Direct GWAS Catalog REST API v2 retrieval:** `fetch_gwas()` now retrieves GWAS associations directly from the EMBL-EBI GWAS Catalog REST API v2.
- **Removed old retrieval dependency:** The former third-party GWAS retrieval dependency was removed from `DESCRIPTION`, package code, tests, and documentation.
- **Updated `fetch_gwas()` output:** `fetch_gwas()` now returns a package-native list containing `associations`, `risk_alleles`, and an internal `cache` table.
- **Updated `run_gwas2crispr()` integration:** `run_gwas2crispr()` now consumes the new package-native `fetch_gwas()` result structure and no longer expects an S4 associations object.
- **Direct endpoint usage:** The package now uses GWAS Catalog REST API v2 endpoints for EFO trait resolution, association retrieval, and SNP metadata retrieval.
- **Optional FASTA export preserved:** CSV and BED outputs are still produced when `out_prefix` is supplied. FASTA output is generated only when `BSgenome.Hsapiens.UCSC.hg38` and `Biostrings` are installed.
- **No default file writing preserved:** `run_gwas2crispr()` still writes no files when `out_prefix = NULL`.
- **Improved written-path handling:** `run_gwas2crispr()` no longer returns `NA` entries inside `written` when FASTA output is unavailable.
- **Improved summary logic:** `SNPs_w_gene` is now calculated as the number of distinct SNPs with gene annotation.
- **Updated tests:** Tests were updated for the new direct REST API output structure, input validation, optional FASTA behaviour, and safe output writing to `tempdir()`.
- **Documentation updated:** README, package metadata, roxygen documentation, and tests were updated to match the direct GWAS Catalog REST API v2 workflow.

### Test environments

Local checks should be run before submission using:

- `devtools::document()`
- `devtools::test()`
- `rcmdcheck::rcmdcheck(args = "--as-cran", error_on = "never")`

### R CMD check results

To be updated after running local checks.

Expected target:

0 errors | 0 warnings | 

### Additional notes

- The package performs network requests to the EMBL-EBI GWAS Catalog REST API v2.
- Network-dependent examples are wrapped in `\donttest{}`.
- Network-dependent tests are skipped on CRAN using `skip_on_cran()` and guarded for offline/API failure conditions.
- FASTA export remains optional and depends on `BSgenome.Hsapiens.UCSC.hg38` and `Biostrings`, both listed in `Suggests`.
- CSV and BED outputs can be produced without optional genome packages.
- The CLI script depends on `optparse`, which is listed in `Suggests`.
- The package prepares computational outputs for downstream CRISPR guide-design workflows. It does not perform therapeutic interpretation, wet-lab validation, or biological efficacy testing.
