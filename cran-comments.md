## Submission (0.1.5)

This submission is a backward-compatible update to trait identifier validation and retrieval robustness.

### Key changes since 0.1.4

- **Selected GWAS Catalog trait identifiers:** `fetch_gwas()` and `run_gwas2crispr()` now accept selected GWAS Catalog trait identifier formats beyond the original EFO-only validation.
- **Backward-compatible interface preserved:** Exported function names are unchanged, the public `efo_id` argument is retained, and the CLI `--efo` option is retained.
- **Identifier normalization:** Underscore and colon identifier formats are normalized internally.
- **Retrieval robustness improved:** Direct identifier-based association retrieval is attempted before label-based retrieval.
- **Coordinate fallback improved:** Missing rsID coordinates can be recovered through GWAS Catalog SNP metadata and an optional non-fatal Ensembl REST fallback.
- **Output contract unchanged:** CSV, BED, and optional FASTA filename patterns are unchanged.
- **Genome policy unchanged:** GRCh38/hg38 remains the only supported genome build.
- **No new required dependencies:** The update uses existing required packages; optional FASTA packages remain in `Suggests`.
- **Tests updated:** Validation, retrieval cascade, response parsing, coordinate fallback, output filename, and returned-object tests were updated.

### Test environments

Local checks should be run before submission using:

- `devtools::document()`
- `devtools::test()`
- `rcmdcheck::rcmdcheck(args = "--as-cran", error_on = "never")`

### R CMD check results

Local `rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"), error_on = "never")` result:

0 errors | 0 warnings | 1 note

The note was:

- unable to verify current time

A full `--as-cran` check was also attempted. Vignettes built successfully after pointing `RSTUDIO_PANDOC` to the bundled RStudio Pandoc, but the indexed PDF manual step could not complete because the local TinyTeX installation does not include `makeindex`. A direct no-index manual PDF build completed successfully.

### Additional notes

- The package performs network requests to the EMBL-EBI GWAS Catalog REST API v2 and, when needed for coordinate fallback, Ensembl REST.
- Network-dependent examples are wrapped in `\donttest{}`.
- Network-dependent tests are skipped on CRAN using `skip_on_cran()` and guarded for offline/API failure conditions.
- FASTA export remains optional and depends on `BSgenome.Hsapiens.UCSC.hg38` and `Biostrings`, both listed in `Suggests`.
- CSV and BED outputs can be produced without optional genome packages.
- The CLI script depends on `optparse`, which is listed in `Suggests`.
- The package prepares computational outputs for downstream CRISPR guide-design workflows. It does not perform biological causality testing, clinical interpretation, therapeutic design, or wet-lab validation.
