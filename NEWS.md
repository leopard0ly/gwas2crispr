# gwas2crispr 0.1.4

- **GWAS Catalog REST API v2 backend:** Replaced the former third-party GWAS retrieval workflow with direct EMBL-EBI GWAS Catalog REST API v2 retrieval.
- **Removed old GWAS retrieval dependency:** Removed `gwasrapidd` from package dependency metadata, code, tests, and user-facing documentation.
- **Updated `fetch_gwas()` output:** `fetch_gwas()` now returns a package-native list with:
  - `associations`
  - `risk_alleles`
  - `cache`
- **Updated `run_gwas2crispr()` integration:** `run_gwas2crispr()` now consumes the new package-native `fetch_gwas()` result and no longer expects an S4 associations object.
- **Direct endpoint workflow:** Added direct handling for GWAS Catalog REST API v2 trait resolution, association retrieval, and SNP metadata retrieval.
- **Optional FASTA behavior improved:** CSV and BED outputs remain available even when FASTA dependencies are unavailable. FASTA is generated only when `BSgenome.Hsapiens.UCSC.hg38` and `Biostrings` are installed.
- **Output path handling improved:** `run_gwas2crispr()` no longer returns `NA` entries inside the `written` vector when FASTA output is skipped.
- **Summary logic improved:** `SNPs_w_gene` is now calculated as the number of distinct SNPs with gene annotation, avoiding inflated counts when variants appear in multiple association rows.
- **Tests updated:** Updated tests for input validation, direct REST API output structure, optional FASTA behavior, safe writing to `tempdir()`, and no-default-write behavior.
- **Documentation updated:** Updated README, `cran-comments.md`, roxygen documentation, and package metadata for the direct GWAS Catalog REST API v2 workflow.
- **Version bump:** Bumped package version to `0.1.4`.

# gwas2crispr 0.1.2

- **CLI improvements:** The command-line script (`inst/scripts/gwas2crispr.R`) no longer forces file output. The `--out`/`-o` option is now truly optional: omit it to run the pipeline in memory without writing any files. A new `--verbose`/`-v` flag prints a concise summary when no outputs are written. The double comma in the option definition has been fixed.
- **Console output:** All side-effects (`print()`/`cat()`) have been replaced with `message()` behind a `verbose` flag. By default, functions run quietly unless `verbose = TRUE` or the new CLI flag is used.
- **Testing:** Added `tests/testthat/test-run_gwas2crispr.R` to ensure that `run_gwas2crispr()` writes only to a provided `out_prefix` and returns objects without writing by default.
- **Documentation:** Updated DESCRIPTION with expanded acronyms and method references; ensured `Language: en-US` and valid Bioconductor `biocViews`. Updated the README and vignette to reflect the new no-default-write behaviour and optional CLI output. Added a simple `inst/CITATION` entry. Bumped version to `0.1.2`.

# gwas2crispr 0.1.1

- Addressed CRAN review comments:
  - Expanded acronyms in DESCRIPTION: GWAS, EFO, SNP, BED, FASTA, and CRISPR.
  - Added method reference: Sudlow et al. (2015) <doi:10.1093/nar/gkv1256>.
  - Replaced `\dontrun{}` with `\donttest{}` where network I/O may occur.
  - Kept examples short.
  - Ensured `run_gwas2crispr(out_prefix = NULL)` returns objects only and does not write files by default.
  - Ensured examples, vignette, and tests write only to `tempdir()` when file output is needed.
  - Made console output suppressible via `verbose` and `message()`.
- Qualified `utils::capture.output` to avoid NOTES.
- Updated vignette with object-only example by default and optional `tempdir()` example for file outputs.
- Fixed CLI script by correcting the `--out` option to use `type = "character"`.
- Fixed CITATION by defining the `meta` object before using `meta$Version`.
- Corrected the `Language` field.
- Expanded DESCRIPTION with acronym explanations and references.

# gwas2crispr 0.1.0

- Initial public release.
- Fetch GWAS associations and export CSV, BED, and FASTA outputs for GRCh38/hg38.
- Added CLI script at `inst/scripts/gwas2crispr.R`.
