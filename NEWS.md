# gwas2crispr 0.1.2

- **CLI improvements:** The command-line script (`inst/scripts/gwas2crispr.R`) no longer forces file output.  The `--out`/`-o` option is now truly optional: omit it to run the pipeline in memory without writing any files.  A new `--verbose`/`-v` flag prints a concise summary when no outputs are written.  The double comma in the option definition has been fixed.
- **Console output:** All side‑effects (`print()`/`cat()`) have been replaced with `message()` behind a `verbose` flag.  By default, functions run quietly unless `verbose = TRUE` or the new CLI flag is used.
- **Testing:** Added `tests/testthat/test-run_gwas2crispr.R` to ensure that `run_gwas2crispr()` writes only to a provided `out_prefix` and returns objects without writing by default.
- **Documentation:** Updated DESCRIPTION with expanded acronyms and method references; ensured `Language: en-US` and valid Bioconductor `biocViews`.  Updated the README and vignette to reflect the new no-default-write behaviour and optional CLI output.  Added a simple `inst/CITATION` entry.  Bumped version to 0.1.2.

# gwas2crispr 0.1.1

- Address CRAN review:
  - Expanded acronyms in DESCRIPTION (GWAS, EFO, SNP, BED, FASTA, CRISPR).
  - Added method reference: Sudlow et al. (2015) <doi:10.1093/nar/gkv1256>.
  - Replaced \dontrun{} with \donttest{} where network I/O may occur; kept other examples short.
  - No default file writing: run_gwas2crispr(out_prefix = NULL) now returns objects only.
  - Examples/vignette/tests write only to tempdir() when needed.
  - Console output is suppressible via `verbose` + `message()`; removed `cat()`/`print()` side-effects.
- Qualified `utils::capture.output` to avoid NOTES.
- Vignette updated: object-only example by default; optional `tempdir()` example for file outputs.
- Fixed CLI script: corrected `--out` option to use `type = "character"`.
- Fixed CITATION: defined `meta` object before using `meta$Version`.
- Corrected `Language` field and expanded Description with acronym explanations and references.
  

# gwas2crispr 0.1.0

- Initial public release.
- Fetch GWAS associations and export CSV/BED/FASTA for GRCh38/hg38.
- CLI script at `inst/scripts/gwas2crispr.R`.
