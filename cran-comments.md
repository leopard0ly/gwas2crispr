## Submission (0.1.2)

This submission builds upon version 0.1.1 to address remaining feedback from CRAN and Bioconductor reviewers.  It further hardens the package for release and ensures that running the pipeline never writes files unless explicitly requested.

### Key changes since 0.1.1

- **Optional CLI output:** The command‑line script (`inst/scripts/gwas2crispr.R`) now treats the `--out`/`-o` option as truly optional.  When omitted, the script runs the pipeline entirely in memory and writes no files.  A new `--verbose`/`-v` flag prints a concise summary when no outputs are written.  A spurious double comma in the `--out` option definition has been removed.
- **Suppressible console output:** All calls to `print()` or `cat()` have been replaced with `message()` behind a `verbose` flag.  Functions run quietly by default; users can enable progress output via `verbose = TRUE` or the CLI flag.
- **Expanded tests:** Added unit tests to verify that `run_gwas2crispr()` writes only when `out_prefix` is supplied and returns objects without writing by default.  Tests avoid network calls on CRAN via `skip_on_cran()`.
- **Metadata and documentation:**  Expanded acronyms (GWAS, EFO, SNP, BED, FASTA, CRISPR) and added method references in CRAN auto‑link format (`<doi:10.1093/nar/gky1120>`, `<doi:10.1093/nar/gkac1010>`, `<doi:10.1126/science.1225829>`, `<https://www.ebi.ac.uk/efo>`).  Corrected the `Language` field to `en-US` and added valid Bioconductor `biocViews` (Software, Genetics, VariantAnnotation, SNP, DataImport).  Added a simple `inst/CITATION` entry.  Bumped version to 0.1.2.  Updated the README and vignette to reflect the new no‑default‑write behaviour.

### Test environments

- Local: Ubuntu 22.04 (Benghazi), R 4.4.1 — 0 errors | 0 warnings | 0 notes.
- win-builder: R-devel — OK.
- R-hub: linux‑x86_64, macOS‑arm64, windows‑x86_64 — OK.

### R CMD check results

0 errors | 0 warnings | 0 notes

### Additional notes

- FASTA export remains optional and depends on `BSgenome.Hsapiens.UCSC.hg38` (in Suggests).  The function checks `requireNamespace()` and skips the FASTA step gracefully if the package is unavailable.
- Examples that perform network I/O are wrapped in `\donttest{}`.  The vignette defaults to non‑evaluated code; an optional example writes to a temporary directory (`tempdir()`).
- The CLI script depends on `optparse`, which is listed in Suggests.  If missing, the script fails with an informative message rather than triggering an error during `R CMD check`.