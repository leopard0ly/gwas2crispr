# gwas2crispr 0.1.1

- Address CRAN review:
  - Expanded acronyms in DESCRIPTION (GWAS, EFO, SNP, BED, FASTA, CRISPR).
  - Added method reference: Sudlow et al. (2015) <doi:10.1093/nar/gkv1256>.
  - Replaced \dontrun{} with \donttest{} where network I/O may occur; kept other examples short.
  - No default file writing: run_gwas2crispr(out_prefix = NULL) now returns objects only.
  - Examples/vignette/tests write only to tempdir() when needed.
  - Console output is suppressible via verbose + message(); removed cat()/print() side-effects.
  - Qualified utils::capture.output to avoid NOTES.
- Vignette updated: object-only example by default; optional tempdir() example for file outputs.

# gwas2crispr 0.1.0

- Initial public release.
- Fetch GWAS associations and export CSV/BED/FASTA for GRCh38/hg38.
- CLI script at `inst/scripts/gwas2crispr.R`.
