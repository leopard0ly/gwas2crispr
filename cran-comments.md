## Resubmission (0.1.1)

This resubmission addresses the reviewer’s comments:
- Expanded acronyms in DESCRIPTION (GWAS, EFO, SNP, BED, FASTA, CRISPR).
- Added an inline reference: Sudlow et al. (2015) <doi:10.1093/nar/gkv1256>.
- Replaced \dontrun{} with \donttest{} where network I/O may occur; kept other examples short and fast.
- No default file writing: `run_gwas2crispr(out_prefix = NULL)` returns objects only.
- Examples, vignette, and tests write only to `tempdir()` when demonstrating file output.
- Replaced unsuppressable console output with `message()` behind a `verbose` flag.
- Qualified `utils::capture.output` to remove a local NOTE.

## Test environments
- Local: Windows 10 x64, R 4.4.1 — 0 errors | 0 warnings | 0 notes.
- win-builder: R-devel — OK.
- R-hub: linux, macOS-arm64, windows — OK.

## R CMD check results
0 errors | 0 warnings | 0 notes

## Additional notes
- **Bioconductor dependency:** FASTA export is optional and uses `BSgenome.Hsapiens.UCSC.hg38` (listed in Suggests). 
  At runtime we check `requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)`; if unavailable, the FASTA step is **skipped gracefully** and the function still succeeds (CSV/BED only). No automatic installation or downloads are performed.
- Examples that require network I/O are wrapped in \donttest{}.
- The vignette defaults to non-evaluated code; an optional example writes to `tempdir()` only.
- The CLI writes only when the user supplies `-o/--out` (no default path).
