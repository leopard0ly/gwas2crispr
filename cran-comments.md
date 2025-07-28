## Test environments
- Local: R 4.x on Windows/Linux (devtools::check)
- CI: (to be added) macOS-latest, windows-latest, ubuntu-latest via r-lib/actions

## R CMD check results
0 errors | 0 warnings | 0 notes

## Notes for CRAN
- Examples are wrapped in \dontrun{} due to network I/O and file writing.
- Package targets GRCh38/hg38 (BSgenome.Hsapiens.UCSC.hg38).
