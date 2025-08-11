## Disable Quarto calls during tests to avoid spurious errors when Quarto
## is unavailable in the testing environment. See zzz-quarto-guard.R for details.
Sys.setenv("_GWAS2CRISPR_DISABLE_QUARTO" = "1")

library(testthat)
library(gwas2crispr)
test_check("gwas2crispr")
