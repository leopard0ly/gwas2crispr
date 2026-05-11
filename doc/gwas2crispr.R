## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE,
  message = FALSE,
  warning = FALSE
)

## -----------------------------------------------------------------------------
# install.packages("gwas2crispr")

## -----------------------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c(
#   "Biostrings",
#   "GenomeInfoDb",
#   "BSgenome.Hsapiens.UCSC.hg38"
# ))

## -----------------------------------------------------------------------------
# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# 
# devtools::install_github("leopard0ly/gwas2crispr")

## -----------------------------------------------------------------------------
# library(gwas2crispr)
# 
# gwas_data <- fetch_gwas(
#   efo_id  = "EFO_0000707",
#   p_cut   = 1e-6,
#   verbose = FALSE
# )
# 
# names(gwas_data)
# head(gwas_data$associations)

## -----------------------------------------------------------------------------
# res <- run_gwas2crispr(
#   efo_id     = "EFO_0000707",
#   p_cut      = 1e-6,
#   flank_bp   = 300,
#   out_prefix = NULL,
#   verbose    = FALSE
# )
# 
# res$summary
# head(res$snps_full)
# head(res$bed)

## -----------------------------------------------------------------------------
# out_prefix <- file.path(tempdir(), "lung")
# 
# res <- run_gwas2crispr(
#   efo_id     = "EFO_0000707",
#   p_cut      = 1e-6,
#   flank_bp   = 300,
#   out_prefix = out_prefix,
#   verbose    = FALSE
# )
# 
# res$written

## -----------------------------------------------------------------------------
# paste0(out_prefix, "_snps_full.csv")
# paste0(out_prefix, "_snps_hg38.bed")
# paste0(out_prefix, "_snps_flank300.fa")

## -----------------------------------------------------------------------------
# names(res)

## -----------------------------------------------------------------------------
# res$summary
# res$snps_full
# res$bed
# res$fasta
# res$written

## ----eval=TRUE----------------------------------------------------------------
sessionInfo()

