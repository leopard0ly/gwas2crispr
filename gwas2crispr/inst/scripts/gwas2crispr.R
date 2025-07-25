#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(gwas2crispr))

opt <- OptionParser() |>
  add_option(c("-e", "--efo"),    type = "character",
             help = "EFO trait ID (مطلوب)") |>
  add_option(c("-p", "--pthresh"), type = "double",  default = 5e-8,
             help = "P‑value cut‑off [default %default]") |>
  add_option(c("-f", "--flank"),   type = "integer", default = 200,
             help = "bp flank for FASTA [default %default]") |>
  add_option(c("-o", "--out"),     type = "character", default = "output",
             help = "Prefix for output files [default %default]") |>
  parse_args()

if (is.null(opt$efo))
  stop("EFO trait ID is required: -e EFO_0001663", call. = FALSE)

gwas2crispr::run_gwas2crispr(
  efo_id     = opt$efo,
  p_cut      = opt$pthresh,
  flank_bp   = opt$flank,
  out_prefix = opt$out
)
