#!/usr/bin/env Rscript

main <- function() {
  # Require optparse without attaching; fail with clear message
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("The 'optparse' package is required. Install it with install.packages('optparse').", call. = FALSE)
  }
  if (!requireNamespace("gwas2crispr", quietly = TRUE)) {
    stop("The 'gwas2crispr' package must be installed.", call. = FALSE)
  }

  option_list <- list(
    optparse::make_option(c("-e", "--efo"), type = "character",
                          help = "EFO trait ID (required), e.g. EFO_0001663"),
    optparse::make_option(c("-p", "--pthresh"), type = "double", default = 5e-8,
                          help = "P-value cut-off [default %default]"),
    optparse::make_option(c("-f", "--flank"), type = "integer", default = 200,
                          help = "bp flank for FASTA [default %default]"),
    # out_prefix is optional; if omitted, no files are written
    optparse::make_option(c("-o", "--out"), type = "character",
                          default = NULL,
                          help = "Prefix for output files (optional; omit to suppress file writing)"),
    # verbose flag to control message output
    optparse::make_option(c("-v", "--verbose"), action = "store_true",
                          default = FALSE,
                          help = "Print progress messages")
  )

  parser <- optparse::OptionParser(
    option_list = option_list,
    description = "gwas2crispr: GWAS-to-CRISPR file generator (hg38 only)"
  )
  opt <- optparse::parse_args(parser)

  if (is.null(opt$efo) || !grepl("^EFO_\\d+$", opt$efo)) {
    optparse::print_help(parser)
    stop("EFO trait ID is required and must look like 'EFO_0001663'.", call. = FALSE)
  }

  res <- gwas2crispr::run_gwas2crispr(
    efo_id     = opt$efo,
    p_cut      = opt$pthresh,
    flank_bp   = opt$flank,
    out_prefix = opt$out,
    verbose    = opt$verbose
  )

  # If no out_prefix was provided and verbose is TRUE, print a concise summary
  if (isTRUE(opt$verbose) && is.null(opt$out)) {
    sum_tbl <- res$summary
    message("Summary:\n", paste(utils::capture.output(print(sum_tbl)), collapse = "\n"))
  }
}

# Execute only when called as a script, not when sourced
if (identical(environment(), globalenv())) {
  tryCatch(main(), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(status = 1L)
  })
}
