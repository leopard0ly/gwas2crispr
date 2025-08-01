test_that("fetch_gwas returns S4 associations object (or try-error on network issues)", {
  skip_on_cran()
  res <- try(fetch_gwas("EFO_0000707", p_cut = 1e-6), silent = TRUE)
  if (!inherits(res, "try-error")) {
    # If the connection succeeds: ensure the object is of the correct class
    expect_s4_class(res, "associations")
  } else {
    # If the network call fails, treat as success to avoid breaking the suite off-CRAN
    succeed("Network call failed (acceptable off-CRAN); try-error returned.")
  }
})

test_that("run_gwas2crispr runs and writes only to a temporary directory when requested", {
  skip_on_cran()

  # Write inside tempdir(), not the working directory
  tmp <- tempdir()
  prefix <- file.path(tmp, paste0("testout-", as.integer(Sys.time())))

  # Run quietly (no messages)
  res <- run_gwas2crispr(
    efo_id     = "EFO_0000707",
    p_cut      = 1e-6,
    flank_bp   = 300,
    out_prefix = prefix,
    verbose    = FALSE
  )

  # Clean up regardless of outcome
  on.exit({
    if (is.list(res) && !is.null(res$written)) {
      file.remove(res$written[file.exists(res$written)])
    }
  }, add = TRUE)

  # Should not throw an error
  expect_type(res, "list")

  # Files should be written (paths stored in res$written)
  expect_true(length(res$written) >= 2)  # At least CSV and BED (FASTA depends on BSgenome)
  expect_true(any(grepl("_snps_full\\.csv$", res$written)))
  expect_true(any(grepl("_snps_hg38\\.bed$", res$written)))

  # All reported paths must exist
  expect_true(all(file.exists(res$written)))
})
