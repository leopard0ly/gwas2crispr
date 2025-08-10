test_that("run_gwas2crispr writes only to tempdir() when requested", {
  skip_on_cran()
  # Skip this test entirely when offline.  If the GWAS REST API is
  # unreachable, run_gwas2crispr() will emit warnings about failed
  # requests that should not cause the test to fail.  testthat provides
  # skip_if_offline() to detect this situation.
  skip_if_offline()
  tmp <- tempdir()
  prefix <- file.path(tmp, paste0("testout-", as.integer(Sys.time())))
  # Wrap the call in suppressWarnings() to avoid printing HTTP 500
  # warnings when the API occasionally returns errors.  The test only
  # checks that the expected files are created when an output prefix is
  # provided.
  res <- suppressWarnings(
    run_gwas2crispr("EFO_0001663", 5e-8, 50, out_prefix = prefix, verbose = FALSE)
  )
  expect_true(all(file.exists(res$written)))
})


 test_that("run_gwas2crispr does not write by default and returns objects", {
  skip_on_cran()
  skip_if_offline()
  # When no output prefix is supplied, run_gwas2crispr() should not write
  # any files.  Wrap the call in suppressWarnings() to silence transient
  # network errors from the underlying GWAS API.
  res <- suppressWarnings(
    run_gwas2crispr("EFO_0001663", 5e-8, 50, out_prefix = NULL, verbose = FALSE)
  )
  expect_type(res, "list")
  expect_equal(length(res$written), 0L)
})
