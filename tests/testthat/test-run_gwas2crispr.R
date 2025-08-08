test_that("run_gwas2crispr writes only to tempdir() when requested", {
  skip_on_cran()
  tmp <- tempdir()
  prefix <- file.path(tmp, paste0("testout-", as.integer(Sys.time())))
  res <- run_gwas2crispr("EFO_0001663", 5e-8, 50, out_prefix = prefix, verbose = FALSE)
  # Expect that returned paths exist on disk
  expect_true(all(file.exists(res$written)))
})

test_that("run_gwas2crispr does not write by default and returns objects", {
  skip_on_cran()
  res <- run_gwas2crispr("EFO_0001663", 5e-8, 50, out_prefix = NULL, verbose = FALSE)
  expect_type(res, "list")
  expect_equal(length(res$written), 0L)
})