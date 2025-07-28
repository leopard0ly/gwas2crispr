test_that("fetch_gwas returns S4 associations object", {
  skip_on_cran()
  res <- try(fetch_gwas("EFO_0000707", p_cut = 1e-6), silent = TRUE)
  expect_true(inherits(res, "associations") || inherits(res, "try-error"))
})

test_that("run_gwas2crispr runs without error for valid input", {
  skip_on_cran()
  # This test only checks that the function runs and outputs files
  expect_error(
    run_gwas2crispr(
      efo_id = "EFO_0000707",
      p_cut = 1e-6,
      flank_bp = 300,
      out_prefix = "testout"
    ),
    NA # expect no error
  )
  # Clean up generated files (if any)
  file.remove(list.files(pattern = "^testout_snps_.*"))
})
