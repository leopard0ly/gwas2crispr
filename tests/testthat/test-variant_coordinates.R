test_that("fill_missing_coordinates merges fallback rows without overwriting existing coordinates", {
  variant_df <- tibble::tibble(
    variant_id = "rs1",
    chromosome_name = "1",
    chromosome_position = 100L
  )
  fallback_df <- tibble::tibble(
    variant_id = c("rs1", "rs2", "rs2"),
    chromosome_name = c("2", "3", "3"),
    chromosome_position = c(200L, 300L, 300L)
  )

  out <- gwas2crispr:::fill_missing_coordinates(variant_df, fallback_df)

  expect_equal(nrow(out), 2L)
  expect_equal(out$chromosome_name[out$variant_id == "rs1"], "1")
  expect_equal(out$chromosome_position[out$variant_id == "rs1"], 100L)
  expect_equal(out$chromosome_name[out$variant_id == "rs2"], "3")
  expect_equal(out$chromosome_position[out$variant_id == "rs2"], 300L)
})

test_that("fill_missing_coordinates handles empty fallback rows", {
  variant_df <- tibble::tibble(
    variant_id = "rs1",
    chromosome_name = "1",
    chromosome_position = 100L
  )

  out <- gwas2crispr:::fill_missing_coordinates(
    variant_df,
    tibble::tibble()
  )

  expect_equal(nrow(out), 1L)
  expect_equal(names(out), c("variant_id", "chromosome_name", "chromosome_position"))
})

test_that("pull_ensembl_variant_details returns empty coordinate table for empty input", {
  out <- gwas2crispr:::pull_ensembl_variant_details(character(), verbose = FALSE)

  expect_equal(nrow(out), 0L)
  expect_equal(names(out), c("variant_id", "chromosome_name", "chromosome_position"))
})
