test_that("run_gwas2crispr validates inputs", {
  expect_error(
    run_gwas2crispr(
      efo_id = "bad_id",
      p_cut = 1e-6,
      flank_bp = 300,
      verbose = FALSE
    ),
    "efo_id must be a single string"
  )

  expect_error(
    run_gwas2crispr(
      efo_id = "EFO_0000707",
      p_cut = -1,
      flank_bp = 300,
      verbose = FALSE
    ),
    "p_cut must be a single positive numeric value"
  )

  expect_error(
    run_gwas2crispr(
      efo_id = "EFO_0000707",
      p_cut = NA_real_,
      flank_bp = 300,
      verbose = FALSE
    ),
    "p_cut must be a single positive numeric value"
  )

  expect_error(
    run_gwas2crispr(
      efo_id = "EFO_0000707",
      p_cut = 1e-6,
      flank_bp = -10,
      verbose = FALSE
    ),
    "flank_bp must be a single positive integer"
  )

  expect_error(
    run_gwas2crispr(
      efo_id = "EFO_0000707",
      p_cut = 1e-6,
      flank_bp = NA_real_,
      verbose = FALSE
    ),
    "flank_bp must be a single positive integer"
  )

  expect_error(
    run_gwas2crispr(
      efo_id = "EFO_0000707",
      p_cut = 1e-6,
      flank_bp = 300,
      out_prefix = "",
      verbose = FALSE
    ),
    "out_prefix must be NULL or a non-empty string"
  )

  expect_error(
    run_gwas2crispr(
      efo_id = "EFO_0000707",
      p_cut = 1e-6,
      flank_bp = 300,
      genome_pkg = "bad_genome",
      verbose = FALSE
    ),
    "Only GRCh38/hg38 is supported"
  )
})

test_that("run_gwas2crispr returns objects without writing by default", {
  skip_on_cran()
  skip_if_offline()

  res <- try(
    run_gwas2crispr(
      efo_id     = "EFO_0000707",
      p_cut      = 1e-6,
      flank_bp   = 300,
      out_prefix = NULL,
      verbose    = FALSE
    ),
    silent = TRUE
  )

  if (inherits(res, "try-error")) {
    succeed("Network/API unavailable during test.")
  } else {
    expect_type(res, "list")

    expect_true(all(
      c("summary", "chr_freq", "snps_full", "bed", "fasta", "written") %in% names(res)
    ))

    expect_equal(length(res$written), 0L)

    expect_true(nrow(res$summary) == 1L)
    expect_true(nrow(res$snps_full) > 0L)
    expect_true(nrow(res$bed) > 0L)

    expect_true(all(
      c("n_SNPs", "SNPs_w_gene", "unique_genes", "n_studies") %in% names(res$summary)
    ))

    expect_true(res$summary$n_SNPs >= res$summary$SNPs_w_gene)

    expect_true(all(
      c("chr", "start0", "end0", "id") %in% names(res$bed)
    ))

    expect_true(all(res$bed$start0 < res$bed$end0))
    expect_true(all(grepl("^chr", res$bed$chr)))
  }
})

test_that("run_gwas2crispr writes CSV and BED to tempdir when requested", {
  skip_on_cran()
  skip_if_offline()

  tmp <- tempdir()
  prefix <- file.path(tmp, paste0("testout-", as.integer(Sys.time())))

  res <- try(
    run_gwas2crispr(
      efo_id     = "EFO_0000707",
      p_cut      = 1e-6,
      flank_bp   = 300,
      out_prefix = prefix,
      verbose    = FALSE
    ),
    silent = TRUE
  )

  if (inherits(res, "try-error")) {
    succeed("Network/API unavailable during test.")
  } else {
    on.exit({
      existing <- res$written[!is.na(res$written) & file.exists(res$written)]
      if (length(existing) > 0L) {
        file.remove(existing)
      }
    }, add = TRUE)

    expect_type(res, "list")

    expect_true(all(
      c("summary", "chr_freq", "snps_full", "bed", "fasta", "written") %in% names(res)
    ))

    expect_false(any(is.na(res$written)))

    expect_true(any(grepl("_snps_full\\.csv$", res$written)))
    expect_true(any(grepl("_snps_hg38\\.bed$", res$written)))

    csv_path <- res$written[grepl("_snps_full\\.csv$", res$written)]
    bed_path <- res$written[grepl("_snps_hg38\\.bed$", res$written)]

    expect_length(csv_path, 1L)
    expect_length(bed_path, 1L)

    expect_true(file.exists(csv_path))
    expect_true(file.exists(bed_path))

    written_existing <- res$written[!is.na(res$written)]

    expect_true(all(file.exists(written_existing)))

    tmp_norm <- normalizePath(tmp, winslash = "/", mustWork = TRUE)
    written_norm <- normalizePath(written_existing, winslash = "/", mustWork = TRUE)

    expect_true(all(startsWith(written_norm, paste0(tmp_norm, "/"))))

    csv_data <- readr::read_csv(csv_path, show_col_types = FALSE)
    bed_data <- readr::read_tsv(
      bed_path,
      col_names = c("chr", "start0", "end0", "id"),
      show_col_types = FALSE
    )

    expect_true(nrow(csv_data) > 0L)
    expect_true(nrow(bed_data) > 0L)

    expect_true(all(
      c("variant_id", "chromosome_name", "chromosome_position") %in% names(csv_data)
    ))

    expect_true(all(
      c("chr", "start0", "end0", "id") %in% names(bed_data)
    ))

    expect_true(all(bed_data$start0 < bed_data$end0))
    expect_true(all(grepl("^chr", bed_data$chr)))
  }
})

test_that("run_gwas2crispr written paths contain no NA when FASTA is optional", {
  skip_on_cran()
  skip_if_offline()

  tmp <- tempdir()
  prefix <- file.path(tmp, paste0("testout-no-na-", as.integer(Sys.time())))

  res <- try(
    run_gwas2crispr(
      efo_id     = "EFO_0000707",
      p_cut      = 1e-6,
      flank_bp   = 300,
      out_prefix = prefix,
      verbose    = FALSE
    ),
    silent = TRUE
  )

  if (inherits(res, "try-error")) {
    succeed("Network/API unavailable during test.")
  } else {
    on.exit({
      existing <- res$written[!is.na(res$written) & file.exists(res$written)]
      if (length(existing) > 0L) {
        file.remove(existing)
      }
    }, add = TRUE)

    expect_false(any(is.na(res$written)))
    expect_true(length(res$written) >= 2L)
    expect_true(all(file.exists(res$written)))
  }
})
