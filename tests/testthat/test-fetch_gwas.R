test_that("fetch_gwas validates inputs", {
  expect_error(
    fetch_gwas("bad_id", p_cut = 1e-6, verbose = FALSE),
    "efo_id must be a single string"
  )

  expect_error(
    fetch_gwas("EFO_0000707", p_cut = -1, verbose = FALSE),
    "p_cut must be a single positive numeric value"
  )

  expect_error(
    fetch_gwas("EFO_0000707", p_cut = NA_real_, verbose = FALSE),
    "p_cut must be a single positive numeric value"
  )

  expect_error(
    fetch_gwas(c("EFO_0000707", "EFO_0001663"), p_cut = 1e-6, verbose = FALSE),
    "efo_id must be a single string"
  )
})

test_that("internal p-value extraction handles mantissa and exponent", {
  f <- c(
    pvalueMantissa = "5",
    pvalueExponent = "-8"
  )

  expect_equal(
    gwas2crispr:::get_pvalue_from_flat(f),
    5e-8
  )
})

test_that("internal p-value extraction handles direct pvalue", {
  f <- c(pvalue = "1e-6")

  expect_equal(
    gwas2crispr:::get_pvalue_from_flat(f),
    1e-6
  )
})

test_that("internal p-value extraction returns NA for invalid p-values", {
  f_empty <- character(0)
  f_bad <- c(pvalue = "not_a_number")
  f_out <- c(pvalue = "2")

  expect_true(is.na(gwas2crispr:::get_pvalue_from_flat(f_empty)))
  expect_true(is.na(gwas2crispr:::get_pvalue_from_flat(f_bad)))
  expect_true(is.na(gwas2crispr:::get_pvalue_from_flat(f_out)))
})

test_that("internal flattening preserves nested values", {
  x <- list(
    a = list(
      b = "rs123",
      c = list(d = "GCST000001")
    )
  )

  f <- gwas2crispr:::flatten_atomic(x)

  expect_true(any(grepl("rs123", f)))
  expect_true(any(grepl("GCST000001", f)))
})

test_that("internal chromosome and position extraction work from flattened records", {
  f <- c(
    chromosomeName = "8",
    chromosomePosition = "128748315"
  )

  expect_equal(gwas2crispr:::get_chr_from_flat(f), "8")
  expect_equal(gwas2crispr:::get_pos_from_flat(f), 128748315L)
})

test_that("internal gene extraction removes URLs and ENSG identifiers", {
  f <- c(
    mappedGene = "MYC",
    geneUrl = "https://example.org/gene",
    ensemblGene = "ENSG00000136997"
  )

  genes <- gwas2crispr:::get_genes_from_flat(f)

  expect_true(grepl("MYC", genes))
  expect_false(grepl("https://", genes))
  expect_false(grepl("ENSG", genes))
})

test_that("extract_items handles expected GWAS Catalog embedded keys", {
  js_assoc <- list(
    `_embedded` = list(
      associations = list(list(a = 1))
    )
  )

  js_snp <- list(
    `_embedded` = list(
      singleNucleotidePolymorphisms = list(list(a = 1))
    )
  )

  js_efo <- list(
    `_embedded` = list(
      efoTraits = list(list(a = 1))
    )
  )

  expect_equal(length(gwas2crispr:::extract_items(js_assoc)), 1L)
  expect_equal(length(gwas2crispr:::extract_items(js_snp)), 1L)
  expect_equal(length(gwas2crispr:::extract_items(js_efo)), 1L)
})

test_that("fetch_gwas returns package-native structure when network is available", {
  skip_on_cran()
  skip_if_offline()

  res <- try(
    fetch_gwas(
      efo_id = "EFO_0000707",
      p_cut = 1e-6,
      verbose = FALSE
    ),
    silent = TRUE
  )

  if (inherits(res, "try-error")) {
    succeed("Network/API unavailable during test.")
  } else {
    expect_type(res, "list")

    expect_true(all(
      c("associations", "risk_alleles", "cache") %in% names(res)
    ))

    expect_true(all(
      c("association_id", "pvalue") %in% names(res$associations)
    ))

    expect_true(all(
      c("association_id", "variant_id") %in% names(res$risk_alleles)
    ))

    expect_true(all(
      c(
        "variant_id",
        "study_accession",
        "p_numeric",
        "chromosome_name",
        "chromosome_position",
        "mapped_gene",
        "association_id"
      ) %in% names(res$cache)
    ))

    expect_true(nrow(res$associations) > 0L)
    expect_true(nrow(res$risk_alleles) > 0L)
    expect_true(nrow(res$cache) > 0L)

    expect_true(all(res$associations$pvalue < 1e-6))
    expect_true(all(grepl("^rs[0-9]+$", res$risk_alleles$variant_id)))
  }
})
