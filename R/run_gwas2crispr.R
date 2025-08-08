#' Run the GWAS to CRISPR export pipeline (hg38)
#'
#' @description End-to-end pipeline: fetch significant associations, annotate,
#' and optionally write CSV/BED/FASTA outputs. By default no files are written;
#' set \code{out_prefix} to write results.
#'
#' @param efo_id   character. Experimental Factor Ontology (EFO) identifier,
#'   e.g., "EFO_0001663".
#' @param p_cut    numeric. P-value threshold for significance (default 5e-8).
#' @param flank_bp integer. Flanking bases for FASTA sequences (default 200).
#' @param out_prefix character or NULL. File prefix (including path) for outputs.
#'   If \code{NULL} (default), nothing is written to disk and a result object is returned.
#'   To write files safely in examples/tests, use \code{file.path(tempdir(), "prefix")}.
#' @param genome_pkg character. BSgenome package to use for FASTA
#'   (default "BSgenome.Hsapiens.UCSC.hg38"); FASTA step is skipped if not installed.
#' @param verbose  logical. If TRUE, emit progress via \code{message()}.
#'
#' @return (Invisibly) a list with elements:
#'   \itemize{
#'     \item \code{summary}: tibble with basic counts.
#'     \item \code{snps_full}: tibble of SNP metadata.
#'     \item \code{bed}: tibble of BED intervals (if computed).
#'     \item \code{fasta}: \code{Biostrings::DNAStringSet} (if computed).
#'     \item \code{written}: character vector of file paths written (possibly empty).
#'   }
#'
#' @details
#' Network I/O may occur when fetching data. Only GRCh38/hg38 is supported.
#'
#' @seealso \code{\link{fetch_gwas}}
#'
#' @examples
#' \donttest{
#'   # Write into a temporary directory so we don't touch the user's filespace:
#'   tmp <- tempdir()
#'   res <- run_gwas2crispr(
#'     efo_id     = "EFO_0001663",
#'     p_cut      = 5e-8,
#'     flank_bp   = 200,
#'     out_prefix = file.path(tmp, "prostate"),
#'     verbose    = FALSE
#'   )
#'
#'   # If you omit 'out_prefix', nothing is written; an object is returned:
#'   res2 <- run_gwas2crispr(
#'     efo_id     = "EFO_0001663",
#'     p_cut      = 5e-8,
#'     flank_bp   = 200,
#'     out_prefix = NULL,
#'     verbose    = FALSE
#'   )
#' }
#'
#' @export
run_gwas2crispr <- function(efo_id,
                            p_cut = 5e-8,
                            flank_bp = 200,
                            out_prefix = NULL,
                            genome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
                            verbose = interactive()) {
  # --- validate inputs -------------------------------------------------------
  if (!is.character(efo_id) || length(efo_id) != 1L || !grepl("^EFO_\\d+$", efo_id)) {
    stop("efo_id must be a single string like 'EFO_0001663'.", call. = FALSE)
  }
  if (!is.numeric(p_cut) || length(p_cut) != 1L || !is.finite(p_cut) || p_cut <= 0) {
    stop("p_cut must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(flank_bp) || length(flank_bp) != 1L || !is.finite(flank_bp) || flank_bp <= 0) {
    stop("flank_bp must be a single positive integer.", call. = FALSE)
  }
  if (!is.null(out_prefix) && (!is.character(out_prefix) || length(out_prefix) != 1L || nchar(out_prefix) == 0L)) {
    stop("out_prefix must be NULL or a non-empty string.", call. = FALSE)
  }
  if (!identical(genome_pkg, "BSgenome.Hsapiens.UCSC.hg38")) {
    stop("Only GRCh38/hg38 is supported: set genome_pkg = 'BSgenome.Hsapiens.UCSC.hg38'.",
         call. = FALSE)
  }

  vmsg <- function(...) if (isTRUE(verbose)) message(...)

  # gwasrapidd is required for variant annotation; fail early if absent
  if (!requireNamespace("gwasrapidd", quietly = TRUE)) {
    stop(
      "The 'gwasrapidd' package is required for variant annotation. Please install it (e.g., via Bioconductor).",
      call. = FALSE
    )
  }

  # --- 1) Fetch associations -------------------------------------------------
  vmsg("Fetching GWAS associations for: ", efo_id)
  cat_assocs <- fetch_gwas(efo_id = efo_id, p_cut = p_cut)

  lead_assoc <- cat_assocs@associations |>
    dplyr::filter(pvalue < p_cut)

  lead_variants <- cat_assocs@risk_alleles |>
    dplyr::filter(association_id %in% lead_assoc$association_id) |>
    dplyr::pull(variant_id) |>
    unique()

  # --- 2) Variant & gene metadata -------------------------------------------
  get_variant_chunks <- function(v, size = 100) {
    split(v, ceiling(seq_along(v) / size)) |>
      purrr::map(~ gwasrapidd::get_variants(variant_id = .x, verbose = FALSE))
  }

  v_chunks   <- get_variant_chunks(lead_variants)
  variant_df <- dplyr::bind_rows(purrr::map(v_chunks, ~ .x@variants))
  context_df <- dplyr::bind_rows(purrr::map(v_chunks, ~ .x@genomic_contexts))

  if (nrow(context_df) == 0L || !all(c("variant_id", "gene_name", "is_mapped_gene") %in% names(context_df))) {
    gene_nest <- tibble::tibble(variant_id = character(), genes = list())
  } else {
    gene_nest <- context_df |>
      dplyr::filter(is_mapped_gene) |>
      dplyr::select(variant_id, gene_name) |>
      dplyr::distinct() |>
      tidyr::nest(genes = c(gene_name))
  }

  variant_tbl <- variant_df |>
    dplyr::left_join(gene_nest, by = "variant_id") |>
    dplyr::mutate(
      genes = purrr::map(genes, ~ if (is.null(.x) || nrow(.x) == 0)
        tibble::tibble(gene_name = character()) else .x)
    ) |>
    dplyr::arrange(chromosome_name, chromosome_position)

  # --- 3) Study & effect metrics --------------------------------------------
  assoc2study <- cat_assocs@risk_alleles |>
    dplyr::transmute(
      association_id,
      study_accession = sub(".*:", "", association_id)
    ) |>
    dplyr::distinct()

  effect_cols <- c("risk_allele", "other_allele", "beta", "odds_ratio",
                   "effect_allele_frequency", "p_value", "pvalue")
  present_cols <- intersect(names(cat_assocs@associations), effect_cols)

  assoc_extra <- cat_assocs@associations |>
    dplyr::select(association_id, dplyr::all_of(present_cols)) |>
    dplyr::distinct()

  risk_map <- cat_assocs@risk_alleles |>
    dplyr::select(variant_id, association_id) |>
    dplyr::distinct()

  study_raw <- tryCatch(cat_assocs@studies, error = function(e) NULL)
  study_df <- if (is.null(study_raw) || nrow(study_raw) == 0) {
    tibble::tibble(study_accession = character(),
                   pubmed_id = character(),
                   publication_date = as.Date(character()))
  } else {
    tibble::as_tibble(study_raw) |>
      dplyr::select(dplyr::any_of(c("study_accession", "pubmed_id", "publication_date", "sample_size")))
  }

  variant_full <- variant_tbl |>
    dplyr::left_join(risk_map,    by = "variant_id") |>
    dplyr::left_join(assoc2study, by = "association_id") |>
    dplyr::left_join(assoc_extra, by = "association_id") |>
    dplyr::left_join(study_df,    by = "study_accession") |>
    dplyr::distinct()

  # --- 4) Summary ------------------------------------------------------------
  summary_tbl <- variant_full |>
    dplyr::mutate(
      has_gene = purrr::map_lgl(genes, ~ nrow(.x) > 0),
      gene_vec = purrr::map(genes, ~ .x$gene_name)
    ) |>
    dplyr::summarise(
      n_SNPs       = dplyr::n_distinct(variant_id),
      SNPs_w_gene  = sum(has_gene),
      unique_genes = dplyr::n_distinct(unlist(gene_vec, use.names = FALSE)),
      n_studies    = dplyr::n_distinct(study_accession, na.rm = TRUE)
    )

  chr_freq <- variant_full |>
    dplyr::count(chromosome_name, name = "SNPs") |>
    dplyr::arrange(dplyr::desc(SNPs)) |>
    dplyr::slice_head(n = 10)

  vmsg("\n===== GWAS -> CRISPR Summary =====")
  vmsg(paste(utils::capture.output(print(summary_tbl)), collapse = "\n"))
  vmsg("\nTop-10 chromosomes by SNP count:")
  vmsg(paste(utils::capture.output(print(chr_freq)), collapse = "\n"))

  # --- 5) Optional output files ---------------------------------------------
  written <- character(0L)
  bed_df <- variant_full |>
    dplyr::filter(!is.na(chromosome_name)) |>
    dplyr::transmute(
      chr    = paste0("chr", chromosome_name),
      start0 = chromosome_position - 1L,
      end0   = chromosome_position,
      id     = variant_id
    )

  fasta_set <- NULL

  if (!is.null(out_prefix)) {
    csv_path <- paste0(out_prefix, "_snps_full.csv")
    bed_path <- paste0(out_prefix, "_snps_hg38.bed")
    fa_path  <- paste0(out_prefix, "_snps_flank", flank_bp, ".fa")

    readr::write_csv(variant_full, csv_path)
    readr::write_tsv(bed_df, bed_path, col_names = FALSE)

    if (requireNamespace(genome_pkg, quietly = TRUE)) {
      genome_obj <- getExportedValue(genome_pkg, "Hsapiens")
      start1 <- pmax(bed_df$start0 - flank_bp + 1L, 1L)
      end1   <- bed_df$end0 + flank_bp
      fasta_set <- Biostrings::getSeq(genome_obj,
                                      names = bed_df$chr,
                                      start = start1,
                                      end   = end1)
      names(fasta_set) <- bed_df$id
      Biostrings::writeXStringSet(fasta_set, fa_path)
    } else {
      vmsg(genome_pkg, " not installed: FASTA step skipped")
    }

    written <- c(
      paste0(out_prefix, "_snps_full.csv"),
      paste0(out_prefix, "_snps_hg38.bed"),
      paste0(out_prefix, "_snps_flank", flank_bp, ".fa")
    )
    vmsg("Done. Output files:\n ", paste(basename(written), collapse = " \n "))
  }

  invisible(list(
    summary   = summary_tbl,
    snps_full = variant_full,
    bed       = bed_df,
    fasta     = fasta_set,
    written   = written
  ))
}
