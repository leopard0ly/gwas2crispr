#' Run the GWAS-to-CRISPR export pipeline using GRCh38/hg38
#'
#' @description
#' Runs the complete computational preparation workflow: retrieves GWAS Catalog
#' associations through \code{\link{fetch_gwas}}, prepares SNP metadata, creates
#' BED intervals, and optionally writes CSV, BED, and FASTA files for downstream
#' CRISPR guide-design preparation.
#'
#' @param efo_id character. EFO trait identifier, such as EFO_0001663.
#' @param p_cut numeric. P-value threshold for significance.
#' @param flank_bp integer. Number of flanking bases for FASTA sequence extraction.
#' @param out_prefix character or \code{NULL}. Prefix for output files. If
#'   \code{NULL}, no files are written.
#' @param genome_pkg character. BSgenome package name used for hg38 FASTA extraction.
#' @param verbose logical. If \code{TRUE}, prints a compact progress line.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item \code{summary}: one-row tibble with basic counts.
#'   \item \code{chr_freq}: chromosome frequency table.
#'   \item \code{snps_full}: harmonized SNP metadata.
#'   \item \code{bed}: BED-style interval table.
#'   \item \code{fasta}: DNAStringSet if FASTA was generated; otherwise \code{NULL}.
#'   \item \code{written}: character vector of written file paths.
#' }
#'
#' @details
#' Only GRCh38/hg38 is supported. CSV and BED outputs can be produced without
#' genome packages. FASTA output is generated only when
#' \pkg{BSgenome.Hsapiens.UCSC.hg38} and \pkg{Biostrings} are installed.
#' If FASTA dependencies are unavailable, the function still writes CSV and BED.
#'
#' @seealso \code{\link{fetch_gwas}}
#'
#' @examples
#' \donttest{
#'   res <- run_gwas2crispr(
#'     efo_id     = "EFO_0000707",
#'     p_cut      = 1e-6,
#'     flank_bp   = 300,
#'     out_prefix = file.path(tempdir(), "lung"),
#'     verbose    = FALSE
#'   )
#'   res$summary
#'   res$written
#' }
#'
#' @export
run_gwas2crispr <- function(efo_id,
                            p_cut = 5e-8,
                            flank_bp = 200,
                            out_prefix = NULL,
                            genome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
                            verbose = interactive()) {
  if (!is.character(efo_id) || length(efo_id) != 1L || !grepl("^EFO_\\d+$", efo_id)) {
    stop("efo_id must be a single string like 'EFO_0001663'.", call. = FALSE)
  }

  if (!is.numeric(p_cut) || length(p_cut) != 1L || !is.finite(p_cut) || p_cut <= 0) {
    stop("p_cut must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(flank_bp) || length(flank_bp) != 1L || !is.finite(flank_bp) || flank_bp <= 0) {
    stop("flank_bp must be a single positive integer.", call. = FALSE)
  }

  if (!is.null(out_prefix) &&
      (!is.character(out_prefix) || length(out_prefix) != 1L || nchar(out_prefix) == 0L)) {
    stop("out_prefix must be NULL or a non-empty string.", call. = FALSE)
  }

  if (!identical(genome_pkg, "BSgenome.Hsapiens.UCSC.hg38")) {
    stop(
      "Only GRCh38/hg38 is supported: set genome_pkg = 'BSgenome.Hsapiens.UCSC.hg38'.",
      call. = FALSE
    )
  }

  cat_assocs <- fetch_gwas(
    efo_id = efo_id,
    p_cut = p_cut,
    verbose = verbose
  )

  lead_assoc <- cat_assocs$associations |>
    dplyr::filter(.data$pvalue < p_cut)

  lead_variants <- cat_assocs$risk_alleles |>
    dplyr::filter(.data$association_id %in% lead_assoc$association_id) |>
    dplyr::pull("variant_id") |>
    unique()

  if (length(lead_variants) == 0L) {
    stop("No lead variants available after filtering.", call. = FALSE)
  }

  cache <- cat_assocs$cache

  variant_df <- cache |>
    dplyr::filter(
      !is.na(.data$chromosome_name),
      !is.na(.data$chromosome_position)
    ) |>
    dplyr::transmute(
      variant_id = .data$variant_id,
      chromosome_name = as.character(.data$chromosome_name),
      chromosome_position = as.integer(.data$chromosome_position)
    ) |>
    dplyr::distinct()

  if (nrow(variant_df) == 0L) {
    variant_df <- pull_v2_snp_details(
      rsids = lead_variants,
      verbose = verbose
    )
  }

  if (nrow(variant_df) == 0L) {
    stop("Variant annotation returned zero rows.", call. = FALSE)
  }

  progress_genomic(85, "building tables", verbose = verbose)

  context_df <- cache |>
    dplyr::filter(
      !is.na(.data$mapped_gene),
      .data$mapped_gene != ""
    ) |>
    tidyr::separate_rows(tidyselect::all_of("mapped_gene"), sep = ",|;") |>
    dplyr::transmute(
      variant_id = .data$variant_id,
      gene_name = trimws(.data$mapped_gene),
      is_mapped_gene = TRUE
    ) |>
    dplyr::filter(.data$gene_name != "") |>
    dplyr::distinct()

  if (nrow(context_df) == 0L ||
      !all(c("variant_id", "gene_name", "is_mapped_gene") %in% names(context_df))) {

    gene_nest <- tibble::tibble(
      variant_id = character(),
      genes = list()
    )

  } else {

    gene_nest <- context_df |>
      dplyr::filter(.data$is_mapped_gene) |>
      dplyr::select(tidyselect::all_of(c("variant_id", "gene_name"))) |>
      dplyr::distinct() |>
      tidyr::nest(genes = tidyselect::all_of("gene_name"))
  }

  variant_tbl <- variant_df |>
    dplyr::left_join(gene_nest, by = "variant_id") |>
    dplyr::mutate(
      genes = purrr::map(
        .data$genes,
        ~ if (is.null(.x) || nrow(.x) == 0) {
          tibble::tibble(gene_name = character())
        } else {
          .x
        }
      )
    ) |>
    dplyr::arrange(.data$chromosome_name, .data$chromosome_position)

  risk_map <- cat_assocs$risk_alleles |>
    dplyr::select(tidyselect::all_of(c("variant_id", "association_id"))) |>
    dplyr::distinct()

  assoc_extra <- cat_assocs$associations |>
    dplyr::select(tidyselect::all_of(c("association_id", "pvalue"))) |>
    dplyr::distinct()

  assoc2study <- risk_map |>
    dplyr::mutate(
      study_accession = sub(".*:", "", .data$association_id)
    ) |>
    dplyr::select(tidyselect::all_of(c("association_id", "study_accession"))) |>
    dplyr::distinct()

  variant_full <- variant_tbl |>
    dplyr::left_join(risk_map, by = "variant_id") |>
    dplyr::left_join(assoc_extra, by = "association_id") |>
    dplyr::left_join(assoc2study, by = "association_id") |>
    dplyr::distinct()

  summary_tbl <- variant_full |>
    dplyr::mutate(
      has_gene = purrr::map_lgl(.data$genes, ~ nrow(.x) > 0),
      gene_vec = purrr::map(.data$genes, ~ .x$gene_name)
    ) |>
    dplyr::summarise(
      n_SNPs       = dplyr::n_distinct(.data$variant_id),
      SNPs_w_gene  = dplyr::n_distinct(.data$variant_id[.data$has_gene]),
      unique_genes = dplyr::n_distinct(unlist(.data$gene_vec, use.names = FALSE)),
      n_studies    = dplyr::n_distinct(.data$study_accession, na.rm = TRUE)
    )

  chr_freq <- variant_full |>
    dplyr::count(.data$chromosome_name, name = "SNPs") |>
    dplyr::arrange(dplyr::desc(.data$SNPs)) |>
    dplyr::slice_head(n = 10)

  bed_df <- variant_full |>
    dplyr::filter(
      !is.na(.data$chromosome_name),
      !is.na(.data$chromosome_position)
    ) |>
    dplyr::mutate(
      chromosome_name = as.character(.data$chromosome_name),
      chromosome_name = dplyr::case_when(
        chromosome_name %in% c("23") ~ "X",
        chromosome_name %in% c("24") ~ "Y",
        chromosome_name %in% c("25", "MT", "M") ~ "M",
        TRUE ~ chromosome_name
      )
    ) |>
    dplyr::transmute(
      chr = ifelse(
        grepl("^chr", .data$chromosome_name),
        .data$chromosome_name,
        paste0("chr", .data$chromosome_name)
      ),
      start0 = as.integer(.data$chromosome_position) - 1L,
      end0   = as.integer(.data$chromosome_position),
      id     = .data$variant_id
    ) |>
    dplyr::filter(
      !is.na(.data$start0),
      !is.na(.data$end0),
      .data$start0 < .data$end0
    ) |>
    dplyr::distinct()

  fasta_set <- NULL
  written <- character(0L)

  progress_genomic(92, "writing outputs", verbose = verbose)

  if (!is.null(out_prefix)) {

    csv_path <- paste0(out_prefix, "_snps_full.csv")
    bed_path <- paste0(out_prefix, "_snps_hg38.bed")
    fa_path  <- paste0(out_prefix, "_snps_flank", flank_bp, ".fa")

    variant_export <- variant_full |>
      dplyr::mutate(
        genes = purrr::map_chr(
          .data$genes,
          ~ if (is.null(.x) || nrow(.x) == 0) {
            ""
          } else {
            paste(unique(.x$gene_name), collapse = ";")
          }
        )
      )

    readr::write_csv(variant_export, csv_path)
    readr::write_tsv(bed_df, bed_path, col_names = FALSE)

    if (requireNamespace(genome_pkg, quietly = TRUE) &&
        requireNamespace("Biostrings", quietly = TRUE)) {

      genome_obj <- getExportedValue(genome_pkg, "Hsapiens")

      if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
        available_chr <- GenomeInfoDb::seqnames(genome_obj)

        bed_fa <- bed_df |>
          dplyr::filter(.data$chr %in% available_chr)
      } else {
        bed_fa <- bed_df
      }

      if (nrow(bed_fa) > 0) {

        start1 <- pmax(bed_fa$start0 - flank_bp + 1L, 1L)
        end1   <- bed_fa$end0 + flank_bp

        fasta_set <- Biostrings::getSeq(
          genome_obj,
          names = bed_fa$chr,
          start = start1,
          end   = end1
        )

        names(fasta_set) <- bed_fa$id
        Biostrings::writeXStringSet(fasta_set, fa_path)

      } else {
        fa_path <- NA_character_
      }

    } else {
      fa_path <- NA_character_
    }

    written <- c(csv_path, bed_path, fa_path)
    written <- written[!is.na(written)]
    written <- as.character(written)
  }

  progress_genomic(100, "completed", verbose = verbose)

  invisible(list(
    summary   = summary_tbl,
    chr_freq  = chr_freq,
    snps_full = variant_full,
    bed       = bed_df,
    fasta     = fasta_set,
    written   = written
  ))
}
