#' GWAS → CRISPR Target Builder (Royal, Battle-Tested)
#'
#' @param efo_id     EFO trait identifier (e.g. "EFO_0001663" for prostate cancer)
#' @param p_cut      P-value threshold for significant SNPs (default 5e-8)
#' @param out_prefix Prefix for output files (default "output")
#' @param flank_bp   Number of bp flanking each SNP for FASTA (default 200)
#' @param genome_pkg Name of BSgenome package (default "BSgenome.Hsapiens.UCSC.hg38")
#' @export
run_gwas2crispr <- function(
    efo_id,
    p_cut,
    flank_bp,
    out_prefix = "output",
    genome_pkg = "BSgenome.Hsapiens.UCSC.hg38"
) {
  if (!identical(genome_pkg, "BSgenome.Hsapiens.UCSC.hg38"))
    stop("Only GRCh38/hg38 is supported. استخدم BSgenome.Hsapiens.UCSC.hg38", call. = FALSE)
  # 1. Fetch GWAS associations
  cat("Fetching GWAS associations for:", efo_id, "\n")
  cat_assocs <- fetch_gwas(efo_id, p_cut)
  lead_assoc <- cat_assocs@associations |>
    dplyr::filter(pvalue < p_cut)
  lead_variants <- cat_assocs@risk_alleles |>
    dplyr::filter(association_id %in% lead_assoc$association_id) |>
    dplyr::pull(variant_id) |> unique()

  # 2. Variant & gene metadata
  get_variant_chunks <- function(v, size = 100) {
    split(v, ceiling(seq_along(v)/size)) |>
      purrr::map(~gwasrapidd::get_variants(variant_id = .x, verbose = FALSE))
  }
  v_chunks   <- get_variant_chunks(lead_variants)
  variant_df <- dplyr::bind_rows(purrr::map(v_chunks, ~.x@variants))
  context_df <- dplyr::bind_rows(purrr::map(v_chunks, ~.x@genomic_contexts))

  gene_nest <- context_df |>
    dplyr::filter(is_mapped_gene) |>
    dplyr::select(variant_id, gene_name) |>
    dplyr::distinct() |>
    tidyr::nest(genes = gene_name)

  variant_tbl <- variant_df |>
    dplyr::left_join(gene_nest, by = "variant_id") |>
    dplyr::mutate(genes = purrr::map(genes, ~if (is.null(.x) || nrow(.x) == 0)
      tibble::tibble(gene_name=character()) else .x)) |>
    dplyr::arrange(chromosome_name, chromosome_position)

  # 3. Study & effect metrics
  assoc2study <- cat_assocs@risk_alleles |>
    dplyr::transmute(variant_id,
                     study_accession = sub(".*:", "", association_id)) |>
    dplyr::distinct()
  effect_cols <- c("risk_allele", "other_allele", "beta", "odds_ratio", "effect_allele_frequency", "p_value")
  present_cols <- intersect(names(cat_assocs@associations), effect_cols)
  assoc_extra <- cat_assocs@associations |>
    dplyr::select(association_id, dplyr::all_of(present_cols)) |>
    dplyr::distinct()
  risk_map <- cat_assocs@risk_alleles |>
    dplyr::select(variant_id, association_id) |> dplyr::distinct()
  study_raw <- tryCatch(cat_assocs@studies, error = function(e) NULL)
  study_df <- if (is.null(study_raw) || nrow(study_raw) == 0) {
    tibble::tibble(study_accession = character(), pubmed_id = character(),
                   publication_date = as.Date(character()))
  } else {
    tibble::as_tibble(study_raw) |>
      dplyr::select(dplyr::any_of(c("study_accession", "pubmed_id",
                                    "publication_date", "sample_size")))
  }
  variant_full <- variant_tbl |>
    dplyr::left_join(risk_map,    by = "variant_id") |>
    dplyr::left_join(assoc2study, by = "variant_id") |>
    dplyr::left_join(assoc_extra, by = "association_id") |>
    dplyr::left_join(study_df,    by = "study_accession") |> dplyr::distinct()

  # 4. Summary
  summary_tbl <- variant_full |>
    dplyr::mutate(has_gene = purrr::map_lgl(genes, ~nrow(.x) > 0),
                  gene_vec = purrr::map(genes, ~.x$gene_name)) |>
    dplyr::summarise(
      n_SNPs         = dplyr::n_distinct(variant_id),
      SNPs_w_gene    = sum(has_gene),
      unique_genes   = dplyr::n_distinct(unlist(gene_vec)),
      n_studies      = dplyr::n_distinct(study_accession, na.rm = TRUE)
    )
  chr_freq <- variant_full |>
    dplyr::count(chromosome_name, name = "SNPs") |>
    dplyr::arrange(desc(SNPs)) |>
    dplyr::slice_head(n = 10)

  cat("\n===== GWAS → CRISPR Summary =====\n")
  print(summary_tbl)
  cat("\nTop-10 chromosomes by SNP count:\n")
  print(chr_freq)

  # 5. Output files
  csv_path <- paste0(out_prefix, "_snps_full.csv")
  bed_path <- paste0(out_prefix, "_snps_hg38.bed")
  fa_path  <- paste0(out_prefix, "_snps_flank", flank_bp, ".fa")

  readr::write_csv(variant_full, csv_path)

  bed_df <- variant_full |>
    dplyr::filter(!is.na(chromosome_name)) |>
    dplyr::transmute(chr = paste0("chr", chromosome_name),
                     start0 = chromosome_position - 1,
                     end0   = chromosome_position,
                     id     = variant_id)
  readr::write_tsv(bed_df, bed_path, col_names = FALSE)

  # FASTA ±flank_bp bp (must have BSgenome installed)
  if (requireNamespace(genome_pkg, quietly = TRUE)) {
    genome_obj <- getExportedValue(genome_pkg, "Hsapiens")
    start1 <- pmax(bed_df$start0 - flank_bp + 1, 1)
    end1   <- bed_df$end0 + flank_bp
    seqs   <- Biostrings::getSeq(genome_obj,
                                 names = bed_df$chr,
                                 start = start1,
                                 end   = end1)
    names(seqs) <- bed_df$id
    Biostrings::writeXStringSet(seqs, fa_path)
  } else {
    cat(genome_pkg, "not installed: FASTA step skipped\n")
  }

  cat("Done. Output files:\n", csv_path, "\n", bed_path, "\n", fa_path, "\n")
  invisible(list(csv = csv_path, bed = bed_path, fasta = fa_path,
                 summary = summary_tbl, chr_freq = chr_freq))
}
