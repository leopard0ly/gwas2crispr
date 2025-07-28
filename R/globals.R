# internal bindings for NSE (dplyr/tidyr) to satisfy R CMD check
utils::globalVariables(c(
  "pvalue",
  "association_id",
  "variant_id",
  "is_mapped_gene",
  "gene_name",
  "genes",
  "chromosome_name",
  "chromosome_position",
  "has_gene",
  "gene_vec",
  "study_accession",
  "SNPs"
))
