normalize_trait_id <- function(x) {
  if (length(x) != 1L) {
    return(NA_character_)
  }

  x <- tryCatch(
    as.character(x),
    error = function(e) NA_character_
  )

  if (length(x) != 1L || is.na(x)) {
    return(NA_character_)
  }

  x <- trimws(x)
  gsub(":", "_", x, fixed = TRUE)
}

trait_id_prefix <- function(x) {
  x <- normalize_trait_id(x)

  if (is.na(x) || !nzchar(x) || !grepl("_", x, fixed = TRUE)) {
    return(NA_character_)
  }

  sub("_.*$", "", x)
}

is_supported_trait_id <- function(x) {
  x <- normalize_trait_id(x)

  if (is.na(x) || !nzchar(x)) {
    return(FALSE)
  }

  grepl(
    "^(EFO|MONDO|HP)_[0-9]+$|^NCIT_C[0-9]+$|^(Orphanet|ORPHA)_[0-9]+$",
    x,
    perl = TRUE
  )
}

validate_trait_id <- function(x, arg = "efo_id") {
  trait_id <- normalize_trait_id(x)

  if (!is.na(trait_id) && grepl("^GO_[0-9]+$", trait_id, perl = TRUE)) {
    stop(
      "GO identifiers are not supported as primary GWAS Catalog trait identifiers in gwas2crispr 0.1.5.",
      call. = FALSE
    )
  }

  if (!is_supported_trait_id(x)) {
    stop(
      sprintf(
        "%s must be a single supported GWAS Catalog trait identifier such as 'EFO_0001663', 'MONDO_0007254', or 'NCIT_C4872'.",
        arg
      ),
      call. = FALSE
    )
  }

  trait_id
}
