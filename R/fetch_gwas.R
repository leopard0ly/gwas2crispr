#' Fetch significant GWAS associations for an EFO trait
#'
#' @description
#' Tries \code{gwasrapidd::get_associations()} first; if it returns no rows
#' or fails, falls back to the EBI GWAS Summary Statistics REST API to
#' retrieve significant associations up to the given p-value threshold.
#'
#' @param efo_id character. Experimental Factor Ontology (EFO) trait identifier
#'   (e.g., "EFO_0001663").
#' @param p_cut numeric. P-value threshold for significance (default 5e-8).
#'
#' @return An S4 object of class \code{"associations"} with slots:
#' \itemize{
#'   \item \code{associations}: data frame with \code{association_id} and \code{pvalue}.
#'   \item \code{risk_alleles}: data frame mapping \code{association_id} to \code{variant_id}.
#' }
#'
#' @details
#' This function performs network calls and may be rate-limited. Column names
#' returned by the REST API may change; defensive checks are applied.
#'
#' @seealso \code{\link{run_gwas2crispr}}
#'
#' @examples
#' \donttest{
#'   # Network call; may be rate-limited, so we mark it as \donttest.
#'   a <- try(fetch_gwas("EFO_0001663", p_cut = 5e-8), silent = TRUE)
#'   if (!inherits(a, "try-error")) {
#'     head(a@associations)
#'   }
#' }
#'
#' @export
fetch_gwas <- function(efo_id = "EFO_0001663", p_cut = 5e-8) {
  # --- validate inputs -------------------------------------------------------
  if (!is.character(efo_id) || length(efo_id) != 1L || !grepl("^EFO_\\d+$", efo_id)) {
    stop("efo_id must be a single string like 'EFO_0001663'.", call. = FALSE)
  }
  if (!is.numeric(p_cut) || length(p_cut) != 1L || !is.finite(p_cut) || p_cut <= 0) {
    stop("p_cut must be a single positive numeric value.", call. = FALSE)
  }

  # --- 1) gwasrapidd first ---------------------------------------------------
  a <- try({
    gwasrapidd::get_associations(efo_id = efo_id, verbose = FALSE)
  }, silent = TRUE)

  if (!inherits(a, "try-error")) {
    if (!is.null(a@associations) && nrow(a@associations) > 0) {
      return(a)
    }
  }

  # --- 2) REST fallback (EBI GWAS REST API) ---------------------------------
  api <- "https://www.ebi.ac.uk/gwas/summary-statistics/api"
  url <- sprintf("%s/traits/%s/associations?p_upper=%g&size=5000", api, efo_id, p_cut)

  pull_all <- function(u) {
    out <- list()
    repeat {
      r  <- httr::GET(u, httr::add_headers(Accept = "application/hal+json, application/json"))
      httr::stop_for_status(r)
      js <- httr::content(r, as = "parsed", type = "application/json")
      emb <- js[["_embedded"]]
      if (is.null(emb) || is.null(emb$associations)) break

      page <- purrr::flatten(emb$associations)

      # normalize p-value naming if needed
      if ("pvalue" %in% names(page) && !"p_value" %in% names(page)) {
        page <- dplyr::rename(page, p_value = pvalue)
      }

      out[[length(out) + 1L]] <- page
      nxt <- js[["_links"]][["next"]][["href"]]
      if (is.null(nxt) || identical(nxt, u)) break
      u <- nxt
    }
    if (length(out) == 0L) return(tibble::tibble())
    dplyr::bind_rows(out)
  }

  ss <- pull_all(url)
  if (nrow(ss) == 0L) {
    stop("No associations returned for the specified trait and threshold.", call. = FALSE)
  }

  # --- 3) return an object analogous to gwasrapidd::get_associations() -------
  methods::new(
    "associations",
    associations = tibble::tibble(
      association_id = paste0(ss$variant_id, ":", ss$study_accession),
      pvalue = suppressWarnings(as.numeric(ss$p_value))
    ),
    risk_alleles = tibble::tibble(
      association_id = paste0(ss$variant_id, ":", ss$study_accession),
      variant_id = ss$variant_id
    )
  )
}
