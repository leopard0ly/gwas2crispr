#' Fetch significant GWAS associations for a trait (EFO ID)
#'
#' @param efo_id     EFO trait identifier (e.g., "EFO_0001663" for prostate cancer)
#' @param p_cut      P-value threshold (default 5e-8)
#' @return           List with associations and risk alleles tables
#' @export
fetch_gwas <- function(efo_id = "EFO_0001663", p_cut = 5e-8) {
  # Try gwasrapidd first
  try({
    a <- gwasrapidd::get_associations(efo_id = efo_id, verbose = FALSE)
    if (nrow(a@associations) > 0) return(a)
  }, silent = TRUE)

  # REST fallback (EBI GWAS REST API)
  api <- "https://www.ebi.ac.uk/gwas/summary-statistics/api"
  url <- sprintf("%s/traits/%s/associations?p_upper=%g&size=5000", api, efo_id, p_cut)

  pull_all <- function(u) {
    out <- list()
    repeat {
      r  <- httr::GET(u, httr::add_headers(Accept = "application/hal+json, application/json"))
      httr::stop_for_status(r)
      js <- httr::content(r, as = "parsed", type = "application/json")
      if (!"associations" %in% names(js$`_embedded`)) break
      page <- purrr::flatten(js$`_embedded`$associations)
      if ("pvalue" %in% names(page) && !"p_value" %in% names(page))
        page <- dplyr::rename(page, p_value = pvalue)
      out[[length(out) + 1]] <- page
      nxt <- js$`_links`$`next`$href
      if (is.null(nxt) || identical(nxt, u)) break
      u <- nxt
    }
    dplyr::bind_rows(out)
  }

  ss <- pull_all(url)
  stopifnot(nrow(ss) > 0)

  # Return an object similar to gwasrapidd's output
  methods::new("associations",
               associations = tibble::tibble(
                 association_id = paste0(ss$variant_id, ":", ss$study_accession),
                 pvalue = as.numeric(ss$p_value)
               ),
               risk_alleles = tibble::tibble(
                 association_id = paste0(ss$variant_id, ":", ss$study_accession),
                 variant_id = ss$variant_id
               )
  )
}
