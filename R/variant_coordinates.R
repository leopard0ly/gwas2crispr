variant_coordinate_template <- function() {
  tibble::tibble(
    variant_id = character(),
    chromosome_name = character(),
    chromosome_position = integer()
  )
}

standardize_variant_coordinates <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(variant_coordinate_template())
  }

  if (!"variant_id" %in% names(df)) {
    df$variant_id <- NA_character_
  }

  if (!"chromosome_name" %in% names(df)) {
    df$chromosome_name <- NA_character_
  }

  if (!"chromosome_position" %in% names(df)) {
    df$chromosome_position <- NA_integer_
  }

  tibble::tibble(
    variant_id = as.character(df$variant_id),
    chromosome_name = as.character(df$chromosome_name),
    chromosome_position = as.integer(df$chromosome_position)
  ) |>
    dplyr::filter(
      !is.na(.data$variant_id),
      .data$variant_id != "",
      !is.na(.data$chromosome_name),
      .data$chromosome_name != "",
      !is.na(.data$chromosome_position)
    ) |>
    dplyr::distinct(.data$variant_id, .keep_all = TRUE)
}

fill_missing_coordinates <- function(variant_df, ensembl_df) {
  existing <- standardize_variant_coordinates(variant_df)
  fallback <- standardize_variant_coordinates(ensembl_df)

  if (nrow(fallback) > 0L && nrow(existing) > 0L) {
    fallback <- fallback |>
      dplyr::filter(!.data$variant_id %in% existing$variant_id)
  }

  dplyr::bind_rows(existing, fallback) |>
    dplyr::distinct(.data$variant_id, .keep_all = TRUE)
}

first_mapping_value <- function(x) {
  vals <- unlist(x, recursive = TRUE, use.names = FALSE)
  vals <- vals[!is.na(vals)]

  if (length(vals) == 0L) {
    return(NA_character_)
  }

  as.character(vals[[1]])
}

pull_ensembl_variant_details <- function(rsids, verbose = interactive()) {
  rsids <- unique(as.character(rsids))
  rsids <- rsids[!is.na(rsids) & grepl("^rs[0-9]+$", rsids)]

  if (length(rsids) == 0L) {
    return(variant_coordinate_template())
  }

  if (isTRUE(verbose)) {
    message("Trying Ensembl REST coordinate fallback for ", length(rsids), " variant(s).")
  }

  primary_chr <- c(as.character(seq_len(22)), "X", "Y", "MT", "M")
  rows <- list()
  parse_one_variant <- function(rs, js) {
    mappings <- js$mappings

    if (is.null(mappings) || !is.list(mappings) || length(mappings) == 0L) {
      return(variant_coordinate_template())
    }

    map_rows <- list()

    for (mapping in mappings) {
      if (is.null(mapping) || !is.list(mapping)) {
        next
      }

      chr <- first_mapping_value(mapping$seq_region_name)
      pos <- suppressWarnings(as.integer(first_mapping_value(mapping$start)))
      assembly <- first_mapping_value(mapping$assembly_name)

      if (is.na(chr) || chr == "" || is.na(pos)) {
        next
      }

      map_rows[[length(map_rows) + 1L]] <- tibble::tibble(
        variant_id = rs,
        chromosome_name = chr,
        chromosome_position = pos,
        assembly_name = assembly
      )
    }

    if (length(map_rows) == 0L) {
      return(variant_coordinate_template())
    }

    candidate_rows <- dplyr::bind_rows(map_rows)
    grch38_rows <- candidate_rows |>
      dplyr::filter(!is.na(.data$assembly_name), .data$assembly_name == "GRCh38")

    if (nrow(grch38_rows) > 0L) {
      candidate_rows <- grch38_rows
    } else if (any(!is.na(candidate_rows$assembly_name) & candidate_rows$assembly_name != "")) {
      return(variant_coordinate_template())
    }

    candidate_rows |>
      dplyr::mutate(
        is_primary = .data$chromosome_name %in% primary_chr,
        chr_rank = match(.data$chromosome_name, primary_chr),
        chr_rank = ifelse(is.na(.data$chr_rank), Inf, .data$chr_rank)
      ) |>
      dplyr::arrange(
        dplyr::desc(.data$is_primary),
        .data$chr_rank,
        .data$chromosome_position
      ) |>
      dplyr::slice_head(n = 1L) |>
      dplyr::select(tidyselect::all_of(c(
        "variant_id",
        "chromosome_name",
        "chromosome_position"
      )))
  }

  chunk_size <- 200L
  chunks <- split(rsids, ceiling(seq_along(rsids) / chunk_size))

  for (chunk in chunks) {
    response <- tryCatch(
      httr::POST(
        "https://rest.ensembl.org/variation/human",
        body = list(ids = unname(chunk)),
        encode = "json",
        httr::add_headers(
          Accept = "application/json",
          `Content-Type` = "application/json"
        ),
        httr::timeout(30)
      ),
      error = function(e) NULL
    )

    if (is.null(response)) {
      next
    }

    status <- httr::status_code(response)

    if (status == 429) {
      Sys.sleep(2)
      next
    }

    if (status >= 400) {
      next
    }

    js <- tryCatch(
      httr::content(response, as = "parsed", type = "application/json"),
      error = function(e) NULL
    )

    if (is.null(js) || !is.list(js)) {
      next
    }

    for (rs in chunk) {
      if (is.null(js[[rs]]) || !is.list(js[[rs]])) {
        next
      }

      candidate_rows <- parse_one_variant(rs, js[[rs]])

      if (nrow(candidate_rows) > 0L) {
        rows[[length(rows) + 1L]] <- candidate_rows
      }
    }

    Sys.sleep(0.2)
  }

  if (length(rows) == 0L) {
    return(variant_coordinate_template())
  }

  dplyr::bind_rows(rows) |>
    dplyr::distinct(.data$variant_id, .keep_all = TRUE)
}
