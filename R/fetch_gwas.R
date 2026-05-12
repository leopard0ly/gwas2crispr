progress_genomic <- function(percent, label = "", verbose = interactive()) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }

  percent <- max(0, min(100, round(percent)))
  icon <- ""

  cat(sprintf(
    "\r%sLoading genomic data for GWAS-to-CRISPR analysis... %3d%%  %s",
    icon,
    percent,
    label
  ))

  flush.console()

  if (percent >= 100) {
    cat("\n")
  }

  invisible(NULL)
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

clean_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

extract_first_rs <- function(x) {
  stringr::str_extract(as.character(x), "rs[0-9]+")
}

find_col <- function(df, candidates) {
  nms <- names(df)
  nms_clean <- clean_name(nms)
  cand_clean <- clean_name(candidates)

  hit <- which(nms_clean %in% cand_clean)

  if (length(hit) > 0) {
    return(nms[hit[1]])
  }

  for (cc in cand_clean) {
    hit2 <- which(grepl(cc, nms_clean, fixed = TRUE))

    if (length(hit2) > 0) {
      return(nms[hit2[1]])
    }
  }

  NA_character_
}

flatten_atomic <- function(x, prefix = "") {
  rec <- function(obj, path) {
    if (is.null(obj)) {
      return(character(0))
    }

    if (is.atomic(obj)) {
      vals <- as.character(obj)
      vals <- vals[!is.na(vals)]

      if (length(vals) == 0) {
        return(character(0))
      }

      nms <- paste0(
        path,
        ifelse(length(vals) > 1, paste0("[", seq_along(vals), "]"), "")
      )

      stats::setNames(vals, nms)

    } else if (is.list(obj)) {
      y <- character(0)
      nm <- names(obj)

      if (is.null(nm)) {
        nm <- as.character(seq_along(obj))
      }

      for (i in seq_along(obj)) {
        child_path <- if (nzchar(path)) {
          paste0(path, ".", nm[i])
        } else {
          nm[i]
        }

        y <- c(y, rec(obj[[i]], child_path))
      }

      y
    } else {
      character(0)
    }
  }

  rec(x, prefix)
}

as_flat_list <- function(x) {
  if (is.data.frame(x)) {
    x <- as.list(x[1, , drop = FALSE])
  }

  if (is.atomic(x) && !is.null(names(x))) {
    x <- as.list(x)
  }

  if (!is.list(x)) {
    return(list())
  }

  x
}

first_flat_value <- function(x, patterns, exclude = character()) {
  x <- as_flat_list(x)
  nms <- names(x)

  if (is.null(nms) || !length(nms)) {
    return(NA_character_)
  }

  for (pat in patterns) {
    idx <- which(grepl(pat, nms, ignore.case = TRUE, perl = TRUE))

    if (length(exclude)) {
      ex <- paste(exclude, collapse = "|")
      idx <- idx[!grepl(ex, nms[idx], ignore.case = TRUE, perl = TRUE)]
    }

    if (!length(idx)) {
      next
    }

    for (i in idx) {
      val <- unlist(x[[i]], recursive = TRUE, use.names = FALSE)
      val <- val[!is.na(val)]

      if (!length(val)) {
        next
      }

      val_chr <- trimws(as.character(val))
      val_chr <- val_chr[nzchar(val_chr)]

      if (length(val_chr)) {
        return(val_chr[[1]])
      }
    }
  }

  NA_character_
}

all_flat_values <- function(x, patterns) {
  x <- as_flat_list(x)
  nms <- names(x)

  if (is.null(nms) || !length(nms)) {
    return(character())
  }

  idx <- integer()

  for (pat in patterns) {
    idx <- c(idx, which(grepl(pat, nms, ignore.case = TRUE, perl = TRUE)))
  }

  idx <- unique(idx)

  if (!length(idx)) {
    return(character())
  }

  vals <- unlist(x[idx], recursive = TRUE, use.names = FALSE)
  vals <- vals[!is.na(vals)]
  vals <- trimws(as.character(vals))
  vals[nzchar(vals)]
}

parse_numeric_value <- function(x) {
  if (length(x) == 0 || is.na(x)) {
    return(NA_real_)
  }

  x <- trimws(as.character(x[[1]]))
  x <- gsub(",", "", x, fixed = TRUE)

  out <- suppressWarnings(as.numeric(x))

  if (!is.na(out)) {
    return(out)
  }

  sci <- regmatches(
    x,
    regexpr(
      "[0-9]+\\.?[0-9]*\\s*([eE]|x|X|x)\\s*10?\\^?\\s*-?[0-9]+|[0-9]+\\.?[0-9]*[eE]-?[0-9]+",
      x,
      perl = TRUE
    )
  )

  if (length(sci) && nzchar(sci)) {
    sci <- gsub("\\s+", "", sci)
    sci <- gsub("x10\\^?", "e", sci)
    sci <- gsub("[xX]10\\^?", "e", sci)
    sci <- gsub("x", "e", sci)
    sci <- gsub("[xX]", "e", sci)

    out <- suppressWarnings(as.numeric(sci))

    if (!is.na(out)) {
      return(out)
    }
  }

  NA_real_
}

get_pvalue_from_flat <- function(f) {
  if (length(f) == 0) {
    return(NA_real_)
  }

  valid_p <- function(z) {
    is.finite(z) && !is.na(z) && z > 0 && z <= 1
  }

  direct <- first_flat_value(
    f,
    patterns = c(
      "(^|[._])p[_\\.]*value($|[._])",
      "(^|[._])pvalue($|[._])",
      "(^|[._])p[_\\.]*val($|[._])",
      "p[_\\.]*value$",
      "pvalue$",
      "pval$"
    ),
    exclude = c("mantissa", "exponent")
  )

  direct_num <- parse_numeric_value(direct)

  if (!is.na(direct_num)) {
    if (valid_p(direct_num)) {
      return(direct_num)
    }
    return(NA_real_)
  }

  mantissa <- first_flat_value(
    f,
    patterns = c(
      "p[_\\.]?value.*mantissa$",
      "pvalue.*mantissa$",
      "p[_\\.]?val.*mantissa$",
      "mantissa$"
    )
  )

  exponent <- first_flat_value(
    f,
    patterns = c(
      "p[_\\.]?value.*exponent$",
      "pvalue.*exponent$",
      "p[_\\.]?val.*exponent$",
      "exponent$"
    )
  )

  mantissa_num <- parse_numeric_value(mantissa)
  exponent_num <- parse_numeric_value(exponent)

  if (!is.na(mantissa_num) && !is.na(exponent_num)) {
    p <- mantissa_num * 10 ^ exponent_num

    if (valid_p(p)) {
      return(p)
    }

    return(NA_real_)
  }

  nm <- clean_name(names(f))

  mant_i <- which(grepl("pvalue", nm) & grepl("mantissa", nm))
  exp_i  <- which(grepl("pvalue", nm) & grepl("exponent", nm))

  if (length(mant_i) > 0 && length(exp_i) > 0) {
    mant <- safe_num(f[mant_i[1]])
    expo <- safe_num(f[exp_i[1]])

    if (is.finite(mant) && is.finite(expo)) {
      p <- mant * 10 ^ expo

      if (valid_p(p)) {
        return(p)
      }

      return(NA_real_)
    }
  }

  p_i <- which(
    grepl("pvalue", nm) |
      grepl("^p$", nm) |
      grepl("pval", nm)
  )

  p_i <- p_i[!grepl("mantissa|exponent", nm[p_i])]

  if (length(p_i) > 0) {
    nums <- safe_num(f[p_i])
    nums <- nums[is.finite(nums) & nums > 0 & nums <= 1]

    if (length(nums) > 0) {
      return(nums[1])
    }
  }

  NA_real_
}

get_study_from_flat <- function(f) {
  vals <- as.character(f)
  hit <- stringr::str_extract(vals, "GCST[0-9]+")
  hit <- unique(na.omit(hit))

  if (length(hit) > 0) {
    return(hit[1])
  }

  NA_character_
}

get_chr_from_flat <- function(f) {
  nm <- clean_name(names(f))

  idx <- which(
    grepl("chromosomename", nm) |
      grepl("chromosome", nm) |
      grepl("chrid", nm) |
      grepl("^chr$", nm)
  )

  if (length(idx) == 0) {
    return(NA_character_)
  }

  vals <- as.character(f[idx])
  vals <- vals[!is.na(vals)]
  vals <- vals[vals != ""]
  vals <- vals[grepl("^[0-9]+$|^X$|^Y$|^MT$|^M$", vals, ignore.case = TRUE)]

  if (length(vals) == 0) {
    return(NA_character_)
  }

  vals[1]
}

get_pos_from_flat <- function(f) {
  if (length(f) == 0) {
    return(NA_integer_)
  }

  pos <- first_flat_value(
    f,
    patterns = c(
      "base[_\\.]?pair.*location$",
      "basepairlocation$",
      "(^|[._])bp[_\\.]?location$",
      "(^|[._])position$",
      "(^|[._])pos$",
      "chromosome.*position$",
      "chromosome[_\\.]?position$",
      "chr.*position$",
      "chr.*pos$"
    )
  )

  pos_num <- parse_numeric_value(pos)

  if (!is.na(pos_num)) {
    return(as.integer(pos_num))
  }

  if (length(pos) && !is.na(pos)) {
    nums <- regmatches(
      as.character(pos),
      gregexpr("[0-9]+", as.character(pos), perl = TRUE)
    )[[1]]

    if (length(nums)) {
      return(as.integer(tail(nums, 1)))
    }
  }

  nm <- clean_name(names(f))

  idx <- which(
    grepl("chromosomeposition", nm) |
      grepl("bplocation", nm) |
      grepl("basepair", nm) |
      grepl("chrpos", nm) |
      grepl("^position$", nm)
  )

  if (length(idx) == 0) {
    return(NA_integer_)
  }

  vals <- safe_num(f[idx])
  vals <- vals[is.finite(vals) & vals > 0]

  if (length(vals) == 0) {
    return(NA_integer_)
  }

  as.integer(vals[1])
}

get_genes_from_flat <- function(f) {
  if (length(f) == 0) {
    return("")
  }

  vals <- all_flat_values(
    f,
    patterns = c(
      "mapped.*gene",
      "reported.*gene",
      "gene.*name",
      "gene[_\\.]?name",
      "(^|[._])gene($|[._])",
      "genes",
      "ensembl.*gene"
    )
  )

  if (!length(vals)) {
    nm <- clean_name(names(f))

    idx <- which(
      grepl("genename", nm) |
        grepl("mappedgene", nm) |
        grepl("reportedgene", nm) |
        grepl("gene", nm)
    )

    if (length(idx) > 0) {
      vals <- as.character(f[idx])
    }
  }

  vals <- vals[!is.na(vals)]
  vals <- vals[vals != ""]

  if (!length(vals)) {
    return("")
  }

  vals <- paste(vals, collapse = " ")

  vals <- gsub("https?://[^[:space:],;|]+", " ", vals, perl = TRUE)
  vals <- gsub("\\bENSG[0-9]+(\\.[0-9]+)?\\b", " ", vals, perl = TRUE)
  vals <- gsub("\\bENST[0-9]+(\\.[0-9]+)?\\b", " ", vals, perl = TRUE)
  vals <- gsub("[\r\n\t]", " ", vals)

  pieces <- unlist(strsplit(vals, "[,;|/\\s]+", perl = TRUE))
  pieces <- trimws(pieces)
  pieces <- gsub("^[()\\[\\]{}]+|[()\\[\\]{}]+$", "", pieces, perl = TRUE)
  pieces <- pieces[nzchar(pieces)]

  if (!length(pieces)) {
    return("")
  }

  drop_words <- c(
    "http", "https", "www", "ensembl", "org", "ebi", "ac", "uk",
    "gene", "genes", "mapped", "reported", "id", "summary",
    "homo", "sapiens", "grch38", "hg38", "api", "rest"
  )

  pieces <- pieces[!tolower(pieces) %in% drop_words]
  pieces <- pieces[!grepl("^ENSG", pieces)]
  pieces <- pieces[!grepl("^ENST", pieces)]
  pieces <- pieces[!grepl("^https?", pieces)]

  pieces <- pieces[grepl("^[A-Za-z][A-Za-z0-9._-]{1,30}$", pieces)]
  pieces <- unique(pieces)

  if (!length(pieces)) {
    return("")
  }

  paste(pieces, collapse = ";")
}

as_record_list <- function(x) {
  if (is.null(x) || !is.list(x) || length(x) == 0L) {
    return(list())
  }

  if (!all(vapply(x, is.list, logical(1)))) {
    return(list())
  }

  nms <- names(x)

  if (!is.null(nms)) {
    meta_names <- c("page", "links", "embedded", "metadata")

    if (all(clean_name(nms) %in% meta_names)) {
      return(list())
    }
  }

  unname(x)
}

extract_items <- function(js) {
  if (is.null(js) || !is.list(js)) {
    return(list())
  }

  emb <- js[["_embedded"]]

  if (!is.null(emb) && is.list(emb)) {
    for (key in c("associations", "singleNucleotidePolymorphisms", "snps", "efoTraits")) {
      items <- as_record_list(emb[[key]])

      if (length(items) > 0L) {
        return(items)
      }
    }

    if (length(emb) == 1L) {
      items <- as_record_list(emb[[1]])

      if (length(items) > 0L) {
        return(items)
      }
    }
  }

  for (key in c("content", "associations", "singleNucleotidePolymorphisms", "snps", "efoTraits")) {
    items <- as_record_list(js[[key]])

    if (length(items) > 0L) {
      return(items)
    }
  }

  as_record_list(js)
}

get_json <- function(url, query = list(), sleep_on_429 = 10, timeout_sec = 120) {
  repeat {
    r <- tryCatch(
      httr::GET(
        url,
        query = query,
        httr::add_headers(
          Accept = "application/json",
          `User-Agent` = "gwas2crispr-direct-v2/0.1.5"
        ),
        httr::timeout(timeout_sec)
      ),
      error = function(e) NULL
    )

    if (is.null(r)) {
      return(NULL)
    }

    status <- httr::status_code(r)

    if (status == 429) {
      Sys.sleep(sleep_on_429)
      next
    }

    if (status >= 400) {
      return(NULL)
    }

    return(tryCatch(
      httr::content(r, as = "parsed", type = "application/json"),
      error = function(e) NULL
    ))
  }
}

trait_id_aliases <- function(trait_id) {
  trait_id <- normalize_trait_id(trait_id)

  if (is.na(trait_id) || !nzchar(trait_id)) {
    return(character(0))
  }

  aliases <- trait_id
  prefix <- trait_id_prefix(trait_id)
  accession <- sub("^[^_]+_", "", trait_id)

  if (identical(prefix, "Orphanet")) {
    aliases <- c(aliases, paste0("ORPHA_", accession))
  } else if (identical(prefix, "ORPHA")) {
    aliases <- c(aliases, paste0("Orphanet_", accession))
  }

  unique(aliases)
}

dedupe_query_plan <- function(plan) {
  if (length(plan) == 0L) {
    return(plan)
  }

  keys <- vapply(
    plan,
    function(item) {
      query <- item$query
      paste(names(query), query, sep = "=", collapse = "&")
    },
    character(1)
  )

  plan[!duplicated(keys)]
}

association_query_plan <- function(trait_id, trait_labels = character()) {
  ids <- trait_id_aliases(trait_id)

  if (length(ids) == 0L) {
    return(list())
  }

  plan <- list(
    list(stage = 1L, route = "direct_id", value = ids[[1]], query = list(efo_id = ids[[1]])),
    list(stage = 2L, route = "trait_id", value = ids[[1]], query = list(efo_trait = ids[[1]]))
  )

  if (length(ids) > 1L) {
    for (alias in ids[-1]) {
      plan[[length(plan) + 1L]] <- list(
        stage = 3L,
        route = "alias_direct_id",
        value = alias,
        query = list(efo_id = alias)
      )
      plan[[length(plan) + 1L]] <- list(
        stage = 3L,
        route = "alias_trait_id",
        value = alias,
        query = list(efo_trait = alias)
      )
    }
  }

  labels <- unique(trimws(as.character(trait_labels)))
  labels <- labels[!is.na(labels) & labels != "" & !grepl("^http", labels)]

  for (label in labels) {
    plan[[length(plan) + 1L]] <- list(
      stage = 4L,
      route = "resolved_label",
      value = label,
      query = list(efo_trait = label)
    )
  }

  dedupe_query_plan(plan)
}

resolve_trait_labels <- function(trait_id, verbose = interactive()) {
  labels <- character(0)
  ids <- trait_id_aliases(trait_id)

  # Always keep identifier candidates themselves. GWAS Catalog v2 can use
  # identifiers through trait filters, which is more stable than relying only
  # on resolved text labels.
  labels <- c(labels, ids)

  known_labels <- list(
    EFO_0001663 = c("prostate cancer", "prostate carcinoma"),
    EFO_0000707 = c("lung cancer", "lung carcinoma", "lung disease"),
    EFO_0001071 = character(0),
    EFO_0001072 = character(0)
  )

  for (id in ids) {
    if (!is.null(known_labels[[id]])) {
      labels <- c(labels, known_labels[[id]])
    }
  }

  extract_trait_labels <- function(js) {
    out <- character(0)

    if (is.null(js)) {
      return(out)
    }

    items <- extract_items(js)

    if (length(items) == 0L) {
      return(out)
    }

    for (it in items) {
      f <- flatten_atomic(it)
      vals <- as.character(f)
      vals <- vals[!is.na(vals)]
      vals <- vals[vals != ""]
      vals <- vals[!grepl("^https?://", vals)]

      nm <- clean_name(names(f))
      keep <- grepl("trait", nm) | grepl("label", nm) | grepl("name", nm)

      out <- c(out, vals[keep])
    }

    out
  }

  query_list <- list()

  for (id in ids) {
    query_list[[length(query_list) + 1L]] <- list(efo_id = id, size = 20, page = 0)
    query_list[[length(query_list) + 1L]] <- list(efo_trait = id, size = 20, page = 0)
    query_list[[length(query_list) + 1L]] <- list(query = id, size = 20, page = 0)
  }

  for (qq in query_list) {
    js <- get_json(
      "https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits",
      query = qq
    )

    labels <- c(labels, extract_trait_labels(js))
  }

  labels <- unique(trimws(labels))
  labels <- labels[labels != ""]
  labels <- labels[!grepl("^http", labels)]

  labels
}

resolve_efo_labels <- function(efo_id, verbose = interactive()) {
  resolve_trait_labels(efo_id, verbose = verbose)
}

pull_v2_associations <- function(efo_id,
                                 trait_labels,
                                 p_cut,
                                 verbose = interactive()) {
  base_url <- "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"

  one_query <- function(plan_item,
                        page_size = 500L,
                        sleep_sec = 0.35,
                        base_progress = 10,
                        max_progress = 70) {
    query <- plan_item$query
    direct_id_route <- plan_item$route %in% c("direct_id", "alias_direct_id")
    rows <- list()
    page <- 0L
    total_pages <- Inf

    repeat {
      js <- get_json(
        base_url,
        query = c(
          query,
          list(
            show_child_traits = "true",
            size = page_size,
            page = page
          )
        )
      )

      if (is.null(js)) {
        break
      }

      if (!is.null(js$page$totalPages)) {
        total_pages <- as.integer(js$page$totalPages)

        if (isTRUE(direct_id_route) &&
            page == 0L &&
            is.finite(total_pages) &&
            total_pages > 50L) {
          break
        }

        if (is.finite(total_pages) && total_pages > 0) {
          p <- base_progress + ((page + 1) / total_pages) *
            (max_progress - base_progress)
          progress_genomic(p, "GWAS Catalog v2", verbose = verbose)
        }
      }

      items <- extract_items(js)

      if (length(items) == 0L) {
        break
      }

      if (isTRUE(direct_id_route)) {
        item_has_id <- vapply(
          items,
          function(it) {
            vals <- as.character(flatten_atomic(it))
            vals <- vals[!is.na(vals)]
            any(
              vals == plan_item$value |
                vals == gsub("_", ":", plan_item$value, fixed = TRUE)
            )
          },
          logical(1)
        )

        if (!any(item_has_id)) {
          break
        }

        items <- items[item_has_id]
      }

      for (it in items) {
        f <- flatten_atomic(it)

        pval <- get_pvalue_from_flat(f)

        if (!is.finite(pval) || is.na(pval) || pval >= p_cut) {
          next
        }

        all_vals <- as.character(f)
        rsids <- unique(na.omit(stringr::str_extract(all_vals, "rs[0-9]+")))
        rsids <- rsids[grepl("^rs[0-9]+$", rsids)]

        if (length(rsids) == 0) {
          next
        }

        study <- get_study_from_flat(f)

        if (is.na(study) || study == "") {
          study <- paste0("study_", length(rows) + 1L)
        }

        chr <- get_chr_from_flat(f)
        pos <- get_pos_from_flat(f)
        genes <- get_genes_from_flat(f)

        for (rs in rsids) {
          rows[[length(rows) + 1L]] <- tibble::tibble(
            variant_id = rs,
            study_accession = study,
            p_numeric = pval,
            chromosome_name = chr,
            chromosome_position = pos,
            mapped_gene = genes
          )
        }
      }

      page <- page + 1L

      if (page >= total_pages) {
        break
      }

      Sys.sleep(sleep_sec)
    }

    if (length(rows) == 0L) {
      return(tibble::tibble())
    }

    dplyr::bind_rows(rows)
  }

  all_rows <- list()
  query_plan <- association_query_plan(efo_id, trait_labels)

  if (length(query_plan) == 0) {
    return(tibble::tibble())
  }

  n_queries <- length(query_plan)

  for (stage in sort(unique(vapply(query_plan, `[[`, integer(1), "stage")))) {
    stage_rows <- list()
    stage_idx <- which(vapply(query_plan, `[[`, integer(1), "stage") == stage)

    for (i in stage_idx) {
      base_p <- 10 + ((i - 1) / n_queries) * 60
      max_p  <- 10 + (i / n_queries) * 60

      stage_rows[[length(stage_rows) + 1L]] <- one_query(
        query_plan[[i]],
        base_progress = base_p,
        max_progress = max_p
      )
    }

    stage_df <- dplyr::bind_rows(stage_rows)

    if (nrow(stage_df) > 0L) {
      all_rows[[length(all_rows) + 1L]] <- stage_df
      break
    }
  }

  ss <- dplyr::bind_rows(all_rows)

  if (nrow(ss) == 0L) {
    return(tibble::tibble())
  }

  ss |>
    dplyr::filter(
      !is.na(.data$p_numeric),
      .data$p_numeric < p_cut,
      !is.na(.data$variant_id),
      .data$variant_id != ""
    ) |>
    dplyr::distinct(.data$variant_id, .data$study_accession, .keep_all = TRUE)
}

pull_v2_snp_details <- function(rsids, verbose = interactive()) {
  base_url <- "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms"

  rows <- list()
  parse_snp_items <- function(items, rs) {
    out <- list()

    if (length(items) == 0L) {
      return(out)
    }

    for (it in items) {
      f <- flatten_atomic(it)

      vals <- as.character(f)
      hit_rs <- unique(na.omit(stringr::str_extract(vals, "rs[0-9]+")))
      hit_rs <- hit_rs[hit_rs == rs]

      if (length(hit_rs) == 0) {
        next
      }

      chr <- get_chr_from_flat(f)
      pos <- get_pos_from_flat(f)
      genes <- get_genes_from_flat(f)

      out[[length(out) + 1L]] <- tibble::tibble(
        variant_id = rs,
        chromosome_name = chr,
        chromosome_position = pos,
        mapped_gene = genes
      )
    }

    out
  }

  if (length(rsids) == 0) {
    return(tibble::tibble())
  }

  for (i in seq_along(rsids)) {
    progress_genomic(
      70 + (i / length(rsids)) * 15,
      "variant metadata",
      verbose = verbose
    )

    rs <- rsids[i]

    js <- get_json(
      paste0(base_url, "/", utils::URLencode(rs, reserved = TRUE)),
      sleep_on_429 = 2,
      timeout_sec = 20
    )

    items <- extract_items(js)

    if (length(items) == 0L && !is.null(js) && is.list(js)) {
      items <- list(js)
    }

    if (length(items) == 0L) {
      js <- get_json(
        base_url,
        query = list(
          rs_id = rs,
          size = 5,
          page = 0
        ),
        sleep_on_429 = 2,
        timeout_sec = 20
      )

      items <- extract_items(js)
    }

    parsed_rows <- parse_snp_items(items, rs)

    if (length(parsed_rows) > 0L) {
      rows <- c(rows, parsed_rows)
    }

    Sys.sleep(0.05)
  }

  if (length(rows) == 0L) {
    return(tibble::tibble())
  }

  dplyr::bind_rows(rows) |>
    dplyr::filter(
      !is.na(.data$variant_id),
      !is.na(.data$chromosome_name),
      !is.na(.data$chromosome_position)
    ) |>
    dplyr::distinct(.data$variant_id, .keep_all = TRUE)
}

#' Fetch significant GWAS associations for a GWAS Catalog trait identifier
#'
#' @description
#' Retrieves significant GWAS Catalog associations directly from the
#' EMBL-EBI GWAS Catalog REST API v2. The function resolves the supplied
#' GWAS Catalog trait identifier to direct identifier queries and trait labels,
#' retrieves paginated association records, filters by p-value, and returns a
#' list used by \code{\link{run_gwas2crispr}}.
#'
#' @param efo_id character. GWAS Catalog trait identifier. The argument name is
#'   retained for backward compatibility. Examples include EFO_0001663,
#'   MONDO_0007254, and NCIT_C4872 when supported by the GWAS Catalog API.
#' @param p_cut numeric. P-value threshold for significance.
#' @param verbose logical. If \code{TRUE}, prints a compact progress line.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{associations}: tibble with \code{association_id} and \code{pvalue}.
#'   \item \code{risk_alleles}: tibble mapping \code{association_id} to \code{variant_id}.
#'   \item \code{cache}: internal tibble with variant metadata used downstream.
#' }
#'
#' @details
#' This function performs network calls to the GWAS Catalog REST API v2 and may
#' be affected by service availability or rate limits. Selected supported
#' disease and cancer trait identifier prefixes include EFO, MONDO, and NCIT.
#' HP, Orphanet, and ORPHA are accepted for compatibility. GO identifiers are
#' not supported as primary GWAS Catalog trait identifiers in gwas2crispr 0.1.5.
#'
#' @seealso \code{\link{run_gwas2crispr}}
#'
#' @examples
#' \donttest{
#'   a <- fetch_gwas("EFO_0000707", p_cut = 1e-6, verbose = FALSE)
#'   head(a$associations)
#' }
#'
#' @export
fetch_gwas <- function(efo_id = "EFO_0001663",
                       p_cut = 5e-8,
                       verbose = interactive()) {
  trait_id <- validate_trait_id(efo_id, arg = "efo_id")

  if (!is.numeric(p_cut) ||
      length(p_cut) != 1L ||
      !is.finite(p_cut) ||
      p_cut <= 0) {
    stop("p_cut must be a single positive numeric value.", call. = FALSE)
  }

  progress_genomic(5, "initializing", verbose = verbose)

  trait_labels <- resolve_trait_labels(trait_id, verbose = verbose)

  if (length(trait_labels) == 0L) {
    stop("Could not resolve trait labels for: ", trait_id, call. = FALSE)
  }

  progress_genomic(10, "trait resolved", verbose = verbose)

  ss <- pull_v2_associations(
    efo_id = trait_id,
    trait_labels = trait_labels,
    p_cut = p_cut,
    verbose = verbose
  )

  if (nrow(ss) == 0L) {
    stop("No associations returned from GWAS Catalog REST API v2.", call. = FALSE)
  }

  ss <- ss |>
    dplyr::filter(
      !is.na(.data$p_numeric),
      .data$p_numeric < p_cut,
      !is.na(.data$variant_id),
      .data$variant_id != "",
      grepl("^rs[0-9]+$", .data$variant_id)
    ) |>
    dplyr::distinct(.data$variant_id, .data$study_accession, .keep_all = TRUE)

  if (nrow(ss) == 0L) {
    stop("No rsID associations passed filtering.", call. = FALSE)
  }

  assoc_id <- paste0(ss$variant_id, ":", ss$study_accession)
  ss$association_id <- assoc_id

  progress_genomic(70, "associations ready", verbose = verbose)

  list(
    associations = tibble::tibble(
      association_id = assoc_id,
      pvalue = ss$p_numeric
    ),
    risk_alleles = tibble::tibble(
      association_id = assoc_id,
      variant_id = ss$variant_id
    ),
    cache = ss
  )
}
