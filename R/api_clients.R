# R/api_clients.R — API client wrappers for OncoKB, PubMed, Scopus, Unpaywall

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(logger)
  library(glue)
  library(dplyr)
  library(purrr)
})

# ══════════════════════════════════════════════════════════════════════════════
# OncoKB API
# Docs: https://api.oncokb.org/
# ══════════════════════════════════════════════════════════════════════════════

#' Create an authenticated OncoKB API request
#' @param endpoint API endpoint path (e.g., "/annotate/mutations/byProteinChange")
#' @return httr2 request object
oncokb_request <- function(endpoint) {
  token <- get_api_key("ONCOKB_API_KEY")
  request("https://www.oncokb.org/api/v1") |>
    req_url_path_append(endpoint) |>
    req_headers(
      Authorization = glue("Bearer {token}"),
      Accept = "application/json"
    ) |>
    req_retry(max_tries = 3, backoff = ~ 2^.x) |>
    req_throttle(rate = 10 / 60)  # 10 requests per minute
}

#' Annotate a somatic mutation via OncoKB
#' @param hugo_symbol Gene symbol (e.g., "BRAF")
#' @param protein_change Protein change (e.g., "V600E")
#' @param tumor_type OncoKB tumor type (e.g., "NSCLC", "MEL")
#' @return List with oncogenic, mutationEffect, treatments, levels
oncokb_annotate_mutation <- function(hugo_symbol, protein_change, tumor_type) {
  log_info("OncoKB: annotating {hugo_symbol} {protein_change} in {tumor_type}")

  resp <- oncokb_request("/annotate/mutations/byProteinChange") |>
    req_url_query(
      hugoSymbol = hugo_symbol,
      alteration = protein_change,
      tumorType = tumor_type
    ) |>
    req_perform()

  result <- resp_body_json(resp)

  list(
    gene = hugo_symbol,
    alteration = protein_change,
    tumor_type = tumor_type,
    oncogenic = result$oncogenic %||% "Unknown",
    mutation_effect = result$mutationEffect$knownEffect %||% "Unknown",
    highest_sensitive_level = result$highestSensitiveLevel %||% NA_character_,
    highest_resistance_level = result$highestResistanceLevel %||% NA_character_,
    treatments = parse_oncokb_treatments(result$treatments %||% list()),
    raw = result
  )
}

#' Annotate a copy number alteration via OncoKB
#' @param hugo_symbol Gene symbol
#' @param cna_type "AMPLIFICATION" or "DELETION"
#' @param tumor_type OncoKB tumor type
#' @return Annotation result list
oncokb_annotate_cna <- function(hugo_symbol, cna_type, tumor_type) {
  log_info("OncoKB: annotating {hugo_symbol} {cna_type} in {tumor_type}")

  valid_types <- c("AMPLIFICATION", "DELETION", "GAIN", "LOSS")
  cna_type <- toupper(cna_type)
  if (!cna_type %in% valid_types) {
    stop(glue("Invalid CNA type: {cna_type}. Must be one of: {paste(valid_types, collapse=', ')}"))
  }

  resp <- oncokb_request("/annotate/copyNumberAlterations") |>
    req_url_query(
      hugoSymbol = hugo_symbol,
      copyNameAlterationType = cna_type,
      tumorType = tumor_type
    ) |>
    req_perform()

  result <- resp_body_json(resp)

  list(
    gene = hugo_symbol,
    alteration = cna_type,
    tumor_type = tumor_type,
    oncogenic = result$oncogenic %||% "Unknown",
    highest_sensitive_level = result$highestSensitiveLevel %||% NA_character_,
    treatments = parse_oncokb_treatments(result$treatments %||% list()),
    raw = result
  )
}

#' Annotate a structural variant / fusion via OncoKB
#' @param gene_a 5' gene symbol
#' @param gene_b 3' gene symbol
#' @param tumor_type OncoKB tumor type
#' @param sv_type Structural variant type (default: "FUSION")
#' @return Annotation result list
oncokb_annotate_fusion <- function(gene_a, gene_b, tumor_type,
                                   sv_type = "FUSION") {
  log_info("OncoKB: annotating {gene_a}::{gene_b} {sv_type} in {tumor_type}")

  resp <- oncokb_request("/annotate/structuralVariants") |>
    req_url_query(
      hugoSymbolA = gene_a,
      hugoSymbolB = gene_b,
      structuralVariantType = sv_type,
      tumorType = tumor_type
    ) |>
    req_perform()

  result <- resp_body_json(resp)

  list(
    gene_a = gene_a,
    gene_b = gene_b,
    alteration = glue("{gene_a}::{gene_b}"),
    tumor_type = tumor_type,
    oncogenic = result$oncogenic %||% "Unknown",
    highest_sensitive_level = result$highestSensitiveLevel %||% NA_character_,
    treatments = parse_oncokb_treatments(result$treatments %||% list()),
    raw = result
  )
}

#' Parse OncoKB treatments list into a tidy tibble
#' @param treatments List of treatment objects from OncoKB response
#' @return Tibble with drug, level, description columns
parse_oncokb_treatments <- function(treatments) {
  if (length(treatments) == 0) {
    return(tibble(
      drugs = character(),
      level = character(),
      description = character(),
      fda_approved = logical()
    ))
  }

  map_dfr(treatments, function(tx) {
    drugs <- map_chr(tx$drugs %||% list(), ~ .x$drugName %||% "Unknown")
    tibble(
      drugs = paste(drugs, collapse = " + "),
      level = tx$level %||% NA_character_,
      description = tx$description %||% "",
      fda_approved = grepl("^LEVEL_[12]$", tx$level %||% "")
    )
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# PubMed (NCBI E-utilities via rentrez)
# ══════════════════════════════════════════════════════════════════════════════

#' Search PubMed for variant-specific literature
#' @param gene Gene symbol
#' @param variant Variant description (e.g., "V600E")
#' @param tumor_type Cancer type for context
#' @param max_results Maximum results to return
#' @return Tibble of PubMed articles
search_pubmed <- function(gene, variant = NULL, tumor_type = NULL,
                          max_results = 10) {
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    stop("Package 'rentrez' required. Install with: install.packages('rentrez')")
  }

  api_key <- tryCatch(get_api_key("PUBMED_API_KEY"), error = function(e) NULL)

  # Build search query
  terms <- gene
  if (!is.null(variant)) terms <- glue("{terms} AND {variant}")
  if (!is.null(tumor_type)) terms <- glue("{terms} AND {tumor_type}")
  terms <- glue("({terms}) AND (cancer OR oncology OR tumor)")

  log_info("PubMed search: {terms}")

  search_result <- rentrez::entrez_search(
    db = "pubmed",
    term = terms,
    retmax = max_results,
    api_key = api_key
  )

  if (length(search_result$ids) == 0) {
    log_info("No PubMed results for: {terms}")
    return(tibble(
      pmid = character(), title = character(), authors = character(),
      journal = character(), year = character(), abstract = character()
    ))
  }

  # Fetch article details
  records <- rentrez::entrez_fetch(
    db = "pubmed",
    id = search_result$ids,
    rettype = "xml",
    api_key = api_key
  )

  parse_pubmed_xml(records)
}

#' Parse PubMed XML response into tidy tibble
#' @param xml_text Raw XML text from entrez_fetch
#' @return Tibble of article metadata
parse_pubmed_xml <- function(xml_text) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Package 'xml2' required for PubMed parsing")
  }

  doc <- xml2::read_xml(xml_text)
  articles <- xml2::xml_find_all(doc, ".//PubmedArticle")

  map_dfr(articles, function(art) {
    pmid <- xml2::xml_text(xml2::xml_find_first(art, ".//PMID"))
    title <- xml2::xml_text(xml2::xml_find_first(art, ".//ArticleTitle"))
    journal <- xml2::xml_text(xml2::xml_find_first(art, ".//Journal/Title"))
    year <- xml2::xml_text(xml2::xml_find_first(art, ".//PubDate/Year"))

    author_nodes <- xml2::xml_find_all(art, ".//Author")
    authors <- paste(map_chr(author_nodes, function(a) {
      last <- xml2::xml_text(xml2::xml_find_first(a, ".//LastName"))
      init <- xml2::xml_text(xml2::xml_find_first(a, ".//Initials"))
      if (is.na(last)) return(NA_character_)
      paste0(last, " ", init)
    }), collapse = ", ")

    abstract_parts <- xml2::xml_find_all(art, ".//AbstractText")
    abstract <- paste(xml2::xml_text(abstract_parts), collapse = " ")

    tibble(
      pmid = pmid %||% NA_character_,
      title = title %||% NA_character_,
      authors = authors,
      journal = journal %||% NA_character_,
      year = year %||% NA_character_,
      abstract = abstract
    )
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# Scopus (Elsevier)
# ══════════════════════════════════════════════════════════════════════════════

#' Search Scopus for variant-specific literature
#' @param gene Gene symbol
#' @param variant Variant description
#' @param tumor_type Cancer type
#' @param max_results Maximum results
#' @return Tibble of Scopus articles
search_scopus <- function(gene, variant = NULL, tumor_type = NULL,
                          max_results = 10) {
  api_key <- get_api_key("SCOPUS_API_KEY")

  # Build search query
  query <- glue('TITLE-ABS-KEY("{gene}")')
  if (!is.null(variant)) query <- glue('{query} AND TITLE-ABS-KEY("{variant}")')
  if (!is.null(tumor_type)) query <- glue('{query} AND TITLE-ABS-KEY("{tumor_type}")')
  query <- glue('{query} AND TITLE-ABS-KEY(cancer OR oncology)')

  log_info("Scopus search: {query}")

  resp <- request("https://api.elsevier.com/content/search/scopus") |>
    req_url_query(
      query = query,
      count = max_results,
      sort = "-citedby-count"
    ) |>
    req_headers(
      `X-ELS-APIKey` = api_key,
      Accept = "application/json"
    ) |>
    req_retry(max_tries = 3) |>
    req_throttle(rate = 6 / 60) |>
    req_perform()

  result <- resp_body_json(resp)
  entries <- result$`search-results`$entry %||% list()

  if (length(entries) == 0) {
    log_info("No Scopus results for: {query}")
    return(tibble(
      scopus_id = character(), title = character(), authors = character(),
      journal = character(), year = character(), doi = character(),
      cited_by = integer()
    ))
  }

  map_dfr(entries, function(e) {
    tibble(
      scopus_id = e$`dc:identifier` %||% NA_character_,
      title = e$`dc:title` %||% NA_character_,
      authors = e$`dc:creator` %||% NA_character_,
      journal = e$`prism:publicationName` %||% NA_character_,
      year = substr(e$`prism:coverDate` %||% "", 1, 4),
      doi = e$`prism:doi` %||% NA_character_,
      cited_by = as.integer(e$`citedby-count` %||% 0L)
    )
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# Unpaywall (open access lookup)
# ══════════════════════════════════════════════════════════════════════════════

#' Look up open access status and PDF URL via Unpaywall
#' @param doi DOI string
#' @return List with is_oa, oa_url, oa_status
unpaywall_lookup <- function(doi) {
  email <- get_api_key("UNPAYWALL_EMAIL")

  resp <- request(glue("https://api.unpaywall.org/v2/{doi}")) |>
    req_url_query(email = email) |>
    req_retry(max_tries = 2) |>
    req_throttle(rate = 100 / 60) |>
    req_error(is_error = ~ FALSE) |>
    req_perform()

  if (resp_status(resp) != 200) {
    return(list(is_oa = FALSE, oa_url = NA_character_, oa_status = "unknown"))
  }

  result <- resp_body_json(resp)
  best_oa <- result$best_oa_location

  list(
    is_oa = result$is_oa %||% FALSE,
    oa_url = best_oa$url_for_pdf %||% best_oa$url %||% NA_character_,
    oa_status = result$oa_status %||% "unknown"
  )
}
