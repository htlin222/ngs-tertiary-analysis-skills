# 07-literature/search_scopus.R — Search Scopus for variant-specific literature with open access lookup

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(glue)
  library(logger)
  library(jsonlite)
  library(tidyr)
})

source("R/api_clients.R")

#' Search Scopus for clinically significant variants
#'
#' @param variants Tibble with gene, hgvsp, consequence, classification columns
#' @param config List with sample$tumor_type and literature$scopus settings
#' @param sample_id Sample identifier for output paths
#'
#' @return Tibble of Scopus articles sorted by citation count with open access information
#' @export
#'
#' @details
#' Searches Scopus for clinically significant variants (pathogenic, likely_pathogenic, or high-impact VUS).
#' Results are sorted by citation count (most cited first).
#' For articles with DOIs, enriches results with open access status and URLs via Unpaywall.
#' Saves results to reports/{sample_id}/07-literature/scopus_hits.json
#'
#' @examples
#' \dontrun{
#'   variants <- tibble(
#'     gene = c("BRAF", "TP53"),
#'     hgvsp = c("p.V600E", "p.R273H"),
#'     consequence = c("missense_variant", "missense_variant"),
#'     classification = c("pathogenic", "pathogenic")
#'   )
#'   config <- list(
#'     sample = list(tumor_type = "MEL"),
#'     literature = list(scopus = list(max_results = 20))
#'   )
#'   results <- search_scopus_for_variants(variants, config, "SAMPLE123")
#' }
search_scopus_for_variants <- function(variants, config, sample_id) {
  log_info("Starting Scopus search for sample: {sample_id}")

  # Validate inputs
  if (!is.data.frame(variants) || nrow(variants) == 0) {
    log_warn("No variants provided for Scopus search")
    return(tibble(
      scopus_id = character(), title = character(), authors = character(),
      journal = character(), year = character(), doi = character(),
      cited_by = integer(), is_oa = logical(), oa_url = character(),
      oa_status = character(), query_gene = character(), query_variant = character()
    ))
  }

  required_cols <- c("gene", "hgvsp", "consequence")
  if (!all(required_cols %in% names(variants))) {
    stop(glue("Variants tibble missing required columns: {paste(required_cols, collapse=', ')}"))
  }

  if (is.null(config$sample$tumor_type)) {
    stop("config$sample$tumor_type is required")
  }

  max_results <- config$literature$scopus$max_results %||% 20
  tumor_type <- config$sample$tumor_type

  # Filter for clinically significant variants
  significant_variants <- variants %>%
    filter(
      classification %in% c("pathogenic", "likely_pathogenic") |
        (classification == "VUS" & consequence %in% c(
          "frameshift_variant", "inframe_deletion", "inframe_insertion",
          "disruptive_inframe_deletion", "disruptive_inframe_insertion",
          "start_lost", "stop_gained", "stop_lost", "splice_acceptor_variant",
          "splice_donor_variant"
        ))
    )

  log_info("Found {nrow(significant_variants)} clinically significant variants for Scopus")

  if (nrow(significant_variants) == 0) {
    log_warn("No significant variants for Scopus search")
    return(tibble(
      scopus_id = character(), title = character(), authors = character(),
      journal = character(), year = character(), doi = character(),
      cited_by = integer(), is_oa = logical(), oa_url = character(),
      oa_status = character(), query_gene = character(), query_variant = character()
    ))
  }

  # Deduplicate: one search per unique gene (gene-level is sufficient for literature)
  unique_queries <- significant_variants |>
    distinct(gene, .keep_all = TRUE)

  log_info("Deduplicated {nrow(significant_variants)} variants to {nrow(unique_queries)} unique gene queries")

  # Search each unique gene
  gene_results <- map_dfr(seq_len(nrow(unique_queries)), function(i) {
    var_row <- unique_queries[i, ]
    gene <- var_row$gene
    variant <- var_row$hgvsp

    log_info("Searching Scopus [{i}/{nrow(unique_queries)}]: {gene} {variant}")

    tryCatch(
      {
        results <- search_scopus(
          gene = gene,
          variant = variant,
          tumor_type = tumor_type,
          max_results = max_results
        )

        if (nrow(results) > 0) {
          results %>%
            mutate(
              query_gene = gene,
              query_variant = variant,
              .before = 1
            )
        } else {
          tibble(
            query_gene = character(), query_variant = character(),
            scopus_id = character(), title = character(), authors = character(),
            journal = character(), year = character(), doi = character(),
            cited_by = integer()
          )
        }
      },
      error = function(e) {
        log_error("Scopus search failed for {gene} {variant}: {conditionMessage(e)}")
        tibble(
          query_gene = character(), query_variant = character(),
          scopus_id = character(), title = character(), authors = character(),
          journal = character(), year = character(), doi = character(),
          cited_by = integer()
        )
      }
    )
  })

  # Join results back to all original variants by gene
  variant_gene_map <- significant_variants |>
    select(gene, hgvsp) |>
    filter(!gene %in% unique_queries$gene | hgvsp != unique_queries$hgvsp[match(gene, unique_queries$gene)])

  all_results <- gene_results
  if (nrow(variant_gene_map) > 0 && nrow(gene_results) > 0) {
    # Expand results to cover all variants that share a gene
    extra_rows <- variant_gene_map |>
      inner_join(gene_results |> select(-query_variant), by = c("gene" = "query_gene"), relationship = "many-to-many") |>
      rename(query_gene = gene, query_variant = hgvsp)
    all_results <- bind_rows(gene_results, extra_rows)
  }

  if (nrow(all_results) == 0) {
    log_warn("No Scopus results found")
    return(tibble(
      scopus_id = character(), title = character(), authors = character(),
      journal = character(), year = character(), doi = character(),
      cited_by = integer(), is_oa = logical(), oa_url = character(),
      oa_status = character(), query_gene = character(), query_variant = character()
    ))
  }

  # Deduplicate by scopus_id
  unique_results <- all_results %>%
    distinct(scopus_id, .keep_all = TRUE) %>%
    arrange(desc(cited_by), desc(year))

  log_info("Retrieved {nrow(unique_results)} unique Scopus articles")

  # Enrich with open access information for articles with DOIs (parallel requests)
  log_info("Enriching with open access status via Unpaywall (parallel)")

  dois <- unique(unique_results$doi[!is.na(unique_results$doi) & unique_results$doi != ""])

  if (length(dois) > 0) {
    unpaywall_email <- tryCatch(get_api_key("UNPAYWALL_EMAIL"), error = function(e) "")

    # Build all Unpaywall requests upfront
    reqs <- lapply(dois, function(doi_val) {
      doi_clean <- sub("^https?://doi\\.org/", "", doi_val)
      httr2::request(glue("https://api.unpaywall.org/v2/{doi_clean}")) |>
        httr2::req_url_query(email = unpaywall_email) |>
        httr2::req_retry(max_tries = 2) |>
        httr2::req_error(is_error = ~ FALSE)
    })

    resps <- httr2::req_perform_parallel(
      reqs,
      pool = curl::new_pool(total_con = 10),
      on_error = "continue"
    )

    # Parse responses
    oa_results <- map2_dfr(dois, resps, function(doi_val, resp) {
      tryCatch(
        {
          if (httr2::resp_status(resp) != 200) {
            return(tibble(doi = doi_val, is_oa = FALSE, oa_url = NA_character_, oa_status = "unknown"))
          }
          body <- httr2::resp_body_json(resp)
          best_oa <- body$best_oa_location
          tibble(
            doi = doi_val,
            is_oa = body$is_oa %||% FALSE,
            oa_url = best_oa$url_for_pdf %||% best_oa$url %||% NA_character_,
            oa_status = body$oa_status %||% "unknown"
          )
        },
        error = function(e) {
          log_debug("Unpaywall parse failed for {doi_val}: {conditionMessage(e)}")
          tibble(doi = doi_val, is_oa = NA, oa_url = NA_character_, oa_status = NA_character_)
        }
      )
    })

    log_info("Unpaywall lookup complete for {length(dois)} DOIs")

    enriched_results <- unique_results |>
      left_join(oa_results, by = "doi")
  } else {
    enriched_results <- unique_results |>
      mutate(is_oa = NA, oa_url = NA_character_, oa_status = NA_character_)
  }

  # Create output directory if needed
  output_dir <- glue("reports/{sample_id}/07-literature")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Save results
  output_path <- file.path(output_dir, "scopus_hits.json")
  tryCatch(
    {
      write_json(
        enriched_results %>%
          select(
            scopus_id, title, authors, journal, year, doi, cited_by,
            is_oa, oa_url, oa_status, query_gene, query_variant
          ),
        output_path,
        pretty = TRUE
      )
      log_info("Saved Scopus results to {output_path}")
    },
    error = function(e) {
      log_error("Failed to save Scopus results: {conditionMessage(e)}")
    }
  )

  enriched_results %>%
    select(
      scopus_id, title, authors, journal, year, doi, cited_by,
      is_oa, oa_url, oa_status, query_gene, query_variant
    )
}
