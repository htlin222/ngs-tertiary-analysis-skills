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

  # Search each variant
  all_results <- map_dfr(seq_len(nrow(significant_variants)), function(i) {
    var_row <- significant_variants[i, ]
    gene <- var_row$gene
    variant <- var_row$hgvsp

    log_info("Searching Scopus [{i}/{nrow(significant_variants)}]: {gene} {variant}")

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

  # Enrich with open access information for articles with DOIs
  log_info("Enriching with open access status via Unpaywall")

  enriched_results <- unique_results %>%
    mutate(
      oa_lookup = map_df(doi, function(doi_val) {
        if (is.na(doi_val) || is.null(doi_val) || doi_val == "") {
          return(tibble(is_oa = NA, oa_url = NA_character_, oa_status = NA_character_))
        }

        tryCatch(
          {
            # Ensure DOI is in standard format
            doi_clean <- sub("^https?://doi\\.org/", "", doi_val)
            result <- unpaywall_lookup(doi_clean)
            tibble(
              is_oa = result$is_oa,
              oa_url = result$oa_url,
              oa_status = result$oa_status
            )
          },
          error = function(e) {
            log_debug("Unpaywall lookup failed for {doi_val}: {conditionMessage(e)}")
            tibble(is_oa = NA, oa_url = NA_character_, oa_status = NA_character_)
          }
        )
      }, .progress = FALSE),
      is_oa = oa_lookup$is_oa,
      oa_url = oa_lookup$oa_url,
      oa_status = oa_lookup$oa_status,
      .keep = "unused"
    )

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
