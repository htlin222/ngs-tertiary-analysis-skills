# 07-literature/search_pubmed.R — Search PubMed for variant-specific literature

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(glue)
  library(logger)
  library(jsonlite)
  library(tidyr)
})

source("R/api_clients.R")

#' Search PubMed for clinically significant variants
#'
#' @param variants Tibble with gene, hgvsp, consequence, classification columns
#' @param config List with sample$tumor_type and literature$pubmed settings
#' @param sample_id Sample identifier for output paths
#'
#' @return Tibble of unique PubMed articles (pmid, title, authors, journal, year, abstract)
#' @export
#'
#' @details
#' Searches for clinically significant variants (pathogenic, likely_pathogenic, or high-impact VUS).
#' Constructs queries as: "{gene} {protein_change} cancer"
#' Deduplicates results by PMID and saves to reports/{sample_id}/07-literature/pubmed_hits.json
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
#'     literature = list(pubmed = list(max_results = 20))
#'   )
#'   results <- search_pubmed_for_variants(variants, config, "SAMPLE123")
#' }
search_pubmed_for_variants <- function(variants, config, sample_id) {
  log_info("Starting PubMed search for sample: {sample_id}")

  # Validate inputs
  if (!is.data.frame(variants) || nrow(variants) == 0) {
    log_warn("No variants provided for PubMed search")
    return(tibble(
      pmid = character(), title = character(), authors = character(),
      journal = character(), year = character(), abstract = character(),
      query_gene = character(), query_variant = character()
    ))
  }

  required_cols <- c("gene", "hgvsp", "consequence")
  if (!all(required_cols %in% names(variants))) {
    stop(glue("Variants tibble missing required columns: {paste(required_cols, collapse=', ')}"))
  }

  if (is.null(config$sample$tumor_type)) {
    stop("config$sample$tumor_type is required")
  }

  max_results <- config$literature$pubmed$max_results %||% 20
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

  log_info("Found {nrow(significant_variants)} clinically significant variants")

  if (nrow(significant_variants) == 0) {
    log_warn("No significant variants for PubMed search")
    return(tibble(
      pmid = character(), title = character(), authors = character(),
      journal = character(), year = character(), abstract = character(),
      query_gene = character(), query_variant = character()
    ))
  }

  # Search each variant
  all_results <- map_dfr(seq_len(nrow(significant_variants)), function(i) {
    var_row <- significant_variants[i, ]
    gene <- var_row$gene
    variant <- var_row$hgvsp

    log_info("Searching PubMed [{i}/{nrow(significant_variants)}]: {gene} {variant}")

    tryCatch(
      {
        results <- search_pubmed(
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
            pmid = character(), title = character(), authors = character(),
            journal = character(), year = character(), abstract = character()
          )
        }
      },
      error = function(e) {
        log_error("PubMed search failed for {gene} {variant}: {conditionMessage(e)}")
        tibble(
          query_gene = character(), query_variant = character(),
          pmid = character(), title = character(), authors = character(),
          journal = character(), year = character(), abstract = character()
        )
      }
    )
  })

  # Add biomarker-specific searches if relevant
  biomarker_searches <- list()

  # TMB search if high mutational burden
  if (nrow(variants) > 500) {
    log_info("High mutation burden detected, searching for TMB-related literature")
    biomarker_searches[["TMB"]] <- tryCatch(
      {
        results <- search_pubmed(
          gene = "tumor mutational burden",
          variant = NULL,
          tumor_type = tumor_type,
          max_results = 10
        )
        if (nrow(results) > 0) {
          results %>%
            mutate(
              query_gene = "TMB",
              query_variant = "biomarker",
              .before = 1
            )
        } else {
          NULL
        }
      },
      error = function(e) {
        log_error("TMB search failed: {conditionMessage(e)}")
        NULL
      }
    )
  }

  # MSI search if high MSI status would be relevant
  if (!is.null(config$sample$msi_status) && config$sample$msi_status == "high") {
    log_info("High MSI status detected, searching for MSI-related literature")
    biomarker_searches[["MSI"]] <- tryCatch(
      {
        results <- search_pubmed(
          gene = "microsatellite instability",
          variant = NULL,
          tumor_type = tumor_type,
          max_results = 10
        )
        if (nrow(results) > 0) {
          results %>%
            mutate(
              query_gene = "MSI",
              query_variant = "biomarker",
              .before = 1
            )
        } else {
          NULL
        }
      },
      error = function(e) {
        log_error("MSI search failed: {conditionMessage(e)}")
        NULL
      }
    )
  }

  # Combine biomarker results
  biomarker_results <- compact(biomarker_searches) %>%
    bind_rows()

  if (nrow(biomarker_results) > 0) {
    all_results <- bind_rows(all_results, biomarker_results)
  }

  # Deduplicate by PMID, keeping first occurrence
  unique_results <- all_results %>%
    distinct(pmid, .keep_all = TRUE) %>%
    arrange(desc(year), query_gene)

  log_info("Retrieved {nrow(unique_results)} unique PubMed articles")

  # Create output directory if needed
  output_dir <- glue("reports/{sample_id}/07-literature")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Save results
  output_path <- file.path(output_dir, "pubmed_hits.json")
  tryCatch(
    {
      write_json(
        unique_results %>%
          select(pmid, title, authors, journal, year, abstract, query_gene, query_variant),
        output_path,
        pretty = TRUE
      )
      log_info("Saved PubMed results to {output_path}")
    },
    error = function(e) {
      log_error("Failed to save PubMed results: {conditionMessage(e)}")
    }
  )

  unique_results
}
