# 06-clinical-annotation/query_oncokb.R — OncoKB API querying for SNVs, CNVs, and fusions

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
  library(glue)
  library(logger)
  library(purrr)
})

# Source required functions
source(here::here("R/api_clients.R"))
source(here::here("R/utils.R"))

#' Check and load cached OncoKB results
#'
#' Attempts to load a previously cached result from disk to avoid redundant API calls.
#'
#' @param cache_key Unique identifier for the alteration (e.g., "BRAF_V600E_NSCLC")
#' @param cache_dir Directory where cache files are stored
#'
#' @return List of cached result, or NULL if not found or cache disabled
#'
#' @keywords internal
load_oncokb_cache <- function(cache_key, cache_dir) {
  if (is.null(cache_dir) || !dir.exists(cache_dir)) {
    return(NULL)
  }

  cache_file <- file.path(cache_dir, glue("{cache_key}.json"))
  if (!file.exists(cache_file)) {
    return(NULL)
  }

  tryCatch(
    {
      cached <- fromJSON(cache_file)
      log_debug("Loaded cached OncoKB result: {cache_key}")
      return(cached)
    },
    error = function(e) {
      log_warn("Failed to load cache for {cache_key}: {e$message}")
      return(NULL)
    }
  )
}

#' Save OncoKB result to cache
#'
#' Writes a result to disk for future use, avoiding redundant API calls.
#'
#' @param result List of OncoKB annotation result
#' @param cache_key Unique identifier for the alteration
#' @param cache_dir Directory where cache files are stored
#'
#' @return TRUE if saved successfully, FALSE otherwise
#'
#' @keywords internal
save_oncokb_cache <- function(result, cache_key, cache_dir) {
  if (is.null(cache_dir)) {
    return(FALSE)
  }

  tryCatch(
    {
      dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
      cache_file <- file.path(cache_dir, glue("{cache_key}.json"))
      write_json(result, cache_file, pretty = TRUE)
      log_debug("Cached OncoKB result: {cache_key}")
      return(TRUE)
    },
    error = function(e) {
      log_warn("Failed to save cache for {cache_key}: {e$message}")
      return(FALSE)
    }
  )
}

#' Annotate mutations via OncoKB with caching
#'
#' Queries OncoKB API for each mutation, with optional disk caching to avoid
#' redundant requests. Gracefully handles errors for individual mutations.
#'
#' @param merged_annotations Tibble with variant annotations (from stage 02)
#'   Must contain: gene, hgvsp columns
#' @param tumor_type OncoKB tumor type (e.g., "NSCLC", "MEL")
#' @param cache_dir Optional directory for caching results
#'
#' @return List of mutation annotation results (one per variant)
#'
#' @keywords internal
query_mutations <- function(merged_annotations, tumor_type, cache_dir = NULL) {
  if (nrow(merged_annotations) == 0) {
    log_info("No mutations to annotate")
    return(list())
  }

  log_info("Querying OncoKB for {nrow(merged_annotations)} mutations")

  results <- map(seq_len(nrow(merged_annotations)), function(i) {
    row <- merged_annotations[i, ]

    gene <- row$gene
    hgvsp <- row$hgvsp

    # Skip missing values
    if (is.na(gene) || is.na(hgvsp)) {
      log_debug("Skipping row {i}: missing gene or hgvsp")
      return(NA)
    }

    # Extract protein change (strip "p." prefix if present)
    protein_change <- str_remove(as.character(hgvsp), "^p\\.")

    # Create cache key
    cache_key <- glue("{gene}_{protein_change}_{tumor_type}")

    # Check cache first
    cached_result <- load_oncokb_cache(cache_key, cache_dir)
    if (!is.null(cached_result)) {
      return(cached_result)
    }

    # Query OncoKB API with error handling
    result <- tryCatch(
      {
        oncokb_annotate_mutation(gene, protein_change, tumor_type)
      },
      error = function(e) {
        log_warn("OncoKB annotation failed for {gene} {protein_change}: {e$message}")
        return(NA)
      }
    )

    # Cache successful result
    if (!is.na(result)) {
      save_oncokb_cache(result, cache_key, cache_dir)
    }

    return(result)
  })

  # Filter out NA results but keep track of failures
  n_failed <- sum(is.na(results))
  if (n_failed > 0) {
    log_warn("Failed to annotate {n_failed} mutations")
  }

  results <- results[!is.na(results)]
  log_info("Successfully annotated {length(results)} mutations")

  results
}

#' Annotate copy number alterations via OncoKB with caching
#'
#' Queries OncoKB for CNAs (amplifications and deletions), with optional caching.
#' Converts local type (AMP/DEL) to OncoKB format (AMPLIFICATION/DELETION).
#'
#' @param parsed_cnv Tibble with parsed CNV data (from stage 03)
#'   Must contain: gene, type columns (AMP/DEL)
#' @param tumor_type OncoKB tumor type
#' @param cache_dir Optional cache directory
#'
#' @return List of CNA annotation results
#'
#' @keywords internal
query_cnas <- function(parsed_cnv, tumor_type, cache_dir = NULL) {
  if (nrow(parsed_cnv) == 0) {
    log_info("No CNVs to annotate")
    return(list())
  }

  # Filter to significant events only (AMP or DEL, skip GAIN/LOSS)
  significant_cnv <- parsed_cnv %>%
    filter(type %in% c("AMP", "DEL"))

  if (nrow(significant_cnv) == 0) {
    log_info("No significant CNVs (AMP/DEL) found to annotate")
    return(list())
  }

  log_info("Querying OncoKB for {nrow(significant_cnv)} CNVs")

  results <- map(seq_len(nrow(significant_cnv)), function(i) {
    row <- significant_cnv[i, ]

    gene <- row$gene
    cna_type <- row$type

    # Skip missing values
    if (is.na(gene) || is.na(cna_type)) {
      log_debug("Skipping row {i}: missing gene or type")
      return(NA)
    }

    # Map local type to OncoKB format
    oncokb_type <- if (cna_type == "AMP") "AMPLIFICATION" else "DELETION"

    # Create cache key
    cache_key <- glue("{gene}_{oncokb_type}_{tumor_type}")

    # Check cache first
    cached_result <- load_oncokb_cache(cache_key, cache_dir)
    if (!is.null(cached_result)) {
      return(cached_result)
    }

    # Query OncoKB API with error handling
    result <- tryCatch(
      {
        oncokb_annotate_cna(gene, oncokb_type, tumor_type)
      },
      error = function(e) {
        log_warn("OncoKB CNA annotation failed for {gene} {oncokb_type}: {e$message}")
        return(NA)
      }
    )

    # Cache successful result
    if (!is.na(result)) {
      save_oncokb_cache(result, cache_key, cache_dir)
    }

    return(result)
  })

  # Filter out NA results
  n_failed <- sum(is.na(results))
  if (n_failed > 0) {
    log_warn("Failed to annotate {n_failed} CNVs")
  }

  results <- results[!is.na(results)]
  log_info("Successfully annotated {length(results)} CNVs")

  results
}

#' Annotate fusions via OncoKB with caching
#'
#' Queries OncoKB for structural variants (fusions), with optional caching.
#'
#' @param parsed_fusions Tibble with parsed fusion data (from stage 04)
#'   Must contain: gene_a, gene_b columns
#' @param tumor_type OncoKB tumor type
#' @param cache_dir Optional cache directory
#'
#' @return List of fusion annotation results
#'
#' @keywords internal
query_fusions <- function(parsed_fusions, tumor_type, cache_dir = NULL) {
  if (nrow(parsed_fusions) == 0) {
    log_info("No fusions to annotate")
    return(list())
  }

  log_info("Querying OncoKB for {nrow(parsed_fusions)} fusions")

  results <- map(seq_len(nrow(parsed_fusions)), function(i) {
    row <- parsed_fusions[i, ]

    gene_a <- row$gene_a
    gene_b <- row$gene_b

    # Skip missing values
    if (is.na(gene_a) || is.na(gene_b)) {
      log_debug("Skipping row {i}: missing gene_a or gene_b")
      return(NA)
    }

    # Create cache key
    cache_key <- glue("{gene_a}_{gene_b}_{tumor_type}")

    # Check cache first
    cached_result <- load_oncokb_cache(cache_key, cache_dir)
    if (!is.null(cached_result)) {
      return(cached_result)
    }

    # Query OncoKB API with error handling
    result <- tryCatch(
      {
        oncokb_annotate_fusion(gene_a, gene_b, tumor_type)
      },
      error = function(e) {
        log_warn("OncoKB fusion annotation failed for {gene_a}::{gene_b}: {e$message}")
        return(NA)
      }
    )

    # Cache successful result
    if (!is.na(result)) {
      save_oncokb_cache(result, cache_key, cache_dir)
    }

    return(result)
  })

  # Filter out NA results
  n_failed <- sum(is.na(results))
  if (n_failed > 0) {
    log_warn("Failed to annotate {n_failed} fusions")
  }

  results <- results[!is.na(results)]
  log_info("Successfully annotated {length(results)} fusions")

  results
}

#' Query OncoKB for all variant types in a sample
#'
#' Main entry point for OncoKB annotation. Annotates SNVs/indels, CNVs, and
#' structural variants using the OncoKB API. Results can be cached to disk to
#' avoid redundant API requests.
#'
#' @param variants Tibble of variant annotations (from stage 02 merge_annotations)
#' @param cnv Tibble of parsed CNV calls (from stage 03 parse_cnv)
#' @param fusions Tibble of parsed fusion calls (from stage 04 parse_fusions)
#' @param config List of pipeline configuration
#'   - config$sample$tumor_type: OncoKB tumor type (required, e.g., "NSCLC")
#'   - config$clinical_annotation$oncokb$cache_dir: Optional cache directory path
#' @param sample_id Sample identifier for logging and output organization
#'
#' @return List with three elements:
#'   - mutations: List of OncoKB mutation annotation results
#'   - cnas: List of OncoKB CNA annotation results
#'   - fusions: List of OncoKB fusion annotation results
#'
#'   Each result contains: gene, alteration, tumor_type, oncogenic, mutation_effect
#'   (for mutations), highest_sensitive_level, highest_resistance_level, treatments
#'
#' @details
#' The function performs the following steps:
#' 1. Validate configuration (tumor_type required)
#' 2. Query mutations from merged variant annotations
#' 3. Query CNAs from parsed CNV data (AMP/DEL only)
#' 4. Query fusions from parsed structural variant data
#' 5. Cache results if configured
#' 6. Save all results to JSON file
#' 7. Return combined results list
#'
#' API throttling and retry logic are handled by the underlying API client
#' functions. Individual annotation failures are logged but do not stop
#' processing of remaining variants.
#'
#' @export
query_oncokb <- function(variants, cnv, fusions, config, sample_id) {
  log_info("Starting OncoKB annotation for sample: {sample_id}")

  # Validate required configuration
  tumor_type <- config$sample$tumor_type %||% NA_character_
  if (is.na(tumor_type)) {
    log_error("tumor_type not configured in config$sample$tumor_type")
    stop("tumor_type is required in configuration")
  }

  log_info("Using tumor type: {tumor_type}")

  # Get cache directory if configured
  cache_dir <- config$clinical_annotation$oncokb$cache_dir %||% NULL

  # Query each variant type
  mutations <- query_mutations(variants, tumor_type, cache_dir)
  cnas <- query_cnas(cnv, tumor_type, cache_dir)
  fusions_results <- query_fusions(fusions, tumor_type, cache_dir)

  # Combine results
  results <- list(
    mutations = mutations,
    cnas = cnas,
    fusions = fusions_results
  )

  # Save to JSON
  output_dir <- stage_output_dir(sample_id, "06-clinical-annotation")
  output_file <- file.path(output_dir, "oncokb_results.json")

  tryCatch(
    {
      write_json(results, output_file, pretty = TRUE, auto_unbox = TRUE)
      log_info("OncoKB results saved to {output_file}")
    },
    error = function(e) {
      log_warn("Failed to save OncoKB results to JSON: {e$message}")
    }
  )

  # Summary logging
  log_info(
    "OncoKB annotation complete: {length(mutations)} mutations, ",
    "{length(cnas)} CNVs, {length(fusions_results)} fusions"
  )

  return(results)
}
