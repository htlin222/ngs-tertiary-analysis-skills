# 06-clinical-annotation/query_civic.R — CiVIC API querying for community variant evidence

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
  library(glue)
  library(logger)
  library(purrr)
})

# Source required functions
source(here::here("R/civic_client.R"))
source(here::here("R/api_clients.R"))
source(here::here("R/utils.R"))

#' Query CiVIC for community variant evidence
#'
#' Builds a combined variant table from SNVs, CNVs, and fusions, then queries
#' the CiVIC GraphQL API for evidence items, gene summaries, and AMP/ASCO/CAP
#' assertions. Results are cached to disk.
#'
#' @param variants Tibble of variant annotations (from stage 02 merge_annotations).
#'   Must contain: gene, hgvsp columns.
#' @param cnv Tibble of parsed CNV calls (from stage 03 parse_cnv).
#'   Must contain: gene, type columns.
#' @param fusions Tibble of parsed fusion calls (from stage 04 parse_fusions).
#'   Must contain: gene_a, gene_b columns.
#' @param config List of pipeline configuration.
#'   - config$clinical_annotation$civic$enabled: Whether CiVIC is enabled.
#'   - config$clinical_annotation$civic$cache_ttl_hours: Cache TTL (default 168).
#' @param sample_id Sample identifier for logging and output organization.
#'
#' @return List with:
#'   - variant_evidence: Tibble of CiVIC evidence items
#'   - gene_summaries: Named list of gene summary info
#'   - assertions: Tibble of AMP/ASCO/CAP assertions
#'
#' @export
query_civic <- function(variants, cnv, fusions, config, sample_id) {
  log_info("Starting CiVIC annotation for sample: {sample_id}")

  # Empty result structure
  empty_result <- list(
    variant_evidence = tibble(
      evidence_id = integer(), evidence_type = character(),
      evidence_level = character(), evidence_direction = character(),
      significance = character(), disease = character(),
      therapies = character(), source_citation = character(),
      source_url = character(), description = character(),
      query_gene = character(), query_variant = character()
    ),
    gene_summaries = list(),
    assertions = tibble(
      assertion_id = integer(), amp_level = character(),
      assertion_type = character(), assertion_direction = character(),
      significance = character(), disease = character(),
      therapies = character(), variant_name = character(),
      gene = character(), description = character(),
      nccn_guideline = character(), regulatory_approval = logical(),
      fda_companion_test = logical()
    )
  )

  # Check if CiVIC is enabled
  civic_enabled <- config$clinical_annotation$civic$enabled %||% FALSE
  if (!isTRUE(civic_enabled)) {
    log_info("CiVIC annotation disabled in config — skipping")
    return(empty_result)
  }

  # Build combined variant table
  combined <- tibble(gene = character(), variant = character())

  # SNVs: gene + protein change (strip p. prefix)
  if (nrow(variants) > 0 && "hgvsp" %in% names(variants)) {
    snv_variants <- variants |>
      filter(!is.na(gene), !is.na(hgvsp)) |>
      transmute(
        gene = gene,
        variant = str_remove(as.character(hgvsp), "^p\\.")
      )
    combined <- bind_rows(combined, snv_variants)
  }

  # CNVs: gene + type
  if (nrow(cnv) > 0 && "type" %in% names(cnv)) {
    cnv_variants <- cnv |>
      filter(!is.na(gene), !is.na(type)) |>
      transmute(
        gene = gene,
        variant = type
      )
    combined <- bind_rows(combined, cnv_variants)
  }

  # Fusions: both gene partners
  if (nrow(fusions) > 0 && all(c("gene_a", "gene_b") %in% names(fusions))) {
    fusion_a <- fusions |>
      filter(!is.na(gene_a), !is.na(gene_b)) |>
      transmute(
        gene = gene_a,
        variant = paste0(gene_a, "::", gene_b)
      )
    fusion_b <- fusions |>
      filter(!is.na(gene_a), !is.na(gene_b)) |>
      transmute(
        gene = gene_b,
        variant = paste0(gene_a, "::", gene_b)
      )
    combined <- bind_rows(combined, fusion_a, fusion_b)
  }

  # Deduplicate and filter
  combined <- combined |>
    filter(!is.na(gene), !is.na(variant)) |>
    distinct()

  if (nrow(combined) == 0) {
    log_info("No variants to query CiVIC — returning empty results")
    return(empty_result)
  }

  log_info("Querying CiVIC for {nrow(combined)} unique variants")

  # Use cached API call
  cache_ttl <- config$clinical_annotation$civic$cache_ttl_hours %||% 168

  results <- cached_api_call(
    cache_key = glue("civic_batch_{sample_id}"),
    fn = function() civic_annotate_variants(combined, sample_id),
    sample_id = sample_id,
    ttl_hours = cache_ttl
  )

  # Save results as RDS
  tryCatch({
    save_stage_result(results, sample_id, "06-clinical-annotation", "civic_results.rds")
    log_info("CiVIC results saved for sample {sample_id}")
  }, error = function(e) {
    log_warn("Failed to save CiVIC results: {e$message}")
  })

  # Summary logging
  log_info(
    "CiVIC annotation complete: {nrow(results$variant_evidence)} evidence items, ",
    "{length(results$gene_summaries)} gene summaries, ",
    "{nrow(results$assertions)} assertions"
  )

  results
}
