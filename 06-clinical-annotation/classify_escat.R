# 06-clinical-annotation/classify_escat.R — ESCAT tier classification with treatment enrichment

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(logger)
  library(purrr)
})

# Source required functions
source(here::here("R/esmo_helpers.R"))
source(here::here("R/utils.R"))

#' Extract treatment information from OncoKB result
#'
#' Parses the treatments list from an OncoKB annotation result into a
#' human-readable summary with drug names, approval status, and evidence.
#'
#' @param oncokb_result List from OncoKB API response (e.g., from mutations/cnas/fusions)
#'
#' @return Tibble with columns: drugs, level, fda_approved, description
#'
#' @keywords internal
extract_treatments_from_result <- function(oncokb_result) {
  if (is.null(oncokb_result$treatments) || nrow(oncokb_result$treatments) == 0) {
    return(tibble(
      drugs = character(),
      level = character(),
      fda_approved = logical(),
      description = character()
    ))
  }

  as_tibble(oncokb_result$treatments)
}

#' Summarize actionable alterations by ESCAT tier
#'
#' Creates a summary table counting alterations in each ESCAT tier,
#' with lists of actionable genes for clinical implementation.
#'
#' @param escat_tibble Tibble from classify_all_escat with ESCAT classifications
#'
#' @return Tibble with columns: tier, count, genes, sample_size
#'
#' @keywords internal
summarize_tiers <- function(escat_tibble) {
  tier_order <- c("I", "I-R", "II", "II-R", "III", "IV", "V", "X")

  summary <- escat_tibble %>%
    group_by(escat_tier) %>%
    summarise(
      count = n(),
      genes = paste(unique(gene), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(
      escat_tier = factor(escat_tier, levels = tier_order),
      tier_description = map_chr(escat_tier, escat_description)
    ) %>%
    arrange(escat_tier) %>%
    select(
      tier = escat_tier,
      count,
      genes,
      description = tier_description
    )

  return(summary)
}

#' Identify actionable alterations for clinical decisions
#'
#' Filters alterations to ESCAT tiers I-III (ready for implementation or
#' strong investigational evidence), excluding resistance biomarkers.
#'
#' @param escat_tibble Tibble from classify_all_escat
#'
#' @return Tibble of actionable alterations with full ESCAT and treatment data
#'
#' @keywords internal
get_actionable_alterations <- function(escat_tibble) {
  escat_tibble %>%
    filter(escat_tier %in% c("I", "II", "III")) %>%
    select(
      gene, alteration, type, oncogenic, escat_tier,
      escat_description, sensitive_level
    ) %>%
    arrange(factor(escat_tier, levels = c("I", "II", "III")), gene)
}

#' Identify potential clinical trial targets
#'
#' Filters alterations to ESCAT tiers II-III (investigational evidence),
#' which may be eligible for clinical trial enrollment.
#'
#' @param escat_tibble Tibble from classify_all_escat
#'
#' @return Tibble of trial-eligible alterations
#'
#' @keywords internal
get_trial_eligible <- function(escat_tibble) {
  escat_tibble %>%
    filter(escat_tier %in% c("II", "III")) %>%
    select(
      gene, alteration, type, escat_tier, escat_description
    ) %>%
    arrange(factor(escat_tier, levels = c("II", "III")), gene)
}

#' Enrich ESCAT classifications with OncoKB treatment information
#'
#' Takes the base ESCAT classification and adds detailed treatment data
#' from OncoKB results, including drug names, approval status, and clinical
#' evidence levels.
#'
#' @param escat_tibble Tibble from classify_all_escat
#' @param oncokb_results List of OncoKB results (mutations, cnas, fusions)
#'
#' @return Enriched tibble with additional treatment columns
#'
#' @keywords internal
enrich_with_treatments <- function(escat_tibble, oncokb_results) {
  # Build lookup of treatments by alteration
  all_results <- c(
    oncokb_results$mutations %||% list(),
    oncokb_results$cnas %||% list(),
    oncokb_results$fusions %||% list()
  )

  treatment_lookup <- map_df(all_results, function(result) {
    treatments <- extract_treatments_from_result(result)
    if (nrow(treatments) == 0) {
      return(tibble(
        alteration = character(),
        drugs = character(),
        level = character(),
        fda_approved = logical()
      ))
    }
    treatments %>%
      mutate(alteration = result$alteration %||% NA_character_, .before = 1)
  })

  # Join with ESCAT tibble and keep all ESCAT data
  if (nrow(treatment_lookup) > 0) {
    enriched <- escat_tibble %>%
      left_join(
        treatment_lookup %>%
          select(alteration, drugs, level, fda_approved) %>%
          distinct(),
        by = "alteration",
        relationship = "many-to-many"
      ) %>%
      # Consolidate treatments for alterations with multiple therapies
      group_by(across(-c(drugs, level, fda_approved))) %>%
      mutate(
        drugs = paste(unique(na.omit(drugs)), collapse = " | "),
        level = paste(unique(na.omit(level)), collapse = " | "),
        .groups = "drop"
      ) %>%
      ungroup() %>%
      distinct()
  } else {
    # No treatments found, add empty columns
    enriched <- escat_tibble %>%
      mutate(
        drugs = NA_character_,
        level = NA_character_,
        fda_approved = NA
      )
  }

  return(enriched)
}

#' Classify all alterations with ESCAT tiers and treatment information
#'
#' Main entry point for ESCAT classification. Calls classify_all_escat from
#' esmo_helpers.R and enriches results with OncoKB treatment data. Generates
#' summary statistics and identifies actionable alterations for clinical
#' decision-making.
#'
#' @param oncokb_results List of OncoKB annotation results
#'   - mutations: List of mutation annotation results
#'   - cnas: List of CNA annotation results
#'   - fusions: List of fusion annotation results
#'   Each result should contain: gene, alteration, oncogenic, highest_sensitive_level,
#'   highest_resistance_level, treatments
#'
#' @param config List of pipeline configuration
#'   - Used by classify_all_escat for sample metadata and output paths
#'
#' @param sample_id Sample identifier for logging and output organization
#'
#' @return List with classified ESCAT tibble plus summary statistics:
#'   - escat_results: Enriched tibble with ESCAT tiers and treatments
#'   - tier_summary: Tibble with counts and genes per ESCAT tier
#'   - actionable: Tibble of tier I-III alterations ready for clinical use
#'   - trial_eligible: Tibble of tier II-III alterations for trial matching
#'
#' @details
#' The function performs the following steps:
#' 1. Call classify_all_escat to assign ESCAT tiers based on OncoKB levels
#' 2. Enrich with OncoKB treatment information (drug names, approval status)
#' 3. Generate tier summary with counts and gene lists
#' 4. Extract actionable alterations (tiers I-III)
#' 5. Extract trial-eligible alterations (tiers II-III)
#' 6. Save enriched ESCAT table to TSV
#' 7. Log clinical actionability summary
#' 8. Return enriched tibble with treatment columns added
#'
#' ESCAT tier meanings:
#' - I: FDA-approved standard of care or strong guideline recommendation
#' - I-R: Resistance biomarker, standard of care
#' - II: Investigational with strong clinical evidence; enables trial enrollment
#' - II-R: Resistance biomarker, investigational
#' - III: Clinical benefit in other tumor types or similar molecular features
#' - IV: Preclinical evidence
#' - V: Co-targeting approaches
#' - X: No evidence or insufficient data
#'
#' @export
classify_escat <- function(oncokb_results, config, sample_id) {
  log_info("Classifying ESCAT tiers for sample: {sample_id}")

  # Step 1: Base ESCAT classification from OncoKB levels
  escat_base <- classify_all_escat(oncokb_results, config, sample_id)
  log_info("Classified {nrow(escat_base)} alterations into ESCAT tiers")

  # Step 2: Enrich with treatment information
  escat_enriched <- enrich_with_treatments(escat_base, oncokb_results)
  log_debug("Enriched ESCAT results with treatment information")

  # Step 3: Generate summary statistics
  tier_summary <- summarize_tiers(escat_enriched)
  log_info("ESCAT tier distribution:\n{paste(capture.output(print(tier_summary)), collapse = '\n')}")

  # Step 4: Extract actionable and trial-eligible subsets
  actionable <- get_actionable_alterations(escat_enriched)
  trial_eligible <- get_trial_eligible(escat_enriched)

  log_info(
    "Clinical actionability: {nrow(actionable)} tier I-III alterations, ",
    "{nrow(trial_eligible)} trial-eligible (tier II-III)"
  )

  # Step 5: Save enriched ESCAT results to TSV
  output_dir <- stage_output_dir(sample_id, "06-clinical-annotation")
  output_file <- file.path(output_dir, "escat_tiers.tsv")

  tryCatch(
    {
      write_tsv(escat_enriched, output_file)
      log_info("ESCAT classifications saved to {output_file}")
    },
    error = function(e) {
      log_warn("Failed to save ESCAT results to TSV: {e$message}")
    }
  )

  # Step 6: Log clinical actionability summary
  if (nrow(actionable) > 0) {
    log_info("Actionable alterations (ESCAT I-III):")
    for (i in seq_len(min(nrow(actionable), 5))) {
      row <- actionable[i, ]
      log_info("  - {row$gene} {row$alteration} ({row$escat_tier}): {row$escat_description}")
    }
    if (nrow(actionable) > 5) {
      log_info("  ... and {nrow(actionable) - 5} more")
    }
  } else {
    log_info("No ESCAT I-III (actionable) alterations identified")
  }

  # Return enriched results with metadata
  results <- list(
    escat_results = escat_enriched,
    tier_summary = tier_summary,
    actionable = actionable,
    trial_eligible = trial_eligible
  )

  log_info("ESCAT classification and enrichment complete")

  return(results)
}
