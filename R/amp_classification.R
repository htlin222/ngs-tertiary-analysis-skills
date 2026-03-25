# R/amp_classification.R — AMP/ASCO/CAP four-tier somatic variant classification
# Per Li et al., J Mol Diagn 2017; 19(1):4-23

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(logger)
})

# ==============================================================================
# Single-Variant Classification
# ==============================================================================

#' Classify a single variant using the AMP/ASCO/CAP four-tier system
#'
#' @param oncokb_oncogenic Character. OncoKB oncogenicity call (e.g., "Oncogenic",
#'   "Likely Oncogenic", "Neutral", "Inconclusive").
#' @param oncokb_level Character. OncoKB therapeutic level (e.g., "LEVEL_1").
#' @param civic_amp_level Character. CiVIC AMP assertion level (e.g., "TIER_I_LEVEL_A").
#' @param vep_impact Character. VEP functional impact (e.g., "HIGH", "MODERATE").
#' @param clinvar Character. ClinVar clinical significance (e.g., "pathogenic",
#'   "benign", "uncertain_significance").
#'
#' @return Named list with amp_tier, amp_level, amp_evidence.
#'
#' @details
#' Classification tiers (Li et al., J Mol Diagn 2017):
#' - Tier I Level A: FDA-approved therapy or professional guidelines
#' - Tier I Level B: Well-powered studies with consensus
#' - Tier II Level C: FDA-approved for different tumor type
#' - Tier II Level D: Preclinical/case studies or oncogenic with high impact
#' - Tier III: Variants of uncertain significance (VUS)
#' - Tier IV: Benign or likely benign
#'
#' @export
classify_amp_tier <- function(oncokb_oncogenic = NA_character_,
                              oncokb_level = NA_character_,
                              civic_amp_level = NA_character_,
                              vep_impact = NA_character_,
                              clinvar = NA_character_) {
  # Normalise inputs
  oncokb_oncogenic <- as.character(oncokb_oncogenic %||% NA_character_)
  oncokb_level     <- as.character(oncokb_level %||% NA_character_)
  civic_amp_level  <- as.character(civic_amp_level %||% NA_character_)
  vep_impact       <- toupper(as.character(vep_impact %||% NA_character_))
  clinvar          <- tolower(as.character(clinvar %||% NA_character_))

  evidence_parts <- character()

  # --- Tier I Level A --------------------------------------------------------
  # FDA-approved therapy or professional guideline recommendation
  is_tier_i_a <- FALSE
  if (!is.na(oncokb_level) && oncokb_level %in% c("LEVEL_1", "LEVEL_2")) {
    is_tier_i_a <- TRUE
    evidence_parts <- c(evidence_parts,
                        paste0("OncoKB ", oncokb_level, " (FDA-approved/guideline)"))
  }
  if (!is.na(civic_amp_level) && grepl("TIER_I_LEVEL_A", civic_amp_level, fixed = TRUE)) {
    is_tier_i_a <- TRUE
    evidence_parts <- c(evidence_parts, "CiVIC Tier I Level A assertion")
  }
  if (is_tier_i_a) {
    return(list(
      amp_tier     = "Tier I",
      amp_level    = "Level A",
      amp_evidence = paste(evidence_parts, collapse = "; ")
    ))
  }

  # --- Tier I Level B --------------------------------------------------------
  # Well-powered studies with expert consensus
  is_tier_i_b <- FALSE
  if (!is.na(oncokb_level) && oncokb_level == "LEVEL_3A") {
    is_tier_i_b <- TRUE
    evidence_parts <- c(evidence_parts,
                        "OncoKB LEVEL_3A (well-powered clinical studies)")
  }
  if (!is.na(civic_amp_level) && grepl("TIER_I_LEVEL_B", civic_amp_level, fixed = TRUE)) {
    is_tier_i_b <- TRUE
    evidence_parts <- c(evidence_parts, "CiVIC Tier I Level B assertion")
  }
  if (is_tier_i_b) {
    return(list(
      amp_tier     = "Tier I",
      amp_level    = "Level B",
      amp_evidence = paste(evidence_parts, collapse = "; ")
    ))
  }

  # --- Tier II Level C -------------------------------------------------------
  # FDA-approved therapy for a different tumor type
  is_tier_ii_c <- FALSE
  if (!is.na(oncokb_level) && oncokb_level == "LEVEL_3B") {
    is_tier_ii_c <- TRUE
    evidence_parts <- c(evidence_parts,
                        "OncoKB LEVEL_3B (approved in different tumor type)")
  }
  if (!is.na(civic_amp_level) && grepl("TIER_II_LEVEL_C", civic_amp_level, fixed = TRUE)) {
    is_tier_ii_c <- TRUE
    evidence_parts <- c(evidence_parts, "CiVIC Tier II Level C assertion")
  }
  if (is_tier_ii_c) {
    return(list(
      amp_tier     = "Tier II",
      amp_level    = "Level C",
      amp_evidence = paste(evidence_parts, collapse = "; ")
    ))
  }

  # --- Tier II Level D -------------------------------------------------------
  # Preclinical/case study evidence, or oncogenic with high functional impact
  is_tier_ii_d <- FALSE
  if (!is.na(oncokb_level) && oncokb_level == "LEVEL_4") {
    is_tier_ii_d <- TRUE
    evidence_parts <- c(evidence_parts,
                        "OncoKB LEVEL_4 (biological evidence)")
  }
  if (!is.na(civic_amp_level) && grepl("TIER_II_LEVEL_D", civic_amp_level, fixed = TRUE)) {
    is_tier_ii_d <- TRUE
    evidence_parts <- c(evidence_parts, "CiVIC Tier II Level D assertion")
  }
  is_oncogenic <- !is.na(oncokb_oncogenic) &&
    oncokb_oncogenic %in% c("Oncogenic", "Likely Oncogenic")
  is_high_impact <- !is.na(vep_impact) && vep_impact == "HIGH"
  if (!is_tier_ii_d && is_oncogenic && is_high_impact) {
    is_tier_ii_d <- TRUE
    evidence_parts <- c(evidence_parts,
                        paste0("Oncogenic (", oncokb_oncogenic, ") + HIGH VEP impact"))
  }
  if (is_tier_ii_d) {
    return(list(
      amp_tier     = "Tier II",
      amp_level    = "Level D",
      amp_evidence = paste(evidence_parts, collapse = "; ")
    ))
  }

  # --- Tier IV Benign --------------------------------------------------------
  # ClinVar benign/likely benign OR OncoKB Neutral/Inconclusive
  is_tier_iv <- FALSE
  if (!is.na(clinvar) && clinvar %in% c("benign", "likely_benign", "likely benign")) {
    is_tier_iv <- TRUE
    evidence_parts <- c(evidence_parts, paste0("ClinVar: ", clinvar))
  }
  if (!is.na(oncokb_oncogenic) &&
      oncokb_oncogenic %in% c("Neutral", "Inconclusive")) {
    is_tier_iv <- TRUE
    evidence_parts <- c(evidence_parts,
                        paste0("OncoKB: ", oncokb_oncogenic))
  }
  if (is_tier_iv) {
    return(list(
      amp_tier     = "Tier IV",
      amp_level    = "Benign",
      amp_evidence = paste(evidence_parts, collapse = "; ")
    ))
  }

  # --- Tier III VUS ----------------------------------------------------------
  # Insufficient evidence — everything else
  evidence_parts <- c(evidence_parts, "Insufficient evidence for classification")
  if (!is.na(oncokb_oncogenic)) {
    evidence_parts <- c(evidence_parts, paste0("OncoKB: ", oncokb_oncogenic))
  }
  if (!is.na(clinvar)) {
    evidence_parts <- c(evidence_parts, paste0("ClinVar: ", clinvar))
  }

  list(
    amp_tier     = "Tier III",
    amp_level    = "VUS",
    amp_evidence = paste(evidence_parts, collapse = "; ")
  )
}

# ==============================================================================
# Batch Classification
# ==============================================================================

#' Classify all variants using the AMP/ASCO/CAP system
#'
#' @param variants Tibble with columns: gene, hgvsp, oncogenic, sensitive_level,
#'   impact, clinvar_significance.
#' @param civic_assertions Optional tibble from CiVIC with columns: gene,
#'   variant_name, amp_level. Default NULL.
#'
#' @return The input tibble with three new columns: amp_tier, amp_level, amp_evidence.
#'
#' @export
classify_all_amp <- function(variants, civic_assertions = NULL) {
  # Handle empty input
  if (is.null(variants) || nrow(variants) == 0) {
    result <- variants
    if (is.null(result)) {
      result <- tibble(
        gene = character(), hgvsp = character(), oncogenic = character(),
        sensitive_level = character(), impact = character(),
        clinvar_significance = character()
      )
    }
    result$amp_tier     <- character()
    result$amp_level    <- character()
    result$amp_evidence <- character()
    return(result)
  }

  # Build CiVIC lookup: gene + stripped variant name -> amp_level
  civic_lookup <- NULL
  if (!is.null(civic_assertions) && nrow(civic_assertions) > 0) {
    civic_lookup <- civic_assertions |>
      filter(!is.na(amp_level), !is.na(gene), !is.na(variant_name)) |>
      mutate(variant_key = str_remove(variant_name, "^p\\.")) |>
      select(gene, variant_key, amp_level) |>
      distinct(gene, variant_key, .keep_all = TRUE)
  }

  # Strip p. prefix from hgvsp for CiVIC matching
  variant_keys <- str_remove(variants$hgvsp %||% rep(NA_character_, nrow(variants)),
                             "^p\\.")

  # Look up CiVIC AMP level for each variant
  civic_levels <- if (!is.null(civic_lookup) && nrow(civic_lookup) > 0) {
    pmap_chr(list(g = variants$gene, vk = variant_keys), function(g, vk) {
      if (is.na(g) || is.na(vk)) return(NA_character_)
      match_row <- civic_lookup |>
        filter(gene == g, variant_key == vk)
      if (nrow(match_row) > 0) match_row$amp_level[1] else NA_character_
    })
  } else {
    rep(NA_character_, nrow(variants))
  }

  # Classify each variant using pmap (vectorised over rows)
  classifications <- pmap(
    list(
      oncokb_oncogenic = variants$oncogenic,
      oncokb_level     = variants$sensitive_level,
      civic_amp_level  = civic_levels,
      vep_impact       = variants$impact,
      clinvar          = variants$clinvar_significance
    ),
    classify_amp_tier
  )

  # Extract results and bind to original tibble
  variants |>
    mutate(
      amp_tier     = map_chr(classifications, "amp_tier"),
      amp_level    = map_chr(classifications, "amp_level"),
      amp_evidence = map_chr(classifications, "amp_evidence")
    )
}
