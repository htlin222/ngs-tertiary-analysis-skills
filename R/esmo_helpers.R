# R/esmo_helpers.R — ESCAT classification and ESMO report formatting helpers

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(logger)
})

# ══════════════════════════════════════════════════════════════════════════════
# ESCAT Tier Classification
# Maps OncoKB therapeutic levels to ESMO Scale for Clinical Actionability
# ══════════════════════════════════════════════════════════════════════════════

#' Map OncoKB level to ESCAT tier
#' @param oncokb_level OncoKB level string (e.g., "LEVEL_1", "LEVEL_3A")
#' @return ESCAT tier string
oncokb_to_escat <- function(oncokb_level) {
  if (is.na(oncokb_level) || is.null(oncokb_level)) return("X")

  mapping <- c(
    "LEVEL_1"  = "I",     # FDA-approved, standard of care
    "LEVEL_2"  = "I",     # Standard of care
    "LEVEL_3A" = "II",    # Investigational, compelling clinical evidence
    "LEVEL_3B" = "III",   # Standard of care in different tumor type
    "LEVEL_4"  = "IV",    # Compelling biological evidence
    "LEVEL_R1" = "I-R",   # Standard of care resistance
    "LEVEL_R2" = "II-R"   # Investigational resistance
  )

  tier <- mapping[oncokb_level]
  if (is.na(tier)) return("X")
  unname(tier)
}

#' Get ESCAT tier description
#' @param tier ESCAT tier string
#' @return Human-readable description
escat_description <- function(tier) {
  descriptions <- c(
    "I"    = "Target ready for implementation in routine clinical decisions",
    "II"   = "Investigational target with strong evidence; enables clinical trial enrollment",
    "III"  = "Clinical benefit demonstrated in other tumor types or similar molecular alterations",
    "IV"   = "Preclinical evidence of actionability",
    "V"    = "Evidence supporting co-targeting approaches",
    "X"    = "No evidence for actionability or lack of data",
    "I-R"  = "Resistance biomarker, standard of care",
    "II-R" = "Resistance biomarker, investigational"
  )
  result <- descriptions[tier]
  if (is.na(result)) return("Unclassified")
  unname(result)
}

#' Classify all variants with ESCAT tiers
#' @param oncokb_results List of OncoKB annotation results
#' @param config Pipeline config
#' @param sample_id Sample identifier
#' @return Tibble with alteration, gene, escat_tier, description, treatments
classify_all_escat <- function(oncokb_results, config, sample_id) {
  log_info("Classifying ESCAT tiers for {sample_id}")

  results <- bind_rows(
    # SNVs/Indels
    map_dfr(oncokb_results$mutations %||% list(), function(m) {
      tier <- oncokb_to_escat(m$highest_sensitive_level)
      tibble(
        gene = m$gene,
        alteration = m$alteration,
        type = "mutation",
        oncogenic = m$oncogenic,
        escat_tier = tier,
        escat_description = escat_description(tier),
        sensitive_level = m$highest_sensitive_level %||% NA_character_,
        resistance_level = m$highest_resistance_level %||% NA_character_
      )
    }),
    # CNAs
    map_dfr(oncokb_results$cnas %||% list(), function(c) {
      tier <- oncokb_to_escat(c$highest_sensitive_level)
      tibble(
        gene = c$gene,
        alteration = c$alteration,
        type = "cna",
        oncogenic = c$oncogenic,
        escat_tier = tier,
        escat_description = escat_description(tier),
        sensitive_level = c$highest_sensitive_level %||% NA_character_,
        resistance_level = NA_character_
      )
    }),
    # Fusions
    map_dfr(oncokb_results$fusions %||% list(), function(f) {
      tier <- oncokb_to_escat(f$highest_sensitive_level)
      tibble(
        gene = f$alteration,
        alteration = f$alteration,
        type = "fusion",
        oncogenic = f$oncogenic,
        escat_tier = tier,
        escat_description = escat_description(tier),
        sensitive_level = f$highest_sensitive_level %||% NA_character_,
        resistance_level = NA_character_
      )
    })
  )

  # Sort by clinical actionability (ESCAT I first)
  tier_order <- c("I", "I-R", "II", "II-R", "III", "IV", "V", "X")
  results <- results |>
    mutate(escat_rank = match(escat_tier, tier_order, nomatch = 99)) |>
    arrange(escat_rank) |>
    select(-escat_rank)

  # Save to reports
  out_dir <- stage_output_dir(sample_id, "06-clinical-annotation")
  write.csv(results, file.path(out_dir, "escat_tiers.csv"), row.names = FALSE)
  log_info("ESCAT classification: {nrow(results)} alterations classified")

  results
}

# ══════════════════════════════════════════════════════════════════════════════
# ESMO Report Formatting Helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Format variant for ESMO report display
#' Produces both HGVS protein and coding notation
#' @param gene Gene symbol
#' @param hgvsp Protein change (HGVS p. notation)
#' @param hgvsc Coding change (HGVS c. notation)
#' @param vaf Variant allele frequency
#' @return Formatted string for report
format_variant_esmo <- function(gene, hgvsp, hgvsc = NA, vaf = NA) {
  parts <- gene
  if (!is.na(hgvsp)) parts <- glue("{parts} {hgvsp}")
  if (!is.na(hgvsc)) parts <- glue("{parts} ({hgvsc})")
  if (!is.na(vaf)) parts <- glue("{parts} [VAF: {round(vaf * 100, 1)}%]")
  as.character(parts)
}

#' Format fusion for ESMO report (gene1::gene2 notation)
#' @param gene_a 5' partner gene
#' @param gene_b 3' partner gene
#' @param exon_a Exon of gene A (optional)
#' @param exon_b Exon of gene B (optional)
#' @return Formatted fusion string
format_fusion_esmo <- function(gene_a, gene_b, exon_a = NA, exon_b = NA) {
  fusion <- glue("{gene_a}::{gene_b}")
  if (!is.na(exon_a) && !is.na(exon_b)) {
    fusion <- glue("{fusion} (exon {exon_a}::exon {exon_b})")
  }
  as.character(fusion)
}

#' Format CNV for ESMO report table
#' @param gene Gene symbol
#' @param type "AMPLIFICATION" or "DELETION"
#' @param log2_ratio Log2 ratio
#' @param copy_number Estimated copy number
#' @return Named list for table row
format_cnv_esmo <- function(gene, type, log2_ratio = NA, copy_number = NA) {
  list(
    gene = gene,
    type = str_to_title(type),
    log2_ratio = if (!is.na(log2_ratio)) round(log2_ratio, 2) else NA,
    copy_number = if (!is.na(copy_number)) round(copy_number, 1) else NA
  )
}

#' Generate ESCAT tier badge HTML for report
#' @param tier ESCAT tier string
#' @return HTML span with colored badge
escat_badge_html <- function(tier) {
  colors <- c(
    "I"    = "#d32f2f",  # Red - highest actionability
    "I-R"  = "#e57373",
    "II"   = "#f57c00",  # Orange
    "II-R" = "#ffb74d",
    "III"  = "#fbc02d",  # Yellow
    "IV"   = "#7cb342",  # Green
    "V"    = "#42a5f5",  # Blue
    "X"    = "#9e9e9e"   # Grey
  )
  color <- colors[tier] %||% "#9e9e9e"
  glue('<span style="background-color:{color};color:white;padding:2px 8px;',
       'border-radius:4px;font-weight:bold;font-size:0.9em;">',
       'ESCAT {tier}</span>')
}

#' Create the clinical actionability summary table data
#' @param escat_results Tibble from classify_all_escat
#' @param oncokb_results Full OncoKB results for treatment details
#' @return Tibble ready for gt table rendering
actionability_table <- function(escat_results, oncokb_results) {
  escat_results |>
    filter(escat_tier != "X") |>
    mutate(
      badge = map_chr(escat_tier, escat_badge_html),
      .before = 1
    ) |>
    select(
      `ESCAT Tier` = badge,
      Gene = gene,
      Alteration = alteration,
      Type = type,
      Oncogenic = oncogenic,
      Description = escat_description
    )
}

#' Identify genes with insufficient coverage for ESMO reporting
#' @param coverage_df Tibble with gene and mean_coverage columns
#' @param min_coverage Minimum required coverage (default: 200)
#' @return Tibble of genes below threshold
coverage_gaps <- function(coverage_df, min_coverage = 200) {
  gaps <- coverage_df |>
    filter(mean_coverage < min_coverage) |>
    arrange(mean_coverage)

  if (nrow(gaps) > 0) {
    log_warn("{nrow(gaps)} genes below {min_coverage}x coverage threshold")
  }
  gaps
}
