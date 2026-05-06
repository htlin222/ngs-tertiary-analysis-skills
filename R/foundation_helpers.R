# R/foundation_helpers.R — pure-data helpers for the actionable report.
# No HTML emission here; just reshape OncoKB / CiVIC / biomarker data into
# tibbles and short blurbs for the section builders to render.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(stringr); library(glue); library(purrr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

ONCOKB_LEVEL_LABEL <- c(
  LEVEL_1   = "Level 1 — FDA-approved",
  LEVEL_2   = "Level 2 — Standard care",
  LEVEL_3A  = "Level 3A — Clinical evidence",
  LEVEL_3B  = "Level 3B — Other tumor type",
  LEVEL_4   = "Level 4 — Preclinical/biological",
  LEVEL_R1  = "R1 — FDA-approved resistance",
  LEVEL_R2  = "R2 — Clinical resistance"
)

ONCOKB_LEVEL_TIER <- c(
  LEVEL_1 = "on-label", LEVEL_2 = "on-label", LEVEL_3A = "on-label",
  LEVEL_3B = "off-label", LEVEL_4 = "off-label",
  LEVEL_R1 = "resistance", LEVEL_R2 = "resistance"
)

#' Flatten OncoKB mutations / cnas / fusions into one-row-per-therapy tibble.
extract_therapies <- function(oncokb) {
  rows <- list()
  push <- function(gene, alt, alt_type, oncogenic, txdf) {
    if (!is.data.frame(txdf) || nrow(txdf) == 0) return(invisible())
    for (i in seq_len(nrow(txdf))) {
      lvl <- txdf$level[i]
      if (is.na(lvl) || !nzchar(lvl)) next
      tier <- ONCOKB_LEVEL_TIER[[lvl]] %||% NA_character_
      if (is.na(tier)) next
      rows[[length(rows) + 1L]] <<- tibble(
        gene = gene, alteration = alt, alt_type = alt_type,
        oncogenic = oncogenic %||% NA_character_,
        drugs = txdf$drugs[i],
        level = lvl, tier = tier,
        evidence = txdf$description[i] %||% ""
      )
    }
  }

  for (m in oncokb$mutations %||% list()) {
    if (is.null(m) || isTRUE(is.na(m))) next
    push(m$gene, m$alteration, "mutation", m$oncogenic, m$treatments)
  }
  for (m in oncokb$cnas %||% list()) {
    if (is.null(m) || isTRUE(is.na(m))) next
    push(m$gene, m$alteration, "cna", m$oncogenic, m$treatments)
  }
  for (m in oncokb$fusions %||% list()) {
    if (is.null(m) || isTRUE(is.na(m))) next
    push(m$alteration %||% paste0(m$gene_a, "::", m$gene_b),
         "Fusion", "fusion", m$oncogenic, m$treatments)
  }

  if (length(rows) == 0) {
    return(tibble(gene = character(), alteration = character(), alt_type = character(),
                  oncogenic = character(), drugs = character(), level = character(),
                  tier = character(), evidence = character()))
  }
  bind_rows(rows) |>
    distinct(gene, alteration, drugs, level, .keep_all = TRUE) |>
    arrange(factor(level, levels = names(ONCOKB_LEVEL_LABEL)), gene)
}

#' One row per actionable / oncogenic mutation, carrying tumorTypeSummary +
#' mutation effect description for the per-variant card.
oncokb_mutation_summary <- function(oncokb) {
  ms <- oncokb$mutations %||% list()
  rows <- list()
  for (m in ms) {
    if (is.null(m) || isTRUE(is.na(m))) next
    lvl_s <- m$highest_sensitive_level
    lvl_r <- m$highest_resistance_level
    has_sig <- (length(lvl_s) > 0 && !is.na(lvl_s) && nzchar(lvl_s)) ||
               (length(lvl_r) > 0 && !is.na(lvl_r) && nzchar(lvl_r))
    onc <- m$oncogenic %||% "Unknown"
    if (!has_sig && onc %in% c("Unknown", "Inconclusive", "Neutral")) next
    rows[[length(rows) + 1L]] <- tibble(
      gene = m$gene,
      alteration = m$alteration,
      oncogenic = onc,
      mutation_effect = m$mutation_effect %||% "Unknown",
      sensitive_level = lvl_s %||% NA_character_,
      resistance_level = lvl_r %||% NA_character_,
      tumor_type_summary = m$raw$tumorTypeSummary %||% "",
      hotspot = isTRUE(m$raw$hotspot),
      mut_eff_desc = m$raw$mutationEffect$description %||% ""
    )
  }
  if (length(rows) == 0) {
    return(tibble(gene = character(), alteration = character(), oncogenic = character(),
                  mutation_effect = character(), sensitive_level = character(),
                  resistance_level = character(), tumor_type_summary = character(),
                  hotspot = logical(), mut_eff_desc = character()))
  }
  bind_rows(rows) |> distinct(gene, alteration, .keep_all = TRUE)
}

biomarker_blurb <- function(kind, score, status) {
  if (length(score) == 0 || is.na(score)) return("Not assessed for this sample.")
  if (kind == "tmb") {
    if (score >= 10) return("TMB-High — pembrolizumab is FDA-approved across solid tumors at TMB ≥10 mut/Mb (NCCN-listed).")
    if (score >= 6)  return("Intermediate TMB — single-agent ICI has limited benefit on this marker alone.")
    return("TMB-Low — ICI unlikely to benefit on this marker alone; consider other biomarkers.")
  }
  if (kind == "msi") {
    if (grepl("MSI-High", status)) return("MSI-High — strong indication for pembrolizumab / dostarlimab; consider germline Lynch syndrome work-up.")
    return("Microsatellite stable — MSI-driven ICI not indicated.")
  }
  if (kind == "hrd") {
    if (score >= 42) return("HRD-Positive (GIS ≥42) — consider PARP inhibitor (olaparib, rucaparib, niraparib).")
    if (score >= 30) return("Borderline HRD — confirm with a dedicated HRD assay before PARPi decision.")
    return("HRD-Negative — PARPi unlikely to benefit on this marker alone.")
  }
  ""
}

# Display label mapping helper (re-exported here for the renderer to use)
oncokb_level_pretty <- function(level) ONCOKB_LEVEL_LABEL[level] %||% level
