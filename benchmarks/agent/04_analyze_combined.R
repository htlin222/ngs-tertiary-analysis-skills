#!/usr/bin/env Rscript
# benchmarks/agent/04_analyze_combined.R
#
# Combined analysis of cohort (real_*) + stress (stress_*) cohorts.
# Reports per-stratum metrics and writes a summary CSV + per-stratum stats.
# This is the headline output that feeds the manuscript.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(jsonlite)
  library(logger)
  library(glue)
  library(fs)
})

source(here::here("R/utils.R"))
source(here::here("agents/io.R"))

setup_logging("INFO")

input_dir  <- here::here("agents/inputs/evidence_reconciler")
output_dir <- here::here("agents/outputs/evidence_reconciler")

# Pick up any case (real_* + stress_*); skip the 5 synthetic smoke cases
input_files <- list.files(input_dir, pattern = "^(real|stress)_.*\\.json$",
                          full.names = FALSE)
case_ids    <- tools::file_path_sans_ext(input_files)
log_info("Found {length(case_ids)} cases ({sum(startsWith(case_ids, 'real_'))} real + {sum(startsWith(case_ids, 'stress_'))} stress)")

rows <- map(case_ids, function(cid) {
  input_path  <- file.path(input_dir, paste0(cid, ".json"))
  output_path <- file.path(output_dir, paste0(cid, ".json"))
  input  <- fromJSON(input_path, simplifyVector = FALSE)
  output_result <- read_agent_output("evidence_reconciler", cid)

  baseline_tier_level <- paste(input$metadata$baseline_tier,
                               input$metadata$baseline_level)
  stratum <- if (startsWith(cid, "real_")) "cohort" else "stress"

  if (!output_result$valid) {
    return(tibble(
      case_id = cid, stratum = stratum,
      gene = input$metadata$gene, variant = input$metadata$variant,
      tumor_type = input$metadata$tumor_type,
      oncokb_level = input$metadata$oncokb_level,
      civic_amp_level = input$metadata$civic_amp_level,
      civic_evidence_in_different_tumor = input$metadata$civic_evidence_in_different_tumor,
      baseline_tier = baseline_tier_level,
      agent_tier = NA_character_, agent_confidence = NA_character_,
      agent_discordance_nature = NA_character_,
      primary_source = NA_character_, dissenting_noted = NA,
      valid = FALSE,
      errors = paste(output_result$errors, collapse = "; ")
    ))
  }
  o <- output_result$output
  tibble(
    case_id = cid, stratum = stratum,
    gene = input$metadata$gene, variant = input$metadata$variant,
    tumor_type = input$metadata$tumor_type,
    oncokb_level = input$metadata$oncokb_level,
    civic_amp_level = input$metadata$civic_amp_level,
    civic_evidence_in_different_tumor = input$metadata$civic_evidence_in_different_tumor,
    baseline_tier = baseline_tier_level,
    agent_tier = o$unified_amp_tier,
    agent_confidence = o$confidence,
    agent_discordance_nature = o$discordance_nature,
    primary_source = o$primary_evidence_source,
    dissenting_noted = o$dissenting_evidence_noted,
    valid = TRUE, errors = NA_character_
  )
})

summary_tbl <- bind_rows(rows) |>
  mutate(agrees_with_baseline = baseline_tier == agent_tier)

# ── Per-stratum metrics ──────────────────────────────────────────────────────

per_stratum <- summary_tbl |>
  group_by(stratum) |>
  summarise(
    n_total                = n(),
    n_valid                = sum(valid),
    n_agree                = sum(agrees_with_baseline, na.rm = TRUE),
    pct_agree              = round(100 * mean(agrees_with_baseline, na.rm = TRUE), 1),
    n_high_conf            = sum(agent_confidence == "high",     na.rm = TRUE),
    n_moderate_conf        = sum(agent_confidence == "moderate", na.rm = TRUE),
    n_low_conf             = sum(agent_confidence == "low",      na.rm = TRUE),
    n_concordant_flagged   = sum(agent_discordance_nature == "concordant",                 na.rm = TRUE),
    n_tumor_specificity    = sum(agent_discordance_nature == "tumor_specificity_mismatch", na.rm = TRUE),
    n_therapeutic_level    = sum(agent_discordance_nature == "therapeutic_level_mismatch", na.rm = TRUE),
    n_oncogenicity         = sum(agent_discordance_nature == "oncogenicity_mismatch",      na.rm = TRUE),
    n_coverage_mismatch    = sum(agent_discordance_nature == "coverage_mismatch",          na.rm = TRUE),
    .groups = "drop"
  )

cat("\n─── Per-stratum metrics ───────────────────────────────────────\n")
print(per_stratum)

# ── Disagreement breakdown ──────────────────────────────────────────────────

disagreements <- summary_tbl |>
  filter(valid, !agrees_with_baseline)

cat(sprintf("\nTotal disagreements: %d / %d (%.1f%%)\n",
            nrow(disagreements), nrow(summary_tbl),
            100 * nrow(disagreements) / nrow(summary_tbl)))

cat("\nDisagreements by tier transition:\n")
print(disagreements |>
        count(baseline_tier, agent_tier, agent_confidence) |>
        arrange(desc(n)),
      n = 30)

# ── Tier-I safety: agent NEVER overrides Tier I-A baseline downward (?)─────

tier_i_baseline <- summary_tbl |>
  filter(valid, grepl("Tier I", baseline_tier))
n_tier_i <- nrow(tier_i_baseline)
n_tier_i_agreed <- sum(tier_i_baseline$agrees_with_baseline)
cat(sprintf("\nTier I baseline cases: %d, agent-agreed: %d (%.1f%%)\n",
            n_tier_i, n_tier_i_agreed,
            if (n_tier_i > 0) 100 * n_tier_i_agreed / n_tier_i else NA_real_))

cat("\nTier I baseline cases where agent disagreed:\n")
print(tier_i_baseline |>
        filter(!agrees_with_baseline) |>
        select(case_id, baseline_tier, agent_tier, agent_confidence,
               agent_discordance_nature),
      n = 20)

# ── Save ─────────────────────────────────────────────────────────────────────

out_dir <- here::here("benchmarks/results")
dir_create(out_dir)
utils::write.csv(summary_tbl,  file.path(out_dir, "agent_combined_summary.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(per_stratum, file.path(out_dir, "agent_combined_per_stratum.csv"),
                 row.names = FALSE, na = "")
log_info("Wrote {file.path(out_dir, 'agent_combined_summary.csv')}")
log_info("Wrote {file.path(out_dir, 'agent_combined_per_stratum.csv')}")
