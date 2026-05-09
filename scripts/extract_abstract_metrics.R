#!/usr/bin/env Rscript
# Print every number the ESMO abstract needs, sourced from existing benchmark CSVs.
# Usage: Rscript scripts/extract_abstract_metrics.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(here)
})

results_dir <- here("benchmarks", "results")

combined         <- read_csv(file.path(results_dir, "agent_combined_summary.csv"), show_col_types = FALSE)
per_stratum      <- read_csv(file.path(results_dir, "agent_combined_per_stratum.csv"), show_col_types = FALSE)
consistency      <- read_csv(file.path(results_dir, "agent_consistency_summary.csv"), show_col_types = FALSE)
cross_model_pw   <- tryCatch(read_csv(file.path(results_dir, "agent_cross_model_pairwise.csv"), show_col_types = FALSE), error = function(e) NULL)
cross_model_sum  <- tryCatch(read_csv(file.path(results_dir, "agent_cross_model_summary.csv"), show_col_types = FALSE), error = function(e) NULL)
ev_sum           <- read_csv(file.path(results_dir, "evidence_concordance_summary.csv"), show_col_types = FALSE)

cat("\n=== ABSTRACT NUMBERS ===\n\n")

n_total <- nrow(combined)
n_valid <- sum(combined$valid, na.rm = TRUE)
n_agree <- sum(combined$agrees_with_baseline, na.rm = TRUE)
pct_agree <- round(100 * n_agree / n_valid, 1)

cat(sprintf("Total benchmark cases       : %d\n", n_total))
cat(sprintf("Schema-valid JSON           : %d / %d (%.0f%%)\n", n_valid, n_total, 100 * n_valid / n_total))
cat(sprintf("Agent-baseline agreement    : %d / %d (%.1f%%)\n", n_agree, n_valid, pct_agree))

cat("\n--- Per stratum ---\n")
print(per_stratum %>% select(stratum, n_total, n_agree, pct_agree, n_high_conf, n_tumor_specificity, n_coverage_mismatch))

cat("\n--- Discordance breakdown (cohort) ---\n")
cohort <- combined %>% filter(stratum == "cohort")
print(cohort %>% count(agent_discordance_nature, sort = TRUE))

cat("\n--- OncoKB-only cases (no CiVIC evidence) — full breakdown ---\n")
oncokb_only <- cohort %>% filter(is.na(civic_amp_level))
n_oo <- nrow(oncokb_only)
n_oo_override <- sum(!oncokb_only$agrees_with_baseline, na.rm = TRUE)
cat(sprintf("OncoKB-only cases (any tier)  : %d\n", n_oo))
cat(sprintf("  Agent overrode baseline     : %d / %d   (mostly Tier III VUS <-> IID boundary)\n",
            n_oo_override, n_oo))

# Safety-critical = OncoKB-only AND baseline tier is clinically actionable (I or II).
# Overrides on Tier III VUS in this subset are the documented Tier III<->IID
# boundary, not a safety risk; the abstract reports those separately.
safety_critical <- oncokb_only %>%
  filter(grepl("^Tier I[ A-]|^Tier II[ A-]", baseline_tier))
n_sc <- nrow(safety_critical)
n_sc_override <- sum(!safety_critical$agrees_with_baseline, na.rm = TRUE)
cat(sprintf("\nSafety-critical (Tier I/II)   : %d\n", n_sc))
cat(sprintf("  Agent overrode baseline     : %d / %d   <- abstract cites this number\n",
            n_sc_override, n_sc))

cat("\n--- Tumor-specificity mismatches (cross-tumor CiVIC evidence) ---\n")
ts <- cohort %>% filter(agent_discordance_nature == "tumor_specificity_mismatch")
print(ts %>% select(case_id, gene, variant, tumor_type, baseline_tier, agent_tier, agrees_with_baseline))

cat("\n--- Tier III VUS <-> Tier II Level D boundary ---\n")
vus_boundary <- combined %>%
  filter(baseline_tier == "Tier III VUS" & agent_tier == "Tier II Level D")
cat(sprintf("Cases at boundary           : %d\n", nrow(vus_boundary)))

cat("\n--- Within-model consistency (3-run Opus) ---\n")
print(consistency)

if (!is.null(cross_model_sum)) {
  cat("\n--- Cross-model reproducibility ---\n")
  print(cross_model_sum)
}

cat("\n--- KB coverage on real cohort ---\n")
print(ev_sum)

cat("\n--- Illustrative case for Results paragraph ---\n")
illus <- cohort %>%
  filter(agent_discordance_nature == "tumor_specificity_mismatch",
         !agrees_with_baseline)
if (nrow(illus) > 0) {
  print(illus %>% select(case_id, gene, variant, tumor_type, baseline_tier, agent_tier, civic_amp_level))
} else {
  cat("(no down-tier tumor-specificity case in this slice — fall back to RET M918T NSCLC if needed)\n")
}

cat("\n=== END ===\n")
