#!/usr/bin/env Rscript
# benchmarks/agent/05_analyze_consistency.R
#
# Phase 3 consistency analysis. For each case_id in the consistency manifest,
# read the 3 run outputs (case_id__run1.json … run3.json), compute:
#   - tier-agreement count (3/3 unanimous, 2/3 majority, 1/3 split)
#   - confidence-agreement count
#   - discordance_nature-agreement count
#   - Fleiss kappa on tier across runs
# and write a per-case summary CSV.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(jsonlite)
  library(logger)
  library(glue)
})

source(here::here("R/utils.R"))
setup_logging("INFO")

manifest_path <- "/tmp/batches/manifest_consistency.txt"
case_ids <- readLines(manifest_path) |> trimws() |> Filter(\(x) nchar(x) > 0, x = _)
runs <- c("run1", "run2", "run3")

output_dir <- here::here("agents/outputs/evidence_reconciler")

read_run <- function(case_id, run) {
  p <- file.path(output_dir, paste0(case_id, "__", run, ".json"))
  if (!file.exists(p)) return(NULL)
  fromJSON(p, simplifyVector = FALSE)
}

rows <- map(case_ids, function(cid) {
  outputs <- map(runs, ~ read_run(cid, .x))
  tiers   <- map_chr(outputs, ~ .x$unified_amp_tier %||% NA_character_)
  conf    <- map_chr(outputs, ~ .x$confidence %||% NA_character_)
  disc    <- map_chr(outputs, ~ .x$discordance_nature %||% NA_character_)
  src     <- map_chr(outputs, ~ .x$primary_evidence_source %||% NA_character_)

  tier_uniq <- length(unique(tiers))
  conf_uniq <- length(unique(conf))
  disc_uniq <- length(unique(disc))
  src_uniq  <- length(unique(src))

  tibble(
    case_id           = cid,
    run1_tier         = tiers[1], run2_tier = tiers[2], run3_tier = tiers[3],
    tier_unique       = tier_uniq,
    tier_unanimous    = tier_uniq == 1,
    run1_conf         = conf[1],  run2_conf = conf[2],  run3_conf = conf[3],
    conf_unique       = conf_uniq,
    conf_unanimous    = conf_uniq == 1,
    run1_disc         = disc[1],  run2_disc = disc[2],  run3_disc = disc[3],
    disc_unique       = disc_uniq,
    disc_unanimous    = disc_uniq == 1,
    primary_src_runs  = paste(src, collapse = "/")
  )
}) |> bind_rows()

cat("\n─── Per-case consistency (3 runs) ─────────────────────────────\n")
print(rows |>
        select(case_id, run1_tier, run2_tier, run3_tier,
               tier_unanimous, conf_unanimous, disc_unanimous),
      n = 20)

n_total      <- nrow(rows)
n_tier_unan  <- sum(rows$tier_unanimous)
n_conf_unan  <- sum(rows$conf_unanimous)
n_disc_unan  <- sum(rows$disc_unanimous)

cat(sprintf("\nTier unanimous (3/3 same): %d / %d (%.1f%%)\n",
            n_tier_unan, n_total, 100 * n_tier_unan / n_total))
cat(sprintf("Confidence unanimous     : %d / %d (%.1f%%)\n",
            n_conf_unan, n_total, 100 * n_conf_unan / n_total))
cat(sprintf("Discordance-type unanimous: %d / %d (%.1f%%)\n",
            n_disc_unan, n_total, 100 * n_disc_unan / n_total))

# ── Fleiss kappa on tier across runs ─────────────────────────────────────────
# 3 raters, 6 categories. Build rating count matrix [N items × K categories].
tier_levels <- c("Tier I Level A", "Tier I Level B",
                 "Tier II Level C", "Tier II Level D",
                 "Tier III VUS", "Tier IV Benign")
counts <- map(seq_len(nrow(rows)), function(i) {
  tiers <- c(rows$run1_tier[i], rows$run2_tier[i], rows$run3_tier[i])
  vapply(tier_levels, \(t) sum(tiers == t, na.rm = TRUE), integer(1))
})
mat <- do.call(rbind, counts)
n_raters <- 3L
N <- nrow(mat)
K <- length(tier_levels)

# Per-item agreement P_i = (sum_j n_ij^2 - n) / (n*(n-1))
P_i <- apply(mat, 1, \(r) (sum(r^2) - n_raters) / (n_raters * (n_raters - 1)))
P_bar <- mean(P_i)
# Per-category proportion p_j
p_j <- colSums(mat) / (N * n_raters)
P_e <- sum(p_j^2)

fleiss_kappa <- (P_bar - P_e) / (1 - P_e)

cat(sprintf("\nFleiss kappa (tier across 3 runs, n=%d): %.3f\n", N, fleiss_kappa))
cat("  Interpretation: ≥0.81 almost perfect; 0.61-0.80 substantial; 0.41-0.60 moderate; ≤0.40 fair-poor.\n")

# ── Save ─────────────────────────────────────────────────────────────────────

out_path <- here::here("benchmarks/results/agent_consistency.csv")
utils::write.csv(rows, out_path, row.names = FALSE, na = "")
log_info("Per-case consistency: {out_path}")

summary_path <- here::here("benchmarks/results/agent_consistency_summary.csv")
utils::write.csv(
  tibble(
    metric = c("n_cases", "tier_unanimous", "tier_unanimous_pct",
               "confidence_unanimous", "discordance_unanimous",
               "fleiss_kappa_tier"),
    value = c(n_total, n_tier_unan,
              round(100 * n_tier_unan / n_total, 1),
              n_conf_unan, n_disc_unan,
              round(fleiss_kappa, 3))
  ),
  summary_path, row.names = FALSE, na = ""
)
log_info("Summary: {summary_path}")
