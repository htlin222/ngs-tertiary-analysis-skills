#!/usr/bin/env Rscript
# benchmarks/agent/06_analyze_cross_model.R
#
# Cross-model robustness check: compare evidence_reconciler outputs across
# Claude Opus 4.7 (3 runs from Phase 3) vs Sonnet 4.6 (1 run) vs Haiku 4.5 (1 run)
# on the same 10-case consistency manifest.
#
# Reports:
#   - per-case tier across 5 model-runs
#   - cross-model unanimity (Opus run1 vs Sonnet vs Haiku)
#   - Fleiss kappa across the 5 raters

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

manifest <- readLines("/tmp/batches/manifest_consistency.txt") |> trimws() |>
  Filter(\(x) nchar(x) > 0, x = _)
output_dir <- here::here("agents/outputs/evidence_reconciler")
runs <- c("run1" = "Opus 4.7 (run 1)",
          "run2" = "Opus 4.7 (run 2)",
          "run3" = "Opus 4.7 (run 3)",
          "sonnet" = "Sonnet 4.6",
          "haiku"  = "Haiku 4.5")

read_run <- function(case_id, run) {
  p <- file.path(output_dir, paste0(case_id, "__", run, ".json"))
  if (!file.exists(p)) return(NULL)
  fromJSON(p, simplifyVector = FALSE)
}

rows <- map(manifest, function(cid) {
  outs <- map(names(runs), \(r) read_run(cid, r))
  tiers <- map_chr(outs, ~ .x$unified_amp_tier %||% NA_character_)
  conf  <- map_chr(outs, ~ .x$confidence %||% NA_character_)
  disc  <- map_chr(outs, ~ .x$discordance_nature %||% NA_character_)
  tibble(
    case_id   = cid,
    opus_r1_tier = tiers[1],
    opus_r2_tier = tiers[2],
    opus_r3_tier = tiers[3],
    sonnet_tier  = tiers[4],
    haiku_tier   = tiers[5],
    cross_model_unanimous = length(unique(c(tiers[1], tiers[4], tiers[5]))) == 1,
    all_5_runs_unanimous  = length(unique(tiers)) == 1
  )
}) |> bind_rows()

cat("\n─── Per-case tier across model-runs ───────────────────────────\n")
print(rows, n = 20)

cat(sprintf("\nUnanimous Opus-r1 / Sonnet / Haiku: %d / %d (%.1f%%)\n",
            sum(rows$cross_model_unanimous), nrow(rows),
            100 * mean(rows$cross_model_unanimous)))
cat(sprintf("Unanimous across all 5 model-runs : %d / %d (%.1f%%)\n",
            sum(rows$all_5_runs_unanimous), nrow(rows),
            100 * mean(rows$all_5_runs_unanimous)))

# ── Pairwise tier-agreement (Opus run 1 vs each other) ───────────────────────

pairwise <- tibble(
  comparison = c("Opus r1 vs Opus r2",
                 "Opus r1 vs Opus r3",
                 "Opus r1 vs Sonnet 4.6",
                 "Opus r1 vs Haiku 4.5"),
  tier_agreement = c(
    mean(rows$opus_r1_tier == rows$opus_r2_tier),
    mean(rows$opus_r1_tier == rows$opus_r3_tier),
    mean(rows$opus_r1_tier == rows$sonnet_tier),
    mean(rows$opus_r1_tier == rows$haiku_tier)
  ),
  cohen_kappa = c(
    psych_kappa <- function(a, b) {
      if (length(unique(c(a, b))) == 1) return(1)
      tab <- table(factor(a, levels = unique(c(a, b))),
                   factor(b, levels = unique(c(a, b))))
      n  <- sum(tab)
      po <- sum(diag(tab)) / n
      pe <- sum(rowSums(tab) * colSums(tab)) / (n * n)
      (po - pe) / (1 - pe)
    },
    psych_kappa(rows$opus_r1_tier, rows$opus_r2_tier),
    psych_kappa(rows$opus_r1_tier, rows$opus_r3_tier),
    psych_kappa(rows$opus_r1_tier, rows$sonnet_tier),
    psych_kappa(rows$opus_r1_tier, rows$haiku_tier)
  )[-1]
)
pairwise$tier_agreement <- round(pairwise$tier_agreement, 3)
pairwise$cohen_kappa    <- round(as.numeric(pairwise$cohen_kappa), 3)

cat("\nPairwise (vs Opus run 1):\n")
print(pairwise)

# ── Fleiss kappa across all 5 raters ─────────────────────────────────────────

tier_levels <- c("Tier I Level A", "Tier I Level B",
                 "Tier II Level C", "Tier II Level D",
                 "Tier III VUS", "Tier IV Benign")
counts <- map(seq_len(nrow(rows)), function(i) {
  ts <- c(rows$opus_r1_tier[i], rows$opus_r2_tier[i], rows$opus_r3_tier[i],
          rows$sonnet_tier[i], rows$haiku_tier[i])
  vapply(tier_levels, \(t) sum(ts == t, na.rm = TRUE), integer(1))
})
mat <- do.call(rbind, counts)
n_raters <- 5L; N <- nrow(mat); K <- length(tier_levels)
P_i <- apply(mat, 1, \(r) (sum(r^2) - n_raters) / (n_raters * (n_raters - 1)))
P_bar <- mean(P_i)
p_j <- colSums(mat) / (N * n_raters)
P_e <- sum(p_j^2)
fleiss_5 <- (P_bar - P_e) / (1 - P_e)

cat(sprintf("\nFleiss kappa across 5 raters (3xOpus + Sonnet + Haiku, n=%d): %.3f\n",
            N, fleiss_5))

# ── Save ─────────────────────────────────────────────────────────────────────

out_path <- here::here("benchmarks/results/agent_cross_model.csv")
utils::write.csv(rows, out_path, row.names = FALSE, na = "")
log_info("Per-case cross-model: {out_path}")

pw_path <- here::here("benchmarks/results/agent_cross_model_pairwise.csv")
utils::write.csv(pairwise, pw_path, row.names = FALSE, na = "")
log_info("Pairwise: {pw_path}")

sm_path <- here::here("benchmarks/results/agent_cross_model_summary.csv")
utils::write.csv(
  tibble(
    metric = c("n_cases", "cross_model_unanimous_opus_sonnet_haiku",
               "all_5_runs_unanimous", "fleiss_kappa_5_raters"),
    value = c(nrow(rows), sum(rows$cross_model_unanimous),
              sum(rows$all_5_runs_unanimous), round(fleiss_5, 3))
  ),
  sm_path, row.names = FALSE, na = ""
)
log_info("Summary: {sm_path}")
