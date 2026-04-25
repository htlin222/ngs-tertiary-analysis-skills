#!/usr/bin/env Rscript
# benchmarks/agent/02_analyze_real_slice.R
#
# Phase 3 of Benchmark A real-data slice for the evidence_reconciler agent.
# Reads agent outputs for cases with case_id starting "real_", validates against
# the JSON schema, compares each agent call to the deterministic baseline that
# was recorded in the input metadata, and writes a summary CSV plus prints
# concordance/discordance breakdown.
#
# Usage: Rscript benchmarks/agent/02_analyze_real_slice.R

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

input_files <- list.files(input_dir, pattern = "^real_.*\\.json$", full.names = FALSE)
case_ids    <- tools::file_path_sans_ext(input_files)

if (length(case_ids) == 0) {
  stop("No real-slice inputs found under agents/inputs/evidence_reconciler/")
}

missing_outputs <- setdiff(
  paste0(case_ids, ".json"),
  list.files(output_dir, pattern = "^real_.*\\.json$", full.names = FALSE)
)
if (length(missing_outputs) > 0) {
  log_warn("Missing {length(missing_outputs)} output(s): {paste(head(missing_outputs, 5), collapse=', ')}{if (length(missing_outputs) > 5) ' ...' else ''}")
}

# ── Build comparison table ───────────────────────────────────────────────────

rows <- map(case_ids, function(cid) {
  input  <- fromJSON(file.path(input_dir, paste0(cid, ".json")), simplifyVector = FALSE)
  output_result <- read_agent_output("evidence_reconciler", cid)

  baseline_tier_level <- paste(input$metadata$baseline_tier,
                               input$metadata$baseline_level)

  if (!output_result$valid) {
    return(tibble(
      case_id            = cid,
      gene               = input$metadata$gene,
      variant            = input$metadata$variant,
      tumor_type         = input$metadata$tumor_type,
      oncokb_level       = input$metadata$oncokb_level,
      civic_amp_level    = input$metadata$civic_amp_level,
      discordance_signal = input$metadata$discordance_signal,
      baseline_tier      = baseline_tier_level,
      agent_tier         = NA_character_,
      agent_confidence   = NA_character_,
      agent_discordance_nature = NA_character_,
      primary_source     = NA_character_,
      dissenting_noted   = NA,
      valid              = FALSE,
      errors             = paste(output_result$errors, collapse = "; ")
    ))
  }

  o <- output_result$output
  tibble(
    case_id            = cid,
    gene               = input$metadata$gene,
    variant            = input$metadata$variant,
    tumor_type         = input$metadata$tumor_type,
    oncokb_level       = input$metadata$oncokb_level,
    civic_amp_level    = input$metadata$civic_amp_level,
    discordance_signal = input$metadata$discordance_signal,
    baseline_tier      = baseline_tier_level,
    agent_tier         = o$unified_amp_tier,
    agent_confidence   = o$confidence,
    agent_discordance_nature = o$discordance_nature,
    primary_source     = o$primary_evidence_source,
    dissenting_noted   = o$dissenting_evidence_noted,
    valid              = TRUE,
    errors             = NA_character_
  )
})

summary_tbl <- bind_rows(rows) |>
  mutate(agrees_with_baseline = baseline_tier == agent_tier)

cat("\n─── evidence_reconciler real-slice summary ──────────────────\n")
print(summary_tbl |>
        select(case_id, baseline_tier, agent_tier, agent_confidence,
               agent_discordance_nature, agrees_with_baseline),
      n = 50)

n_total <- nrow(summary_tbl)
n_valid <- sum(summary_tbl$valid, na.rm = TRUE)
n_agree <- sum(summary_tbl$agrees_with_baseline, na.rm = TRUE)

cat(sprintf("\nSchema-valid outputs : %d / %d\n", n_valid, n_total))
cat(sprintf("Agree with baseline  : %d / %d (%.1f%%)\n",
            n_agree, n_total, 100 * n_agree / n_total))

cat("\nAgent confidence distribution:\n")
print(summary_tbl |> filter(valid) |> count(agent_confidence))

cat("\nAgent-detected discordance distribution:\n")
print(summary_tbl |> filter(valid) |> count(agent_discordance_nature))

cat("\nDisagreements (where agent overruled baseline):\n")
disagreements <- summary_tbl |> filter(valid, !agrees_with_baseline)
if (nrow(disagreements) > 0) {
  print(disagreements |>
          select(case_id, baseline_tier, agent_tier, agent_confidence,
                 agent_discordance_nature, primary_source),
        n = 50)
} else {
  cat("  (none — agent always agreed with baseline on this slice)\n")
}

# ── Save ─────────────────────────────────────────────────────────────────────

out_path <- here::here("benchmarks/results/agent_evidence_reconciler_real_slice.csv")
dir_create(dirname(out_path))
utils::write.csv(summary_tbl, out_path, row.names = FALSE, na = "")
log_info("Summary CSV written: {out_path}")
