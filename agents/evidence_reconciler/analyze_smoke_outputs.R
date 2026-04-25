# agents/evidence_reconciler/analyze_smoke_outputs.R
# Phase 3 of the evidence_reconciler smoke test.
#
# Reads the 5 output JSONs produced by the operator, validates against the
# agent's JSON schema, and compares to the deterministic baseline that was
# recorded in each input's metadata block. Writes a summary CSV.
#
# Usage (from repo root):
#   Rscript agents/evidence_reconciler/analyze_smoke_outputs.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(jsonlite)
  library(logger)
})

source(here::here("R/utils.R"))
source(here::here("agents/io.R"))

setup_logging("INFO")

# Discover all cases that have both an input and an output
input_dir  <- here::here("agents/inputs/evidence_reconciler")
output_dir <- here::here("agents/outputs/evidence_reconciler")

input_files  <- list.files(input_dir,  pattern = "\\.json$", full.names = FALSE)
output_files <- list.files(output_dir, pattern = "\\.json$", full.names = FALSE)

case_ids <- tools::file_path_sans_ext(input_files)
missing_outputs <- setdiff(paste0(case_ids, ".json"), output_files)
if (length(missing_outputs) > 0) {
  log_warn("Missing {length(missing_outputs)} output file(s): {paste(missing_outputs, collapse=', ')}")
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
      baseline_tier      = baseline_tier_level,
      agent_tier         = NA_character_,
      agent_confidence   = NA_character_,
      discordance_nature = NA_character_,
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
    baseline_tier      = baseline_tier_level,
    agent_tier         = o$unified_amp_tier,
    agent_confidence   = o$confidence,
    discordance_nature = o$discordance_nature,
    primary_source     = o$primary_evidence_source,
    dissenting_noted   = o$dissenting_evidence_noted,
    valid              = TRUE,
    errors             = NA_character_
  )
})

summary_tbl <- bind_rows(rows) |>
  mutate(agrees_with_baseline = baseline_tier == agent_tier)

# ── Report ───────────────────────────────────────────────────────────────────

cat("\n─── evidence_reconciler smoke-test summary ──────────────────\n")
print(summary_tbl)

n_valid <- sum(summary_tbl$valid, na.rm = TRUE)
n_total <- nrow(summary_tbl)
n_agree <- sum(summary_tbl$agrees_with_baseline, na.rm = TRUE)

cat(sprintf("\nSchema-valid outputs : %d / %d\n", n_valid, n_total))
cat(sprintf("Agree with baseline  : %d / %d (expected: concordant case agrees; others test discordance)\n",
            n_agree, n_total))

# Breakdown by discordance_nature
cat("\nDiscordance distribution (agent):\n")
print(summary_tbl |>
      filter(valid) |>
      count(discordance_nature))

out_path <- here::here("agents/logs/smoke_test_evidence_reconciler.csv")
dir_create(dirname(out_path))
utils::write.csv(summary_tbl, out_path, row.names = FALSE, na = "")
log_info("Summary written to {out_path}")
