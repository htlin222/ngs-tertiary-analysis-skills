#!/usr/bin/env Rscript
# benchmarks/scripts/02_run_seqc2_benchmark.R
# Run the pipeline on SEQC2 HCC1395 data and record stage timings.

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(here)
  library(dplyr)
  library(tibble)
})

source(here::here("R/utils.R"))

log_info("=== SEQC2 Pipeline Benchmark ===")

# --- Configuration -----------------------------------------------------------

sample_id  <- "SEQC2_HCC1395"
tumor_type <- "Breast Cancer"
results_dir <- here::here("benchmarks/results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# --- Check for SEQC2 BAM -----------------------------------------------------

# The full BAM is expected from testdata/download_testdata.sh --full
bam_candidates <- c(

  here::here("testdata", "HCC1395.bam"),
  here::here("testdata", "SEQC2_HCC1395.bam"),
  Sys.getenv("SEQC2_BAM_PATH", unset = "")
)

bam_path <- NULL
for (candidate in bam_candidates) {
  if (nchar(candidate) > 0 && file.exists(candidate)) {
    bam_path <- candidate
    break
  }
}

if (is.null(bam_path)) {
  log_warn("SEQC2 BAM not found in testdata/")
  log_warn("Expected locations: {paste(bam_candidates[1:2], collapse = ', ')}")
  log_info("You can download it with: bash testdata/download_testdata.sh --full")
  log_info("Or set SEQC2_BAM_PATH environment variable")

  # Offer to download
  if (interactive()) {
    answer <- readline("Download SEQC2 BAM now? (y/N): ")
    if (tolower(answer) == "y") {
      download_script <- here::here("testdata/download_testdata.sh")
      if (file.exists(download_script)) {
        log_info("Running download script...")
        system2("bash", args = c(download_script, "--full"))
        # Re-check
        for (candidate in bam_candidates[1:2]) {
          if (file.exists(candidate)) {
            bam_path <- candidate
            break
          }
        }
      } else {
        log_error("Download script not found: {download_script}")
      }
    }
  }

  if (is.null(bam_path)) {
    stop("SEQC2 BAM file not available. Cannot proceed with benchmark.")
  }
}

log_info("Using BAM: {bam_path}")

# --- Set environment for pipeline run -----------------------------------------

Sys.setenv(
  BAM_PATH  = bam_path,
  SAMPLE_ID = sample_id
)

# --- Define pipeline stages to time ------------------------------------------

# We time each stage by running the targets pipeline and extracting
# timestamps from tar_meta() after completion.

pipeline_stages <- c(
  "00-qc",
  "01-variant-calling",
  "02-annotation",
  "03-cnv",
  "04-fusions",
  "05-biomarkers",
  "06-clinical-annotation",
  "07-literature",
  "08-report"
)

# --- Run pipeline with timing ------------------------------------------------

log_info("Starting pipeline run for {sample_id} ({tumor_type})...")

# Load pipeline config with overrides for SEQC2
config_overrides <- list(
  "sample.id"         = sample_id,
  "sample.tumor_type" = tumor_type
)

# Time the entire pipeline
total_time <- system.time({
  tryCatch({
    # Use targets if _targets.R exists
    if (requireNamespace("targets", quietly = TRUE) &&
        file.exists(here::here("_targets.R"))) {
      log_info("Running targets pipeline...")
      targets::tar_make()
    } else {
      log_warn("targets package or _targets.R not found")
      log_info("Running pipeline stages manually...")

      # Fallback: source stage scripts sequentially
      stage_dirs <- list.dirs(here::here(), recursive = FALSE)
      stage_dirs <- sort(stage_dirs[grepl("^\\d{2}-", basename(stage_dirs))])
      for (stage_dir in stage_dirs) {
        run_script <- file.path(stage_dir, "run.R")
        if (file.exists(run_script)) {
          log_info("Running: {basename(stage_dir)}")
          source(run_script, local = new.env(parent = globalenv()))
        }
      }
    }
  }, error = function(e) {
    log_error("Pipeline failed: {e$message}")
    stop(e)
  })
})

log_info("Total pipeline time: {round(total_time['elapsed'], 1)}s")

# --- Extract per-stage timing from targets metadata ---------------------------

timing_records <- list()

if (requireNamespace("targets", quietly = TRUE)) {
  tryCatch({
    meta <- targets::tar_meta()

    if (nrow(meta) > 0 && "seconds" %in% names(meta)) {
      # Map target names to stages
      for (i in seq_len(nrow(meta))) {
        target_name <- meta$name[i]
        seconds     <- meta$seconds[i]
        if (!is.na(seconds)) {
          timing_records[[length(timing_records) + 1]] <- tibble(
            sample_id = sample_id,
            target    = target_name,
            seconds   = seconds,
            timestamp = Sys.time()
          )
        }
      }
    }
  }, error = function(e) {
    log_warn("Could not extract targets metadata: {e$message}")
  })
}

# Add total time record
timing_records[[length(timing_records) + 1]] <- tibble(
  sample_id = sample_id,
  target    = "TOTAL",
  seconds   = as.numeric(total_time["elapsed"]),
  timestamp = Sys.time()
)

# --- Save timing results ------------------------------------------------------

timing_df <- bind_rows(timing_records)
output_path <- file.path(results_dir, "seqc2_timing.csv")
write.csv(timing_df, output_path, row.names = FALSE)

log_info("Timing results saved to: {output_path}")
log_info("Stages timed: {nrow(timing_df)}")

# Print summary
cat("\n--- SEQC2 Benchmark Timing Summary ---\n")
print(timing_df |> select(target, seconds) |> arrange(desc(seconds)), n = 50)

log_info("=== SEQC2 Benchmark Complete ===")
