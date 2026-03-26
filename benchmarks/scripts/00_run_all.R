#!/usr/bin/env Rscript
# benchmarks/scripts/00_run_all.R
# Master script: run all benchmark scripts in sequence.

suppressPackageStartupMessages({
  library(logger)
  library(here)
})

log_info("========================================")
log_info("NGS Tertiary Analysis Pipeline Benchmarks")
log_info("========================================")
log_info("Start time: {Sys.time()}")

benchmark_start <- Sys.time()

scripts <- c(
  "01_download_seqc2.R",
  "02_run_seqc2_benchmark.R",
  "03_variant_accuracy.R",
  "04_clinvar_concordance.R",
  "05_evidence_concordance.R",
  "06_runtime_benchmark.R",
  "07_generate_figures.R",
  "08_generate_tables.R"
)

results <- list()

for (script in scripts) {
  script_path <- here::here("benchmarks/scripts", script)

  if (!file.exists(script_path)) {
    log_error("Script not found: {script}")
    results[[script]] <- list(status = "MISSING", time = NA_real_)
    next
  }

  log_info("")
  log_info(">>> Running: {script}")
  log_info(paste(rep("-", 60), collapse = ""))

  script_time <- system.time({
    status <- tryCatch({
      source(script_path, local = new.env(parent = globalenv()))
      "SUCCESS"
    }, error = function(e) {
      log_error("Script {script} failed: {e$message}")
      paste0("FAILED: ", e$message)
    })
  })

  results[[script]] <- list(
    status = status,
    time   = round(script_time["elapsed"], 1)
  )

  log_info("<<< {script}: {status} ({round(script_time['elapsed'], 1)}s)")
}

# --- Summary ------------------------------------------------------------------

benchmark_end <- Sys.time()
total_time <- as.numeric(difftime(benchmark_end, benchmark_start, units = "secs"))

cat("\n")
cat("========================================\n")
cat("       BENCHMARK SUMMARY\n")
cat("========================================\n")
cat(sprintf("Total time: %.1f seconds (%.1f minutes)\n\n",
            total_time, total_time / 60))

for (script in scripts) {
  r <- results[[script]]
  icon <- if (grepl("SUCCESS", r$status)) "[OK]" else "[!!]"
  time_str <- if (!is.na(r$time)) sprintf("%.1fs", r$time) else "N/A"
  cat(sprintf("  %s %-35s %s  %s\n", icon, script, time_str, r$status))
}

cat("\nAll benchmarks complete! Results in benchmarks/results/\n")
cat(sprintf("  Figures: benchmarks/results/figures/\n"))
cat(sprintf("  Tables:  benchmarks/results/tables/\n"))

log_info("=== All Benchmarks Complete ===")
