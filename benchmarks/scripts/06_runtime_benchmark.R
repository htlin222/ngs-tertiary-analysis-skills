#!/usr/bin/env Rscript
# benchmarks/scripts/06_runtime_benchmark.R
# Time pipeline stages on synthetic VCFs of varying sizes.

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
})

source(here::here("R/utils.R"))

log_info("=== Runtime Benchmark ===")

# --- Configuration -----------------------------------------------------------

results_dir <- here::here("benchmarks/results")
synthetic_dir <- here::here("benchmarks/datasets/synthetic")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(synthetic_dir)) dir.create(synthetic_dir, recursive = TRUE)

variant_counts <- c(10, 50, 100, 500)
panel_bed <- here::here("config/tso500_panel.bed")

# --- Read panel BED for coordinate sampling -----------------------------------

read_panel_regions <- function(bed_path) {
  if (!file.exists(bed_path)) {
    log_error("Panel BED not found: {bed_path}")
    stop("Missing panel BED file")
  }
  bed <- read.table(bed_path, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE,
                    col.names = c("chr", "start", "end", "name")[
                      seq_len(min(4, ncol(read.table(bed_path, nrows = 1))))
                    ])
  if (ncol(bed) < 3) stop("BED file must have at least 3 columns")
  names(bed)[1:3] <- c("chr", "start", "end")
  bed
}

panel_regions <- read_panel_regions(panel_bed)
log_info("Loaded {nrow(panel_regions)} panel regions")

# --- Generate synthetic VCFs --------------------------------------------------

generate_synthetic_vcf <- function(n_variants, output_path, regions) {
  set.seed(42 + n_variants)  # reproducible per size

  # Sample random positions from panel regions
  region_lengths <- regions$end - regions$start
  region_probs <- region_lengths / sum(region_lengths)

  sampled_regions <- sample(seq_len(nrow(regions)), n_variants,
                            replace = TRUE, prob = region_probs)

  bases <- c("A", "C", "G", "T")

  variants <- tibble(
    CHROM = regions$chr[sampled_regions],
    POS   = mapply(function(s, e) sample(s:e, 1),
                   regions$start[sampled_regions],
                   regions$end[sampled_regions]),
    ID    = paste0("synth_", seq_len(n_variants)),
    REF   = sample(bases, n_variants, replace = TRUE),
    ALT   = NA_character_,
    QUAL  = round(runif(n_variants, 30, 200), 1),
    FILTER = "PASS",
    INFO  = paste0(
      "DP=", sample(50:500, n_variants, replace = TRUE),
      ";AF=", round(runif(n_variants, 0.05, 0.95), 3)
    )
  )

  # Ensure ALT differs from REF
  variants$ALT <- sapply(variants$REF, function(ref) {
    alts <- setdiff(bases, ref)
    sample(alts, 1)
  })

  # Add some indels (~20%)
  n_indels <- round(n_variants * 0.2)
  if (n_indels > 0) {
    indel_idx <- sample(seq_len(n_variants), n_indels)
    for (idx in indel_idx) {
      if (runif(1) > 0.5) {
        # Insertion
        variants$ALT[idx] <- paste0(variants$REF[idx],
                                     paste0(sample(bases, sample(1:5, 1), replace = TRUE),
                                            collapse = ""))
      } else {
        # Deletion
        del_len <- sample(1:5, 1)
        variants$REF[idx] <- paste0(variants$REF[idx],
                                     paste0(sample(bases, del_len, replace = TRUE),
                                            collapse = ""))
      }
    }
  }

  # Sort by chromosome and position
  variants <- variants |> arrange(CHROM, POS)

  # Write VCF
  header <- c(
    "##fileformat=VCFv4.2",
    glue("##source=ngs-tertiary-benchmark-synthetic-{n_variants}"),
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depth\">",
    "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC"
  )

  # Create sample genotype column
  gt_col <- sapply(seq_len(n_variants), function(i) {
    af <- as.numeric(str_extract(variants$INFO[i], "AF=([0-9.]+)", group = 1))
    dp <- as.integer(str_extract(variants$INFO[i], "DP=([0-9]+)", group = 1))
    alt_reads <- round(dp * af)
    ref_reads <- dp - alt_reads
    paste0("0/1:", ref_reads, ",", alt_reads, ":", round(af, 3))
  })

  data_lines <- paste(
    variants$CHROM, variants$POS, variants$ID, variants$REF, variants$ALT,
    variants$QUAL, variants$FILTER, variants$INFO,
    "GT:AD:AF", gt_col,
    sep = "\t"
  )

  writeLines(c(header, data_lines), output_path)
  log_info("Generated synthetic VCF: {basename(output_path)} ({n_variants} variants)")
  output_path
}

# Generate VCFs for each size
synthetic_vcfs <- map_chr(variant_counts, function(n) {
  vcf_path <- file.path(synthetic_dir, glue("synthetic_{n}_variants.vcf"))
  generate_synthetic_vcf(n, vcf_path, panel_regions)
  vcf_path
})

names(synthetic_vcfs) <- as.character(variant_counts)

# --- Run pipeline on each synthetic VCF and time stages -----------------------

log_info("Running pipeline on synthetic VCFs...")

has_targets <- requireNamespace("targets", quietly = TRUE) &&
  file.exists(here::here("_targets.R"))

timing_results <- list()

for (i in seq_along(variant_counts)) {
  n <- variant_counts[i]
  vcf_path <- synthetic_vcfs[i]
  sample_id <- glue("BENCHMARK_{n}_variants")

  log_info("--- Running pipeline for {n} variants ---")

  # Set environment for VCF input mode
  Sys.setenv(
    VCF_PATH  = vcf_path,
    SAMPLE_ID = sample_id,
    INPUT_MODE = "vcf"
  )

  if (has_targets) {
    # Time via targets
    total_elapsed <- system.time({
      tryCatch({
        targets::tar_make()
      }, error = function(e) {
        log_warn("Pipeline failed for {n} variants: {e$message}")
      })
    })["elapsed"]

    # Extract per-stage timing from targets metadata
    tryCatch({
      meta <- targets::tar_meta()
      if (nrow(meta) > 0 && "seconds" %in% names(meta)) {
        for (j in seq_len(nrow(meta))) {
          if (!is.na(meta$seconds[j])) {
            timing_results[[length(timing_results) + 1]] <- tibble(
              n_variants = n,
              stage      = meta$name[j],
              seconds    = meta$seconds[j]
            )
          }
        }
      }
    }, error = function(e) {
      log_warn("Could not extract targets metadata: {e$message}")
    })

    timing_results[[length(timing_results) + 1]] <- tibble(
      n_variants = n,
      stage      = "TOTAL",
      seconds    = as.numeric(total_elapsed)
    )
  } else {
    # Fallback: run stage scripts and time each
    stage_dirs <- list.dirs(here::here(), recursive = FALSE)
    stage_dirs <- sort(stage_dirs[grepl("^\\d{2}-", basename(stage_dirs))])

    for (stage_dir in stage_dirs) {
      run_script <- file.path(stage_dir, "run.R")
      if (file.exists(run_script)) {
        stage_name <- basename(stage_dir)
        stage_time <- system.time({
          tryCatch({
            source(run_script, local = new.env(parent = globalenv()))
          }, error = function(e) {
            log_warn("Stage {stage_name} failed: {e$message}")
          })
        })
        timing_results[[length(timing_results) + 1]] <- tibble(
          n_variants = n,
          stage      = stage_name,
          seconds    = as.numeric(stage_time["elapsed"])
        )
      }
    }
  }

  elapsed_val <- if (exists("total_elapsed")) round(as.numeric(total_elapsed), 1) else NA
  log_info("Completed pipeline for {n} variants in {elapsed_val}s")
}

# --- Save results -------------------------------------------------------------

timing_df <- bind_rows(timing_results)
output_path <- file.path(results_dir, "runtime_benchmark.csv")
write.csv(timing_df, output_path, row.names = FALSE)

log_info("Runtime benchmark saved to: {output_path}")

cat("\n--- Runtime Benchmark Summary ---\n")
timing_df |>
  filter(stage == "TOTAL") |>
  select(n_variants, seconds) |>
  print()

cat("\n--- Per-Stage Timing ---\n")
timing_df |>
  filter(stage != "TOTAL") |>
  tidyr::pivot_wider(names_from = n_variants, values_from = seconds,
                     names_prefix = "n_") |>
  print(n = 30)

log_info("=== Runtime Benchmark Complete ===")
