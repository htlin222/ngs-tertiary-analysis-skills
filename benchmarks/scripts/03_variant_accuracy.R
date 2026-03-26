#!/usr/bin/env Rscript
# benchmarks/scripts/03_variant_accuracy.R
# Compare pipeline variants against SEQC2 truth set.
# Calculate sensitivity, precision, F1 stratified by variant type and VAF.

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

source(here::here("R/utils.R"))

log_info("=== Variant Accuracy Benchmark ===")

# --- Configuration -----------------------------------------------------------

sample_id    <- "SEQC2_HCC1395"
results_dir  <- here::here("benchmarks/results")
figures_dir  <- here::here("benchmarks/results/figures")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

truth_vcf_path   <- here::here("benchmarks/datasets/seqc2/seqc2_truth_tso500_filtered.vcf.gz")
pipeline_vcf_dir <- here::here("reports", sample_id, "01-variant-calling")

# --- Helper: Read VCF into a variant tibble -----------------------------------

read_vcf_variants <- function(vcf_path) {
  if (!file.exists(vcf_path)) {
    log_error("VCF not found: {vcf_path}")
    stop(glue("VCF file not found: {vcf_path}"))
  }

  if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
    log_info("Reading VCF with VariantAnnotation: {basename(vcf_path)}")
    vcf <- VariantAnnotation::readVcf(vcf_path, genome = "GRCh38")
    rd  <- SummarizedExperiment::rowRanges(vcf)

    variants <- tibble(
      chr = as.character(GenomicRanges::seqnames(rd)),
      pos = GenomicRanges::start(rd),
      ref = as.character(rd$REF),
      alt = as.character(unlist(rd$ALT))
    )

    # Extract VAF from genotype fields if available
    geno <- VariantAnnotation::geno(vcf)
    if ("AF" %in% names(geno)) {
      af_vals <- geno$AF
      if (is.matrix(af_vals) || is.list(af_vals)) {
        variants$vaf <- as.numeric(sapply(
          if (is.matrix(af_vals)) af_vals[, 1] else af_vals,
          function(x) x[1]
        ))
      }
    } else if ("AD" %in% names(geno)) {
      # Calculate VAF from allelic depths
      ad_vals <- geno$AD
      variants$vaf <- sapply(seq_len(nrow(variants)), function(i) {
        ad <- if (is.matrix(ad_vals)) ad_vals[i, 1] else ad_vals[[i]]
        if (length(ad) >= 2 && sum(ad) > 0) ad[2] / sum(ad) else NA_real_
      })
    } else {
      variants$vaf <- NA_real_
    }

    # Classify variant type
    variants <- variants |>
      mutate(
        var_type = case_when(
          nchar(ref) == 1 & nchar(alt) == 1 ~ "SNV",
          nchar(ref) > nchar(alt) ~ "DEL",
          nchar(ref) < nchar(alt) ~ "INS",
          TRUE ~ "COMPLEX"
        ),
        indel = var_type %in% c("DEL", "INS", "COMPLEX"),
        key = paste(chr, pos, ref, alt, sep = ":")
      )

    return(variants)
  }

  # Fallback: parse VCF with base R
  log_info("Reading VCF with base R parser: {basename(vcf_path)}")
  if (grepl("\\.gz$", vcf_path)) {
    con <- gzfile(vcf_path, "r")
  } else {
    con <- file(vcf_path, "r")
  }
  on.exit(close(con))

  lines <- readLines(con)
  data_lines <- lines[!grepl("^#", lines)]

  if (length(data_lines) == 0) {
    return(tibble(chr = character(), pos = integer(), ref = character(),
                  alt = character(), vaf = double(), var_type = character(),
                  indel = logical(), key = character()))
  }

  fields <- strsplit(data_lines, "\t")
  variants <- tibble(
    chr = sapply(fields, `[`, 1),
    pos = as.integer(sapply(fields, `[`, 2)),
    ref = sapply(fields, `[`, 4),
    alt = sapply(fields, `[`, 5),
    vaf = NA_real_
  ) |>
    mutate(
      var_type = case_when(
        nchar(ref) == 1 & nchar(alt) == 1 ~ "SNV",
        nchar(ref) > nchar(alt) ~ "DEL",
        nchar(ref) < nchar(alt) ~ "INS",
        TRUE ~ "COMPLEX"
      ),
      indel = var_type %in% c("DEL", "INS", "COMPLEX"),
      key = paste(chr, pos, ref, alt, sep = ":")
    )

  variants
}

# --- Read truth VCF -----------------------------------------------------------

if (!file.exists(truth_vcf_path)) {
  # Fall back to unfiltered truth
  truth_vcf_path <- here::here(
    "benchmarks/datasets/seqc2",
    "high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz"
  )
}

log_info("Reading truth VCF: {truth_vcf_path}")
truth <- read_vcf_variants(truth_vcf_path)
log_info("Truth set: {nrow(truth)} variants ({sum(truth$var_type == 'SNV')} SNVs, {sum(truth$indel)} indels)")

# --- Read pipeline output VCF -------------------------------------------------

pipeline_vcf_candidates <- c(
  list.files(pipeline_vcf_dir, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE),
  here::here("reports", sample_id, "02-annotation",
             list.files(here::here("reports", sample_id, "02-annotation"),
                        pattern = "\\.vcf(\\.gz)?$"))
)

if (length(pipeline_vcf_candidates) == 0) {
  log_error("No pipeline output VCF found in: {pipeline_vcf_dir}")
  stop("Run 02_run_seqc2_benchmark.R first to generate pipeline output")
}

pipeline_vcf_path <- pipeline_vcf_candidates[1]
log_info("Reading pipeline VCF: {pipeline_vcf_path}")
pipeline <- read_vcf_variants(pipeline_vcf_path)
log_info("Pipeline output: {nrow(pipeline)} variants")

# --- Match variants by chr:pos:ref:alt ----------------------------------------

log_info("Matching variants by chr:pos:ref:alt...")

tp_keys <- intersect(truth$key, pipeline$key)
fp_keys <- setdiff(pipeline$key, truth$key)
fn_keys <- setdiff(truth$key, pipeline$key)

tp <- length(tp_keys)
fp <- length(fp_keys)
fn <- length(fn_keys)

sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
f1_score    <- if (!is.na(sensitivity) && !is.na(precision) && (sensitivity + precision) > 0) {
  2 * sensitivity * precision / (sensitivity + precision)
} else {
  NA_real_
}

log_info("Overall: TP={tp}, FP={fp}, FN={fn}")
log_info("Sensitivity={round(sensitivity, 4)}, Precision={round(precision, 4)}, F1={round(f1_score, 4)}")

# --- Stratify by variant type -------------------------------------------------

calc_metrics <- function(truth_subset, pipeline_subset) {
  tp <- length(intersect(truth_subset$key, pipeline_subset$key))
  fp <- length(setdiff(pipeline_subset$key, truth_subset$key))
  fn <- length(setdiff(truth_subset$key, pipeline_subset$key))
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1   <- if (!is.na(sens) && !is.na(prec) && (sens + prec) > 0) {
    2 * sens * prec / (sens + prec)
  } else {
    NA_real_
  }
  tibble(TP = tp, FP = fp, FN = fn,
         sensitivity = sens, precision = prec, F1 = f1)
}

type_results <- bind_rows(
  calc_metrics(
    truth |> filter(var_type == "SNV"),
    pipeline |> filter(var_type == "SNV")
  ) |> mutate(stratum = "SNV", .before = 1),
  calc_metrics(
    truth |> filter(indel),
    pipeline |> filter(indel)
  ) |> mutate(stratum = "Indel", .before = 1),
  calc_metrics(truth, pipeline) |>
    mutate(stratum = "All", .before = 1)
)

# --- Stratify by VAF bucket ---------------------------------------------------

# Use pipeline VAFs for VAF stratification of true positives
pipeline_tp <- pipeline |> filter(key %in% tp_keys)

vaf_buckets <- c(0, 0.05, 0.10, 0.20, 1.0)
vaf_labels  <- c("<5%", "5-10%", "10-20%", ">20%")

if (any(!is.na(pipeline_tp$vaf))) {
  pipeline_tp <- pipeline_tp |>
    mutate(vaf_bucket = cut(vaf, breaks = vaf_buckets, labels = vaf_labels,
                            include.lowest = TRUE, right = FALSE))

  # For VAF stratification, compute sensitivity at each VAF threshold
  vaf_thresholds <- seq(0.01, 0.50, by = 0.01)
  vaf_sensitivity <- map_dfr(vaf_thresholds, function(thresh) {
    detected <- pipeline |> filter(!is.na(vaf) & vaf >= thresh)
    detected_tp <- sum(detected$key %in% truth$key)
    total_truth <- nrow(truth)
    tibble(
      vaf_threshold = thresh,
      sensitivity = if (total_truth > 0) detected_tp / total_truth else NA_real_,
      n_detected = nrow(detected),
      n_tp = detected_tp
    )
  })

  # VAF bucket summary
  vaf_bucket_results <- pipeline_tp |>
    filter(!is.na(vaf_bucket)) |>
    group_by(vaf_bucket) |>
    summarise(n_tp = n(), .groups = "drop") |>
    mutate(stratum = paste0("VAF_", vaf_bucket))

  type_results <- bind_rows(type_results, vaf_bucket_results |>
                              select(stratum, TP = n_tp))
} else {
  vaf_sensitivity <- tibble(vaf_threshold = double(), sensitivity = double(),
                            n_detected = integer(), n_tp = integer())
  log_warn("No VAF data available in pipeline VCF; skipping VAF stratification")
}

# --- Save results -------------------------------------------------------------

output_path <- file.path(results_dir, "variant_accuracy.csv")
write.csv(type_results, output_path, row.names = FALSE)
log_info("Accuracy results saved to: {output_path}")

if (nrow(vaf_sensitivity) > 0) {
  vaf_output <- file.path(results_dir, "variant_accuracy_by_vaf.csv")
  write.csv(vaf_sensitivity, vaf_output, row.names = FALSE)
  log_info("VAF sensitivity saved to: {vaf_output}")
}

# --- Generate accuracy plot ---------------------------------------------------

if (nrow(vaf_sensitivity) > 0) {
  p <- ggplot(vaf_sensitivity, aes(x = vaf_threshold, y = sensitivity)) +
    geom_line(color = "#2c3e50", linewidth = 1) +
    geom_point(color = "#e74c3c", size = 1.5) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.45, y = 0.955, label = "95% sensitivity",
             color = "grey50", size = 3) +
    scale_x_continuous(labels = scales::percent_format(),
                       breaks = seq(0, 0.5, 0.05)) +
    scale_y_continuous(labels = scales::percent_format(),
                       limits = c(0, 1)) +
    labs(
      title = "Variant Detection Sensitivity vs VAF Threshold",
      subtitle = "SEQC2 HCC1395 truth set, TSO500 panel regions",
      x = "VAF Threshold",
      y = "Sensitivity"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(figures_dir, "sensitivity_vs_vaf.pdf"), p,
         width = 8, height = 5, dpi = 300)
  ggsave(file.path(figures_dir, "sensitivity_vs_vaf.png"), p,
         width = 8, height = 5, dpi = 300)
  log_info("Sensitivity vs VAF plot saved")
}

# --- Print summary ------------------------------------------------------------

cat("\n--- Variant Accuracy Summary ---\n")
print(type_results, n = 20)

log_info("=== Variant Accuracy Benchmark Complete ===")
