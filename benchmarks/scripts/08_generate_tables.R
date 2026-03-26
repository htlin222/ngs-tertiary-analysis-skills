#!/usr/bin/env Rscript
# benchmarks/scripts/08_generate_tables.R
# Generate publication tables as CSV and gt HTML.

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

source(here::here("R/utils.R"))

log_info("=== Generate Publication Tables ===")

# --- Configuration -----------------------------------------------------------

results_dir <- here::here("benchmarks/results")
tables_dir  <- here::here("benchmarks/results/tables")
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

has_gt <- requireNamespace("gt", quietly = TRUE)
if (has_gt) {
  library(gt)
} else {
  log_warn("gt package not available; HTML tables will be skipped")
}

save_table <- function(df, name, title = NULL, subtitle = NULL) {
  # CSV
  csv_path <- file.path(tables_dir, paste0(name, ".csv"))
  write.csv(df, csv_path, row.names = FALSE)
  log_info("Saved table CSV: {name}.csv")

  # gt HTML
  if (has_gt) {
    gt_tbl <- gt(df)
    if (!is.null(title)) {
      gt_tbl <- gt_tbl |> tab_header(title = title, subtitle = subtitle)
    }
    gt_tbl <- gt_tbl |>
      tab_options(
        table.font.size = px(12),
        heading.title.font.size = px(16),
        heading.subtitle.font.size = px(13),
        column_labels.font.weight = "bold"
      )
    html_path <- file.path(tables_dir, paste0(name, ".html"))
    gtsave(gt_tbl, html_path)
    log_info("Saved table HTML: {name}.html")
  }
}

# ==============================================================================
# Table 1: Pipeline Stages Summary
# ==============================================================================

log_info("Generating Table 1: Pipeline stages...")

table1 <- tribble(
  ~Stage, ~`Stage Name`, ~Tool, ~Input, ~Output,
  "00", "Quality Control",
  "Rsamtools, samtools",
  "BAM file",
  "QC metrics, coverage report",

  "01", "Variant Calling",
  "GATK Mutect2",
  "BAM, reference FASTA, panel BED",
  "Somatic VCF (SNVs + indels)",

  "02", "Annotation",
  "VEP, ClinVar, COSMIC",
  "Raw VCF",
  "Annotated VCF with functional impact",

  "03", "Copy Number Variation",
  "CNVkit (hybrid mode)",
  "BAM, reference CNN",
  "CNV segments, gene-level calls",

  "04", "Fusion Detection",
  "Manta",
  "BAM",
  "Structural variant VCF, fusion candidates",

  "05", "Biomarkers",
  "Custom R (TMB, MSI, HRD)",
  "Annotated VCF, MSI sites BED",
  "TMB score, MSI status, HRD score",

  "06", "Clinical Annotation",
  "OncoKB API, CiVIC API",
  "Annotated variants, tumor type",
  "AMP/ASCO/CAP tiered variants, therapies",

  "07", "Literature Search",
  "PubMed, Scopus, Unpaywall APIs",
  "Classified variants",
  "Relevant publications with open-access links",

  "08", "Report Generation",
  "Quarto, gt, ggplot2",
  "All upstream results",
  "ESMO 2024-compliant HTML clinical report"
)

save_table(table1, "table1_pipeline_stages",
           title = "Table 1. Pipeline Stages and Components",
           subtitle = "NGS Tertiary Analysis Pipeline architecture")

# ==============================================================================
# Table 2: Variant Detection Performance
# ==============================================================================

log_info("Generating Table 2: Variant detection performance...")

accuracy_path <- file.path(results_dir, "variant_accuracy.csv")
if (file.exists(accuracy_path)) {
  accuracy <- read.csv(accuracy_path)

  table2 <- accuracy |>
    filter(stratum %in% c("All", "SNV", "Indel")) |>
    select(
      `Variant Type` = stratum,
      TP, FP, FN,
      Sensitivity = sensitivity,
      Precision = precision,
      `F1 Score` = F1
    ) |>
    mutate(
      across(c(Sensitivity, Precision, `F1 Score`), ~ round(., 4))
    )

  save_table(table2, "table2_variant_performance",
             title = "Table 2. Variant Detection Performance",
             subtitle = "SEQC2 HCC1395 truth set, TSO500 panel regions")
} else {
  log_warn("Skipping Table 2: variant_accuracy.csv not found")
  table2 <- tibble(Note = "Run 03_variant_accuracy.R to generate this table")
  save_table(table2, "table2_variant_performance",
             title = "Table 2. Variant Detection Performance (placeholder)")
}

# ==============================================================================
# Table 3: AMP Classification Concordance Matrix
# ==============================================================================

log_info("Generating Table 3: AMP concordance matrix...")

concordance_path <- file.path(results_dir, "clinvar_concordance.csv")
if (file.exists(concordance_path)) {
  concordance <- read.csv(concordance_path)

  # Build concordance matrix as a table
  conf_matrix <- concordance |>
    filter(!is.na(expected_tier) & !is.na(amp_tier_group)) |>
    count(expected_tier, amp_tier_group) |>
    pivot_wider(names_from = amp_tier_group, values_from = n, values_fill = 0) |>
    rename(`ClinVar Expected` = expected_tier)

  # Add row totals
  numeric_cols <- setdiff(names(conf_matrix), "ClinVar Expected")
  conf_matrix <- conf_matrix |>
    mutate(Total = rowSums(across(all_of(numeric_cols))))

  # Append summary row
  summary_row <- tibble(`ClinVar Expected` = "Total") |>
    bind_cols(as_tibble(t(colSums(conf_matrix |> select(all_of(c(numeric_cols, "Total")))))))
  conf_matrix <- bind_rows(conf_matrix, summary_row)

  # Add kappa from summary
  summary_path <- file.path(results_dir, "clinvar_concordance_summary.csv")
  if (file.exists(summary_path)) {
    summary_data <- read.csv(summary_path)
    kappa_val <- summary_data |> filter(metric == "cohens_kappa") |> pull(value)
    concordance_rate <- summary_data |> filter(metric == "concordance_rate") |> pull(value)

    table3_note <- glue("Cohen's kappa = {kappa_val}; Overall concordance = {round(as.numeric(concordance_rate) * 100, 1)}%")
  } else {
    table3_note <- NULL
  }

  save_table(conf_matrix, "table3_amp_concordance",
             title = "Table 3. AMP Classification Concordance Matrix",
             subtitle = table3_note)
} else {
  log_warn("Skipping Table 3: clinvar_concordance.csv not found")
  table3 <- tibble(Note = "Run 04_clinvar_concordance.R to generate this table")
  save_table(table3, "table3_amp_concordance",
             title = "Table 3. AMP Concordance Matrix (placeholder)")
}

# ==============================================================================
# Table 4: OncoKB vs CiVIC Evidence Summary
# ==============================================================================

log_info("Generating Table 4: Evidence summary...")

evidence_path <- file.path(results_dir, "evidence_concordance_by_tumor.csv")
if (file.exists(evidence_path)) {
  evidence_summary <- read.csv(evidence_path)

  table4 <- evidence_summary |>
    select(
      `Tumor Type` = tumor_type,
      `N Variants` = n_variants,
      `OncoKB` = n_oncokb_evidence,
      `CiVIC` = n_civic_evidence,
      `Both` = n_both,
      `OncoKB Only` = n_oncokb_only,
      `CiVIC Only` = n_civic_only,
      Neither = n_neither,
      `OncoKB %` = oncokb_coverage_pct,
      `CiVIC %` = civic_coverage_pct
    )

  # Add totals row
  overall_path <- file.path(results_dir, "evidence_concordance_summary.csv")
  if (file.exists(overall_path)) {
    overall <- read.csv(overall_path)
    total_row <- tibble(
      `Tumor Type` = "TOTAL",
      `N Variants` = as.integer(overall$value[overall$metric == "total_variants"]),
      `OncoKB` = NA_integer_,
      `CiVIC` = NA_integer_,
      `Both` = as.integer(overall$value[overall$metric == "both_oncokb_civic"]),
      `OncoKB Only` = as.integer(overall$value[overall$metric == "oncokb_only"]),
      `CiVIC Only` = as.integer(overall$value[overall$metric == "civic_only"]),
      Neither = as.integer(overall$value[overall$metric == "neither"]),
      `OncoKB %` = as.numeric(overall$value[overall$metric == "oncokb_coverage_pct"]),
      `CiVIC %` = as.numeric(overall$value[overall$metric == "civic_coverage_pct"])
    )
    table4 <- bind_rows(table4, total_row)
  }

  save_table(table4, "table4_evidence_summary",
             title = "Table 4. OncoKB vs CiVIC Evidence Coverage",
             subtitle = "Evidence availability for top 100 COSMIC variants by tumor type")
} else {
  log_warn("Skipping Table 4: evidence_concordance_by_tumor.csv not found")
  table4 <- tibble(Note = "Run 05_evidence_concordance.R to generate this table")
  save_table(table4, "table4_evidence_summary",
             title = "Table 4. Evidence Summary (placeholder)")
}

# ==============================================================================
# Table 5: Runtime Benchmarks
# ==============================================================================

log_info("Generating Table 5: Runtime benchmarks...")

runtime_path <- file.path(results_dir, "runtime_benchmark.csv")
if (file.exists(runtime_path)) {
  runtime <- read.csv(runtime_path)

  # Pivot: stages as rows, variant counts as columns
  table5 <- runtime |>
    select(stage, n_variants, seconds) |>
    mutate(seconds = round(seconds, 2)) |>
    pivot_wider(
      names_from = n_variants,
      values_from = seconds,
      names_prefix = "n="
    ) |>
    rename(Stage = stage)

  save_table(table5, "table5_runtime_benchmarks",
             title = "Table 5. Runtime Benchmarks (seconds)",
             subtitle = "Pipeline execution time by stage and input variant count")
} else {
  log_warn("Skipping Table 5: runtime_benchmark.csv not found")
  table5 <- tibble(Note = "Run 06_runtime_benchmark.R to generate this table")
  save_table(table5, "table5_runtime_benchmarks",
             title = "Table 5. Runtime Benchmarks (placeholder)")
}

# --- Summary ------------------------------------------------------------------

tables_created <- list.files(tables_dir, pattern = "\\.(csv|html)$")
log_info("Created {length(tables_created)} table files in {tables_dir}")

log_info("=== Table Generation Complete ===")
