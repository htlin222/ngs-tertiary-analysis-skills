#!/usr/bin/env Rscript
# benchmarks/scripts/07_generate_figures.R
# Generate publication-quality figures for the methods paper.

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

source(here::here("R/utils.R"))

log_info("=== Generate Publication Figures ===")

# --- Configuration -----------------------------------------------------------

results_dir <- here::here("benchmarks/results")
figures_dir <- here::here("benchmarks/results/figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Publication theme
theme_publication <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2,
                                      hjust = 0, margin = margin(b = 8)),
      plot.subtitle    = element_text(size = base_size, color = "grey40",
                                      margin = margin(b = 12)),
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      strip.text       = element_text(face = "bold", size = base_size),
      plot.margin      = margin(12, 12, 12, 12)
    )
}

# Color palette
palette_main <- c(
  "#2c3e50", "#e74c3c", "#3498db", "#2ecc71", "#f39c12",
  "#9b59b6", "#1abc9c", "#e67e22", "#34495e", "#16a085"
)

save_figure <- function(plot, name, width = 8, height = 5) {
  pdf_path <- file.path(figures_dir, paste0(name, ".pdf"))
  png_path <- file.path(figures_dir, paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, dpi = 300)
  ggsave(png_path, plot, width = width, height = height, dpi = 300)
  log_info("Saved figure: {name} (PDF + PNG)")
}

# ==============================================================================
# Figure 1: Pipeline Overview Diagram
# ==============================================================================

log_info("Generating Figure 1: Pipeline overview...")

pipeline_stages <- tibble(
  stage = factor(c(
    "00 QC", "01 Variant\nCalling", "02 Annotation",
    "03 CNV", "04 Fusions", "05 Biomarkers",
    "06 Clinical\nAnnotation", "07 Literature", "08 Report"
  ), levels = c(
    "00 QC", "01 Variant\nCalling", "02 Annotation",
    "03 CNV", "04 Fusions", "05 Biomarkers",
    "06 Clinical\nAnnotation", "07 Literature", "08 Report"
  )),
  tool = c(
    "Rsamtools\nsamtools", "Mutect2\nGATK", "VEP\nClinVar",
    "CNVkit\nhybrid", "Manta\nSVs", "TMB\nMSI/HRD",
    "OncoKB\nCiVIC", "PubMed\nScopus", "Quarto\nESMO"
  ),
  category = c(
    "QC", "Variant", "Variant", "Structural",
    "Structural", "Biomarker", "Clinical", "Clinical", "Report"
  ),
  x = seq_len(9),
  y = 1
)

category_colors <- c(
  "QC" = "#3498db", "Variant" = "#e74c3c", "Structural" = "#f39c12",
  "Biomarker" = "#2ecc71", "Clinical" = "#9b59b6", "Report" = "#2c3e50"
)

fig1 <- ggplot(pipeline_stages, aes(x = x, y = y)) +
  geom_tile(aes(fill = category), width = 0.9, height = 0.6,
            color = "white", linewidth = 1) +
  geom_text(aes(label = stage), vjust = -0.3, fontface = "bold", size = 3) +
  geom_text(aes(label = tool), vjust = 1.5, size = 2.5, color = "grey30") +
  # Arrows between stages
  geom_segment(
    data = pipeline_stages |> filter(x < max(x)),
    aes(x = x + 0.5, xend = x + 0.5, y = y, yend = y),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
    color = "grey50", linewidth = 0.5
  ) +
  scale_fill_manual(values = category_colors, name = "Category") +
  labs(
    title = "NGS Tertiary Analysis Pipeline Architecture",
    subtitle = "BAM/VCF input through ESMO 2024-compliant clinical report generation"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5,
                                 margin = margin(b = 10)),
    legend.position = "bottom",
    plot.margin   = margin(15, 15, 15, 15)
  ) +
  coord_cartesian(ylim = c(0.5, 1.5))

save_figure(fig1, "fig1_pipeline_overview", width = 12, height = 4)

# ==============================================================================
# Figure 2: Variant Detection Sensitivity vs VAF
# ==============================================================================

log_info("Generating Figure 2: Sensitivity vs VAF...")

vaf_data_path <- file.path(results_dir, "variant_accuracy_by_vaf.csv")
if (file.exists(vaf_data_path)) {
  vaf_data <- read.csv(vaf_data_path)

  fig2 <- ggplot(vaf_data, aes(x = vaf_threshold, y = sensitivity)) +
    geom_ribbon(aes(ymin = pmax(sensitivity - 0.02, 0),
                    ymax = pmin(sensitivity + 0.02, 1)),
                fill = "#3498db", alpha = 0.15) +
    geom_line(color = "#2c3e50", linewidth = 1) +
    geom_point(color = "#e74c3c", size = 1.2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.42, y = 0.96, label = "95% target",
             color = "grey50", size = 3.5, fontface = "italic") +
    scale_x_continuous(labels = percent_format(), breaks = seq(0, 0.5, 0.05)) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1),
                       breaks = seq(0, 1, 0.1)) +
    labs(
      title = "Variant Detection Sensitivity by VAF Threshold",
      subtitle = "SEQC2 HCC1395 truth set filtered to TSO500 panel regions",
      x = "Variant Allele Frequency (VAF) Threshold",
      y = "Sensitivity"
    ) +
    theme_publication()

  save_figure(fig2, "fig2_sensitivity_vs_vaf")
} else {
  log_warn("Skipping Figure 2: variant_accuracy_by_vaf.csv not found")
  log_warn("Run 03_variant_accuracy.R first")
}

# ==============================================================================
# Figure 3: AMP Classification Concordance Heatmap
# ==============================================================================

log_info("Generating Figure 3: AMP concordance heatmap...")

concordance_path <- file.path(results_dir, "clinvar_concordance.csv")
if (file.exists(concordance_path)) {
  concordance <- read.csv(concordance_path)

  # Build confusion matrix
  conf_matrix <- concordance |>
    filter(!is.na(expected_tier) & !is.na(amp_tier_group)) |>
    count(expected_tier, amp_tier_group) |>
    complete(expected_tier, amp_tier_group, fill = list(n = 0))

  # Add percentage labels
  conf_matrix <- conf_matrix |>
    group_by(expected_tier) |>
    mutate(
      total = sum(n),
      pct = ifelse(total > 0, n / total, 0),
      label = ifelse(n > 0, paste0(n, "\n(", round(pct * 100), "%)"), "0")
    ) |>
    ungroup()

  fig3 <- ggplot(conf_matrix, aes(x = amp_tier_group, y = expected_tier,
                                   fill = pct)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = label), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "white", mid = "#3498db", high = "#2c3e50",
                         midpoint = 0.5, name = "Proportion",
                         labels = percent_format()) +
    labs(
      title = "AMP Classification vs ClinVar Concordance",
      subtitle = "Confusion matrix of pipeline AMP tier vs expected ClinVar classification",
      x = "Pipeline AMP Classification",
      y = "Expected (from ClinVar)"
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid = element_blank()
    )

  save_figure(fig3, "fig3_amp_concordance_heatmap", width = 7, height = 5)
} else {
  log_warn("Skipping Figure 3: clinvar_concordance.csv not found")
  log_warn("Run 04_clinvar_concordance.R first")
}

# ==============================================================================
# Figure 4: OncoKB vs CiVIC Venn Diagram (as Euler/Proportional Bar)
# ==============================================================================

log_info("Generating Figure 4: OncoKB vs CiVIC evidence overlap...")

evidence_path <- file.path(results_dir, "evidence_concordance.csv")
if (file.exists(evidence_path)) {
  evidence <- read.csv(evidence_path)

  # Create overlap categories
  overlap <- evidence |>
    mutate(
      category = case_when(
        has_oncokb_evidence & has_civic_evidence ~ "Both",
        has_oncokb_evidence & !has_civic_evidence ~ "OncoKB Only",
        !has_oncokb_evidence & has_civic_evidence ~ "CiVIC Only",
        TRUE ~ "Neither"
      ),
      category = factor(category, levels = c("Both", "OncoKB Only",
                                              "CiVIC Only", "Neither"))
    )

  # Overall counts
  overlap_counts <- overlap |>
    count(category) |>
    mutate(pct = n / sum(n))

  # Stacked bar by tumor type
  overlap_by_tumor <- overlap |>
    count(tumor_type, category) |>
    group_by(tumor_type) |>
    mutate(pct = n / sum(n)) |>
    ungroup()

  fig4 <- ggplot(overlap_by_tumor,
                 aes(x = tumor_type, y = n, fill = category)) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_stack(vjust = 0.5),
              size = 3, fontface = "bold", color = "white") +
    scale_fill_manual(
      values = c("Both" = "#2ecc71", "OncoKB Only" = "#3498db",
                 "CiVIC Only" = "#e74c3c", "Neither" = "#bdc3c7"),
      name = "Evidence Source"
    ) +
    labs(
      title = "OncoKB vs CiVIC Evidence Coverage",
      subtitle = "Evidence availability for top COSMIC variants by tumor type",
      x = "Tumor Type",
      y = "Number of Variants"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  save_figure(fig4, "fig4_oncokb_vs_civic", width = 9, height = 5)

  # Also create a summary Venn-style waffle
  fig4b <- ggplot(overlap_counts, aes(x = "", y = n, fill = category)) +
    geom_col(width = 1) +
    coord_polar("y") +
    geom_text(aes(label = paste0(category, "\n", n, " (", round(pct * 100), "%)")),
              position = position_stack(vjust = 0.5),
              size = 3, fontface = "bold", color = "white") +
    scale_fill_manual(
      values = c("Both" = "#2ecc71", "OncoKB Only" = "#3498db",
                 "CiVIC Only" = "#e74c3c", "Neither" = "#bdc3c7")
    ) +
    labs(title = "Overall Evidence Overlap") +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  save_figure(fig4b, "fig4b_evidence_overlap_pie", width = 6, height = 6)
} else {
  log_warn("Skipping Figure 4: evidence_concordance.csv not found")
  log_warn("Run 05_evidence_concordance.R first")
}

# ==============================================================================
# Figure 5: Runtime by Stage and Variant Count
# ==============================================================================

log_info("Generating Figure 5: Runtime benchmarks...")

runtime_path <- file.path(results_dir, "runtime_benchmark.csv")
if (file.exists(runtime_path)) {
  runtime <- read.csv(runtime_path)

  # Exclude TOTAL row for per-stage plot
  runtime_stages <- runtime |> filter(stage != "TOTAL")
  runtime_total  <- runtime |> filter(stage == "TOTAL")

  if (nrow(runtime_stages) > 0) {
    fig5a <- ggplot(runtime_stages,
                    aes(x = n_variants, y = seconds, color = stage)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      scale_x_continuous(breaks = unique(runtime_stages$n_variants)) +
      scale_color_manual(values = palette_main, name = "Pipeline Stage") +
      labs(
        title = "Pipeline Runtime by Stage",
        subtitle = "Execution time scaling with variant count (VCF input mode)",
        x = "Number of Input Variants",
        y = "Time (seconds)"
      ) +
      theme_publication() +
      theme(legend.position = "right")

    save_figure(fig5a, "fig5a_runtime_by_stage", width = 10, height = 6)
  }

  if (nrow(runtime_total) > 0) {
    fig5b <- ggplot(runtime_total, aes(x = n_variants, y = seconds)) +
      geom_line(color = "#2c3e50", linewidth = 1) +
      geom_point(color = "#e74c3c", size = 3) +
      geom_text(aes(label = paste0(round(seconds, 1), "s")),
                vjust = -1, size = 3.5) +
      scale_x_continuous(breaks = unique(runtime_total$n_variants)) +
      labs(
        title = "Total Pipeline Runtime vs Variant Count",
        subtitle = "End-to-end execution time (VCF input mode, skipping BAM stages)",
        x = "Number of Input Variants",
        y = "Total Time (seconds)"
      ) +
      theme_publication()

    save_figure(fig5b, "fig5b_total_runtime", width = 8, height = 5)
  }

  # Stacked area: proportion of time per stage
  if (nrow(runtime_stages) > 0) {
    runtime_prop <- runtime_stages |>
      group_by(n_variants) |>
      mutate(pct = seconds / sum(seconds)) |>
      ungroup()

    fig5c <- ggplot(runtime_prop,
                    aes(x = n_variants, y = pct, fill = stage)) +
      geom_area(alpha = 0.8, position = "stack") +
      scale_x_continuous(breaks = unique(runtime_prop$n_variants)) +
      scale_y_continuous(labels = percent_format()) +
      scale_fill_manual(values = palette_main, name = "Pipeline Stage") +
      labs(
        title = "Proportion of Runtime by Stage",
        subtitle = "Relative time distribution across pipeline stages",
        x = "Number of Input Variants",
        y = "Proportion of Total Time"
      ) +
      theme_publication()

    save_figure(fig5c, "fig5c_runtime_proportions", width = 10, height = 6)
  }
} else {
  log_warn("Skipping Figure 5: runtime_benchmark.csv not found")
  log_warn("Run 06_runtime_benchmark.R first")
}

# --- Summary ------------------------------------------------------------------

figures_created <- list.files(figures_dir, pattern = "\\.(pdf|png)$")
log_info("Created {length(figures_created)} figure files in {figures_dir}")

log_info("=== Figure Generation Complete ===")
