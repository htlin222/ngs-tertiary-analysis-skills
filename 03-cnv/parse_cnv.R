# 03-cnv/parse_cnv.R — Parse and visualize CNVkit results
# Processes CNVkit gene-level copy number calls and generates clinical summaries

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(logger)
  library(ggplot2)
})

# Load shared utilities
source(here::here("R/utils.R"))

#' Parse and Classify CNVkit Copy Number Calls
#'
#' Reads CNVkit gene-level copy number output (.call.cns), classifies variants
#' by amplification/gain/deletion status, and filters to clinically significant
#' changes. Generates a genome-wide visualization and saves parsed results.
#'
#' @param cnvkit_results List. Output from run_cnvkit() containing:
#'   - cnr: Path to bin-level coverage file
#'   - cns: Path to segmented copy number file
#'   - call_cns: Path to gene-level copy number calls
#' @param config List. Pipeline configuration (from load_config())
#' @param sample_id Character. Sample identifier used for output organization
#'
#' @return Tibble with significant CNAs:
#'   - chromosome, start, end, gene: Genomic location and gene annotation
#'   - log2_ratio: Log2 copy number ratio
#'   - copy_number: Absolute copy number estimate
#'   - length: Segment length in base pairs
#'   - type: CNA classification (AMPLIFICATION, GAIN, DELETION, NEUTRAL)
#'
#' @details
#' Classification thresholds (from config):
#' - AMPLIFICATION: log2 >= config$cnv$amplification_threshold (typically ~1.5)
#' - GAIN: log2 >= config$cnv$min_log2_gain (typically ~0.32, ~3 copies)
#' - DELETION: log2 <= config$cnv$min_log2_loss (typically ~-0.42, ~1 copy)
#' - NEUTRAL: otherwise (within diploid range)
#'
#' Outputs:
#' - cnvkit_segments.tsv: Tab-separated table of all significant CNAs
#' - cnv_plot.png: Genome-wide visualization with chromosomes colored by type
#'
#' @examples
#' \dontrun{
#'   config <- load_config("config/default.yaml")
#'   results <- run_cnvkit("samples/tumor.bam", config, "SAMPLE_001")
#'   cnas <- parse_cnv(results, config, "SAMPLE_001")
#'   nrow(cnas)  # Count of significant CNAs
#' }
parse_cnv <- function(cnvkit_results, config, sample_id) {
  log_info("Starting CNV parsing and classification for {sample_id}")

  # ── Input validation ──────────────────────────────────────────────────────
  if (is.null(cnvkit_results$call_cns) || !file.exists(cnvkit_results$call_cns)) {
    log_error("Gene-level calls file not found: {cnvkit_results$call_cns}")
    stop("CNVkit gene-level calls file (.call.cns) is required")
  }

  # Verify thresholds are configured
  if (is.null(config$cnv$amplification_threshold)) {
    log_warn("amplification_threshold not configured, using default 1.5")
    config$cnv$amplification_threshold <- 1.5
  }
  if (is.null(config$cnv$min_log2_gain)) {
    log_warn("min_log2_gain not configured, using default 0.32 (~3 copies)")
    config$cnv$min_log2_gain <- 0.32
  }
  if (is.null(config$cnv$min_log2_loss)) {
    log_warn("min_log2_loss not configured, using default -0.42 (~1 copy)")
    config$cnv$min_log2_loss <- -0.42
  }

  log_info("CNV thresholds - Amplification: {config$cnv$amplification_threshold}, "
    "Gain: {config$cnv$min_log2_gain}, Loss: {config$cnv$min_log2_loss}")

  # ── Read CNVkit gene-level calls ──────────────────────────────────────────
  log_info("Reading CNVkit gene-level calls: {cnvkit_results$call_cns}")

  # CNVkit .call.cns format is tab-separated with columns:
  # chromosome, start, end, gene, log2, cn (plus optional cn.se, cn.mean)
  cna_data <- readr::read_tsv(
    cnvkit_results$call_cns,
    col_types = readr::cols(
      chromosome = readr::col_character(),
      start = readr::col_integer(),
      end = readr::col_integer(),
      gene = readr::col_character(),
      log2 = readr::col_double(),
      cn = readr::col_double(),
      .default = readr::col_skip()
    ),
    show_col_types = FALSE
  ) %>%
    as_tibble()

  log_info("Read {nrow(cna_data)} gene-level CNV entries")

  # ── Parse and classify CNVs ───────────────────────────────────────────────
  cna_classified <- cna_data %>%
    mutate(
      # Calculate segment length
      length = end - start,
      # Rename columns for clarity
      log2_ratio = log2,
      copy_number = cn,
      # Classify CNA type
      type = case_when(
        log2_ratio >= config$cnv$amplification_threshold ~ "AMPLIFICATION",
        log2_ratio >= config$cnv$min_log2_gain ~ "GAIN",
        log2_ratio <= config$cnv$min_log2_loss ~ "DELETION",
        TRUE ~ "NEUTRAL"
      )
    ) %>%
    select(
      chromosome, start, end, gene, log2_ratio, copy_number, length, type
    )

  log_info("Classified CNVs: {sum(cna_classified$type != 'NEUTRAL')} significant, "
    "{sum(cna_classified$type == 'NEUTRAL')} neutral")

  # ── Filter to clinically significant CNAs ─────────────────────────────────
  cna_significant <- cna_classified %>%
    filter(type != "NEUTRAL")

  log_info("Retained {nrow(cna_significant)} clinically significant CNAs")

  # Log summary by type
  summary_by_type <- cna_significant %>%
    count(type, name = "count") %>%
    mutate(label = glue("{type}: {count}")) %>%
    pull(label)

  if (length(summary_by_type) > 0) {
    log_info("Summary: {paste(summary_by_type, collapse = ', ')}")
  } else {
    log_info("No clinically significant CNAs detected")
  }

  # ── Prepare output directory ──────────────────────────────────────────────
  output_dir <- stage_output_dir(sample_id, "03-cnv")

  # ── Save parsed results to TSV ────────────────────────────────────────────
  output_tsv <- path(output_dir, "cnvkit_segments.tsv")
  log_info("Writing CNV results to: {output_tsv}")

  readr::write_tsv(cna_significant, output_tsv)
  log_info("Saved {nrow(cna_significant)} significant CNVs to {output_tsv}")

  # ── Generate genome-wide visualization ────────────────────────────────────
  log_info("Generating genome-wide CNV plot")

  # Prepare data for plotting
  # Convert chromosome to factor with natural order (1-22, X, Y, MT)
  chrom_order <- c(as.character(1:22), "X", "Y", "MT")
  cna_plot_data <- cna_classified %>%
    filter(type != "NEUTRAL") %>%
    mutate(
      chromosome = factor(chromosome, levels = chrom_order),
      type_color = case_when(
        type == "AMPLIFICATION" ~ "#d62728",  # Red
        type == "GAIN" ~ "#ff7f0e",           # Orange
        type == "DELETION" ~ "#1f77b4",       # Blue
        TRUE ~ "#7f7f7f"                      # Gray
      )
    ) %>%
    arrange(chromosome)

  # Create genome-wide plot
  p <- ggplot(cna_plot_data, aes(
    x = chromosome,
    y = log2_ratio,
    fill = type,
    color = type
  )) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    geom_hline(
      yintercept = config$cnv$amplification_threshold,
      linetype = "dashed",
      color = "#d62728",
      alpha = 0.5,
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = config$cnv$min_log2_gain,
      linetype = "dashed",
      color = "#ff7f0e",
      alpha = 0.5,
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = config$cnv$min_log2_loss,
      linetype = "dashed",
      color = "#1f77b4",
      alpha = 0.5,
      linewidth = 0.8
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
    scale_fill_manual(
      name = "Type",
      values = c(
        AMPLIFICATION = "#d62728",
        GAIN = "#ff7f0e",
        DELETION = "#1f77b4"
      )
    ) +
    scale_color_manual(
      name = "Type",
      values = c(
        AMPLIFICATION = "#d62728",
        GAIN = "#ff7f0e",
        DELETION = "#1f77b4"
      )
    ) +
    labs(
      title = glue("Genome-Wide Copy Number Profile - {sample_id}"),
      x = "Chromosome",
      y = "Log2 Copy Number Ratio",
      caption = glue("Significant CNAs (n={nrow(cna_plot_data)}) | "
        "Amplification >= {config$cnv$amplification_threshold}, "
        "Deletion <= {config$cnv$min_log2_loss}")
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "top"
    )

  # Save plot to PNG
  output_plot <- path(output_dir, "cnv_plot.png")
  log_info("Saving CNV plot to: {output_plot}")

  ggsave(
    output_plot,
    plot = p,
    width = 14,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  log_info("CNV plot saved: {output_plot}")

  # ── Log completion ────────────────────────────────────────────────────────
  log_info("CNV parsing and visualization completed for {sample_id}")
  log_info("Results saved to:")
  log_info("  TSV: {output_tsv}")
  log_info("  Plot: {output_plot}")

  # Return the parsed and classified CNA data
  invisible(cna_significant)
}
