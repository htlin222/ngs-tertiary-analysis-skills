# 00-qc/run_qc.R — Quality control metrics for NGS tertiary analysis
# Stage 0: BAM file QC and per-gene coverage analysis

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(logger)
  library(glue)
})

# ── Helper functions ────────────────────────────────────────────────────────

#' Ensure BAM file has an index
#' @param bam_path Path to BAM file
#' @details Creates .bai index via samtools if missing
ensure_bam_index <- function(bam_path) {
  index_path <- paste0(bam_path, ".bai")

  if (!file.exists(index_path)) {
    log_info("BAM index not found, creating: {index_path}")

    if (!check_tool("samtools")) {
      stop("samtools not found in PATH")
    }

    result <- system2("samtools", c("index", bam_path),
                     stdout = TRUE, stderr = TRUE)
    exit_code <- attr(result, "status") %||% 0L

    if (exit_code != 0) {
      log_error("Failed to create BAM index")
      stop(glue("samtools index failed: {exit_code}"))
    }

    log_info("BAM index created successfully")
  }

  invisible(index_path)
}

#' Load panel BED file and convert to GRanges
#' @param bed_path Path to panel BED file
#' @return GRanges object with genes as metadata
load_panel_bed <- function(bed_path) {
  if (!file.exists(bed_path)) {
    stop(glue("Panel BED file not found: {bed_path}"))
  }

  log_info("Loading panel BED: {bed_path}")

  # Read BED file (chrom, start, end, [name], [score], [strand], ...)
  bed <- read.delim(bed_path, header = FALSE, stringsAsFactors = FALSE,
                    colClasses = c("character", "integer", "integer",
                                   "character", "numeric", "character",
                                   rep("character", 3)))

  colnames(bed)[1:6] <- c("seqname", "start", "end", "name", "score", "strand")

  # Convert to 1-based (BED is 0-based for start)
  bed$start <- bed$start + 1

  # Create GRanges
  granges <- GRanges(
    seqnames = bed$seqname,
    ranges = IRanges(start = bed$start, end = bed$end),
    strand = bed$strand,
    gene = bed$name
  )

  log_info("Loaded {length(granges)} target regions from panel")
  granges
}

#' Calculate per-gene coverage statistics
#' @param bam_path Path to BAM file
#' @param panel_granges GRanges object with genes
#' @param min_mapq Minimum mapping quality
#' @return Tibble with per-gene coverage stats
calculate_per_gene_coverage <- function(bam_path, panel_granges, min_mapq = 30) {
  log_info("Calculating per-gene coverage")

  # Load BAM with appropriate parameters
  param <- ScanBamParam(
    which = panel_granges,
    what = c("rname", "pos", "qwidth", "mapq"),
    flag = scanBamFlag(isUnmappedQuery = FALSE)
  )

  bam <- readGAlignments(bam_path, param = param, use.names = FALSE)

  if (length(bam) == 0) {
    log_warn("No reads found in BAM file within panel regions")
    return(tibble(
      gene = character(),
      mean_coverage = numeric(),
      min_coverage = numeric(),
      pct_above_threshold = numeric()
    ))
  }

  # Filter by mapping quality
  bam <- bam[mcols(bam)$mapq >= min_mapq]

  if (length(bam) == 0) {
    log_warn("No reads remain after mapping quality filter (mapq >= {min_mapq})")
    return(tibble(
      gene = character(),
      mean_coverage = numeric(),
      min_coverage = numeric(),
      pct_above_threshold = numeric()
    ))
  }

  # Calculate coverage per gene
  coverage_results <- list()

  for (gene in unique(mcols(panel_granges)$gene)) {
    # Get regions for this gene
    gene_regions <- panel_granges[mcols(panel_granges)$gene == gene]

    # Count reads overlapping this gene's regions
    reads_in_gene <- bam[overlapsAny(bam, gene_regions)]

    if (length(reads_in_gene) == 0) {
      coverage_results[[gene]] <- tibble(
        gene = gene,
        mean_coverage = 0,
        min_coverage = 0,
        pct_above_threshold = 0
      )
      next
    }

    # Calculate position-level coverage using coverage function
    cov <- coverage(reads_in_gene, weight = rep(1, length(reads_in_gene)))

    # Get coverage values for all positions in gene regions
    cov_values <- integer()
    for (region in gene_regions) {
      chr <- as.character(seqnames(region))
      if (chr %in% names(cov)) {
        region_cov <- as.numeric(cov[[chr]][start(region):end(region)])
        cov_values <- c(cov_values, region_cov)
      }
    }

    if (length(cov_values) == 0) {
      coverage_results[[gene]] <- tibble(
        gene = gene,
        mean_coverage = 0,
        min_coverage = 0,
        pct_above_threshold = 0
      )
    } else {
      mean_cov <- mean(cov_values)
      min_cov <- min(cov_values)
      # Percentage of bases above 20x coverage (standard threshold)
      pct_above <- 100 * sum(cov_values >= 20) / length(cov_values)

      coverage_results[[gene]] <- tibble(
        gene = gene,
        mean_coverage = mean_cov,
        min_coverage = min_cov,
        pct_above_threshold = pct_above
      )
    }
  }

  per_gene_cov <- bind_rows(coverage_results) %>%
    arrange(desc(mean_coverage))

  log_info("Calculated coverage for {nrow(per_gene_cov)} genes")
  per_gene_cov
}

#' Create coverage distribution plot
#' @param bam_path Path to BAM file
#' @param panel_granges GRanges object with panel regions
#' @param output_path Path to save PNG
#' @param min_mapq Minimum mapping quality
create_coverage_plot <- function(bam_path, panel_granges, output_path,
                               min_mapq = 30) {
  log_info("Creating coverage distribution plot")

  # Load BAM with region filter
  param <- ScanBamParam(
    which = panel_granges,
    what = c("rname", "pos", "qwidth", "mapq"),
    flag = scanBamFlag(isUnmappedQuery = FALSE)
  )

  bam <- readGAlignments(bam_path, param = param, use.names = FALSE)

  if (length(bam) == 0) {
    log_warn("No reads found for coverage plot")
    # Create empty plot
    p <- ggplot() +
      theme_minimal() +
      ggtitle("Coverage Distribution (No Data)")

    ggsave(output_path, p, width = 10, height = 6, dpi = 100)
    return(invisible(NULL))
  }

  # Filter by mapping quality
  bam <- bam[mcols(bam)$mapq >= min_mapq]

  # Calculate position-level coverage across all panel regions
  coverage_list <- coverage(bam)

  # Extract coverage values from panel regions
  cov_data <- integer()

  for (region in panel_granges) {
    chr <- as.character(seqnames(region))
    if (chr %in% names(coverage_list)) {
      region_cov <- as.numeric(coverage_list[[chr]][start(region):end(region)])
      cov_data <- c(cov_data, region_cov)
    }
  }

  # Create histogram
  cov_df <- tibble(coverage = cov_data)

  # Calculate summary statistics
  mean_cov <- mean(cov_data)
  median_cov <- median(cov_data)

  p <- ggplot(cov_df, aes(x = coverage)) +
    geom_histogram(bins = 50, fill = "#2E86AB", alpha = 0.7, color = "black") +
    geom_vline(aes(xintercept = mean_cov, color = "Mean"),
               linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = median_cov, color = "Median"),
               linetype = "dotted", size = 1) +
    scale_color_manual(values = c("Mean" = "#A23B72", "Median" = "#F18F01")) +
    labs(
      title = "Coverage Distribution Across Panel Targets",
      x = "Coverage Depth (reads)",
      y = "Number of Positions",
      color = "Statistic"
    ) +
    annotate("text", x = Inf, y = Inf,
             label = glue("Mean: {round(mean_cov, 1)}x\nMedian: {round(median_cov, 1)}x"),
             hjust = 1.1, vjust = 1.1, size = 4, family = "mono") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      legend.position = "top"
    )

  # Ensure output directory exists
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

  ggsave(output_path, p, width = 10, height = 6, dpi = 100)
  log_info("Coverage plot saved to {output_path}")

  invisible(NULL)
}

# ── Main QC function ────────────────────────────────────────────────────────

#' Run quality control metrics on a BAM file
#'
#' @param bam_path Path to aligned BAM file
#' @param config Configuration list (must contain reference and qc sections)
#' @param sample_id Sample identifier for output organization
#'
#' @return List with:
#'   - `summary`: Named list of overall QC metrics
#'   - `per_gene_coverage`: Tibble of per-gene coverage stats
#'   - `pass`: Boolean indicating if sample passes QC thresholds
#'
#' @details
#' Generates:
#'   - reports/{sample_id}/00-qc/qc_summary.tsv
#'   - reports/{sample_id}/00-qc/per_gene_coverage.tsv
#'   - reports/{sample_id}/00-qc/coverage_plot.png
#'
#' @examples
#' \dontrun{
#'   config <- load_config("config/default.yaml")
#'   qc_results <- run_qc("sample.bam", config, "SAMPLE_001")
#'   print(qc_results$summary)
#' }
#'
#' @export
run_qc <- function(bam_path, config, sample_id) {

  log_info("Starting QC analysis for {sample_id}")

  # ── Validation ──────────────────────────────────────────────────────────

  if (!file.exists(bam_path)) {
    stop(glue("BAM file not found: {bam_path}"))
  }

  if (!is.list(config)) {
    stop("config must be a list")
  }

  required_config <- c("reference", "qc")
  if (!all(required_config %in% names(config))) {
    stop(glue("config must contain: {paste(required_config, collapse = ', ')}"))
  }

  # ── Ensure BAM index ────────────────────────────────────────────────────

  ensure_bam_index(bam_path)

  # ── Load panel regions ──────────────────────────────────────────────────

  panel_bed <- config$reference$panel_bed
  panel_granges <- load_panel_bed(panel_bed)

  # ── Get basic BAM statistics ────────────────────────────────────────────

  log_info("Reading basic BAM statistics")

  # Use idxstats for read counts
  idxstats_output <- system2("samtools", c("idxstats", bam_path),
                            stdout = TRUE, stderr = FALSE)

  idxstats_df <- read.delim(text = paste(idxstats_output, collapse = "\n"),
                            header = FALSE, stringsAsFactors = FALSE)
  colnames(idxstats_df) <- c("chrom", "length", "mapped", "unmapped")

  total_reads_mapped <- sum(idxstats_df$mapped, na.rm = TRUE)
  total_reads_unmapped <- sum(idxstats_df$unmapped, na.rm = TRUE)
  total_reads <- total_reads_mapped + total_reads_unmapped

  mapping_rate <- if (total_reads > 0) {
    total_reads_mapped / total_reads
  } else {
    0
  }

  log_info("Total reads: {total_reads}, Mapped: {total_reads_mapped}")

  # ── Get detailed BAM statistics (including duplicates) ──────────────────

  log_info("Computing detailed BAM flag statistics")

  flag_summary <- system2("samtools", c("flagstat", bam_path),
                         stdout = TRUE, stderr = FALSE)

  # Parse flagstat output
  flagstat_lines <- trimws(flag_summary)

  # Extract duplicate count (typically last meaningful entry)
  duplicates <- 0
  for (line in flagstat_lines) {
    if (grepl("duplicates", line, ignore.case = TRUE)) {
      # Extract number from "X + Y duplicates"
      match <- regexpr("^[0-9]+", line)
      if (match > -1) {
        duplicates <- as.numeric(substr(line, match,
                                        match + attr(match, "match.length") - 1))
      }
    }
  }

  duplicate_rate <- if (total_reads > 0) {
    duplicates / total_reads
  } else {
    0
  }

  log_info("Duplicate reads: {duplicates} ({round(100*duplicate_rate, 2)}%)")

  # ── Calculate per-gene coverage ─────────────────────────────────────────

  per_gene_coverage <- calculate_per_gene_coverage(
    bam_path,
    panel_granges,
    min_mapq = config$qc$min_mapping_quality %||% 30
  )

  mean_coverage <- if (nrow(per_gene_coverage) > 0) {
    mean(per_gene_coverage$mean_coverage)
  } else {
    0
  }

  log_info("Mean coverage across genes: {round(mean_coverage, 2)}x")

  # ── Calculate on-target rate ────────────────────────────────────────────

  log_info("Computing on-target rate")

  # Load all reads (not just in panel)
  all_bam <- readGAlignments(bam_path, use.names = FALSE)

  # Count reads overlapping panel regions
  reads_on_target <- sum(overlapsAny(all_bam, panel_granges))
  on_target_rate <- if (length(all_bam) > 0) {
    reads_on_target / length(all_bam)
  } else {
    0
  }

  log_info("On-target reads: {reads_on_target} / {length(all_bam)} ({round(100*on_target_rate, 2)}%)")

  # ── Tumor purity (placeholder) ──────────────────────────────────────────

  log_info("Tumor purity estimation requires external tool (FACETS, sequenza, etc.)")
  tumor_purity <- NA_real_
  purity_note <- "External tool (FACETS/sequenza) needed for accurate estimation"

  # ── Create summary list ─────────────────────────────────────────────────

  summary <- list(
    sample_id = sample_id,
    bam_file = basename(bam_path),
    total_reads = total_reads,
    mapped_reads = total_reads_mapped,
    unmapped_reads = total_reads_unmapped,
    mapping_rate = mapping_rate,
    duplicate_reads = duplicates,
    duplicate_rate = duplicate_rate,
    on_target_reads = reads_on_target,
    on_target_rate = on_target_rate,
    mean_coverage = mean_coverage,
    median_coverage = if (nrow(per_gene_coverage) > 0) {
      median(per_gene_coverage$mean_coverage)
    } else {
      0
    },
    num_genes_covered = nrow(per_gene_coverage),
    tumor_purity = tumor_purity,
    purity_note = purity_note
  )

  # ── QC pass/fail evaluation ─────────────────────────────────────────────

  min_mean_coverage <- config$qc$min_mean_coverage %||% 200
  min_on_target_rate <- config$qc$min_on_target_rate %||% 0.80
  fail_on_low_coverage <- config$qc$fail_on_low_coverage %||% FALSE

  qc_pass <- TRUE
  fail_reasons <- character()

  if (mean_coverage < min_mean_coverage) {
    if (fail_on_low_coverage) {
      qc_pass <- FALSE
      fail_reasons <- c(fail_reasons,
        glue("Mean coverage {round(mean_coverage, 2)}x < {min_mean_coverage}x"))
    } else {
      log_warn("Mean coverage below threshold (but not failing)")
    }
  }

  if (on_target_rate < min_on_target_rate) {
    qc_pass <- FALSE
    fail_reasons <- c(fail_reasons,
      glue("On-target rate {round(100*on_target_rate, 2)}% < {round(100*min_on_target_rate, 2)}%"))
  }

  if (!qc_pass) {
    log_warn("QC FAILED: {paste(fail_reasons, collapse = '; ')}")
  } else {
    log_info("QC PASSED")
  }

  summary$qc_pass <- qc_pass
  summary$qc_notes <- if (qc_pass) "PASS" else paste(fail_reasons, collapse = "; ")

  # ── Create output directory ─────────────────────────────────────────────

  out_dir <- stage_output_dir(sample_id, "00-qc")

  # ── Save summary TSV ────────────────────────────────────────────────────

  summary_path <- file.path(out_dir, "qc_summary.tsv")

  summary_df <- tibble(
    metric = names(summary),
    value = as.character(unlist(summary))
  )

  write.table(summary_df, summary_path, sep = "\t", row.names = FALSE,
              quote = FALSE)

  log_info("Saved QC summary to {summary_path}")

  # ── Save per-gene coverage TSV ──────────────────────────────────────────

  per_gene_path <- file.path(out_dir, "per_gene_coverage.tsv")

  write.table(per_gene_coverage, per_gene_path, sep = "\t", row.names = FALSE,
              quote = FALSE)

  log_info("Saved per-gene coverage to {per_gene_path}")

  # ── Create coverage plot ────────────────────────────────────────────────

  plot_path <- file.path(out_dir, "coverage_plot.png")

  create_coverage_plot(bam_path, panel_granges, plot_path,
                       min_mapq = config$qc$min_mapping_quality %||% 30)

  # ── Return results ──────────────────────────────────────────────────────

  log_info("QC analysis complete for {sample_id}")

  result <- list(
    summary = summary,
    per_gene_coverage = per_gene_coverage,
    pass = qc_pass
  )

  class(result) <- c("qc_result", "list")
  result
}

# ── S3 print method for QC results ──────────────────────────────────────────

#' @export
print.qc_result <- function(x, ...) {
  cat("QC Analysis Results\n")
  cat("===================\n\n")

  cat(glue("Sample: {x$summary$sample_id}\n"))
  cat(glue("QC Status: {if (x$pass) 'PASS' else 'FAIL'}\n\n"))

  cat("Key Metrics:\n")
  cat(glue("  Total Reads:        {format(x$summary$total_reads, big.mark=',')}\n"))
  cat(glue("  Mapped Reads:       {format(x$summary$mapped_reads, big.mark=',')} ({round(100*x$summary$mapping_rate, 2)}%)\n"))
  cat(glue("  Duplicate Rate:     {round(100*x$summary$duplicate_rate, 2)}%\n"))
  cat(glue("  On-Target Rate:     {round(100*x$summary$on_target_rate, 2)}%\n"))
  cat(glue("  Mean Coverage:      {round(x$summary$mean_coverage, 2)}x\n"))
  cat(glue("  Genes Covered:      {x$summary$num_genes_covered}\n"))

  if (!is.na(x$summary$tumor_purity)) {
    cat(glue("  Tumor Purity:       {round(100*x$summary$tumor_purity, 2)}%\n"))
  } else {
    cat(glue("  Tumor Purity:       Not determined ({x$summary$purity_note})\n"))
  }

  if (!x$pass) {
    cat(glue("\nQC Notes: {x$summary$qc_notes}\n"))
  }

  invisible(x)
}
