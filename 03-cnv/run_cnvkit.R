# 03-cnv/run_cnvkit.R — CNVkit copy number variant detection
# Performs copy number variant detection using CNVkit in hybrid-capture mode

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(fs)
})

# Load shared utilities
source(here::here("R/utils.R"))

#' Run CNVkit for Copy Number Variant Detection
#'
#' Performs copy number variant (CNV) detection on a tumor BAM file using CNVkit
#' in hybrid-capture mode. Includes bin-level read depth analysis, segmentation,
#' and gene-level copy number calls.
#'
#' @param bam_path Character. Path to tumor BAM file (indexed)
#' @param config List. Pipeline configuration (from load_config())
#' @param sample_id Character. Sample identifier used for output naming
#'
#' @return List with elements:
#'   - cnr: Character. Path to bin-level coverage file (.cnr)
#'   - cns: Character. Path to segmented copy number file (.cns)
#'   - call_cns: Character. Path to gene-level copy number calls (.call.cns)
#'
#' @details
#' The function:
#' 1. Checks CNVkit availability via PATH
#' 2. Runs cnvkit.py batch in hybrid-capture mode:
#'    - Targets panel regions (--targets)
#'    - Reference genome (--fasta)
#'    - Optional reference CNN model for tumor/normal matching
#'    - Output directory and sample prefix
#' 3. Performs segmentation on bin-level coverage
#' 4. Generates gene-level copy number calls
#' 5. Returns paths to all output files for downstream analysis
#'
#' @examples
#' \dontrun{
#'   config <- load_config("config/default.yaml")
#'   results <- run_cnvkit(
#'     "samples/tumor.bam",
#'     config,
#'     "SAMPLE_001"
#'   )
#'   # results$cnr, results$cns, results$call_cns available for parse_cnv()
#' }
run_cnvkit <- function(bam_path, config, sample_id) {
  log_info("Starting CNVkit analysis for {sample_id}")

  # ── Input validation ──────────────────────────────────────────────────────
  if (!file.exists(bam_path)) {
    log_error("BAM file not found: {bam_path}")
    stop(glue("BAM file not found: {bam_path}"))
  }

  if (is.null(config$reference$panel_bed) || !file.exists(config$reference$panel_bed)) {
    log_error("Panel BED file not found: {config$reference$panel_bed}")
    stop("Panel BED file is required for CNVkit hybrid-capture mode")
  }

  if (is.null(config$reference$fasta) || !file.exists(config$reference$fasta)) {
    log_error("Reference FASTA not found: {config$reference$fasta}")
    stop("Reference FASTA is required for CNVkit")
  }

  # ── Check CNVkit availability ─────────────────────────────────────────────
  cnvkit_cmd <- "cnvkit.py"

  if (nchar(Sys.which(cnvkit_cmd)) == 0) {
    log_error("CNVkit not found in PATH")
    stop("CNVkit is required but not available in PATH")
  }

  log_debug("Using CNVkit from PATH: {cnvkit_cmd}")

  # ── Prepare output directory and paths ─────────────────────────────────────
  output_dir <- stage_output_dir(sample_id, "03-cnv")
  sample_prefix <- path(output_dir, sample_id)

  log_info("Output directory: {output_dir}")
  log_info("Sample prefix: {sample_prefix}")

  # Define output files
  cnr_file <- glue("{sample_prefix}.cnr")
  cns_file <- glue("{sample_prefix}.cns")
  call_cns_file <- glue("{sample_prefix}.call.cns")

  # ── Build CNVkit batch command ────────────────────────────────────────────
  batch_args <- c("batch")
  batch_args <- c(batch_args, bam_path)
  batch_args <- c(batch_args, "--method", "hybrid")
  batch_args <- c(batch_args, "--targets", config$reference$panel_bed)
  batch_args <- c(batch_args, "--fasta", config$reference$fasta)
  batch_args <- c(batch_args, "--output-dir", output_dir)

  # Optional: Reference CNN model
  # If available, use it for tumor/normal matching; otherwise use flat reference
  if (!is.null(config$cnv$reference_cnn) && file.exists(config$cnv$reference_cnn)) {
    batch_args <- c(batch_args, "--reference", config$cnv$reference_cnn)
    log_info("Using reference CNN model: {config$cnv$reference_cnn}")
  } else {
    batch_args <- c(batch_args, "--normal")
    log_info("Using flat reference (tumor-only mode, no reference CNN)")
  }

  # ── Run CNVkit batch (coverage + reference) ───────────────────────────────
  log_info("Running CNVkit batch: coverage and reference generation")
  run_tool(
    cnvkit_cmd,
    batch_args,
    description = glue("CNVkit batch for {sample_id}")
  )

  # Verify coverage file (.cnr) was created
  if (!file.exists(cnr_file)) {
    log_error("CNVkit batch did not produce coverage file: {cnr_file}")
    stop("CNVkit batch failed to generate .cnr file")
  }

  log_info("Coverage file generated: {cnr_file}")

  # ── Run CNVkit segmentation ───────────────────────────────────────────────
  log_info("Running CNVkit segmentation on coverage data")

  segment_args <- c("segment", cnr_file, "-o", cns_file)

  run_tool(
    cnvkit_cmd,
    segment_args,
    description = glue("CNVkit segmentation for {sample_id}")
  )

  # Verify segmentation output
  if (!file.exists(cns_file)) {
    log_error("CNVkit segmentation did not produce file: {cns_file}")
    stop("CNVkit segmentation failed to generate .cns file")
  }

  log_info("Segmentation file generated: {cns_file}")

  # ── Run CNVkit gene-level calls ───────────────────────────────────────────
  log_info("Running CNVkit gene-level copy number calls")

  call_args <- c("call", cns_file, "-o", call_cns_file)

  run_tool(
    cnvkit_cmd,
    call_args,
    description = glue("CNVkit gene-level calls for {sample_id}")
  )

  # Verify gene-level calls output
  if (!file.exists(call_cns_file)) {
    log_error("CNVkit call did not produce file: {call_cns_file}")
    stop("CNVkit gene-level calls failed to generate .call.cns file")
  }

  log_info("Gene-level calls file generated: {call_cns_file}")

  # ── Log completion ────────────────────────────────────────────────────────
  log_info("CNVkit analysis completed for {sample_id}")
  log_info("CNR file (bin-level): {cnr_file}")
  log_info("CNS file (segments): {cns_file}")
  log_info("Call.CNS file (gene-level): {call_cns_file}")

  # Return list of output file paths
  invisible(list(
    cnr = as.character(cnr_file),
    cns = as.character(cns_file),
    call_cns = as.character(call_cns_file)
  ))
}
