# 01-variant-calling/run_mutect2.R — GATK Mutect2 variant calling
# Performs somatic variant calling on tumor BAM files using GATK Mutect2

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(fs)
})

# Load shared utilities
source(here::here("R/utils.R"))

#' Run GATK Mutect2 for somatic variant calling
#'
#' Performs somatic variant calling on a tumor BAM file using GATK Mutect2.
#' Integrates optional germline and panel-of-normals resources when available.
#'
#' @param bam_path Character. Path to tumor BAM file (indexed)
#' @param config List. Pipeline configuration (from load_config())
#' @param sample_id Character. Sample identifier used for output naming
#'
#' @return Character. Path to output VCF file (gzip compressed)
#'
#' @details
#' The function:
#' 1. Checks GATK availability (via PATH or tools/gatk-4.6.1.0/)
#' 2. Constructs command with:
#'    - BAM input and reference genome
#'    - Interval restriction to panel regions
#'    - Optional germline and PoN resources
#'    - Extra arguments from config
#'    - Tumor sample identifier
#' 3. Runs command via run_tool() with error checking
#' 4. Returns path to compressed output VCF
#'
#' @examples
#' \dontrun{
#'   config <- load_config("config/default.yaml")
#'   vcf <- run_mutect2(
#'     "samples/tumor.bam",
#'     config,
#'     "SAMPLE_001"
#'   )
#' }
run_mutect2 <- function(bam_path, config, sample_id) {
  log_info("Starting Mutect2 variant calling for {sample_id}")

  # ── Input validation ──────────────────────────────────────────────────────
  if (!file.exists(bam_path)) {
    log_error("BAM file not found: {bam_path}")
    stop(glue("BAM file not found: {bam_path}"))
  }

  if (is.null(config$reference$fasta) || !file.exists(config$reference$fasta)) {
    log_error("Reference FASTA not found: {config$reference$fasta}")
    stop("Reference FASTA is required and must exist")
  }

  # ── Check GATK availability ───────────────────────────────────────────────
  gatk_cmd <- "gatk"

  # Try PATH first
  if (nchar(Sys.which(gatk_cmd)) == 0) {
    # Try bundled location
    gatk_bundled <- here::here("tools/gatk-4.6.1.0/gatk")
    if (file.exists(gatk_bundled)) {
      gatk_cmd <- gatk_bundled
      log_info("Using bundled GATK: {gatk_cmd}")
    } else {
      log_error("GATK not found in PATH or {gatk_bundled}")
      stop("GATK is required but not available")
    }
  } else {
    log_debug("Using GATK from PATH")
  }

  # ── Prepare output directory and paths ─────────────────────────────────────
  output_dir <- stage_output_dir(sample_id, "01-variant-calling")
  output_vcf <- path(output_dir, "mutect2_raw.vcf.gz")

  log_info("Output VCF: {output_vcf}")

  # ── Build Mutect2 command ─────────────────────────────────────────────────
  args <- c("Mutect2")

  # Required arguments
  args <- c(args, "--input", bam_path)
  args <- c(args, "--reference", config$reference$fasta)

  # Panel of regions (restrict to sequenced regions for efficiency)
  if (!is.null(config$reference$panel_bed) &&
    file.exists(config$reference$panel_bed)) {
    args <- c(args, "--intervals", config$reference$panel_bed)
    log_info("Restricting to panel regions: {config$reference$panel_bed}")
  }

  # Output
  args <- c(args, "--output", output_vcf)

  # Optional: Germline resource (population allele frequency)
  if (!is.null(config$variant_calling$gnomad_resource) &&
    file.exists(config$variant_calling$gnomad_resource)) {
    args <- c(args, "--germline-resource", config$variant_calling$gnomad_resource)
    log_info("Using gnomAD resource: {config$variant_calling$gnomad_resource}")
  }

  # Optional: Panel of normals
  if (!is.null(config$variant_calling$panel_of_normals) &&
    file.exists(config$variant_calling$panel_of_normals)) {
    args <- c(args, "--panel-of-normals", config$variant_calling$panel_of_normals)
    log_info("Using panel of normals: {config$variant_calling$panel_of_normals}")
  }

  # Tumor sample name (required for proper annotation)
  args <- c(args, "--tumor-sample", sample_id)

  # Extra arguments from config
  if (!is.null(config$variant_calling$mutect2_extra_args)) {
    extra_args <- strsplit(config$variant_calling$mutect2_extra_args, "\\s+")[[1]]
    args <- c(args, extra_args)
    log_debug("Added extra arguments: {paste(extra_args, collapse = ' ')}")
  }

  # ── PairHMM threading ──────────────────────────────────────────────────
  hmm_threads <- config$variant_calling$mutect2_threads %||% 4
  args <- c(args, "--native-pair-hmm-threads", as.character(hmm_threads))
  log_info("Using {hmm_threads} PairHMM threads")

  # ── Run Mutect2 ───────────────────────────────────────────────────────────
  run_tool(
    gatk_cmd,
    args,
    description = glue("GATK Mutect2 for {sample_id}")
  )

  # ── Verify output ─────────────────────────────────────────────────────────
  if (!file.exists(output_vcf)) {
    log_error("Output VCF was not created: {output_vcf}")
    stop("Mutect2 failed to produce output VCF")
  }

  log_info("Mutect2 variant calling completed for {sample_id}")
  log_info("Output VCF: {output_vcf}")

  invisible(as.character(output_vcf))
}
