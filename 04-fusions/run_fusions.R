suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(logger)
})

#' Run Manta Structural Variant Caller for Gene Fusion Detection
#'
#' Executes Manta SV caller on a tumor BAM file to detect structural variants
#' that may represent gene fusions. Falls back to TSO500 Local App output if
#' Manta is unavailable.
#'
#' @param bam_path Character. Path to input BAM file (must be indexed).
#' @param config List. Configuration object containing reference genome paths
#'   and analysis parameters.
#' @param sample_id Character. Sample identifier for output organization.
#'
#' @return Character. Path to the output VCF file containing candidate SVs
#'   (either candidateSV.vcf.gz from Manta or tumorSV.vcf.gz).
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input BAM file exists and is indexed
#' 2. Prepares panel BED file (bgzip and tabix if needed)
#' 3. Runs configManta.py to configure the analysis
#' 4. Executes runWorkflow.py locally
#' 5. Falls back to TSO500 fusion output if Manta fails
#' 6. Validates and returns the output VCF path
#'
#' @keywords internal
run_fusions <- function(bam_path, config, sample_id) {
  log_info("Stage 4: Gene Fusion Detection - Starting Manta SV calling")

  # Validate input BAM
  if (!file.exists(bam_path)) {
    log_error("BAM file not found: {bam_path}")
    stop("Input BAM file does not exist")
  }

  bam_index <- paste0(bam_path, ".bai")
  if (!file.exists(bam_index)) {
    log_error("BAM index not found: {bam_index}")
    stop("BAM file must be indexed (.bai)")
  }

  log_debug("Input BAM validated: {bam_path}")

  # Setup output directory
  output_dir <- stage_output_dir(sample_id, "04-fusions")
  run_dir <- file.path(output_dir, "manta_run")
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  log_info("Output directory: {output_dir}")

  # Prepare reference FASTA
  ref_fasta <- config$reference$fasta
  if (!file.exists(ref_fasta)) {
    log_error("Reference FASTA not found: {ref_fasta}")
    stop("Reference FASTA file does not exist")
  }

  # Prepare panel BED file - bgzip and tabix if needed
  panel_bed <- config$fusions$panel_bed
  if (!is.null(panel_bed) && file.exists(panel_bed)) {
    panel_bed_gz <- paste0(panel_bed, ".gz")
    panel_bed_tbi <- paste0(panel_bed_gz, ".tbi")

    if (!file.exists(panel_bed_gz)) {
      log_info("Compressing panel BED file with bgzip")
      cmd_bgzip <- glue("bgzip -c {panel_bed} > {panel_bed_gz}")
      result <- run_tool(cmd_bgzip, tool_name = "bgzip", sample_id = sample_id)
      if (!result$success) {
        log_warn("Failed to bgzip panel BED, attempting to run Manta without --callRegions")
        panel_bed_gz <- NULL
      }
    }

    if (!is.null(panel_bed_gz) && file.exists(panel_bed_gz) && !file.exists(panel_bed_tbi)) {
      log_info("Indexing bgzipped BED file with tabix")
      cmd_tabix <- glue("tabix -p bed {panel_bed_gz}")
      result <- run_tool(cmd_tabix, tool_name = "tabix", sample_id = sample_id)
      if (!result$success) {
        log_warn("Failed to tabix panel BED")
        panel_bed_gz <- NULL
      }
    }
  } else {
    panel_bed_gz <- NULL
    log_info("No panel BED file configured, running Manta without --callRegions")
  }

  # Configure Manta
  log_info("Configuring Manta")
  config_cmd <- glue(
    "configManta.py \\
      --tumorBam {bam_path} \\
      --referenceFasta {ref_fasta} \\
      --runDir {run_dir}"
  )

  if (!is.null(panel_bed_gz) && file.exists(panel_bed_gz)) {
    config_cmd <- glue("{config_cmd} --callRegions {panel_bed_gz}")
  }

  config_result <- run_tool(config_cmd, tool_name = "configManta.py", sample_id = sample_id)
  if (!config_result$success) {
    log_error("Manta configuration failed")
    return(.fallback_tso500_fusions(output_dir, sample_id))
  }

  # Execute Manta workflow
  log_info("Executing Manta workflow")
  workflow_script <- file.path(run_dir, "runWorkflow.py")
  if (!file.exists(workflow_script)) {
    log_error("Manta workflow script not found: {workflow_script}")
    return(.fallback_tso500_fusions(output_dir, sample_id))
  }

  workflow_cmd <- glue("{workflow_script} -m local")
  workflow_result <- run_tool(workflow_cmd, tool_name = "runWorkflow.py", sample_id = sample_id)

  if (!workflow_result$success) {
    log_error("Manta workflow execution failed")
    return(.fallback_tso500_fusions(output_dir, sample_id))
  }

  # Locate output VCF
  results_dir <- file.path(run_dir, "results", "variants")
  candidate_sv <- file.path(results_dir, "candidateSV.vcf.gz")
  tumor_sv <- file.path(results_dir, "tumorSV.vcf.gz")

  output_vcf <- NULL
  if (file.exists(candidate_sv)) {
    output_vcf <- candidate_sv
    log_info("Found Manta candidateSV output: {candidate_sv}")
  } else if (file.exists(tumor_sv)) {
    output_vcf <- tumor_sv
    log_info("Found Manta tumorSV output: {tumor_sv}")
  } else {
    log_error("No Manta VCF output found in {results_dir}")
    return(.fallback_tso500_fusions(output_dir, sample_id))
  }

  # Validate output VCF
  if (!file.exists(output_vcf)) {
    log_error("Output VCF file not found or inaccessible: {output_vcf}")
    return(.fallback_tso500_fusions(output_dir, sample_id))
  }

  log_success("Gene Fusion Detection (Manta) completed successfully")
  log_info("Output VCF: {output_vcf}")

  return(output_vcf)
}

#' Fallback to TSO500 Fusion Output
#'
#' If Manta is unavailable or fails, attempt to use TSO500 Local App
#' fusion output CSV if present.
#'
#' @param output_dir Character. Output directory path.
#' @param sample_id Character. Sample identifier.
#'
#' @return Character. Path to TSO500 fusion CSV or NULL if not found.
#'
#' @keywords internal
.fallback_tso500_fusions <- function(output_dir, sample_id) {
  log_info("Attempting fallback to TSO500 Local App fusion output")

  # Common TSO500 output paths
  tso500_paths <- c(
    file.path(output_dir, "fusions.csv"),
    file.path(dirname(output_dir), "tso500_fusions.csv"),
    file.path(dirname(dirname(output_dir)), "TSO500", "fusions.csv")
  )

  for (path in tso500_paths) {
    if (file.exists(path)) {
      log_success("Found TSO500 fusion output: {path}")
      return(path)
    }
  }

  log_error("No Manta output and no TSO500 fallback available")
  return(NULL)
}
