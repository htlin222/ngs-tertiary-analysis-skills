suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(logger)
})

#' Run Ensembl VEP on Filtered VCF
#'
#' Annotates variants using Ensembl Variant Effect Predictor (VEP) with
#' comprehensive functional predictions, regulatory information, and allele
#' frequencies.
#'
#' @param vcf_path Character. Path to filtered VCF file from stage 1.
#' @param config List. Configuration with annotation parameters:
#'   - annotation$vep_cache_dir: Path to VEP cache directory
#'   - annotation$vep_assembly: Assembly (default: GRCh38)
#'   - annotation$vep_species: Species (default: homo_sapiens)
#' @param sample_id Character. Sample identifier for output organization.
#'
#' @return Character. Path to VEP-annotated VCF file (gzipped).
#'
#' @details
#' VEP parameters:
#' - --everything: Includes SIFT, PolyPhen, regulatory, conservation scores
#' - --pick: One consequence per variant (first in canonical order)
#' - --hgvs: HGVS nomenclature (c. and p. notation)
#' - --symbol: Gene symbol annotation
#' - --canonical: Mark canonical transcript
#' - --biotype: Transcript biotype
#' - --af_gnomad: Add gnomAD allele frequencies
#' - --offline: Use local cache (requires pre-downloaded data)
#' - --force_overwrite: Overwrite existing output
#'
#' @export
run_annotation <- function(vcf_path, config, sample_id) {
  log_info("Starting VEP annotation for sample: {sample_id}")

  # Validate input
  if (!file.exists(vcf_path)) {
    log_error("Input VCF not found: {vcf_path}")
    stop("Input VCF not found: ", vcf_path)
  }

  # Create output directory
  output_dir <- stage_output_dir(sample_id, "02-annotation")
  fs::dir_create(output_dir, recurse = TRUE)

  # Define output path
  output_vcf <- file.path(output_dir, "vep_annotated.vcf")
  output_vcf_gz <- paste0(output_vcf, ".gz")

  # Construct VEP command
  vep_cmd <- glue(
    "vep \\",
    "  --input_file {vcf_path} \\",
    "  --output_file {output_vcf} \\",
    "  --format vcf \\",
    "  --vcf \\",
    "  --cache \\",
    "  --dir_cache {config$annotation$vep_cache_dir} \\",
    "  --assembly {config$annotation$vep_assembly} \\",
    "  --species {config$annotation$vep_species} \\",
    "  --offline \\",
    "  --force_overwrite \\",
    "  --everything \\",
    "  --pick \\",
    "  --hgvs \\",
    "  --symbol \\",
    "  --canonical \\",
    "  --biotype \\",
    "  --fork {config$annotation$vep_fork %||% 4} \\",
    "  --af_gnomad",
    .sep = "\n"
  )

  log_info("Running VEP with command:\n{vep_cmd}")

  # Execute VEP
  result <- run_tool(
    command = "bash",
    args = c("-c", vep_cmd),
    description = glue("VEP annotation for {sample_id}")
  )

  if (result$status != 0) {
    log_error("VEP annotation failed with status {result$status}")
    log_error("stderr: {result$stderr}")
    stop("VEP annotation failed")
  }

  log_info("VEP annotation completed")

  # Compress output VCF
  log_info("Compressing annotated VCF")
  compress_cmd <- glue("bgzip -f {output_vcf}")
  compress_result <- run_tool(
    command = "bash",
    args = c("-c", compress_cmd),
    description = glue("Compressing VEP output for {sample_id}")
  )

  if (compress_result$status != 0) {
    log_error("VCF compression failed")
    stop("VCF compression failed")
  }

  # Index compressed VCF
  log_info("Indexing annotated VCF")
  index_cmd <- glue("tabix -p vcf {output_vcf_gz}")
  index_result <- run_tool(
    command = "bash",
    args = c("-c", index_cmd),
    description = glue("Indexing VEP output for {sample_id}")
  )

  if (index_result$status != 0) {
    log_warning("VCF indexing failed, continuing without index")
  }

  # Verify output
  if (!file.exists(output_vcf_gz)) {
    log_error("Output file not created: {output_vcf_gz}")
    stop("Output file not created")
  }

  log_info("VEP annotation complete: {output_vcf_gz}")
  return(output_vcf_gz)
}
