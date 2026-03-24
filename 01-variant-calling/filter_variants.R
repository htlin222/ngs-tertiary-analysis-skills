# 01-variant-calling/filter_variants.R — Variant filtering and refinement
# Applies GATK FilterMutectCalls followed by R-based filtering on VAF and allele counts

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(fs)
  library(VariantAnnotation)
  library(dplyr)
})

# Load shared utilities
source(here::here("R/utils.R"))

#' Filter VCF variants using GATK FilterMutectCalls and R-based criteria
#'
#' Applies two-stage filtering:
#' 1. GATK FilterMutectCalls to mark low-confidence calls
#' 2. R-based filtering by VAF and alternative allele read count
#'
#' @param vcf_path Character. Path to raw VCF file from Mutect2
#' @param config List. Pipeline configuration (from load_config())
#' @param sample_id Character. Sample identifier used for output naming
#'
#' @return Character. Path to filtered and compressed VCF file
#'
#' @details
#' Stage 1 (GATK FilterMutectCalls):
#'   - Flags variants as PASS or REJECT based on artifact detection
#'   - Considers strand bias, orientation bias, map quality
#'
#' Stage 2 (R-based filtering):
#'   - Retains only variants with FILTER=PASS
#'   - Minimum VAF >= config$variant_calling$min_vaf
#'   - Minimum alt reads >= config$variant_calling$min_alt_reads
#'   - Logs variant counts at each filtering step
#'
#' @examples
#' \dontrun{
#'   config <- load_config("config/default.yaml")
#'   filtered_vcf <- filter_variants(
#'     "reports/SAMPLE_001/01-variant-calling/mutect2_raw.vcf.gz",
#'     config,
#'     "SAMPLE_001"
#'   )
#' }
filter_variants <- function(vcf_path, config, sample_id) {
  log_info("Starting variant filtering for {sample_id}")

  # ── Input validation ──────────────────────────────────────────────────────
  if (!file.exists(vcf_path)) {
    log_error("Input VCF not found: {vcf_path}")
    stop(glue("Input VCF not found: {vcf_path}"))
  }

  if (is.null(config$reference$fasta) || !file.exists(config$reference$fasta)) {
    log_error("Reference FASTA not found: {config$reference$fasta}")
    stop("Reference FASTA is required")
  }

  # ── Check GATK availability ───────────────────────────────────────────────
  gatk_cmd <- "gatk"

  if (nchar(Sys.which(gatk_cmd)) == 0) {
    gatk_bundled <- here::here("tools/gatk-4.6.1.0/gatk")
    if (file.exists(gatk_bundled)) {
      gatk_cmd <- gatk_bundled
      log_info("Using bundled GATK: {gatk_cmd}")
    } else {
      log_error("GATK not found in PATH or {gatk_bundled}")
      stop("GATK is required but not available")
    }
  }

  # ── Prepare output paths ──────────────────────────────────────────────────
  output_dir <- stage_output_dir(sample_id, "01-variant-calling")
  gatk_filtered_vcf <- path(output_dir, "gatk_filtered.vcf.gz")
  final_vcf <- path(output_dir, "filtered.vcf.gz")

  log_info("GATK filtered VCF: {gatk_filtered_vcf}")
  log_info("Final filtered VCF: {final_vcf}")

  # ── Stage 1: Run GATK FilterMutectCalls ────────────────────────────────────
  log_info("Stage 1: Running GATK FilterMutectCalls")

  args <- c(
    "FilterMutectCalls",
    "--variant", vcf_path,
    "--reference", config$reference$fasta,
    "--output", gatk_filtered_vcf
  )

  run_tool(
    gatk_cmd,
    args,
    description = glue("GATK FilterMutectCalls for {sample_id}")
  )

  if (!file.exists(gatk_filtered_vcf)) {
    log_error("GATK FilterMutectCalls failed to produce output: {gatk_filtered_vcf}")
    stop("GATK FilterMutectCalls did not produce output")
  }

  log_info("GATK FilterMutectCalls completed")

  # ── Stage 2: R-based filtering ──────────────────────────────────────────────
  log_info("Stage 2: Applying R-based filtering criteria")

  # Read the GATK-filtered VCF
  vcf <- readVcf(gatk_filtered_vcf, "hg38")
  n_variants_input <- nrow(vcf)
  log_info("Input variants from GATK filtering: {n_variants_input}")

  # Extract FILTER field
  vcf_filters <- fixed(vcf)$FILTER
  log_debug("Filter status counts: {paste(table(vcf_filters), collapse = ', ')}")

  # Keep only PASS variants
  vcf_pass <- vcf[elementNROWS(vcf_filters) == 0 |
    sapply(vcf_filters, function(x) all(x == "PASS"))]
  n_variants_pass <- nrow(vcf_pass)
  log_info("Variants with FILTER=PASS: {n_variants_pass}")

  if (n_variants_pass == 0) {
    log_warn("No PASS variants found after GATK filtering")
    # Create empty VCF with same structure
    vcf_final <- vcf_pass
  } else {
    # Extract genotype information for VAF and alt allele count filtering
    gt_data <- geno(vcf_pass)

    # For Mutect2, we need to compute VAF from AD (allelic depth)
    # AD is a matrix with ref (column 1) and alt (column 2) alleles
    if (!is.null(gt_data$AD)) {
      ad <- gt_data$AD
      # Extract alt allele depth (second column for each sample)
      # Assuming single sample (tumor only)
      if (is.array(ad)) {
        alt_reads <- ad[, 2, 1, drop = TRUE]  # [variant, allele, sample]
      } else {
        # Fallback if AD is a list
        alt_reads <- sapply(ad, function(x) {
          if (is.matrix(x)) x[1, 2] else x[2]
        })
      }

      # Total depth
      if (!is.null(gt_data$DP)) {
        total_depth <- gt_data$DP[, 1, drop = TRUE]  # [variant, sample]
        # Compute VAF
        vaf <- ifelse(total_depth > 0, alt_reads / total_depth, 0)
      } else {
        # No DP field, use AD sum
        ref_reads <- ifelse(is.array(ad), ad[, 1, 1, drop = TRUE],
          sapply(ad, function(x) {
            if (is.matrix(x)) x[1, 1] else x[1]
          }))
        total_depth <- ref_reads + alt_reads
        vaf <- ifelse(total_depth > 0, alt_reads / total_depth, 0)
      }
    } else {
      log_warn("No AD (allelic depth) field found; using AF from INFO if available")
      info_data <- info(vcf_pass)
      if (!is.null(info_data$AF)) {
        vaf <- as.numeric(info_data$AF)
        alt_reads <- NA  # Cannot compute without read counts
      } else {
        log_error("Neither AD nor AF found in VCF")
        stop("Cannot compute VAF: no AD or AF field in VCF")
      }
    }

    # Apply VAF filter
    min_vaf <- config$variant_calling$min_vaf %||% 0.05
    mask_vaf <- vaf >= min_vaf
    n_vaf_pass <- sum(mask_vaf, na.rm = TRUE)
    log_info("Variants with VAF >= {min_vaf}: {n_vaf_pass}")

    # Apply minimum alt reads filter (if alt_reads is available)
    if (!all(is.na(alt_reads))) {
      min_alt_reads <- config$variant_calling$min_alt_reads %||% 5
      mask_alt <- alt_reads >= min_alt_reads
      mask_alt[is.na(mask_alt)] <- FALSE
      n_alt_pass <- sum(mask_alt, na.rm = TRUE)
      log_info("Variants with alt reads >= {min_alt_reads}: {n_alt_pass}")

      # Combine filters
      mask_final <- mask_vaf & mask_alt
    } else {
      log_warn("Alt reads not available; using only VAF filter")
      mask_final <- mask_vaf
    }

    n_variants_final <- sum(mask_final, na.rm = TRUE)
    log_info("Variants passing all filters: {n_variants_final}")
    log_info("Filtering summary: {n_variants_input} -> {n_variants_pass} (PASS) -> {n_variants_final} (final)")

    vcf_final <- vcf_pass[mask_final]
  }

  # ── Write final filtered VCF ─────────────────────────────────────────────
  writeVcf(vcf_final, final_vcf, index = TRUE)

  if (!file.exists(final_vcf)) {
    log_error("Failed to write final VCF: {final_vcf}")
    stop("Failed to write final filtered VCF")
  }

  log_info("Variant filtering completed for {sample_id}")
  log_info("Final VCF: {final_vcf}")

  invisible(as.character(final_vcf))
}
