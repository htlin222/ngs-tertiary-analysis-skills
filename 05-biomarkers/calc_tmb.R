suppressPackageStartupMessages({
  library(tidyverse)
  library(VariantAnnotation)
  library(glue)
  library(logger)
  library(jsonlite)
})

#' Calculate Tumor Mutational Burden (TMB)
#'
#' Calculates TMB score from a filtered VCF file, counting somatic mutations
#' eligible for TMB calculation based on variant type and VAF thresholds.
#'
#' @param vcf_path Character. Path to filtered VCF file.
#' @param config List. Configuration object containing biomarker settings.
#' @param sample_id Character. Sample identifier.
#'
#' @return List containing:
#'   - tmb_score: Numeric, mutations per megabase
#'   - tmb_class: Character, classification (TMB-High, TMB-Intermediate, TMB-Low)
#'   - variant_count: Integer, count of eligible variants
#'   - panel_size_mb: Numeric, panel coding region size in Mb
#'
#' @export
calc_tmb <- function(vcf_path, config, sample_id) {
  log_info("Starting TMB calculation for {sample_id}")

  # Read VCF file
  if (!file.exists(vcf_path)) {
    log_error("VCF file not found: {vcf_path}")
    stop("VCF file not found")
  }

  vcf <- readVcf(vcf_path, genome = NA)
  log_info("Loaded VCF with {length(vcf)} variants")

  # Extract configuration parameters
  exclude_synonymous <- config$biomarkers$tmb$exclude_synonymous %||% FALSE
  exclude_known_drivers <- config$biomarkers$tmb$exclude_known_drivers %||% FALSE
  min_vaf_for_tmb <- config$biomarkers$tmb$min_vaf_for_tmb %||% 0.01
  panel_coding_size_mb <- config$biomarkers$tmb$panel_coding_size_mb %||% 1.94
  known_drivers <- config$biomarkers$tmb$known_driver_genes %||% c()

  # Initialize eligible variant count
  eligible_count <- 0

  # Process each variant
  for (i in seq_len(length(vcf))) {
    var <- vcf[i]

    # Get VAF if available in sample data
    vaf <- NA
    if ("AF" %in% names(info(var))) {
      vaf <- info(var)$AF[1]
    } else if ("AD" %in% names(geno(var))) {
      ad <- geno(var)$AD[1, 1, ]
      if (length(ad) >= 2 && sum(ad) > 0) {
        vaf <- ad[2] / sum(ad)
      }
    }

    # Check VAF threshold
    if (!is.na(vaf) && vaf < min_vaf_for_tmb) {
      next
    }

    # Check if non-synonymous (if configured)
    is_nonsynonymous <- TRUE
    if (exclude_synonymous) {
      csq_field <- NULL
      if ("CSQ" %in% names(info(var))) {
        csq_field <- info(var)$CSQ
      } else if ("ANN" %in% names(info(var))) {
        csq_field <- info(var)$ANN
      }

      if (!is.null(csq_field)) {
        csq_str <- as.character(csq_field)[1]
        is_nonsynonymous <- !grepl("synonymous", csq_str, ignore.case = TRUE)
      }
    }

    if (!is_nonsynonymous) {
      next
    }

    # Check if driver gene (if configured)
    is_driver <- FALSE
    if (exclude_known_drivers && length(known_drivers) > 0) {
      gene_field <- NULL
      if ("SYMBOL" %in% names(info(var))) {
        gene <- info(var)$SYMBOL[1]
        if (!is.na(gene) && gene %in% known_drivers) {
          is_driver <- TRUE
        }
      }
    }

    if (is_driver) {
      next
    }

    eligible_count <- eligible_count + 1
  }

  # Calculate TMB score
  tmb_score <- eligible_count / panel_coding_size_mb

  # Classify TMB
  if (tmb_score >= 10) {
    tmb_class <- "TMB-High"
  } else if (tmb_score >= 6) {
    tmb_class <- "TMB-Intermediate"
  } else {
    tmb_class <- "TMB-Low"
  }

  log_info("TMB calculation complete: {eligible_count} eligible variants")
  log_info("TMB score: {round(tmb_score, 2)} mut/Mb ({tmb_class})")

  # Prepare output
  result <- list(
    tmb_score = round(tmb_score, 2),
    tmb_class = tmb_class,
    variant_count = eligible_count,
    panel_size_mb = panel_coding_size_mb
  )

  # Save result
  output_dir <- stage_output_dir(sample_id, "05-biomarkers")
  output_file <- file.path(output_dir, "tmb_result.json")

  write_json(result, output_file, pretty = TRUE)
  log_info("Saved TMB result to {output_file}")

  return(result)
}
