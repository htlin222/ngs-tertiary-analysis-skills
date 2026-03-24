suppressPackageStartupMessages({
  library(tidyverse)
  library(VariantAnnotation)
  library(glue)
  library(logger)
})

#' Parse VEP CSQ Field
#'
#' Parses the VEP consequence field (CSQ) into a tidy data frame with
#' individual consequence records per variant-transcript pair.
#'
#' @param csq Character vector. VEP CSQ field values (pipe-delimited).
#'
#' @return Tibble with parsed consequences and annotations.
#'
#' @keywords internal
parse_vep_csq <- function(csq) {
  if (is.na(csq) || csq == "") {
    return(tibble())
  }

  # Split by comma to get individual consequences
  consequences <- str_split(csq, ",")[[1]]

  # Parse each consequence entry
  parsed <- map_df(consequences, function(cons) {
    fields <- str_split(cons, "\\|")[[1]]

    # VEP CSQ header order (standard with --everything flag):
    # Consequence, SYMBOL, HGVSC, HGVSP, IMPACT, SIFT, PolyPhen,
    # gnomAD_AF, ClinVar, COSMIC, and others...
    tibble(
      consequence = fields[1] %||% NA_character_,
      gene = fields[2] %||% NA_character_,
      hgvsc = fields[3] %||% NA_character_,
      hgvsp = fields[4] %||% NA_character_,
      impact = fields[5] %||% NA_character_,
      sift_prediction = fields[6] %||% NA_character_,
      polyphen_prediction = fields[7] %||% NA_character_,
      gnomad_af = fields[8] %||% NA_character_,
      clinvar_significance = fields[9] %||% NA_character_,
      cosmic_id = fields[10] %||% NA_character_
    )
  })

  return(parsed)
}

#' Classify Variant Pathogenicity
#'
#' Assigns pathogenicity classification based on ClinVar significance,
#' functional predictions (SIFT, PolyPhen), and impact level.
#'
#' @param clinvar_sig Character. ClinVar significance.
#' @param impact Character. VEP impact level.
#' @param sift Character. SIFT prediction.
#' @param polyphen Character. PolyPhen prediction.
#'
#' @return Character. Classification: pathogenic, likely_pathogenic, vus,
#'   likely_benign, or benign.
#'
#' @keywords internal
classify_pathogenicity <- function(clinvar_sig, impact, sift, polyphen) {
  # ClinVar-based classification (highest priority)
  if (!is.na(clinvar_sig)) {
    if (str_detect(tolower(clinvar_sig), "pathogenic")) {
      return("pathogenic")
    }
    if (str_detect(tolower(clinvar_sig), "likely_pathogenic")) {
      return("likely_pathogenic")
    }
    if (str_detect(tolower(clinvar_sig), "benign")) {
      return("benign")
    }
    if (str_detect(tolower(clinvar_sig), "likely_benign")) {
      return("likely_benign")
    }
  }

  # Functional prediction-based classification
  high_impact <- !is.na(impact) && impact == "HIGH"
  sift_deleterious <- !is.na(sift) && str_detect(tolower(sift), "deleterious")
  polyphen_probably_damaging <- !is.na(polyphen) &&
    str_detect(tolower(polyphen), "probably_damaging")

  # High confidence pathogenic
  if (high_impact && (sift_deleterious || polyphen_probably_damaging)) {
    return("likely_pathogenic")
  }

  # Moderate predictions
  moderate_impact <- !is.na(impact) && impact == "MODERATE"
  if (moderate_impact && (sift_deleterious || polyphen_probably_damaging)) {
    return("vus")
  }

  # Default to VUS if any functional prediction exists
  if (sift_deleterious || polyphen_probably_damaging) {
    return("vus")
  }

  # High impact without functional evidence
  if (high_impact) {
    return("likely_pathogenic")
  }

  # Default to VUS
  return("vus")
}

#' Merge VEP Annotations
#'
#' Reads VEP-annotated VCF and merges with sample-specific information
#' (VAF, depths) to create a comprehensive variant annotation table.
#'
#' @param annotated_vcf Character. Path to VEP-annotated VCF file.
#' @param config List. Configuration with filtering parameters:
#'   - annotation$gnomad_af_threshold: Maximum gnomAD AF to retain (default: 0.01)
#' @param sample_id Character. Sample identifier for output organization.
#'
#' @return Tibble with merged annotations including:
#'   - Variant identifiers (CHROM, POS, REF, ALT)
#'   - VEP predictions (gene, impact, consequence, HGVS notation)
#'   - Functional scores (SIFT, PolyPhen, gnomAD AF)
#'   - ClinVar and COSMIC information
#'   - Pathogenicity classification
#'   - Sample-specific metrics (VAF, depth)
#'
#' @details
#' Filtering steps:
#' 1. Parse VEP CSQ field into individual consequence records
#' 2. Select canonical/highest impact consequence per variant
#' 3. Extract VAF and depth from FORMAT fields (GT, AD, DP)
#' 4. Filter common germline variants (gnomAD AF > threshold)
#' 5. Classify pathogenicity
#' 6. Save to TSV and return tibble
#'
#' @export
merge_annotations <- function(annotated_vcf, config, sample_id) {
  log_info("Merging VEP annotations for sample: {sample_id}")

  # Validate input
  if (!file.exists(annotated_vcf)) {
    log_error("Annotated VCF not found: {annotated_vcf}")
    stop("Annotated VCF not found: ", annotated_vcf)
  }

  # Read annotated VCF
  log_info("Reading VEP-annotated VCF: {annotated_vcf}")
  vcf <- readVcf(annotated_vcf, genome = "GRCh38")

  # Extract basic variant information
  variants <- tibble(
    chrom = as.character(seqnames(vcf)),
    pos = start(vcf),
    ref = as.character(ref(vcf)),
    alt = as.character(unlist(alt(vcf))),
    qual = qual(vcf),
    filter = as.character(filt(vcf))
  )

  # Extract VEP CSQ field
  info_data <- info(vcf)
  csq_field <- if ("CSQ" %in% names(info_data)) {
    info_data$CSQ
  } else {
    log_warning("CSQ field not found in VCF info")
    rep(list(character()), nrow(vcf))
  }

  # Parse CSQ for each variant
  log_info("Parsing VEP consequence field")
  annotations <- map_df(seq_along(csq_field), function(i) {
    csq <- csq_field[[i]]
    if (length(csq) == 0 || is.na(csq)) {
      return(tibble(
        variant_idx = i,
        consequence = NA_character_,
        gene = NA_character_,
        hgvsc = NA_character_,
        hgvsp = NA_character_,
        impact = NA_character_,
        sift_prediction = NA_character_,
        polyphen_prediction = NA_character_,
        gnomad_af = NA_character_,
        clinvar_significance = NA_character_,
        cosmic_id = NA_character_
      ))
    }

    # Parse each consequence
    parsed <- parse_vep_csq(csq[1])
    if (nrow(parsed) == 0) {
      return(tibble(
        variant_idx = i,
        consequence = NA_character_,
        gene = NA_character_,
        hgvsc = NA_character_,
        hgvsp = NA_character_,
        impact = NA_character_,
        sift_prediction = NA_character_,
        polyphen_prediction = NA_character_,
        gnomad_af = NA_character_,
        clinvar_significance = NA_character_,
        cosmic_id = NA_character_
      ))
    }

    # Select canonical/highest impact consequence
    parsed <- parsed %>%
      arrange(
        desc(impact %in% c("HIGH", "MODERATE")),
        impact
      ) %>%
      slice(1) %>%
      mutate(variant_idx = i)

    return(parsed)
  })

  # Extract genotype information (VAF, depths)
  log_info("Extracting genotype information")
  genotype_info <- map_df(seq_len(nrow(vcf)), function(i) {
    gt_matrix <- geno(vcf)$GT[i, , drop = FALSE]
    ad_matrix <- geno(vcf)$AD[i, , drop = FALSE]
    dp_matrix <- geno(vcf)$DP[i, , drop = FALSE]

    # Use first sample
    gt <- gt_matrix[1, 1]
    ad <- ad_matrix[1, 1]
    dp <- dp_matrix[1, 1]

    # Calculate VAF from allele depth
    vaf <- NA_real_
    alt_depth <- NA_integer_
    total_depth <- NA_integer_

    if (!is.na(dp) && !is.null(dp)) {
      total_depth <- as.integer(dp)
    }

    if (!is.null(ad) && length(ad) >= 2) {
      ref_depth <- ad[1]
      alt_depth <- ad[2]
      if (!is.na(ref_depth) && !is.na(alt_depth)) {
        total <- ref_depth + alt_depth
        if (total > 0) {
          vaf <- alt_depth / total
        }
      }
    }

    tibble(
      variant_idx = i,
      vaf = vaf,
      alt_depth = alt_depth,
      total_depth = total_depth
    )
  })

  # Combine all data
  merged <- bind_cols(
    variants %>% mutate(variant_idx = row_number()),
    annotations %>% select(-variant_idx),
    genotype_info %>% select(-variant_idx)
  ) %>%
    select(-variant_idx)

  # Parse gnomAD AF (extract numeric value from CSQ format)
  log_info("Processing allele frequencies")
  merged <- merged %>%
    mutate(
      gnomad_af_numeric = if_else(
        is.na(gnomad_af) | gnomad_af == "",
        NA_real_,
        as.numeric(gnomad_af)
      )
    )

  # Filter common germline variants
  gnomad_threshold <- config$annotation$gnomad_af_threshold %||% 0.01
  log_info(
    "Filtering variants with gnomAD AF > {gnomad_threshold}: ",
    "keeping {nrow(merged)} variants"
  )

  merged <- merged %>%
    filter(
      is.na(gnomad_af_numeric) | gnomad_af_numeric <= gnomad_threshold
    )

  log_info(
    "After gnomAD filtering: {nrow(merged)} variants retained"
  )

  # Classify pathogenicity
  log_info("Classifying variant pathogenicity")
  merged <- merged %>%
    mutate(
      pathogenic_class = pmap_chr(
        list(
          clinvar_significance,
          impact,
          sift_prediction,
          polyphen_prediction
        ),
        classify_pathogenicity
      )
    )

  # Clean and reorder columns
  merged <- merged %>%
    select(
      chrom, pos, ref, alt, qual, filter,
      gene, hgvsc, hgvsp, consequence, impact,
      sift_prediction, polyphen_prediction,
      gnomad_af, gnomad_af_numeric,
      clinvar_significance, cosmic_id,
      pathogenic_class,
      vaf, alt_depth, total_depth
    ) %>%
    arrange(chrom, pos)

  # Create output directory
  output_dir <- stage_output_dir(sample_id, "02-annotation")
  fs::dir_create(output_dir, recurse = TRUE)

  # Save merged annotations
  output_file <- file.path(output_dir, "merged_annotations.tsv")
  log_info("Saving merged annotations to {output_file}")

  write_tsv(merged, output_file)

  log_info(
    "Annotation merging complete: {nrow(merged)} variants, ",
    "saved to {output_file}"
  )

  return(merged)
}
