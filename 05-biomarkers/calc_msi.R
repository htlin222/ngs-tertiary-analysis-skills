suppressPackageStartupMessages({
  library(tidyverse)
  library(Rsamtools)
  library(glue)
  library(logger)
  library(jsonlite)
})

#' Default Bethesda Microsatellite Markers
#'
#' Returns a data frame with coordinates of commonly used microsatellite loci
#' for MSI detection (Bethesda panel).
#'
#' @return Data frame with columns: chr, start, end, marker_name
#'
#' @keywords internal
get_default_bethesda_markers <- function() {
  tibble(
    chr = c("2", "2", "2", "17", "17"),
    start = c(47641559, 21067052, 95744450, 41196311, 41197819),
    end = c(47641590, 21067083, 95744481, 41196349, 41197854),
    marker_name = c("BAT25", "BAT26", "D5S346", "D17S250", "D17S261")
  )
}

#' Calculate Microsatellite Instability (MSI)
#'
#' Detects microsatellite instability by analyzing repeat length variation
#' at defined microsatellite loci from BAM alignment data.
#'
#' @param bam_path Character. Path to BAM file.
#' @param config List. Configuration object containing biomarker settings.
#' @param sample_id Character. Sample identifier.
#'
#' @return List containing:
#'   - msi_score: Numeric, percentage of unstable sites
#'   - msi_status: Character, classification (MSI-H, MSI-L, MSS)
#'   - unstable_sites: Integer, count of unstable microsatellites
#'   - total_sites: Integer, count of analyzed sites
#'
#' @export
calc_msi <- function(bam_path, config, sample_id) {
  log_info("Starting MSI calculation for {sample_id}")

  # Check BAM file exists
  if (!file.exists(bam_path)) {
    log_error("BAM file not found: {bam_path}")
    stop("BAM file not found")
  }

  # Load microsatellite sites
  sites_bed <- config$biomarkers$msi$sites_bed
  if (!is.null(sites_bed) && file.exists(sites_bed)) {
    log_info("Loading MSI sites from {sites_bed}")
    msi_sites <- read_tsv(sites_bed, col_names = c("chr", "start", "end", "marker_name"),
                          col_types = "ciic")
  } else {
    log_info("Using default Bethesda microsatellite markers")
    msi_sites <- get_default_bethesda_markers()
  }

  # Extract configuration parameters
  msi_h_threshold <- config$biomarkers$msi$msi_h_threshold %||% 20
  min_depth <- config$biomarkers$msi$min_depth %||% 30

  log_info("Analyzing {nrow(msi_sites)} microsatellite sites")

  # Initialize BAM file connection
  bam_file <- BamFile(bam_path, index = paste0(bam_path, ".bai"))

  # Analyze each site
  unstable_count <- 0
  analyzed_count <- 0

  for (i in seq_len(nrow(msi_sites))) {
    site <- msi_sites[i, ]
    marker_name <- site$marker_name

    tryCatch({
      # Extract reads at this locus
      param <- ScanBamParam(
        which = GRanges(
          seqnames = site$chr,
          ranges = IRanges(start = site$start, end = site$end)
        ),
        what = c("seq", "cigar")
      )

      bam_data <- scanBam(bam_file, param = param)[[1]]

      if (length(bam_data$seq) == 0) {
        log_warn("No reads found for {marker_name} ({site$chr}:{site$start}-{site$end})")
        next
      }

      # Check depth
      depth <- length(bam_data$seq)
      if (depth < min_depth) {
        log_warn("Insufficient depth at {marker_name} ({depth}x, min: {min_depth}x)")
        next
      }

      # Analyze repeat length variation
      # Extract lengths from CIGAR strings and sequences
      cigar_lengths <- nchar(bam_data$seq)
      expected_length <- site$end - site$start + 1

      # Calculate length distribution
      length_table <- table(cigar_lengths)
      length_distribution <- as.numeric(length_table) / sum(length_table)

      # Detect instability: if multiple different lengths present with >10% frequency
      major_lengths <- names(length_table)[length_distribution >= 0.1]

      is_unstable <- FALSE
      if (length(major_lengths) >= 2) {
        length_diffs <- as.numeric(major_lengths)
        max_diff <- max(length_diffs) - min(length_diffs)
        if (max_diff > 0) {
          is_unstable <- TRUE
        }
      }

      if (is_unstable) {
        unstable_count <- unstable_count + 1
        log_debug("{marker_name} is unstable (repeat length variation detected)")
      } else {
        log_debug("{marker_name} is stable")
      }

      analyzed_count <- analyzed_count + 1
    },
    error = function(e) {
      log_warn("Error analyzing {marker_name}: {e$message}")
    })
  }

  close(bam_file)

  # Calculate MSI score and classify
  if (analyzed_count == 0) {
    log_warn("No microsatellite sites could be analyzed")
    msi_score <- NA_real_
    msi_status <- "Undetermined"
  } else {
    msi_score <- (unstable_count / analyzed_count) * 100

    if (msi_score >= msi_h_threshold) {
      msi_status <- "MSI-H"
    } else if (msi_score > 0) {
      msi_status <- "MSI-L"
    } else {
      msi_status <- "MSS"
    }
  }

  log_info("MSI calculation complete: {unstable_count}/{analyzed_count} unstable sites")
  log_info("MSI score: {round(msi_score, 2)}% ({msi_status})")

  # Prepare output
  result <- list(
    msi_score = round(msi_score, 2),
    msi_status = msi_status,
    unstable_sites = unstable_count,
    total_sites = analyzed_count
  )

  # Save result
  output_dir <- stage_output_dir(sample_id, "05-biomarkers")
  output_file <- file.path(output_dir, "msi_result.json")

  write_json(result, output_file, pretty = TRUE)
  log_info("Saved MSI result to {output_file}")

  return(result)
}
