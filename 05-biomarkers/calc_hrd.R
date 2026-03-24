suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(logger)
  library(jsonlite)
})

#' Calculate Homologous Recombination Deficiency (HRD)
#'
#' Estimates HRD score from CNV segment data by calculating three component scores:
#' Loss of Heterozygosity (LOH), Telomeric Allelic Imbalance (TAI), and
#' Large-scale State Transitions (LST). Note: This is a simplified estimation
#' from panel CNV data and has limitations compared to whole-genome approaches.
#'
#' @param cnv_results List or character. Either a list with CNV segment data
#'   (with elements 'segments' containing a data frame) or path to CNV results file.
#' @param config List. Configuration object containing biomarker settings.
#' @param sample_id Character. Sample identifier.
#'
#' @return List containing:
#'   - hrd_score: Numeric, total HRD score (LOH + TAI + LST)
#'   - hrd_status: Character, classification (HRD-positive or HRD-negative)
#'   - loh_score: Integer, LOH component score
#'   - tai_score: Integer, TAI component score
#'   - lst_score: Integer, LST component score
#'   - is_reliable: Logical, FALSE indicating this is panel-based estimation
#'
#' @details
#' Component scores are calculated as:
#' - LOH: Count of segments with loss of heterozygosity >15 Mb
#' - TAI: Count of segments with allelic imbalance extending to telomere
#' - LST: Count of breakpoints between segments with CN state changes and size >10 Mb
#'
#' @export
calc_hrd <- function(cnv_results, config, sample_id) {
  log_info("Starting HRD calculation for {sample_id}")

  # Load CNV segment data
  if (is.character(cnv_results)) {
    if (!file.exists(cnv_results)) {
      log_error("CNV results file not found: {cnv_results}")
      stop("CNV results file not found")
    }
    log_info("Loading CNV results from {cnv_results}")
    cnv_data <- readRDS(cnv_results)
  } else if (is.list(cnv_results)) {
    cnv_data <- cnv_results
  } else {
    log_error("cnv_results must be a list or file path")
    stop("Invalid cnv_results format")
  }

  # Extract segments
  if (!("segments" %in% names(cnv_data))) {
    log_error("CNV data missing 'segments' element")
    stop("CNV data missing 'segments' element")
  }

  segments <- as.data.frame(cnv_data$segments)

  # Ensure required columns exist
  required_cols <- c("chrom", "start", "end", "copy_number")
  if (!all(required_cols %in% names(segments))) {
    log_error("Segments missing required columns: {paste(required_cols, collapse = ', ')}")
    stop("Segments missing required columns")
  }

  # Extract configuration parameters
  loh_size_threshold_mb <- config$biomarkers$hrd$loh_size_threshold_mb %||% 15
  lst_size_threshold_mb <- config$biomarkers$hrd$lst_size_threshold_mb %||% 10
  hrd_score_threshold <- config$biomarkers$hrd$score_threshold %||% 42

  log_info("Calculating HRD components from {nrow(segments)} CNV segments")

  # Convert to proper data frame and ensure numeric columns
  segments <- segments %>%
    mutate(
      chrom = as.character(chrom),
      start = as.numeric(start),
      end = as.numeric(end),
      copy_number = as.numeric(copy_number),
      segment_size_mb = (end - start) / 1e6
    )

  # Calculate LOH (Loss of Heterozygosity) score
  # Count segments with LOH >15 Mb
  loh_score <- segments %>%
    filter(segment_size_mb >= loh_size_threshold_mb) %>%
    nrow()

  log_info("LOH score: {loh_score} segments >{loh_size_threshold_mb}Mb")

  # Calculate TAI (Telomeric Allelic Imbalance) score
  # Simplified: count segments with allelic imbalance (copy_number != 2) near telomeres
  # For panel data, use segments at chromosome ends
  telomere_proximity_mb <- 5  # Consider within 5Mb of telomere

  # Define approximate chromosome sizes for hg38 (in bp)
  chr_sizes <- list(
    "1" = 248e6, "2" = 242e6, "3" = 198e6, "4" = 190e6, "5" = 181e6,
    "6" = 171e6, "7" = 160e6, "8" = 146e6, "9" = 141e6, "10" = 135e6,
    "11" = 134e6, "12" = 133e6, "13" = 114e6, "14" = 107e6, "15" = 101e6,
    "16" = 98e6, "17" = 92e6, "18" = 90e6, "19" = 61e6, "20" = 64e6,
    "21" = 48e6, "22" = 50e6, "X" = 155e6, "Y" = 59e6
  )

  tai_score <- segments %>%
    filter(copy_number != 2) %>%
    mutate(
      is_near_telomere = (start < (telomere_proximity_mb * 1e6)) |
        (end > (chr_sizes[[chrom]] %||% Inf) - (telomere_proximity_mb * 1e6))
    ) %>%
    filter(is_near_telomere) %>%
    nrow()

  log_info("TAI score: {tai_score} allelic imbalance segments near telomeres")

  # Calculate LST (Large-scale State Transitions) score
  # Count breakpoints between segments with different CN states and >10 Mb
  lst_score <- 0

  for (i in seq_len(nrow(segments) - 1)) {
    curr_seg <- segments[i, ]
    next_seg <- segments[i + 1, ]

    # Only count transitions if on same chromosome
    if (curr_seg$chrom == next_seg$chrom) {
      # Check if CN state changes
      if (curr_seg$copy_number != next_seg$copy_number) {
        # Check if both segments >10 Mb
        if (curr_seg$segment_size_mb >= lst_size_threshold_mb &&
            next_seg$segment_size_mb >= lst_size_threshold_mb) {
          lst_score <- lst_score + 1
        }
      }
    }
  }

  log_info("LST score: {lst_score} large-scale state transitions (>{lst_size_threshold_mb}Mb)")

  # Calculate total HRD score
  hrd_score <- loh_score + tai_score + lst_score

  # Classify HRD status
  if (hrd_score >= hrd_score_threshold) {
    hrd_status <- "HRD-positive"
  } else {
    hrd_status <- "HRD-negative"
  }

  log_info("HRD calculation complete: total score {hrd_score} ({hrd_status})")
  log_info("Note: Panel-based HRD estimation has limited reliability compared to WGS")

  # Prepare output
  result <- list(
    hrd_score = hrd_score,
    hrd_status = hrd_status,
    loh_score = loh_score,
    tai_score = tai_score,
    lst_score = lst_score,
    is_reliable = FALSE,
    reliability_note = "This is a simplified panel-based HRD estimation. Accuracy is limited compared to whole-genome sequencing approaches. Results should be interpreted with caution and validated with orthogonal methods if critical for clinical decision-making."
  )

  # Save result
  output_dir <- stage_output_dir(sample_id, "05-biomarkers")
  output_file <- file.path(output_dir, "hrd_result.json")

  write_json(result, output_file, pretty = TRUE)
  log_info("Saved HRD result to {output_file}")

  return(result)
}
