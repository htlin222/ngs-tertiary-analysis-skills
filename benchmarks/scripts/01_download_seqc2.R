#!/usr/bin/env Rscript
# benchmarks/scripts/01_download_seqc2.R
# Download SEQC2 high-confidence somatic variant truth set for HCC1395
# and filter to TSO500 panel regions.

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(here)
})

source(here::here("R/utils.R"))

log_info("=== SEQC2 Truth Set Download ===")

# --- Configuration -----------------------------------------------------------

seqc2_base_url <- paste0(

  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/",
  "Somatic_Mutation_WG/release/latest/"
)
truth_vcf_name <- "high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz"
truth_vcf_url  <- paste0(seqc2_base_url, truth_vcf_name)

output_dir    <- here::here("benchmarks/datasets/seqc2")
panel_bed     <- here::here("config/tso500_panel.bed")
raw_vcf_path  <- file.path(output_dir, truth_vcf_name)
filtered_vcf  <- file.path(output_dir, "seqc2_truth_tso500_filtered.vcf.gz")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Check dependencies ------------------------------------------------------

bcftools_available <- nchar(Sys.which("bcftools")) > 0
if (!bcftools_available) {
  log_warn("bcftools not found in PATH — VCF filtering will be skipped")
  log_warn("Install with: conda install -c bioconda bcftools")
}

if (!file.exists(panel_bed)) {
  log_error("Panel BED file not found: {panel_bed}")
  stop("Missing panel BED file at config/tso500_panel.bed")
}

# --- Download truth VCF -------------------------------------------------------

if (file.exists(raw_vcf_path)) {
  log_info("Truth VCF already downloaded: {raw_vcf_path}")
} else {
  log_info("Downloading SEQC2 truth VCF from NCBI...")
  log_info("URL: {truth_vcf_url}")

  download_result <- tryCatch({
    download.file(
      url      = truth_vcf_url,
      destfile = raw_vcf_path,
      mode     = "wb",
      quiet    = FALSE
    )
    TRUE
  }, error = function(e) {
    log_error("Download failed: {e$message}")
    FALSE
  })

  if (!download_result) {
    stop("Failed to download SEQC2 truth VCF")
  }

  # Download index if available
  tbi_url  <- paste0(truth_vcf_url, ".tbi")
  tbi_path <- paste0(raw_vcf_path, ".tbi")
  tryCatch({
    download.file(url = tbi_url, destfile = tbi_path, mode = "wb", quiet = TRUE)
    log_info("Downloaded VCF index")
  }, error = function(e) {
    log_debug("No .tbi index available; will index locally if needed")
  })

  log_info("Download complete: {raw_vcf_path}")
}

# --- Count total variants in truth set ---------------------------------------

if (bcftools_available) {
  total_count_cmd <- glue("bcftools view -H {raw_vcf_path} | wc -l")
  total_variants <- as.integer(trimws(system(total_count_cmd, intern = TRUE)))
  log_info("Total variants in SEQC2 truth set: {total_variants}")
} else {
  total_variants <- NA_integer_
  log_info("Skipping variant count (bcftools not available)")
}

# --- Filter to TSO500 panel regions -------------------------------------------

if (bcftools_available) {
  if (file.exists(filtered_vcf)) {
    log_info("Filtered VCF already exists: {filtered_vcf}")
  } else {
    log_info("Filtering truth VCF to TSO500 panel regions...")

    # Index the VCF if no .tbi exists
    tbi_path <- paste0(raw_vcf_path, ".tbi")
    if (!file.exists(tbi_path)) {
      log_info("Indexing truth VCF with bcftools...")
      index_cmd <- glue("bcftools index -t {raw_vcf_path}")
      index_result <- system(index_cmd)
      if (index_result != 0) {
        log_warn("bcftools index failed; trying tabix...")
        tabix_cmd <- glue("tabix -p vcf {raw_vcf_path}")
        system(tabix_cmd)
      }
    }

    filter_cmd <- glue(
      "bcftools view -R {panel_bed} {raw_vcf_path} ",
      "-Oz -o {filtered_vcf}"
    )
    filter_result <- system(filter_cmd)

    if (filter_result != 0) {
      log_error("bcftools filtering failed (exit code: {filter_result})")
      stop("VCF filtering failed")
    }

    # Index filtered VCF
    system(glue("bcftools index -t {filtered_vcf}"))
    log_info("Filtered VCF saved: {filtered_vcf}")
  }

  # Count filtered variants
  filtered_count_cmd <- glue("bcftools view -H {filtered_vcf} | wc -l")
  filtered_variants <- as.integer(trimws(system(filtered_count_cmd, intern = TRUE)))
  log_info("Variants in TSO500 panel regions: {filtered_variants}")
} else {
  filtered_variants <- NA_integer_
  log_warn("Skipping VCF filtering (bcftools not available)")
  log_warn("Install bcftools and re-run this script to filter the truth set")
}

# --- Summary ------------------------------------------------------------------

summary_df <- data.frame(
  metric = c("total_truth_variants", "panel_filtered_variants", "panel_bed"),
  value  = c(
    ifelse(is.na(total_variants), "N/A", as.character(total_variants)),
    ifelse(is.na(filtered_variants), "N/A", as.character(filtered_variants)),
    panel_bed
  ),
  stringsAsFactors = FALSE
)

summary_path <- file.path(output_dir, "download_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

log_info("Summary saved to: {summary_path}")
log_info("=== SEQC2 Download Complete ===")
