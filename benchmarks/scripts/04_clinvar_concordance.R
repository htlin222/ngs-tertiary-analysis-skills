#!/usr/bin/env Rscript
# benchmarks/scripts/04_clinvar_concordance.R
# Compare AMP classification against ClinVar annotations.
# Calculate Cohen's kappa for concordance.

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
})

source(here::here("R/utils.R"))
source(here::here("R/amp_classification.R"))
source(here::here("R/api_clients.R"))

# Load .env for OncoKB API key
load_env()

log_info("=== ClinVar Concordance Benchmark ===")

# --- Configuration -----------------------------------------------------------

results_dir <- here::here("benchmarks/results")
clinvar_dir <- here::here("benchmarks/results/clinvar")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(clinvar_dir)) dir.create(clinvar_dir, recursive = TRUE)

clinvar_url <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
clinvar_vcf_path <- file.path(clinvar_dir, "clinvar.vcf.gz")
panel_bed <- here::here("config/tso500_panel.bed")

# --- TSO500 panel genes (commonly included on TSO500) -------------------------

tso500_genes <- c(
  "ABL1", "AKT1", "ALK", "APC", "AR", "ARAF", "ATM", "ATR", "ATRX",
  "BAP1", "BRAF", "BRCA1", "BRCA2", "CDH1", "CDK4", "CDK6", "CDKN2A",
  "CHEK2", "CTNNB1", "DDR2", "EGFR", "ERBB2", "ERBB3", "ERBB4", "ESR1",
  "EZH2", "FGFR1", "FGFR2", "FGFR3", "FLT3", "FOXL2", "GNA11", "GNAQ",
  "GNAS", "HNF1A", "HRAS", "IDH1", "IDH2", "JAK2", "JAK3", "KIT",
  "KRAS", "MAP2K1", "MAP2K2", "MET", "MLH1", "MPL", "MSH2", "MSH6",
  "MTOR", "MYC", "MYCN", "NF1", "NF2", "NOTCH1", "NPM1", "NRAS",
  "NTRK1", "NTRK2", "NTRK3", "PALB2", "PDGFRA", "PIK3CA", "PMS2",
  "POLE", "PTEN", "PTCH1", "RAF1", "RB1", "RET", "ROS1", "SMAD4",
  "SMO", "STK11", "TERT", "TP53", "TSC1", "TSC2", "VHL"
)

# --- Download ClinVar VCF ----------------------------------------------------

bcftools_available <- nchar(Sys.which("bcftools")) > 0

if (file.exists(clinvar_vcf_path)) {
  log_info("ClinVar VCF already downloaded: {clinvar_vcf_path}")
} else {
  log_info("Downloading ClinVar VCF...")
  tryCatch({
    download.file(url = clinvar_url, destfile = clinvar_vcf_path,
                  mode = "wb", quiet = FALSE)
    # Download index
    download.file(url = paste0(clinvar_url, ".tbi"),
                  destfile = paste0(clinvar_vcf_path, ".tbi"),
                  mode = "wb", quiet = TRUE)
    log_info("ClinVar VCF downloaded")
  }, error = function(e) {
    log_error("ClinVar download failed: {e$message}")
    stop("Failed to download ClinVar VCF")
  })
}

# --- Parse ClinVar VCF -------------------------------------------------------

log_info("Parsing ClinVar VCF for TSO500 panel genes...")

parse_clinvar <- function(vcf_path, gene_filter = NULL) {
  # Use bcftools to extract relevant fields efficiently
  if (bcftools_available) {
    log_info("Using bcftools to extract ClinVar variants")

    # Extract: CHROM, POS, REF, ALT, CLNSIG, CLNVC, GENEINFO, MC
    cmd <- glue(
      "bcftools query -f ",
      "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CLNSIG\\t%INFO/CLNVC\\t%INFO/GENEINFO\\t%INFO/MC\\n' ",
      "{vcf_path}"
    )

    lines <- tryCatch(
      system(cmd, intern = TRUE),
      error = function(e) {
        log_warn("bcftools query failed: {e$message}")
        return(character(0))
      }
    )

    if (length(lines) == 0) {
      log_warn("No variants extracted from ClinVar VCF")
      return(tibble())
    }

    fields <- strsplit(lines, "\t")
    variants <- tibble(
      chr      = sapply(fields, `[`, 1),
      pos      = as.integer(sapply(fields, `[`, 2)),
      ref      = sapply(fields, `[`, 3),
      alt      = sapply(fields, `[`, 4),
      clnsig   = sapply(fields, `[`, 5),
      clnvc    = sapply(fields, `[`, 6),
      geneinfo = sapply(fields, `[`, 7),
      mc       = sapply(fields, `[`, 8)
    )
  } else {
    # Fallback: read VCF with VariantAnnotation
    log_info("Using VariantAnnotation to parse ClinVar")
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
      stop("Either bcftools or VariantAnnotation is required to parse ClinVar")
    }

    vcf <- VariantAnnotation::readVcf(vcf_path, genome = "GRCh38")
    info_df <- VariantAnnotation::info(vcf)
    rd <- SummarizedExperiment::rowRanges(vcf)

    variants <- tibble(
      chr      = as.character(GenomicRanges::seqnames(rd)),
      pos      = GenomicRanges::start(rd),
      ref      = as.character(rd$REF),
      alt      = as.character(unlist(rd$ALT)),
      clnsig   = as.character(info_df$CLNSIG),
      clnvc    = as.character(info_df$CLNVC %||% NA),
      geneinfo = as.character(info_df$GENEINFO %||% NA),
      mc       = as.character(info_df$MC %||% NA)
    )
  }

  # Extract gene symbol from GENEINFO field (e.g., "BRAF:673")
  variants <- variants |>
    mutate(gene = str_extract(geneinfo, "^[^:]+"))

  # Filter to panel genes
  if (!is.null(gene_filter)) {
    variants <- variants |> filter(gene %in% gene_filter)
  }

  # Normalize clinical significance
  variants <- variants |>
    mutate(
      clnsig_norm = case_when(
        str_detect(tolower(clnsig), "pathogenic") &
          !str_detect(tolower(clnsig), "likely") &
          !str_detect(tolower(clnsig), "conflicting") ~ "pathogenic",
        str_detect(tolower(clnsig), "likely_pathogenic") ~ "likely_pathogenic",
        str_detect(tolower(clnsig), "benign") &
          !str_detect(tolower(clnsig), "likely") &
          !str_detect(tolower(clnsig), "conflicting") ~ "benign",
        str_detect(tolower(clnsig), "likely_benign") ~ "likely_benign",
        str_detect(tolower(clnsig), "uncertain_significance") ~ "uncertain_significance",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(clnsig_norm))

  # Map MC (molecular consequence) to VEP-like impact
  variants <- variants |>
    mutate(
      vep_impact = case_when(
        str_detect(tolower(mc), "nonsense|frameshift|splice_donor|splice_acceptor") ~ "HIGH",
        str_detect(tolower(mc), "missense") ~ "MODERATE",
        str_detect(tolower(mc), "synonymous|intron") ~ "LOW",
        TRUE ~ "MODIFIER"
      )
    )

  # Extract protein change where possible (e.g., from MC or variant name)
  variants <- variants |>
    mutate(
      protein_change = str_extract(mc, "p\\.[A-Za-z0-9]+"),
      key = paste(chr, pos, ref, alt, sep = ":")
    )

  log_info("Parsed {nrow(variants)} ClinVar variants in panel genes")
  variants
}

clinvar_variants <- parse_clinvar(clinvar_vcf_path, gene_filter = tso500_genes)

log_info("ClinVar variants with definitive classifications: {nrow(clinvar_variants)}")
log_info("Breakdown: {paste(names(table(clinvar_variants$clnsig_norm)), table(clinvar_variants$clnsig_norm), sep='=', collapse=', ')}")

# --- Run AMP classification on ClinVar variants -------------------------------

log_info("Running AMP classification on ClinVar variants...")

# Limit to a manageable subset for OncoKB queries (rate-limited)
max_oncokb_queries <- as.integer(Sys.getenv("BENCHMARK_MAX_ONCOKB", unset = "500"))
clinvar_subset <- clinvar_variants

if (nrow(clinvar_subset) > max_oncokb_queries) {
  log_info("Sampling {max_oncokb_queries} variants for OncoKB annotation (set BENCHMARK_MAX_ONCOKB to change)")
  set.seed(42)
  per_group <- max_oncokb_queries %/% 5
  clinvar_subset <- clinvar_subset |>
    group_by(clnsig_norm) |>
    slice_head(n = per_group) |>
    ungroup()
}

classify_clinvar_variant <- function(row) {
  gene           <- row$gene
  protein_change <- row$protein_change
  clinvar_sig    <- row$clnsig_norm
  vep_impact     <- row$vep_impact

  # Query OncoKB if protein change is available
  oncokb_oncogenic <- NA_character_
  oncokb_level     <- NA_character_

  if (!is.na(protein_change) && !is.na(gene)) {
    tryCatch({
      oncokb <- oncokb_annotate_mutation(
        hugo_symbol    = gene,
        protein_change = protein_change,
        tumor_type     = "CANCER"  # generic tumor type
      )
      oncokb_oncogenic <- oncokb$oncogenic
      oncokb_level     <- oncokb$highest_sensitive_level
    }, error = function(e) {
      log_debug("OncoKB annotation failed for {gene} {protein_change}: {e$message}")
    })
  }

  # Run AMP classification
  amp <- classify_amp_tier(
    oncokb_oncogenic = oncokb_oncogenic,
    oncokb_level     = oncokb_level,
    civic_amp_level  = NA_character_,
    vep_impact       = vep_impact,
    clinvar          = clinvar_sig
  )

  tibble(
    gene             = gene,
    protein_change   = protein_change,
    clinvar_sig      = clinvar_sig,
    vep_impact       = vep_impact,
    oncokb_oncogenic = oncokb_oncogenic,
    oncokb_level     = oncokb_level,
    amp_tier         = amp$amp_tier,
    amp_level        = amp$amp_level,
    amp_evidence     = amp$amp_evidence
  )
}

classification_results <- map_dfr(
  seq_len(nrow(clinvar_subset)),
  function(i) {
    if (i %% 50 == 0) log_info("Classified {i}/{nrow(clinvar_subset)} variants")
    classify_clinvar_variant(clinvar_subset[i, ])
  }
)

log_info("AMP classification complete for {nrow(classification_results)} variants")

# --- Expected mapping: ClinVar -> AMP tier ------------------------------------

# ClinVar pathogenic -> Tier I or Tier II
# ClinVar benign/likely_benign -> Tier IV
# ClinVar VUS -> Tier III

classification_results <- classification_results |>
  mutate(
    expected_tier = case_when(
      clinvar_sig %in% c("pathogenic", "likely_pathogenic") ~ "Tier_I_or_II",
      clinvar_sig %in% c("benign", "likely_benign") ~ "Tier_IV",
      clinvar_sig == "uncertain_significance" ~ "Tier_III",
      TRUE ~ NA_character_
    ),
    amp_tier_group = case_when(
      amp_tier %in% c("Tier I", "Tier II") ~ "Tier_I_or_II",
      amp_tier == "Tier III" ~ "Tier_III",
      amp_tier == "Tier IV" ~ "Tier_IV",
      TRUE ~ "Unclassified"
    ),
    concordant = (expected_tier == amp_tier_group)
  )

# --- Calculate Cohen's Kappa -------------------------------------------------

calc_cohens_kappa <- function(observed, expected) {
  # Create confusion matrix
  levels <- sort(unique(c(observed, expected)))
  conf_mat <- table(factor(observed, levels = levels),
                    factor(expected, levels = levels))

  n <- sum(conf_mat)
  if (n == 0) return(NA_real_)

  # Observed agreement
  p_o <- sum(diag(conf_mat)) / n

  # Expected agreement by chance
  row_sums <- rowSums(conf_mat)
  col_sums <- colSums(conf_mat)
  p_e <- sum(row_sums * col_sums) / (n^2)

  # Kappa
  if (p_e == 1) return(1.0)
  kappa <- (p_o - p_e) / (1 - p_e)
  kappa
}

kappa <- calc_cohens_kappa(
  classification_results$amp_tier_group,
  classification_results$expected_tier
)

concordance_rate <- mean(classification_results$concordant, na.rm = TRUE)

log_info("Concordance rate: {round(concordance_rate * 100, 1)}%")
log_info("Cohen's kappa: {round(kappa, 3)}")

# --- Concordance matrix -------------------------------------------------------

concordance_matrix <- table(
  AMP_Tier = classification_results$amp_tier_group,
  ClinVar  = classification_results$expected_tier
)

cat("\n--- Concordance Matrix ---\n")
print(concordance_matrix)

# --- Save results -------------------------------------------------------------

output_path <- file.path(results_dir, "clinvar_concordance.csv")
write.csv(classification_results, output_path, row.names = FALSE)
log_info("Concordance results saved to: {output_path}")

# Save summary
summary_df <- tibble(
  metric = c("total_variants", "concordance_rate", "cohens_kappa",
             "pathogenic_as_tier_i_ii", "benign_as_tier_iv", "vus_as_tier_iii"),
  value = c(
    nrow(classification_results),
    round(concordance_rate, 4),
    round(kappa, 4),
    sum(classification_results$clinvar_sig %in% c("pathogenic", "likely_pathogenic") &
          classification_results$amp_tier_group == "Tier_I_or_II", na.rm = TRUE) /
      max(1, sum(classification_results$clinvar_sig %in% c("pathogenic", "likely_pathogenic"))),
    sum(classification_results$clinvar_sig %in% c("benign", "likely_benign") &
          classification_results$amp_tier_group == "Tier_IV", na.rm = TRUE) /
      max(1, sum(classification_results$clinvar_sig %in% c("benign", "likely_benign"))),
    sum(classification_results$clinvar_sig == "uncertain_significance" &
          classification_results$amp_tier_group == "Tier_III", na.rm = TRUE) /
      max(1, sum(classification_results$clinvar_sig == "uncertain_significance"))
  )
)

summary_path <- file.path(results_dir, "clinvar_concordance_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)
log_info("Summary saved to: {summary_path}")

log_info("=== ClinVar Concordance Benchmark Complete ===")
