#!/usr/bin/env Rscript
# benchmarks/scripts/05_evidence_concordance.R
# Compare OncoKB vs CiVIC evidence for top COSMIC variants across tumor types.

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

source(here::here("R/utils.R"))
source(here::here("R/api_clients.R"))
source(here::here("R/civic_client.R"))

# Load .env for OncoKB API key
load_env()

log_info("=== OncoKB vs CiVIC Evidence Concordance ===")

# --- Configuration -----------------------------------------------------------

results_dir <- here::here("benchmarks/results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# --- Top 100 COSMIC variants curated across 5 tumor types --------------------
# 20 well-known oncogenic variants per tumor type

cosmic_variants <- tribble(
  ~gene,    ~variant,   ~tumor_type,
  # NSCLC (Non-Small Cell Lung Cancer) — 20 variants
  "EGFR",   "L858R",    "NSCLC",
  "EGFR",   "T790M",    "NSCLC",
  "EGFR",   "C797S",    "NSCLC",
  "EGFR",   "G719S",    "NSCLC",
  "EGFR",   "L861Q",    "NSCLC",
  "EGFR",   "S768I",    "NSCLC",
  "KRAS",   "G12C",     "NSCLC",
  "KRAS",   "G12V",     "NSCLC",
  "KRAS",   "G12D",     "NSCLC",
  "ALK",    "F1174L",   "NSCLC",
  "ROS1",   "G2032R",   "NSCLC",
  "BRAF",   "V600E",    "NSCLC",
  "MET",    "D1228N",   "NSCLC",
  "MET",    "Y1003F",   "NSCLC",
  "ERBB2",  "S310F",    "NSCLC",
  "PIK3CA", "E545K",    "NSCLC",
  "STK11",  "Q37*",     "NSCLC",
  "TP53",   "R175H",    "NSCLC",
  "TP53",   "R248W",    "NSCLC",
  "RET",    "M918T",    "NSCLC",

  # Colorectal Cancer (COAD) — 20 variants
  "KRAS",   "G12D",     "COAD",
  "KRAS",   "G12V",     "COAD",
  "KRAS",   "G13D",     "COAD",
  "KRAS",   "A146T",    "COAD",
  "KRAS",   "Q61H",     "COAD",
  "NRAS",   "Q61K",     "COAD",
  "NRAS",   "Q61R",     "COAD",
  "BRAF",   "V600E",    "COAD",
  "PIK3CA", "E545K",    "COAD",
  "PIK3CA", "H1047R",   "COAD",
  "APC",    "R1450*",   "COAD",
  "TP53",   "R175H",    "COAD",
  "TP53",   "R273H",    "COAD",
  "SMAD4",  "R361H",    "COAD",
  "PTEN",   "R130Q",    "COAD",
  "FBXW7",  "R465C",    "COAD",
  "ERBB2",  "S310F",    "COAD",
  "MET",    "T1010I",   "COAD",
  "POLE",   "P286R",    "COAD",
  "MSH2",   "A636P",    "COAD",

  # Melanoma — 20 variants
  "BRAF",   "V600E",    "MEL",
  "BRAF",   "V600K",    "MEL",
  "BRAF",   "V600R",    "MEL",
  "BRAF",   "K601E",    "MEL",
  "NRAS",   "Q61R",     "MEL",
  "NRAS",   "Q61K",     "MEL",
  "NRAS",   "Q61L",     "MEL",
  "KIT",    "D816V",    "MEL",
  "KIT",    "L576P",    "MEL",
  "KIT",    "V559D",    "MEL",
  "MAP2K1", "P124S",    "MEL",
  "MAP2K1", "Q56P",     "MEL",
  "CTNNB1", "S37F",     "MEL",
  "GNA11",  "Q209L",    "MEL",
  "GNAQ",   "Q209P",    "MEL",
  "IDH1",   "R132H",    "MEL",
  "TP53",   "R248W",    "MEL",
  "CDKN2A", "R58*",     "MEL",
  "NF1",    "R1947*",   "MEL",
  "PTEN",   "R130*",    "MEL",

  # Breast Cancer — 20 variants
  "PIK3CA", "H1047R",   "BRCA",
  "PIK3CA", "E545K",    "BRCA",
  "PIK3CA", "E542K",    "BRCA",
  "PIK3CA", "H1047L",   "BRCA",
  "PIK3CA", "N345K",    "BRCA",
  "AKT1",   "E17K",     "BRCA",
  "ERBB2",  "S310F",    "BRCA",
  "ERBB2",  "L755S",    "BRCA",
  "ERBB2",  "V777L",    "BRCA",
  "ESR1",   "D538G",    "BRCA",
  "ESR1",   "Y537S",    "BRCA",
  "TP53",   "R175H",    "BRCA",
  "TP53",   "Y220C",    "BRCA",
  "BRCA1",  "C61G",     "BRCA",
  "BRCA2",  "E3002K",   "BRCA",
  "CDH1",   "A617T",    "BRCA",
  "GATA3",  "R330fs",   "BRCA",
  "MAP3K1", "E303*",    "BRCA",
  "PTEN",   "R130Q",    "BRCA",
  "KRAS",   "G12D",     "BRCA",

  # Ovarian Cancer — 20 variants
  "BRCA1",  "C61G",     "OVT",
  "BRCA1",  "E23fs",    "OVT",
  "BRCA2",  "E3002K",   "OVT",
  "BRCA2",  "S1982fs",  "OVT",
  "TP53",   "R175H",    "OVT",
  "TP53",   "R248W",    "OVT",
  "TP53",   "R273H",    "OVT",
  "TP53",   "Y220C",    "OVT",
  "PIK3CA", "H1047R",   "OVT",
  "PIK3CA", "E545K",    "OVT",
  "KRAS",   "G12D",     "OVT",
  "KRAS",   "G12V",     "OVT",
  "PTEN",   "R130Q",    "OVT",
  "ARID1A", "Q563*",    "OVT",
  "CTNNB1", "S37F",     "OVT",
  "BRAF",   "V600E",    "OVT",
  "ERBB2",  "S310F",    "OVT",
  "NF1",    "R1947*",   "OVT",
  "CDK12",  "K756fs",   "OVT",
  "RB1",    "R445*",    "OVT"
)

log_info("Curated list: {nrow(cosmic_variants)} variants across {length(unique(cosmic_variants$tumor_type))} tumor types")

# --- Query OncoKB for each variant --------------------------------------------

log_info("Querying OncoKB for {nrow(cosmic_variants)} variants...")

query_oncokb_safe <- function(gene, variant, tumor_type) {
  tryCatch({
    result <- oncokb_annotate_mutation(
      hugo_symbol    = gene,
      protein_change = variant,
      tumor_type     = tumor_type
    )
    tibble(
      has_oncokb_evidence = !is.na(result$oncogenic) && result$oncogenic != "Unknown",
      oncokb_oncogenic    = result$oncogenic %||% NA_character_,
      oncokb_level        = result$highest_sensitive_level %||% NA_character_,
      oncokb_mutation_effect = result$mutation_effect %||% NA_character_
    )
  }, error = function(e) {
    log_debug("OncoKB failed for {gene} {variant}: {e$message}")
    tibble(
      has_oncokb_evidence = FALSE,
      oncokb_oncogenic    = NA_character_,
      oncokb_level        = NA_character_,
      oncokb_mutation_effect = NA_character_
    )
  })
}

oncokb_results <- map_dfr(seq_len(nrow(cosmic_variants)), function(i) {
  row <- cosmic_variants[i, ]
  if (i %% 10 == 0) log_info("OncoKB: queried {i}/{nrow(cosmic_variants)}")
  query_oncokb_safe(row$gene, row$variant, row$tumor_type)
})

# --- Query CiVIC for each variant ---------------------------------------------

log_info("Querying CiVIC for {nrow(cosmic_variants)} variants...")

query_civic_safe <- function(gene, variant) {
  tryCatch({
    search <- civic_search_variant(gene, variant)
    if (nrow(search) == 0) {
      return(tibble(
        has_civic_evidence = FALSE,
        civic_evidence_count = 0L,
        civic_variant_id = NA_integer_
      ))
    }
    best <- search[1, ]
    evidence_count <- best$evidence_count %||% 0L

    # Get detailed evidence if molecular profile is available
    if (!is.na(best$molecular_profile_id)) {
      evidence <- tryCatch(
        civic_get_evidence(best$molecular_profile_id),
        error = function(e) tibble()
      )
      evidence_count <- max(evidence_count, nrow(evidence))
    }

    tibble(
      has_civic_evidence = evidence_count > 0,
      civic_evidence_count = as.integer(evidence_count),
      civic_variant_id = as.integer(best$variant_id)
    )
  }, error = function(e) {
    log_debug("CiVIC failed for {gene} {variant}: {e$message}")
    tibble(
      has_civic_evidence = FALSE,
      civic_evidence_count = 0L,
      civic_variant_id = NA_integer_
    )
  })
}

civic_results <- map_dfr(seq_len(nrow(cosmic_variants)), function(i) {
  row <- cosmic_variants[i, ]
  if (i %% 10 == 0) log_info("CiVIC: queried {i}/{nrow(cosmic_variants)}")
  query_civic_safe(row$gene, row$variant)
})

# --- Combine results ----------------------------------------------------------

combined <- bind_cols(cosmic_variants, oncokb_results, civic_results)

log_info("Combined results: {nrow(combined)} variants")

# --- Overlap statistics -------------------------------------------------------

both_have      <- sum(combined$has_oncokb_evidence & combined$has_civic_evidence)
oncokb_only    <- sum(combined$has_oncokb_evidence & !combined$has_civic_evidence)
civic_only     <- sum(!combined$has_oncokb_evidence & combined$has_civic_evidence)
neither        <- sum(!combined$has_oncokb_evidence & !combined$has_civic_evidence)

log_info("Evidence overlap:")
log_info("  Both OncoKB + CiVIC: {both_have}")
log_info("  OncoKB only: {oncokb_only}")
log_info("  CiVIC only: {civic_only}")
log_info("  Neither: {neither}")

# Jaccard index
jaccard <- if ((both_have + oncokb_only + civic_only) > 0) {
  both_have / (both_have + oncokb_only + civic_only)
} else {
  NA_real_
}
log_info("Jaccard index: {round(jaccard, 3)}")

# --- CiVIC added value -------------------------------------------------------

civic_adds_value <- combined |>
  filter(has_civic_evidence & !has_oncokb_evidence)

if (nrow(civic_adds_value) > 0) {
  log_info("CiVIC adds unique evidence for {nrow(civic_adds_value)} variants:")
  for (i in seq_len(nrow(civic_adds_value))) {
    row <- civic_adds_value[i, ]
    log_info("  {row$gene} {row$variant} ({row$tumor_type}): {row$civic_evidence_count} CiVIC evidence items")
  }
}

# --- Per-tumor-type summary ---------------------------------------------------

tumor_summary <- combined |>
  group_by(tumor_type) |>
  summarise(
    n_variants           = n(),
    n_oncokb_evidence    = sum(has_oncokb_evidence),
    n_civic_evidence     = sum(has_civic_evidence),
    n_both               = sum(has_oncokb_evidence & has_civic_evidence),
    n_oncokb_only        = sum(has_oncokb_evidence & !has_civic_evidence),
    n_civic_only         = sum(!has_oncokb_evidence & has_civic_evidence),
    n_neither            = sum(!has_oncokb_evidence & !has_civic_evidence),
    oncokb_coverage_pct  = round(mean(has_oncokb_evidence) * 100, 1),
    civic_coverage_pct   = round(mean(has_civic_evidence) * 100, 1),
    .groups = "drop"
  )

cat("\n--- Per-Tumor-Type Summary ---\n")
print(tumor_summary, n = 10)

# --- Save results -------------------------------------------------------------

output_path <- file.path(results_dir, "evidence_concordance.csv")
write.csv(combined, output_path, row.names = FALSE)
log_info("Evidence concordance saved to: {output_path}")

summary_path <- file.path(results_dir, "evidence_concordance_summary.csv")
summary_overall <- tibble(
  metric = c("total_variants", "both_oncokb_civic", "oncokb_only",
             "civic_only", "neither", "jaccard_index",
             "oncokb_coverage_pct", "civic_coverage_pct"),
  value = c(
    nrow(combined), both_have, oncokb_only, civic_only, neither,
    round(jaccard, 4),
    round(mean(combined$has_oncokb_evidence) * 100, 1),
    round(mean(combined$has_civic_evidence) * 100, 1)
  )
)
write.csv(summary_overall, summary_path, row.names = FALSE)
log_info("Summary saved to: {summary_path}")

tumor_summary_path <- file.path(results_dir, "evidence_concordance_by_tumor.csv")
write.csv(tumor_summary, tumor_summary_path, row.names = FALSE)

log_info("=== Evidence Concordance Benchmark Complete ===")
