#!/usr/bin/env Rscript
# run_e2e_test.R — End-to-end pipeline test with ovarian cancer mock data
#
# Simulates an ovarian cancer TSO500 case (matching EGA EGAD50000001451 cohort)
# Uses real API calls to OncoKB, PubMed, Scopus
# Produces a final ESMO-compliant HTML report
#
# Usage: Rscript run_e2e_test.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(glue)
  library(fs)
  library(yaml)
  library(jsonlite)
  library(httr2)
  library(logger)
  library(here)
})

# Source pipeline utilities
source(here("R/utils.R"))
source(here("R/api_clients.R"))
source(here("R/esmo_helpers.R"))

# ── Setup ────────────────────────────────────────────────────────────────────

setup_logging(level = "INFO")
load_env()
config <- load_config()

SAMPLE_ID <- "EGA_OV_TSO500_001"
TUMOR_TYPE <- "HGSOC"  # High-grade serous ovarian carcinoma

# Override config for this test
config$sample$id <- SAMPLE_ID
config$sample$tumor_type <- TUMOR_TYPE

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  NGS Tertiary Analysis Pipeline — End-to-End Test          ║\n")
cat("║  Sample: EGA_OV_TSO500_001 (Ovarian Cancer, TSO500)       ║\n")
cat("║  Cohort: EGAD50000001451 (8 ovarian tumor samples)        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 0: Mock QC Results
# ══════════════════════════════════════════════════════════════════════════════

log_info("=== Stage 0: QC (mock) ===")
out_qc <- stage_output_dir(SAMPLE_ID, "00-qc")

qc_results <- list(
  summary = list(
    total_reads = 45200000L,
    mapped_reads = 44100000L,
    mapping_rate = 0.976,
    mean_coverage = 487.3,
    median_coverage = 512.0,
    on_target_rate = 0.92,
    duplicate_rate = 0.15,
    tumor_purity = 0.65,
    insert_size_median = 168L
  ),
  per_gene_coverage = tibble(
    gene = c("TP53", "BRCA1", "BRCA2", "PIK3CA", "KRAS", "NRAS", "BRAF",
             "PTEN", "CDK12", "RB1", "NF1", "CCNE1", "MYC", "EGFR",
             "ERBB2", "ATM", "PALB2", "RAD51C", "RAD51D", "ARID1A",
             "MLH1", "MSH2", "MSH6", "PMS2", "POLE", "TERT", "STK11",
             "APC", "SMAD4", "VHL", "KIT", "PDGFRA", "IDH1", "IDH2",
             "FGFR1", "FGFR2", "FGFR3", "RET", "ALK", "ROS1",
             "NTRK1", "NTRK2", "NTRK3", "MET", "CDK4", "CDK6",
             "CDKN2A", "MDM2", "MYCN", "MAP2K1"),
    mean_coverage = c(523, 498, 412, 501, 545, 510, 489, 476, 445, 502,
                      388, 534, 467, 512, 498, 421, 456, 389, 401, 467,
                      505, 489, 478, 445, 312, 534, 498, 467, 445, 512,
                      489, 476, 501, 498, 445, 423, 456, 478, 489, 467,
                      445, 401, 389, 512, 498, 478, 456, 423, 401, 445),
    min_coverage = round(mean_coverage * 0.6),
    pct_above_200x = c(99.2, 98.8, 95.1, 99.0, 99.5, 98.9, 98.5, 97.8, 96.2, 99.1,
                        92.3, 99.3, 97.1, 99.0, 98.7, 95.5, 96.8, 92.0, 93.5, 97.1,
                        99.0, 98.5, 97.8, 96.2, 78.4, 99.3, 98.7, 97.1, 96.2, 99.0,
                        98.5, 97.8, 99.0, 98.7, 96.2, 95.0, 96.8, 97.8, 98.5, 97.1,
                        96.2, 93.5, 92.0, 99.0, 98.7, 97.8, 96.8, 95.0, 93.5, 96.2)
  ),
  pass = TRUE
)

# Save QC
write.csv(
  tibble(metric = names(qc_results$summary), value = unlist(qc_results$summary)),
  file.path(out_qc, "qc_summary.tsv"), row.names = FALSE
)
write.csv(qc_results$per_gene_coverage, file.path(out_qc, "per_gene_coverage.tsv"),
          row.names = FALSE)
log_info("QC results saved to {out_qc}")

# ══════════════════════════════════════════════════════════════════════════════
# STAGES 1-2: Mock Somatic Variants (realistic ovarian cancer mutations)
# Based on known HGSOC molecular landscape from TCGA + literature
# ══════════════════════════════════════════════════════════════════════════════

log_info("=== Stages 1-2: Variants (mock from known HGSOC landscape) ===")
out_ann <- stage_output_dir(SAMPLE_ID, "02-annotation")

merged_annotations <- tibble(
  gene = c("TP53", "BRCA1", "PIK3CA", "NF1", "RB1", "CDK12", "PTEN",
           "ARID1A", "KRAS", "CCNE1"),
  # Ensembl canonical transcript IDs (GRCh38)
  transcript_id = c("ENST00000269305.9", "ENST00000357654.9", "ENST00000263967.4",
                     "ENST00000358273.9", "ENST00000267163.5", "ENST00000447079.7",
                     "ENST00000371953.8", "ENST00000324856.14", "ENST00000256078.10", NA),
  hgvsp = c("p.R248W", "p.E1836fs", "p.H1047R", "p.R1534*", "p.E137*",
            "p.K765fs", "p.R130*", "p.Q1334*", "p.G12V", "p.E1"),
  hgvsc = c("c.742C>T", "c.5506del", "c.3140A>G", "c.4600C>T", "c.409G>T",
            "c.2293del", "c.388C>T", "c.4000C>T", "c.35G>T", NA),
  consequence = c("missense_variant", "frameshift_variant", "missense_variant",
                  "stop_gained", "stop_gained", "frameshift_variant",
                  "stop_gained", "stop_gained", "missense_variant", "amplification"),
  impact = c("HIGH", "HIGH", "MODERATE", "HIGH", "HIGH", "HIGH",
             "HIGH", "HIGH", "MODERATE", "HIGH"),
  vaf = c(0.68, 0.42, 0.23, 0.31, 0.15, 0.38, 0.45, 0.12, 0.08, NA),
  alt_depth = c(156L, 89L, 52L, 67L, 34L, 82L, 98L, 28L, 18L, NA_integer_),
  total_depth = c(229L, 212L, 226L, 216L, 227L, 216L, 218L, 233L, 225L, NA_integer_),
  gnomad_af = c(0.00002, 0, 0.00001, 0, 0, 0, 0, 0, 0, NA),
  clinvar_significance = c("Pathogenic", "Pathogenic", "Pathogenic",
                           "Pathogenic", "Likely_pathogenic", "Pathogenic",
                           "Pathogenic", "Likely_pathogenic", "Pathogenic", NA),
  cosmic_id = c("COSM10656", "COSM19602", "COSM775", "COSM1485887",
                NA, NA, "COSM5154", NA, "COSM520", NA),
  sift_prediction = c("deleterious", NA, "deleterious", NA, NA, NA,
                      NA, NA, "deleterious", NA),
  polyphen_prediction = c("probably_damaging", NA, "probably_damaging", NA, NA, NA,
                          NA, NA, "probably_damaging", NA),
  pathogenicity = c("pathogenic", "pathogenic", "pathogenic", "pathogenic",
                    "likely_pathogenic", "pathogenic", "pathogenic",
                    "likely_pathogenic", "pathogenic", "pathogenic")
)

# Mutation confidence ranking per Table 3 guidelines
# Based on: mapping quality, VAF vs detection threshold, reproducibility (ClinVar/COSMIC)
detection_threshold <- config$variant_calling$min_vaf  # 0.05

merged_annotations <- merged_annotations |>
  mutate(
    confidence = case_when(
      # High confidence: high VAF + known in ClinVar/COSMIC + high depth
      vaf >= detection_threshold * 2 &
        !is.na(clinvar_significance) &
        alt_depth >= 20 ~ "High confidence",
      # Questionable but potentially impactful: VAF near/below threshold but high impact
      !is.na(vaf) & vaf < detection_threshold * 2 &
        impact == "HIGH" ~ "Questionable but potentially impactful",
      # Uncertain: borderline quality or ambiguous implications
      !is.na(vaf) ~ "Uncertain",
      # CNAs don't have VAF-based confidence
      TRUE ~ "N/A"
    ),
    confidence_reason = case_when(
      confidence == "High confidence" ~
        paste0("VAF ", round(vaf*100,1), "% (>", round(detection_threshold*200,0), "% threshold), ",
               "depth ", alt_depth, "/", total_depth, ", ",
               ifelse(!is.na(clinvar_significance), clinvar_significance, "no ClinVar")),
      confidence == "Questionable but potentially impactful" ~
        paste0("VAF ", round(vaf*100,1), "% near detection limit (", round(detection_threshold*100,0), "%), ",
               "but ", impact, " impact -- recommend orthogonal validation"),
      confidence == "Uncertain" ~
        paste0("Borderline metrics: VAF ", round(vaf*100,1), "%, depth ", alt_depth, "/", total_depth),
      TRUE ~ "Copy number alteration -- confidence based on log2 ratio"
    )
  )

write.csv(merged_annotations, file.path(out_ann, "merged_annotations.tsv"),
          row.names = FALSE)
log_info("Somatic variants saved: {nrow(merged_annotations)} variants")

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 3: Mock CNV Results
# ══════════════════════════════════════════════════════════════════════════════

log_info("=== Stage 3: CNV (mock) ===")
out_cnv <- stage_output_dir(SAMPLE_ID, "03-cnv")

parsed_cnv <- tibble(
  chromosome = c("chr19", "chr10", "chr17", "chr8", "chr17"),
  start = c(30200000L, 87860000L, 41200000L, 128700000L, 37800000L),
  end   = c(30250000L, 87970000L, 41280000L, 128770000L, 37900000L),
  gene = c("CCNE1", "PTEN", "BRCA1", "MYC", "ERBB2"),
  log2_ratio = c(2.8, -2.1, -1.5, 1.9, 0.4),
  copy_number = c(12, 0, 1, 8, 3),
  length = c(50000L, 110000L, 80000L, 70000L, 100000L),
  type = c("AMPLIFICATION", "DELETION", "DELETION", "AMPLIFICATION", "NEUTRAL")
)

# Filter significant
parsed_cnv_sig <- parsed_cnv |> filter(type != "NEUTRAL")

write.csv(parsed_cnv_sig, file.path(out_cnv, "cnvkit_segments.tsv"), row.names = FALSE)
log_info("CNV results saved: {nrow(parsed_cnv_sig)} significant CNAs")

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 4: Mock Fusion Results
# ══════════════════════════════════════════════════════════════════════════════

log_info("=== Stage 4: Fusions (mock) ===")
out_fus <- stage_output_dir(SAMPLE_ID, "04-fusions")

parsed_fusions <- tibble(
  gene_a = c("EWSR1"),
  gene_b = c("FLI1"),
  fusion_name = c("EWSR1::FLI1"),
  exon_a = c(7L),
  exon_b = c(6L),
  chr_a = c("chr22"),
  pos_a = c(29683123L),
  chr_b = c("chr11"),
  pos_b = c(128675261L),
  supporting_reads = c(45L),
  fusion_type = c("in-frame"),
  known_fusion = c(TRUE)
)

write.csv(parsed_fusions, file.path(out_fus, "fusions.tsv"), row.names = FALSE)
log_info("Fusion results saved: {nrow(parsed_fusions)} fusions")

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 5: Biomarkers
# ══════════════════════════════════════════════════════════════════════════════

log_info("=== Stage 5: Biomarkers ===")
out_bio <- stage_output_dir(SAMPLE_ID, "05-biomarkers")

# TMB calculation from mock variants
coding_variants <- merged_annotations |>
  filter(!is.na(vaf), consequence != "amplification")
tmb_score <- nrow(coding_variants) / config$biomarkers$tmb$panel_coding_size_mb

tmb_result <- list(
  tmb_score = round(tmb_score, 1),
  tmb_class = if (tmb_score >= 10) "TMB-High" else if (tmb_score >= 6) "TMB-Intermediate" else "TMB-Low",
  variant_count = nrow(coding_variants),
  panel_size_mb = config$biomarkers$tmb$panel_coding_size_mb
)

# MSI — ovarian cancers are typically MSS unless MMR-deficient
msi_result <- list(
  msi_score = 5.2,
  msi_status = "MSS",
  unstable_sites = 2L,
  total_sites = 38L
)

# HRD — high-grade serous ovarian with BRCA1 loss → expect HRD+
hrd_result <- list(
  hrd_score = 58,
  hrd_status = "HRD-positive",
  loh_score = 22,
  tai_score = 18,
  lst_score = 18,
  is_reliable = FALSE  # panel-based estimation
)

write_json(tmb_result, file.path(out_bio, "tmb_result.json"), auto_unbox = TRUE, pretty = TRUE)
write_json(msi_result, file.path(out_bio, "msi_result.json"), auto_unbox = TRUE, pretty = TRUE)
write_json(hrd_result, file.path(out_bio, "hrd_result.json"), auto_unbox = TRUE, pretty = TRUE)

log_info("TMB: {tmb_result$tmb_score} mut/Mb ({tmb_result$tmb_class})")
log_info("MSI: {msi_result$msi_status} ({msi_result$msi_score}%)")
log_info("HRD: {hrd_result$hrd_status} (score={hrd_result$hrd_score})")

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 6: OncoKB Clinical Annotation (REAL API CALLS)
# ══════════════════════════════════════════════════════════════════════════════

log_info("=== Stage 6: OncoKB Clinical Annotation (REAL API) ===")
out_clin <- stage_output_dir(SAMPLE_ID, "06-clinical-annotation")

oncokb_tumor_type <- "HGSOC"

# -- Annotate mutations --
log_info("Querying OncoKB for {nrow(merged_annotations)} variants...")

mutation_annotations <- list()
for (i in seq_len(nrow(merged_annotations))) {
  row <- merged_annotations[i, ]
  if (row$consequence == "amplification") next  # CNAs handled separately

  protein_change <- str_remove(row$hgvsp, "^p\\.")

  result <- tryCatch({
    oncokb_annotate_mutation(row$gene, protein_change, oncokb_tumor_type)
  }, error = function(e) {
    log_warn("OncoKB mutation error for {row$gene} {protein_change}: {e$message}")
    list(
      gene = row$gene, alteration = protein_change, tumor_type = oncokb_tumor_type,
      oncogenic = "Unknown", mutation_effect = "Unknown",
      highest_sensitive_level = NA_character_,
      highest_resistance_level = NA_character_,
      treatments = tibble(drugs = character(), level = character(),
                          description = character(), fda_approved = logical())
    )
  })

  mutation_annotations[[length(mutation_annotations) + 1]] <- result
  log_info("  {row$gene} {protein_change}: {result$oncogenic} | Level: {result$highest_sensitive_level %||% 'none'}")

  Sys.sleep(0.5)  # Rate limiting courtesy
}

# -- Annotate CNAs --
log_info("Querying OncoKB for {nrow(parsed_cnv_sig)} CNAs...")

cna_annotations <- list()
for (i in seq_len(nrow(parsed_cnv_sig))) {
  row <- parsed_cnv_sig[i, ]
  cna_type <- toupper(row$type)
  if (!cna_type %in% c("AMPLIFICATION", "DELETION")) next

  result <- tryCatch({
    oncokb_annotate_cna(row$gene, cna_type, oncokb_tumor_type)
  }, error = function(e) {
    log_warn("OncoKB CNA error for {row$gene} {cna_type}: {e$message}")
    list(
      gene = row$gene, alteration = cna_type, tumor_type = oncokb_tumor_type,
      oncogenic = "Unknown", highest_sensitive_level = NA_character_,
      treatments = tibble(drugs = character(), level = character(),
                          description = character(), fda_approved = logical())
    )
  })

  cna_annotations[[length(cna_annotations) + 1]] <- result
  log_info("  {row$gene} {cna_type}: {result$oncogenic} | Level: {result$highest_sensitive_level %||% 'none'}")

  Sys.sleep(0.5)
}

# -- Annotate fusions --
log_info("Querying OncoKB for {nrow(parsed_fusions)} fusions...")

fusion_annotations <- list()
for (i in seq_len(nrow(parsed_fusions))) {
  row <- parsed_fusions[i, ]

  result <- tryCatch({
    oncokb_annotate_fusion(row$gene_a, row$gene_b, oncokb_tumor_type)
  }, error = function(e) {
    log_warn("OncoKB fusion error for {row$fusion_name}: {e$message}")
    list(
      gene_a = row$gene_a, gene_b = row$gene_b,
      alteration = row$fusion_name, tumor_type = oncokb_tumor_type,
      oncogenic = "Unknown", highest_sensitive_level = NA_character_,
      treatments = tibble(drugs = character(), level = character(),
                          description = character(), fda_approved = logical())
    )
  })

  fusion_annotations[[length(fusion_annotations) + 1]] <- result
  log_info("  {row$fusion_name}: {result$oncogenic} | Level: {result$highest_sensitive_level %||% 'none'}")
}

oncokb_results <- list(
  mutations = mutation_annotations,
  cnas = cna_annotations,
  fusions = fusion_annotations
)

# Save OncoKB results (strip raw to avoid serialization issues)
oncokb_save <- list(
  mutations = map(mutation_annotations, ~ .x[names(.x) != "raw"]),
  cnas = map(cna_annotations, ~ .x[names(.x) != "raw"]),
  fusions = map(fusion_annotations, ~ .x[names(.x) != "raw"])
)
write_json(oncokb_save, file.path(out_clin, "oncokb_results.json"),
           auto_unbox = TRUE, pretty = TRUE)

# -- ESCAT Classification --
log_info("Classifying ESCAT tiers...")
escat_tiers <- classify_all_escat(oncokb_results, config, SAMPLE_ID)
log_info("ESCAT classification complete:")
tier_counts <- escat_tiers |> count(escat_tier)
for (i in seq_len(nrow(tier_counts))) {
  log_info("  Tier {tier_counts$escat_tier[i]}: {tier_counts$n[i]} alterations")
}

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 7: Literature Search (REAL API CALLS)
# ══════════════════════════════════════════════════════════════════════════════

log_info("=== Stage 7: Literature Search (REAL API) ===")
out_lit <- stage_output_dir(SAMPLE_ID, "07-literature")

# Search PubMed for key variants
key_genes <- c("TP53", "BRCA1", "PIK3CA", "CCNE1", "PTEN")
key_alterations <- c("R248W", "E1836fs", "H1047R", "amplification", "deletion")

pubmed_results <- tibble()
for (i in seq_along(key_genes)) {
  log_info("PubMed search: {key_genes[i]} {key_alterations[i]} ovarian cancer")
  hits <- tryCatch({
    search_pubmed(
      gene = key_genes[i],
      variant = key_alterations[i],
      tumor_type = "ovarian cancer",
      max_results = 5
    )
  }, error = function(e) {
    log_warn("PubMed error for {key_genes[i]}: {e$message}")
    tibble(pmid = character(), title = character(), authors = character(),
           journal = character(), year = character(), abstract = character())
  })

  if (nrow(hits) > 0) {
    hits$query_gene <- key_genes[i]
    hits$query_alteration <- key_alterations[i]
    pubmed_results <- bind_rows(pubmed_results, hits)
    log_info("  Found {nrow(hits)} articles")
  }
  Sys.sleep(0.5)
}

pubmed_results <- distinct(pubmed_results, pmid, .keep_all = TRUE)
write_json(pubmed_results, file.path(out_lit, "pubmed_hits.json"),
           auto_unbox = TRUE, pretty = TRUE)
log_info("PubMed: {nrow(pubmed_results)} unique articles found")

# Search Scopus for key variants
scopus_results <- tibble()
for (i in seq_along(key_genes)) {
  log_info("Scopus search: {key_genes[i]} {key_alterations[i]} ovarian cancer")
  hits <- tryCatch({
    search_scopus(
      gene = key_genes[i],
      variant = key_alterations[i],
      tumor_type = "ovarian cancer",
      max_results = 5
    )
  }, error = function(e) {
    log_warn("Scopus error for {key_genes[i]}: {e$message}")
    tibble(scopus_id = character(), title = character(), authors = character(),
           journal = character(), year = character(), doi = character(),
           cited_by = integer())
  })

  if (nrow(hits) > 0) {
    hits$query_gene <- key_genes[i]
    hits$query_alteration <- key_alterations[i]
    scopus_results <- bind_rows(scopus_results, hits)
    log_info("  Found {nrow(hits)} articles")
  }
  Sys.sleep(0.5)
}

scopus_results <- distinct(scopus_results, scopus_id, .keep_all = TRUE)
write_json(scopus_results, file.path(out_lit, "scopus_hits.json"),
           auto_unbox = TRUE, pretty = TRUE)
log_info("Scopus: {nrow(scopus_results)} unique articles found")

# Generate TREATMENT-FOCUSED narratives with AMA-style references
narratives <- list()
for (i in seq_along(key_genes)) {
  gene <- key_genes[i]
  alt <- key_alterations[i]

  # Find OncoKB annotation for this gene
  oncokb_info <- NULL
  for (m in mutation_annotations) {
    if (m$gene == gene) { oncokb_info <- m; break }
  }
  for (c_ann in cna_annotations) {
    if (c_ann$gene == gene) { oncokb_info <- c_ann; break }
  }

  # Get relevant literature
  gene_pubs <- pubmed_results |> filter(query_gene == gene)

  # Build AMA-style references (Author et al. *Journal*. Year;PMID)
  ama_refs <- ""
  if (nrow(gene_pubs) > 0) {
    top_refs <- head(gene_pubs, 3)
    ref_nums <- seq_len(nrow(top_refs))
    ama_strings <- map_chr(ref_nums, function(j) {
      first_author <- str_extract(top_refs$authors[j], "^[^,]+")
      if (is.na(first_author)) first_author <- "Unknown"
      author_count <- str_count(top_refs$authors[j], ",") + 1
      auth <- if (author_count > 3) paste0(first_author, " et al") else top_refs$authors[j]
      paste0(auth, ". ", top_refs$title[j], ". *", top_refs$journal[j], "*. ",
             top_refs$year[j], ". PMID: ", top_refs$pmid[j])
    })
    ama_refs <- paste0("\n\n**References:** ", paste(paste0(ref_nums, ". ", ama_strings), collapse = " "))
  }

  # Oncogenic classification
  oncogenic_text <- if (!is.null(oncokb_info)) oncokb_info$oncogenic else "Unknown"

  # TREATMENT-FOCUSED narrative
  if (!is.null(oncokb_info) && nrow(oncokb_info$treatments) > 0) {
    tx <- oncokb_info$treatments
    fda_tx <- tx |> filter(fda_approved == TRUE)
    invest_tx <- tx |> filter(fda_approved == FALSE)

    # Standard of care treatments
    soc_text <- ""
    if (nrow(fda_tx) > 0) {
      soc_text <- paste0(
        "**Recommended treatment:** Based on OncoKB Level ",
        paste(unique(fda_tx$level), collapse = "/"),
        " evidence, ",
        paste(fda_tx$drugs, collapse = " or "),
        " is/are indicated as standard-of-care targeted therapy for patients with this alteration. ",
        "These agents have demonstrated clinical benefit in prospective trials and are FDA-approved ",
        "for this indication. "
      )
    }

    # Investigational treatments
    invest_text <- ""
    if (nrow(invest_tx) > 0) {
      invest_text <- paste0(
        "**Investigational options:** ",
        paste(invest_tx$drugs, collapse = ", "),
        " (Evidence Level: ", paste(invest_tx$level, collapse = ", "),
        ") represent investigational strategies with emerging evidence. ",
        "Enrollment in clinical trials targeting this alteration should be considered. "
      )
    }

    narrative <- paste0(
      "*", gene, "* ", alt, " is classified as **", oncogenic_text,
      "** in high-grade serous ovarian carcinoma. ",
      soc_text, invest_text,
      "Treatment decisions should integrate this molecular finding with clinical context, ",
      "prior therapy lines, and patient performance status.",
      ama_refs
    )
  } else {
    # No targeted therapy available
    narrative <- paste0(
      "*", gene, "* ", alt, " is classified as **", oncogenic_text,
      "** in high-grade serous ovarian carcinoma. ",
      "No directly targeted therapies are currently approved for this specific alteration ",
      "in ovarian cancer. However, this alteration may inform prognosis and response to ",
      "standard platinum-based chemotherapy. Consider enrollment in basket trials or ",
      "molecular tumor boards for emerging therapeutic strategies.",
      ama_refs
    )
  }

  narratives[[paste0(gene, "_", alt)]] <- narrative
  log_info("Treatment narrative generated for {gene} {alt}")
}

write_json(narratives, file.path(out_lit, "narratives.json"),
           auto_unbox = TRUE, pretty = TRUE)

literature_results <- list(
  pubmed = pubmed_results,
  scopus = scopus_results,
  narratives = narratives
)

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 8: Report Generation
# ══════════════════════════════════════════════════════════════════════════════

log_info("=== Stage 8: ESMO Report Generation ===")
out_report <- stage_output_dir(SAMPLE_ID, "08-report")

# Save all data as RDS for the Quarto document
data_dir <- file.path(out_report, "data")
dir_create(data_dir)

saveRDS(qc_results, file.path(data_dir, "qc_results.rds"))
saveRDS(merged_annotations, file.path(data_dir, "merged_annotations.rds"))
saveRDS(parsed_cnv_sig, file.path(data_dir, "parsed_cnv.rds"))
saveRDS(parsed_fusions, file.path(data_dir, "parsed_fusions.rds"))
saveRDS(tmb_result, file.path(data_dir, "tmb_result.rds"))
saveRDS(msi_result, file.path(data_dir, "msi_result.rds"))
saveRDS(hrd_result, file.path(data_dir, "hrd_result.rds"))
saveRDS(oncokb_results, file.path(data_dir, "oncokb_results.rds"))
saveRDS(escat_tiers, file.path(data_dir, "escat_tiers.rds"))
saveRDS(literature_results, file.path(data_dir, "literature_results.rds"))
saveRDS(config, file.path(data_dir, "config.rds"))

log_info("All stage data saved to {data_dir}")

# Render the Quarto report
report_qmd <- here("08-report/clinical_report.qmd")
output_html <- file.path(out_report, "clinical_report.html")

log_info("Rendering Quarto report...")
quarto_result <- tryCatch({
  system2("quarto", args = c(
    "render", report_qmd,
    "--output-dir", out_report,
    "-P", glue("sample_id:{SAMPLE_ID}"),
    "-P", glue("reports_dir:{here('reports')}"),
    "-P", glue("config_path:{here('config/default.yaml')}")
  ), stdout = TRUE, stderr = TRUE)
}, error = function(e) {
  log_error("Quarto render failed: {e$message}")
  NULL
})

exit_code <- attr(quarto_result, "status") %||% 0L

if (exit_code == 0 && file.exists(output_html)) {
  log_info("Report rendered successfully: {output_html}")
} else {
  log_warn("Quarto render had issues (exit code: {exit_code})")
  log_warn("Falling back to standalone report generation...")

  # Fallback: generate report directly in R
  source(here("run_standalone_report.R"))
  generate_standalone_report(
    sample_id = SAMPLE_ID,
    config = config,
    qc = qc_results,
    variants = merged_annotations,
    cnv = parsed_cnv_sig,
    fusions = parsed_fusions,
    tmb = tmb_result,
    msi = msi_result,
    hrd = hrd_result,
    escat = escat_tiers,
    literature = literature_results,
    output_path = output_html
  )
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  Pipeline Complete!                                        ║\n")
cat(glue("║  Report: {output_html}"), "\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")
