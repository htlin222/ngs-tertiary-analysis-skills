#!/usr/bin/env Rscript
# scripts/extract_cohort_metrics.R
#
# Aggregate the 62 batch-7/all-batches patient artefacts into cohort-level
# numbers suitable for the second (cohort-focused) ESMO 2026 abstract.
# Reads only the .rds files already on disk — no API calls.

.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library",
            .libPaths()))
setwd("/Users/htlin/ngs-tertiary-analysis-skills")

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(stringr); library(purrr); library(tidyr)
  library(readr); library(here); library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

manifest <- read_tsv(here("inputs/TSO500-HRD/_manifest.tsv"),
                     show_col_types = FALSE)
cat(sprintf("\n=== Cohort manifest: %d patients ===\n", nrow(manifest)))

per_patient <- map_dfr(seq_len(nrow(manifest)), function(i) {
  sid <- manifest$sample_id[i]
  tt  <- manifest$oncokb_code[i]
  d   <- file.path("reports", sid, "08-report", "data")
  if (!file.exists(file.path(d, "tmb.rds"))) {
    return(tibble(sid = sid, tt = tt, has_data = FALSE))
  }
  tmb_d <- readRDS(file.path(d, "tmb.rds"))
  msi_d <- readRDS(file.path(d, "msi.rds"))
  hrd_d <- readRDS(file.path(d, "hrd.rds"))
  amp <- tryCatch(readRDS(file.path(d, "amp.rds")), error = function(e) NULL)
  oncokb <- tryCatch(readRDS(file.path(d, "oncokb.rds")), error = function(e) NULL)
  parsed_cnv <- tryCatch(readRDS(file.path(d, "cnv.rds")), error = function(e) NULL)

  # Therapy-tier counts derived from oncokb mutations + cnas + fusions
  on_label <- off_label <- resistance <- 0L
  oncogenic_genes <- character()
  if (!is.null(oncokb)) {
    for (slot in c("mutations", "cnas", "fusions")) {
      for (m in oncokb[[slot]] %||% list()) {
        if (is.null(m) || isTRUE(is.na(m))) next
        lvl <- m$highest_sensitive_level %||% NA_character_
        rl  <- m$highest_resistance_level %||% NA_character_
        onc <- m$oncogenic %||% "Unknown"
        if (length(lvl) && !is.na(lvl) && nzchar(lvl)) {
          if (lvl %in% c("LEVEL_1", "LEVEL_2", "LEVEL_3A")) on_label  <- on_label + 1L
          if (lvl %in% c("LEVEL_3B", "LEVEL_4"))           off_label <- off_label + 1L
        }
        if (length(rl) && !is.na(rl) && nzchar(rl) && grepl("R[12]", rl))
          resistance <- resistance + 1L
        if (onc %in% c("Oncogenic", "Likely Oncogenic"))
          oncogenic_genes <- c(oncogenic_genes, m$gene %||% NA_character_)
      }
    }
  }
  oncogenic_genes <- unique(stats::na.omit(oncogenic_genes))

  amp_tier_i_ii <- 0L
  if (is.data.frame(amp) && nrow(amp) > 0) {
    amp_tier_i_ii <- amp |>
      filter(amp_tier %in% c("Tier I", "Tier II")) |>
      distinct(gene, hgvsp) |> nrow()
  }

  n_cnv_signif <- if (is.data.frame(parsed_cnv))
    sum(parsed_cnv$type %in% c("AMPLIFICATION", "DELETION"), na.rm = TRUE) else 0L

  tibble(
    sid = sid, tt = tt, has_data = TRUE,
    tmb = tmb_d$tmb_score %||% NA_real_,
    tmb_class = tmb_d$tmb_class %||% "Unknown",
    msi = msi_d$msi_score %||% NA_real_,
    msi_status = msi_d$msi_status %||% "Unknown",
    gis = hrd_d$hrd_score %||% NA_real_,
    n_oncogenic_genes = length(oncogenic_genes),
    n_amp_tier_i_ii = amp_tier_i_ii,
    n_on_label = on_label,
    n_off_label = off_label,
    n_resistance = resistance,
    n_cnv_signif = n_cnv_signif,
    oncogenic_genes_csv = paste(oncogenic_genes, collapse = ",")
  )
})

# Save per-patient table for inspection
write_tsv(per_patient, here("reports", "_cohort_per_patient.tsv"))

ok <- per_patient |> filter(has_data)
n_total <- nrow(ok)
cat(sprintf("Patients with .rds artefacts: %d / %d\n", n_total, nrow(manifest)))

# ── Headline cohort metrics ────────────────────────────────────────────────

cat("\n=== Cancer-type distribution ===\n")
ct <- ok |> count(tt, sort = TRUE)
print(ct)

cat(sprintf("\nDistinct OncoKB tumor types (excl. generic CANCER): %d\n",
            length(setdiff(unique(ok$tt), "CANCER"))))

cat("\n=== Biomarker rates ===\n")
tmb_h    <- sum(ok$tmb >= 10, na.rm = TRUE)
tmb_int  <- sum(ok$tmb >= 6 & ok$tmb < 10, na.rm = TRUE)
msi_h    <- sum(grepl("MSI-High", ok$msi_status))
hrd_pos  <- sum(ok$gis >= 42, na.rm = TRUE)
hrd_int  <- sum(ok$gis >= 30 & ok$gis < 42, na.rm = TRUE)

cat(sprintf("TMB-High (>=10 mut/Mb)     : %d / %d (%.1f%%)\n",
            tmb_h, n_total, 100 * tmb_h / n_total))
cat(sprintf("TMB-Intermediate (6-<10)   : %d / %d (%.1f%%)\n",
            tmb_int, n_total, 100 * tmb_int / n_total))
cat(sprintf("MSI-High                   : %d / %d (%.1f%%)\n",
            msi_h, n_total, 100 * msi_h / n_total))
cat(sprintf("HRD-positive (GIS >=42)    : %d / %d (%.1f%%)\n",
            hrd_pos, n_total, 100 * hrd_pos / n_total))
cat(sprintf("HRD-intermediate (30-<42)  : %d / %d (%.1f%%)\n",
            hrd_int, n_total, 100 * hrd_int / n_total))

cat("\n=== Actionability rates ===\n")
any_amp     <- sum(ok$n_amp_tier_i_ii > 0)
any_on      <- sum(ok$n_on_label > 0)
any_off     <- sum(ok$n_off_label > 0)
any_res     <- sum(ok$n_resistance > 0)
any_therapy <- sum((ok$n_on_label + ok$n_off_label) > 0)
all_with_oncogenic <- sum(ok$n_oncogenic_genes > 0)

cat(sprintf("Patients with >=1 AMP Tier I/II variant          : %d / %d (%.1f%%)\n",
            any_amp, n_total, 100 * any_amp / n_total))
cat(sprintf("Patients with >=1 OncoKB on-label (Level 1/2/3A) : %d / %d (%.1f%%)\n",
            any_on, n_total, 100 * any_on / n_total))
cat(sprintf("Patients with >=1 OncoKB off-label (Level 3B/4)  : %d / %d (%.1f%%)\n",
            any_off, n_total, 100 * any_off / n_total))
cat(sprintf("Patients with >=1 OncoKB resistance (R1/R2)      : %d / %d (%.1f%%)\n",
            any_res, n_total, 100 * any_res / n_total))
cat(sprintf("Patients with >=1 therapy match (on+off-label)   : %d / %d (%.1f%%)\n",
            any_therapy, n_total, 100 * any_therapy / n_total))
cat(sprintf("Patients with >=1 oncogenic / likely-oncogenic   : %d / %d (%.1f%%)\n",
            all_with_oncogenic, n_total, 100 * all_with_oncogenic / n_total))

cat("\n=== Distribution of OncoKB Level-1/2/3A therapy matches per patient ===\n")
print(summary(ok$n_on_label))
cat(sprintf("Median (IQR): %.0f (%.0f-%.0f)\n",
            median(ok$n_on_label), quantile(ok$n_on_label, .25),
            quantile(ok$n_on_label, .75)))

cat("\n=== Top recurrently altered genes (>=2 patients) ===\n")
top_genes <- ok |>
  mutate(genes = str_split(oncogenic_genes_csv, ",")) |>
  unnest(genes) |>
  filter(genes != "") |>
  count(genes, sort = TRUE) |>
  filter(n >= 2)
print(top_genes, n = 30)

cat("\n=== Per-cancer-type actionability ===\n")
per_ct <- ok |>
  group_by(tt) |>
  summarise(
    n = n(),
    tmb_h_n = sum(tmb >= 10, na.rm = TRUE),
    msi_h_n = sum(grepl("MSI-High", msi_status)),
    hrd_pos_n = sum(gis >= 42, na.rm = TRUE),
    any_on_label_n = sum(n_on_label > 0),
    any_therapy_n = sum((n_on_label + n_off_label) > 0),
    .groups = "drop"
  ) |>
  arrange(desc(n))
print(per_ct)

# ── Save summary as JSON for the abstract draft ───────────────────────────

summary_json <- list(
  n_total = n_total,
  n_cancer_types = length(setdiff(unique(ok$tt), "CANCER")),
  collection_period_months = 8,
  collection_period_batches = 8,
  rates = list(
    tmb_high = list(n = tmb_h, pct = round(100 * tmb_h / n_total, 1)),
    tmb_intermediate = list(n = tmb_int, pct = round(100 * tmb_int / n_total, 1)),
    msi_high = list(n = msi_h, pct = round(100 * msi_h / n_total, 1)),
    hrd_positive = list(n = hrd_pos, pct = round(100 * hrd_pos / n_total, 1)),
    hrd_intermediate = list(n = hrd_int, pct = round(100 * hrd_int / n_total, 1)),
    any_amp_tier_i_ii = list(n = any_amp, pct = round(100 * any_amp / n_total, 1)),
    any_on_label = list(n = any_on, pct = round(100 * any_on / n_total, 1)),
    any_off_label = list(n = any_off, pct = round(100 * any_off / n_total, 1)),
    any_resistance = list(n = any_res, pct = round(100 * any_res / n_total, 1)),
    any_therapy_match = list(n = any_therapy, pct = round(100 * any_therapy / n_total, 1))
  ),
  top_genes_ge2 = top_genes |>
    mutate(gene = genes) |> select(gene, n) |>
    head(20) |> as.data.frame()
)
write_json(summary_json, here("reports", "_cohort_summary.json"),
           pretty = TRUE, auto_unbox = TRUE)
cat("\nWrote reports/_cohort_summary.json and reports/_cohort_per_patient.tsv\n")
