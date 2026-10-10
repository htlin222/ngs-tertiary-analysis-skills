#!/usr/bin/env Rscript
# ESMO2026-poster/poster/extract_cohort.R
#
# De-identified per-patient table for the #6004 cohort poster figures.
# Joins reports/_cohort_per_patient.tsv (scripts/extract_cohort_metrics.R) with
# the OncoKB therapy levels in each patient's oncokb.rds. No sample IDs are
# written; output stays under reports/ (gitignored). No API calls.
#
#   Rscript --vanilla ESMO2026-poster/poster/extract_cohort.R

suppressPackageStartupMessages({
  library(readr); library(jsonlite); library(here)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

per_patient <- read_tsv(here("reports", "_cohort_per_patient.tsv"),
                        show_col_types = FALSE)
stopifnot(all(per_patient$has_data))

on_label_levels <- c("LEVEL_1", "LEVEL_2", "LEVEL_3A")

patients <- lapply(seq_len(nrow(per_patient)), function(i) {
  p <- per_patient[i, ]
  oncokb <- readRDS(here("reports", p$sid, "08-report", "data", "oncokb.rds"))
  on_label <- list(); resistance <- list()
  for (slot in c("mutations", "cnas", "fusions")) {
    for (m in oncokb[[slot]] %||% list()) {
      if (is.null(m) || isTRUE(is.na(m))) next
      alt <- list(gene = m$gene %||% NA_character_,
                  alteration = m$alteration %||% m$hgvsp %||% NA_character_)
      lvl <- m$highest_sensitive_level %||% NA_character_
      rl  <- m$highest_resistance_level %||% NA_character_
      if (!is.na(lvl) && lvl %in% on_label_levels)
        on_label[[length(on_label) + 1]] <- c(alt, level = lvl)
      if (!is.na(rl) && grepl("R[12]", rl))
        resistance[[length(resistance) + 1]] <- c(alt, level = rl)
    }
  }
  # Must agree with the per-patient counts the abstract numbers were built on
  stopifnot(length(on_label) == p$n_on_label,
            length(resistance) == p$n_resistance)
  list(
    tumor = p$tt,
    tmb = p$tmb,
    msi_high = grepl("MSI-High", p$msi_status),
    gis = p$gis,
    n_on_label = p$n_on_label,
    n_off_label = p$n_off_label,
    oncogenic_genes = I(strsplit(p$oncogenic_genes_csv %||% "", ",")[[1]]),
    on_label = on_label,
    resistance = resistance
  )
})

out <- here("reports", "_poster_cohort.json")
write_json(patients, out, auto_unbox = TRUE, pretty = TRUE, na = "null")
cat(sprintf("Wrote %d patients -> %s\n", length(patients), out))
