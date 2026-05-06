#!/usr/bin/env Rscript
# scripts/openevidence_queue.R — Build a deduplicated list of (gene, alt, tt)
# tuples that need an OpenEvidence literature lookup. Emits a TSV that
# Claude/MCP can iterate over and fill into reports/.openevidence_cache/.

.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library",
            .libPaths()))
setwd("/Users/htlin/ngs-tertiary-analysis-skills")

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(stringr); library(here)
})
source("R/oncokb_helpers.R")

samples <- tribble(
  ~sid,         ~tt,    ~tumor_label,
  "M26-0165R",  "UCEC", "endometrial carcinoma",
  "M26-0234R",  "BRCA", "breast cancer",
  "M26-0284R",  "COAD", "colorectal cancer",
  "M26-0320R",  "BRCA", "breast cancer",
  "M26-0400R",  "UCEC", "endometrial carcinoma",
  "M26-0401R",  "COAD", "colorectal cancer",
  "M26-0410R",  "UCEC", "endometrial carcinoma",
  "M26-0414R",  "UCEC", "endometrial carcinoma"
)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

rows <- list()
for (i in seq_len(nrow(samples))) {
  sid <- samples$sid[i]; tt <- samples$tt[i]; tl <- samples$tumor_label[i]
  d <- file.path("reports", sid, "08-report", "data")
  amp <- readRDS(file.path(d, "amp.rds"))
  oncokb <- readRDS(file.path(d, "oncokb.rds"))

  # Short variants — Tier I/II
  snv <- amp |> filter(amp_tier %in% c("Tier I", "Tier II")) |>
            distinct(gene, hgvsp)
  for (j in seq_len(nrow(snv))) {
    rows[[length(rows)+1]] <- tibble(sid = sid, tt = tt, tumor_label = tl,
      kind = "snv", gene = snv$gene[j],
      alt_3letter = snv$hgvsp[j],
      alt_short = hgvsp_to_short(as.character(snv$hgvsp[j])),
      alt_display = paste0(snv$gene[j], " ", snv$hgvsp[j])
    )
  }

  # CNVs — only those OncoKB flagged oncogenic / sensitivity-leveled
  for (m in oncokb$cnas %||% list()) {
    if (is.null(m) || isTRUE(is.na(m))) next
    lvl <- m$highest_sensitive_level
    onc <- m$oncogenic %||% "Unknown"
    if (!((length(lvl) > 0 && !is.na(lvl) && nzchar(lvl)) ||
          onc %in% c("Oncogenic", "Likely Oncogenic"))) next
    rows[[length(rows)+1]] <- tibble(sid = sid, tt = tt, tumor_label = tl,
      kind = "cnv", gene = m$gene,
      alt_3letter = m$alteration, alt_short = m$alteration,
      alt_display = paste0(m$gene, " ", tolower(m$alteration))
    )
  }
}

queue <- bind_rows(rows)
unique_queue <- queue |>
  distinct(kind, gene, alt_short, tt, tumor_label) |>
  mutate(cache_key = paste0(gene, "_", gsub("[^A-Za-z0-9]+", "-", alt_short), "_", tt))

dir.create("reports/.openevidence_cache", showWarnings = FALSE, recursive = TRUE)
readr::write_tsv(unique_queue, "reports/.openevidence_cache/_queue.tsv")
readr::write_tsv(queue,        "reports/.openevidence_cache/_per_sample.tsv")

cat(sprintf("Unique queries: %d\nPer-sample expansions: %d\n",
            nrow(unique_queue), nrow(queue)))
print(unique_queue, n = Inf)
