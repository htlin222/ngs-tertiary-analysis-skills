#!/usr/bin/env Rscript
# scripts/audit_batch7.R — Batch 7 reconciliation:
# compare each report's findings against (a) the lab's curated 重點摘要 and
# (b) the per-sample TSO500 deliverable.

.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library",
            .libPaths()))
setwd("/Users/htlin/ngs-tertiary-analysis-skills")

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tibble); library(purrr); library(jsonlite); library(here)
})
source("R/foundation_helpers.R")

# ── Curated reference (from 重點摘要.xlsx, single source of truth) ─────────
expected <- tribble(
  ~sid,         ~tt,    ~tmb,  ~tmb_class, ~msi,    ~gis,
  "M26-0165R",  "UCEC", 10.2,  "TMB-High",  0.84,  36L,
  "M26-0234R",  "BRCA",  5.2,  "TMB-Low",   1.67,  57L,
  "M26-0284R",  "COAD",  5.5,  "TMB-Low",   0.87,   7L,
  "M26-0320R",  "BRCA", 10.5,  "TMB-High",  3.30,  37L,
  "M26-0400R",  "UCEC",  6.3,  "TMB-Low",   2.88,  20L,
  "M26-0401R",  "COAD",  9.7,  "TMB-Low", NA_real_, 16L,
  "M26-0410R",  "UCEC",  3.9,  "TMB-Low",   1.57,  30L,
  "M26-0414R",  "UCEC",  3.9,  "TMB-Low",   1.71,  30L
)

curated_alterations <- list(
  "M26-0165R" = c("ARID1A:p.A139Cfs*89:LEVEL_4",
                  "ATM:p.E2668*:NA",
                  "ATM:c.1236-2A>T:NA",
                  "KRAS:p.G12S:LEVEL_4",
                  "RNF43:p.R397Gfs*22:NA",
                  "BRCA1:LR_LOSS_FC0.74:NA"),
  "M26-0234R" = c("BRCA1:LR_LOSS_FC0.15:LEVEL_1",
                  "TP53:p.Q331Dfs*5:NA"),
  "M26-0284R" = c("NRAS:p.G13D:LEVEL_R1",
                  "TP53:p.R248Q:NA"),
  "M26-0320R" = c("PIK3CA:p.H1047R:LEVEL_1",
                  "ERBB2:Amp:LEVEL_1",
                  "CCNE1:Amp:LEVEL_4",
                  "TP53:p.L308Cfs*37:NA",
                  "TET2:p.E1151*:NA",
                  "SMARCA4:c.3168+2dup:NA",
                  "MYC:Amp:NA"),
  "M26-0400R" = c("ARID1A:p.V1817Yfs*66:LEVEL_4",
                  "PIK3CA:p.Q546R:LEVEL_4",
                  "PTEN:p.R130G:LEVEL_4",
                  "PTEN:p.R130Q:LEVEL_4",
                  "TP53:p.R175H:NA",
                  "MYC:Amp:NA"),
  "M26-0401R" = c("APC:p.S1356*:NA",
                  "KRAS:p.G12D:LEVEL_R1",
                  "FANCA:p.D1325Ifs*38:NA",
                  "TP53:p.R175H:NA",
                  "SOX9:p.A372Hfs*11:NA"),
  "M26-0410R" = c("NRAS:p.Q61L:NA",
                  "TP53:p.Y220S:NA"),
  "M26-0414R" = c("NRAS:p.Q61L:NA",
                  "TP53:p.Y220S:NA")
)

dr <- function(d) here("reports", d, "08-report", "data")

audit_one <- function(row) {
  sid <- row$sid
  d <- dr(sid)

  parsed_tmb <- readRDS(file.path(d, "tmb.rds"))
  parsed_msi <- readRDS(file.path(d, "msi.rds"))
  parsed_hrd <- readRDS(file.path(d, "hrd.rds"))
  oncokb     <- readRDS(file.path(d, "oncokb.rds"))
  amp        <- readRDS(file.path(d, "amp.rds"))

  cat(sprintf("\n=================== %s (%s) ===================\n", sid, row$tt))

  # Biomarkers
  bm_ok <- function(name, got, want, tol = 0.05) {
    flag <- if (is.na(got) && is.na(want)) "OK"
            else if (is.na(got) || is.na(want)) "MISMATCH"
            else if (abs(got - want) <= tol) "OK"
            else sprintf("DELTA=%.2f", abs(got - want))
    cat(sprintf("  %-9s want=%-7s got=%-7s [%s]\n",
                name, format(want), format(got), flag))
  }
  bm_ok("TMB",       parsed_tmb$tmb_score, row$tmb)
  bm_ok("MSI",       parsed_msi$msi_score, row$msi)
  bm_ok("GIS",       parsed_hrd$hrd_score, as.numeric(row$gis))
  cat(sprintf("  %-9s want=%-12s got=%-12s\n", "TMB-class",
              row$tmb_class, parsed_tmb$tmb_class))

  # Alterations
  oncokb_muts <- oncokb$mutations %||% list()
  oncokb_cnas <- oncokb$cnas %||% list()
  cna_genes <- unique(unlist(lapply(oncokb_cnas, function(m)
                       if (!is.null(m$gene)) m$gene else character())))
  all_oncokb_genes <- unique(c(
    unlist(lapply(oncokb_muts, function(m) if (!is.null(m$gene)) m$gene else character())),
    cna_genes
  ))
  exp_strs <- curated_alterations[[sid]]
  exp_genes <- unique(map_chr(exp_strs, ~ str_split_1(.x, ":")[1]))

  cat("\n  Curated alterations vs report:\n")
  for (g in exp_strs) {
    parts <- str_split_1(g, ":")
    gene <- parts[1]; alt <- parts[2]; lvl_want <- parts[3]

    if (grepl("Amp|LR_", alt)) {
      hit <- gene %in% cna_genes
      flag <- if (hit) "OK" else "MISSING"
      cat(sprintf("    %-8s %-25s want=%-9s  in CNV oncokb=%-5s [%s]\n",
                  gene, alt, lvl_want, hit, flag))
    } else if (startsWith(alt, "p.")) {
      proper <- sub("^p\\.", "", alt)
      lvl_got <- NA_character_
      for (m in oncokb_muts) {
        if (!is.null(m$gene) && m$gene == gene &&
            !is.null(m$alteration) && m$alteration == proper) {
          lvl_got <- if (length(m$highest_sensitive_level) > 0 &&
                         !is.na(m$highest_sensitive_level) &&
                         nzchar(m$highest_sensitive_level))
                       m$highest_sensitive_level
                     else if (length(m$highest_resistance_level) > 0 &&
                              !is.na(m$highest_resistance_level) &&
                              nzchar(m$highest_resistance_level))
                       m$highest_resistance_level
                     else NA_character_
          break
        }
      }
      flag <- if (is.na(lvl_want) || lvl_want == "NA") {
                if (is.na(lvl_got)) "OK" else paste0("got_", lvl_got)
              } else {
                want_short <- gsub("LEVEL_", "", lvl_want)
                got_short  <- if (!is.na(lvl_got)) gsub("LEVEL_", "", lvl_got) else "—"
                if (want_short == got_short) "OK"
                else sprintf("WANT=%s GOT=%s", want_short, got_short)
              }
      cat(sprintf("    %-8s %-25s want=%-9s  got=%-9s [%s]\n",
                  gene, alt, lvl_want, lvl_got %||% "—", flag))
    } else {
      cat(sprintf("    %-8s %-25s want=%-9s  (splice / non-protein notation — not OncoKB-queryable)\n",
                  gene, alt, lvl_want))
    }
  }

  amp_sig <- if (is.data.frame(amp))
               amp |> filter(amp_tier %in% c("Tier I","Tier II")) |>
               distinct(gene, hgvsp) |> arrange(gene)
             else tibble()
  cat(sprintf("\n  AMP Tier I/II short variants (%d): %s\n",
              nrow(amp_sig),
              paste(paste0(amp_sig$gene, " ", amp_sig$hgvsp), collapse = "; ")))

  # Surplus actionable in OncoKB not in curated list
  surplus_genes <- setdiff(all_oncokb_genes, exp_genes)
  surplus_actionable <- character()
  for (m in c(oncokb_muts, oncokb_cnas)) {
    if (is.null(m) || is.null(m$gene) || !(m$gene %in% surplus_genes)) next
    lvl <- m$highest_sensitive_level
    onc <- m$oncogenic %||% "Unknown"
    if ((length(lvl) > 0 && !is.na(lvl) && nzchar(lvl)) ||
        onc %in% c("Oncogenic","Likely Oncogenic")) {
      surplus_actionable <- unique(c(surplus_actionable, m$gene))
    }
  }
  if (length(surplus_actionable) > 0) {
    cat(sprintf("  EXTRA actionable genes (in OncoKB, not in curated): %s\n",
                paste(surplus_actionable, collapse=", ")))
  }
  invisible(NULL)
}

for (i in seq_len(nrow(expected))) audit_one(expected[i, ])
