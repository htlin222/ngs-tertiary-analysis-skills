#!/usr/bin/env Rscript
# scripts/build_handoff_index.R — Rewrite the Batch 7 handoff folder's
# index.html with a Foundation One–style table linking every patient's
# actionable HTML and PDF, plus a quick "actionable findings" preview.

.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library",
            .libPaths()))
setwd("/Users/htlin/ngs-tertiary-analysis-skills")

suppressPackageStartupMessages({
  library(htmltools); library(dplyr); library(stringr); library(glue); library(tibble)
  library(here); library(fs); library(purrr)
})
source(here("R/foundation_helpers.R"))

handoff <- here("reports", "_handoff_Batch_7_20260505")
samples <- tribble(
  ~sid,         ~tt,    ~zh,
  "M26-0165R",  "UCEC", "子宮內膜癌",
  "M26-0234R",  "BRCA", "乳癌",
  "M26-0284R",  "COAD", "大腸癌",
  "M26-0320R",  "BRCA", "乳癌",
  "M26-0400R",  "UCEC", "子宮內膜癌",
  "M26-0401R",  "COAD", "乙狀結腸",
  "M26-0410R",  "UCEC", "子宮內膜癌",
  "M26-0414R",  "UCEC", "子宮內膜癌"
)

# Pull a per-sample summary (TMB / MSI / GIS + #actionable + top therapy)
load_one <- function(sid, tt) {
  d <- here("reports", sid, "08-report", "data")
  if (!file.exists(file.path(d, "oncokb.rds"))) return(NULL)
  oncokb <- readRDS(file.path(d, "oncokb.rds"))
  amp    <- readRDS(file.path(d, "amp.rds"))
  parsed_summary <- list(
    tmb = readRDS(file.path(d, "tmb.rds")),
    msi = readRDS(file.path(d, "msi.rds")),
    hrd = readRDS(file.path(d, "hrd.rds"))
  )
  tx <- extract_therapies(oncokb)
  on_label <- sum(tx$tier == "on-label")
  off_label <- sum(tx$tier == "off-label")
  resistance <- sum(tx$tier == "resistance")
  amp_tier12 <- if (is.data.frame(amp)) sum(amp$amp_tier %in% c("Tier I", "Tier II")) else 0L
  top_drug <- if (nrow(tx) > 0) head(tx$drugs, 1) else "—"
  tibble(
    sid = sid, tt = tt,
    tmb = parsed_summary$tmb$tmb_score, tmb_class = parsed_summary$tmb$tmb_class,
    msi = parsed_summary$msi$msi_score, msi_status = parsed_summary$msi$msi_status,
    gis = parsed_summary$hrd$hrd_score, hrd_status = parsed_summary$hrd$hrd_status,
    on_label = on_label, off_label = off_label, resistance = resistance,
    amp_actionable = amp_tier12, top_drug = top_drug
  )
}

rows <- purrr::map2_dfr(samples$sid, samples$tt, load_one)
samples <- samples |> left_join(rows, by = c("sid", "tt"))

fmt_n <- function(x) if (is.na(x)) "—" else format(x, nsmall = 1)
biomarker_pill <- function(value, status, kind) {
  band <- if (kind == "tmb") {
            if (!is.na(value) && value >= 10) "high"
            else if (!is.na(value) && value >= 6) "med" else "low"
          } else if (kind == "msi") {
            if (grepl("MSI-High", status %||% "")) "high" else "low"
          } else {
            if (!is.na(value) && value >= 42) "high"
            else if (!is.na(value) && value >= 30) "med" else "low"
          }
  tags$span(class = paste("pill", band), value)
}

render_row <- function(r) {
  has_html <- file_exists(file.path(handoff, glue("{r$sid}_{r$tt}_actionable.html")))
  has_pdf  <- file_exists(file.path(handoff, glue("{r$sid}_{r$tt}_actionable.pdf")))
  tags$tr(
    tags$td(class = "sid", r$sid),
    tags$td(tags$span(class = paste("tt-badge", tolower(r$tt)), r$tt),
            tags$br(), tags$small(r$zh)),
    tags$td(biomarker_pill(fmt_n(r$tmb), r$tmb_class, "tmb")),
    tags$td(biomarker_pill(fmt_n(r$msi), r$msi_status, "msi")),
    tags$td(biomarker_pill(round(r$gis %||% NA, 0), r$hrd_status, "hrd")),
    tags$td(class = "actions",
            if (r$on_label > 0)  tags$span(class = "tag on",  paste(r$on_label,  "on-label")) else NULL,
            if (r$off_label > 0) tags$span(class = "tag off", paste(r$off_label, "off-label")) else NULL,
            if (r$resistance > 0) tags$span(class = "tag res", paste(r$resistance, "R")) else NULL,
            if (!is.na(r$amp_actionable) && r$amp_actionable > 0)
              tags$span(class = "tag amp", paste(r$amp_actionable, "AMP I/II")) else NULL),
    tags$td(class = "top-drug",
            if (!is.na(r$top_drug)) r$top_drug else "—"),
    tags$td(class = "links",
            if (has_html) tags$a(href = glue("{r$sid}_{r$tt}_actionable.html"), "HTML") else NULL,
            if (has_html && has_pdf) " · " else NULL,
            if (has_pdf)  tags$a(href = glue("{r$sid}_{r$tt}_actionable.pdf"),  "PDF") else NULL)
  )
}

doc <- tags$html(lang = "en",
  tags$head(
    tags$meta(charset = "utf-8"),
    tags$title("TSO500-HRD Batch 7 — Clinician Handoff"),
    tags$style(HTML("
      body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
           max-width:1100px;margin:2em auto;padding:0 1.5em;color:#222;line-height:1.5}
      h1{color:#0d2c54;border-bottom:3px solid #0d2c54;padding-bottom:.4em}
      h2{color:#1a3e7c;margin-top:2em}
      table{width:100%;border-collapse:collapse;margin:1em 0;font-size:.95em}
      th{background:#0d2c54;color:#fff;text-align:left;padding:.6em .8em}
      td{padding:.6em .8em;border-bottom:1px solid #e0e6ef;vertical-align:top}
      tr:nth-child(even) td{background:#f6f9fc}
      .sid{font-weight:600;color:#0d2c54}
      .tt-badge{display:inline-block;padding:1px 8px;border-radius:3px;
                font-size:.78em;color:#fff;font-weight:600}
      .tt-badge.ucec{background:#7b1fa2}.tt-badge.brca{background:#c2185b}
      .tt-badge.coad{background:#00796b}
      small{color:#666;font-size:.85em}
      .pill{display:inline-block;padding:1px 10px;border-radius:3px;
            font-weight:600;font-size:.85em}
      .pill.high{background:#ffefef;color:#b71c1c;border:1px solid #ffcdd2}
      .pill.med {background:#fff7e6;color:#e65100;border:1px solid #ffe0b2}
      .pill.low {background:#f1f8e9;color:#33691e;border:1px solid #dcedc8}
      .tag{display:inline-block;padding:1px 6px;border-radius:3px;font-size:.8em;margin-right:3px;margin-bottom:2px}
      .tag.on{background:#c62828;color:#fff}
      .tag.off{background:#1565c0;color:#fff}
      .tag.res{background:#6a1b9a;color:#fff}
      .tag.amp{background:#0d2c54;color:#fff}
      .top-drug{font-size:.85em;font-style:italic;color:#444;max-width:240px}
      .links a{color:#1565c0;text-decoration:none;font-weight:600}
      .links a:hover{text-decoration:underline}
      .note{background:#fff8e1;border-left:4px solid #ffa000;padding:.8em 1em;
            margin:1em 0;font-size:.92em}
    "))
  ),
  tags$body(
    tags$h1("TSO500-HRD Batch 7 — Clinician Handoff"),
    tags$p(
      tags$b("Run date:"), " ", format(Sys.time(), "%Y-%m-%d"), " · ",
      tags$b("Batch date:"), " 2026-04-28 · ",
      tags$b("Samples:"), " 8 · ",
      tags$b("Panel:"), " Illumina TruSight Oncology 500 (523 genes, ~1.94 Mb)"
    ),
    tags$div(class = "note",
      "Reports are Foundation One–style: clinician-actionable, no interactive ",
      "widgets, designed for plain PDF print. Each ", tags$b("HTML"),
      " is self-contained (~300 KB); ", tags$b("PDF"), " is the print-ready ",
      "headless-Chrome rendering (~1 MB)."),
    tags$h2("Patients"),
    tags$table(
      tags$thead(tags$tr(
        tags$th("Sample"), tags$th("Tumor"),
        tags$th("TMB (mut/Mb)"), tags$th("MSI (% unstable)"), tags$th("GIS"),
        tags$th("Findings"), tags$th("Top therapy"), tags$th("Report"))),
      tags$tbody(lapply(seq_len(nrow(samples)), function(i) render_row(samples[i, ])))
    ),
    tags$h2("Legend"),
    tags$ul(
      tags$li(tags$span(class = "tag on",  "on-label"), " — OncoKB Level 1, 2, 3A in this tumor type"),
      tags$li(tags$span(class = "tag off", "off-label"), " — OncoKB Level 3B / 4 (other tumor types or preclinical)"),
      tags$li(tags$span(class = "tag res", "R"), " — Predicted resistance (Level R1/R2)"),
      tags$li(tags$span(class = "tag amp", "AMP I/II"), " — AMP/ASCO/CAP Tier I or II per Li et al. 2017")
    ),
    tags$p(class = "note",
      tags$b("Caveats: "),
      "BAM-level QC (coverage/duplicates) is unavailable — input was the ",
      "DRAGEN VCF + CombinedVariantOutput.tsv. Tumor fraction, ploidy, ",
      "absolute copy number, and LoH are beta features per Illumina TSO500 v2.6.0.")
  )
)

out <- file.path(handoff, "index.html")
htmltools::save_html(doc, out, libdir = NULL)
cat("Wrote", out, "\n")
