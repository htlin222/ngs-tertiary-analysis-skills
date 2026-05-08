#!/usr/bin/env Rscript
# scripts/build_unified_index.R — Stage every produced clinical_actionable
# report (across all batches) into reports/_handoff_all_batches/ and write a
# single grouped index.html linking each patient.

.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library",
            .libPaths()))
setwd("/Users/htlin/ngs-tertiary-analysis-skills")

suppressPackageStartupMessages({
  library(htmltools); library(dplyr); library(stringr); library(glue); library(tibble)
  library(here); library(fs); library(purrr); library(readr); library(jsonlite)
})
source(here("R/foundation_helpers.R"))

manifest <- read_tsv(here("inputs/TSO500-HRD/_manifest.tsv"),
                     show_col_types = FALSE)

handoff <- here("reports", "_handoff_all_batches")
dir_create(handoff, recurse = TRUE)

# Collect per-sample data, copy HTML/PDF to handoff with friendly names
load_one <- function(sid, batch, tt, zh) {
  data_dir <- here("reports", sid, "08-report")
  src_html <- file.path(data_dir, "clinical_actionable.html")
  src_pdf  <- file.path(data_dir, "clinical_actionable.pdf")
  if (!file.exists(src_html)) return(NULL)

  base <- glue("{sid}_{tt}_actionable")
  file_copy(src_html, file.path(handoff, paste0(base, ".html")), overwrite = TRUE)
  has_pdf <- file.exists(src_pdf)
  if (has_pdf)
    file_copy(src_pdf, file.path(handoff, paste0(base, ".pdf")), overwrite = TRUE)

  d <- file.path(data_dir, "data")
  oncokb_d <- tryCatch(readRDS(file.path(d, "oncokb.rds")), error = function(e) NULL)
  amp_d    <- tryCatch(readRDS(file.path(d, "amp.rds")),    error = function(e) NULL)
  tmb_d    <- tryCatch(readRDS(file.path(d, "tmb.rds")),    error = function(e) NULL)
  msi_d    <- tryCatch(readRDS(file.path(d, "msi.rds")),    error = function(e) NULL)
  hrd_d    <- tryCatch(readRDS(file.path(d, "hrd.rds")),    error = function(e) NULL)

  tx <- if (!is.null(oncokb_d)) extract_therapies(oncokb_d) else tibble()
  on  <- sum(tx$tier == "on-label")
  off <- sum(tx$tier == "off-label")
  res <- sum(tx$tier == "resistance")
  amp_n <- if (is.data.frame(amp_d))
             sum(amp_d$amp_tier %in% c("Tier I", "Tier II")) else 0L
  top_drug <- if (nrow(tx) > 0) head(tx$drugs, 1) else "—"

  tibble(
    sid = sid, batch = batch, tt = tt, zh = zh,
    base = base, has_pdf = has_pdf,
    tmb_val = tmb_d$tmb_score %||% NA_real_,
    tmb_class = tmb_d$tmb_class %||% "Unknown",
    msi_val = msi_d$msi_score %||% NA_real_,
    msi_status = msi_d$msi_status %||% "Unknown",
    gis = hrd_d$hrd_score %||% NA_real_,
    on_label = on, off_label = off, resistance = res,
    amp_actionable = amp_n, top_drug = top_drug
  )
}

rows <- pmap_dfr(list(manifest$sample_id, manifest$batch,
                       manifest$oncokb_code, manifest$tumor_zh), load_one)

cat(glue("Staged {nrow(rows)} reports into {handoff}\n\n"))

fmt_n <- function(x) if (is.na(x)) "—" else formatC(x, digits = 1, format = "f")
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

tt_color <- function(tt) {
  switch(tt,
    UCEC="#7b1fa2", BRCA="#c2185b", COAD="#00796b", READ="#00796b",
    NSCLC="#0d47a1", LUAD="#0d47a1", OV="#6a1b9a", HGSOC="#6a1b9a",
    PAAD="#bf360c", STAD="#5d4037", LIHC="#3e2723", CESC="#880e4f",
    SARC="#33691e", CHOL="#827717", PRAD="#1a237e", "#424242")
}

render_row <- function(r) {
  href_html <- paste0(r$base, ".html")
  href_pdf  <- paste0(r$base, ".pdf")
  badge <- tags$span(class = "tt-badge",
                     style = paste0("background:", tt_color(r$tt)),
                     r$tt)
  tags$tr(
    tags$td(class = "sid", r$sid),
    tags$td(badge, tags$br(), tags$small(r$zh)),
    tags$td(biomarker_pill(fmt_n(r$tmb_val), r$tmb_val_class, "tmb")),
    tags$td(biomarker_pill(fmt_n(r$msi_val), r$msi_status, "msi")),
    tags$td(biomarker_pill(round(r$gis %||% NA, 0), "", "hrd")),
    tags$td(class = "actions",
       if (r$on_label > 0)  tags$span(class="tag on",  paste(r$on_label, "on-label")) else NULL,
       if (r$off_label > 0) tags$span(class="tag off", paste(r$off_label, "off-label")) else NULL,
       if (r$resistance > 0) tags$span(class="tag res", paste(r$resistance, "R")) else NULL,
       if (r$amp_actionable > 0) tags$span(class="tag amp", paste(r$amp_actionable, "AMP I/II")) else NULL),
    tags$td(class = "top-drug", r$top_drug),
    tags$td(class = "links",
            tags$a(href = href_html, "HTML"),
            if (r$has_pdf) tagList(" · ", tags$a(href = href_pdf, "PDF")) else NULL)
  )
}

batch_section <- function(batch, sub) {
  date_str <- gsub("Batch_(\\d+)_(\\d{4})(\\d{2})(\\d{2})",
                   "Batch \\1 — \\2-\\3-\\4", batch)
  tags$details(open = NA,
    tags$summary(
      tags$h2(date_str, " ",
              tags$small(glue("({nrow(sub)} samples)")))),
    tags$table(
      tags$thead(tags$tr(
        tags$th("Sample"), tags$th("Tumor"),
        tags$th("TMB (mut/Mb)"), tags$th("MSI (%)"), tags$th("GIS"),
        tags$th("Findings"), tags$th("Top therapy"), tags$th("Report"))),
      tags$tbody(lapply(seq_len(nrow(sub)), function(i) render_row(sub[i, ])))
    )
  )
}

batches <- split(rows, rows$batch)
batches <- batches[order(names(batches), decreasing = TRUE)]  # newest first

doc <- tags$html(lang = "en",
  tags$head(
    tags$meta(charset = "utf-8"),
    tags$title("KFSYSCC — TSO500-HRD Cancer Genomic Profiles"),
    tags$style(HTML("
      body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
           max-width:1200px;margin:1.5em auto;padding:0 1.5em;color:#222;line-height:1.45}
      h1{color:#0d2c54;border-bottom:3px solid #0d2c54;padding-bottom:.4em;margin-bottom:.3em}
      h1 .sub{font-size:.45em;color:#666;font-weight:400;margin-left:.5em}
      summary{cursor:pointer;list-style:none;padding:.4em .6em;background:#eef3fb;border-radius:5px;margin:1em 0 .4em}
      summary h2{display:inline;color:#1a3e7c;font-size:1.05em;margin:0}
      summary small{color:#666;font-weight:400;margin-left:.5em;font-size:.85em}
      details{margin-bottom:.6em}
      table{width:100%;border-collapse:collapse;margin:.5em 0;font-size:.9em}
      th{background:#0d2c54;color:#fff;text-align:left;padding:.45em .65em;font-weight:600}
      td{padding:.45em .65em;border-bottom:1px solid #e0e6ef;vertical-align:top}
      tr:nth-child(even) td{background:#f6f9fc}
      .sid{font-weight:600;color:#0d2c54;font-family:ui-monospace,Menlo,monospace}
      .tt-badge{display:inline-block;padding:1px 8px;border-radius:3px;font-size:.78em;color:#fff;font-weight:600}
      small{color:#666;font-size:.82em}
      .pill{display:inline-block;padding:1px 9px;border-radius:3px;font-weight:600;font-size:.82em}
      .pill.high{background:#ffefef;color:#b71c1c;border:1px solid #ffcdd2}
      .pill.med {background:#fff7e6;color:#e65100;border:1px solid #ffe0b2}
      .pill.low {background:#f1f8e9;color:#33691e;border:1px solid #dcedc8}
      .tag{display:inline-block;padding:1px 6px;border-radius:3px;font-size:.78em;margin-right:3px;margin-bottom:2px}
      .tag.on{background:#c62828;color:#fff} .tag.off{background:#1565c0;color:#fff}
      .tag.res{background:#6a1b9a;color:#fff} .tag.amp{background:#0d2c54;color:#fff}
      .top-drug{font-size:.82em;font-style:italic;color:#444;max-width:240px}
      .links a{color:#1565c0;text-decoration:none;font-weight:600}
      .links a:hover{text-decoration:underline}
      .note{background:#fff8e1;border-left:4px solid #ffa000;padding:.7em 1em;margin:1em 0;font-size:.88em}
      .summary{display:flex;gap:1em;margin:1em 0;flex-wrap:wrap}
      .summary .stat{background:#f5f7fb;border-radius:5px;padding:.5em .9em;border-left:3px solid #1a3e7c}
      .summary .stat strong{display:block;font-size:1.4em;color:#0d2c54}
    "))
  ),
  tags$body(
    tags$h1("KFSYSCC — TSO500-HRD Cancer Genomic Profiles",
            tags$span(class = "sub",
                      paste0("8 batches · ", nrow(rows), " patients"))),
    tags$p(tags$b("Report generated:"), " ", format(Sys.time(), "%Y-%m-%d"), " · ",
           tags$b("Panel:"), " Illumina TruSight Oncology 500 (523 genes, ~1.94 Mb) · ",
           tags$b("Pipeline:"), " DRAGEN TSO500 v2.6 → Foundation One–style synthesis"),
    tags$div(class = "summary",
      tags$div(class = "stat", tags$strong(nrow(rows)), "patients"),
      tags$div(class = "stat", tags$strong(sum(rows$on_label > 0, na.rm = TRUE)),
               "with on-label therapies"),
      tags$div(class = "stat", tags$strong(sum(rows$gis >= 42, na.rm = TRUE)),
               "HRD-positive (GIS ≥ 42)"),
      tags$div(class = "stat", tags$strong(sum(rows$tmb_val >= 10, na.rm = TRUE)),
               "TMB-High (≥10)"),
      tags$div(class = "stat", tags$strong(sum(grepl("MSI-High", rows$msi_status))),
               "MSI-High")
    ),
    tags$div(class = "note",
      "Foundation One–style HTML reports, self-contained for offline browsing. ",
      "Patient names are NOT in the report bodies — only de-identifiable sample IDs. ",
      "PDFs are headless-Chrome prints of the HTML; HTML is the source of truth."),
    do.call(tagList,
            lapply(names(batches), function(b) batch_section(b, batches[[b]])))
  )
)

out <- file.path(handoff, "index.html")
save_html(doc, out, libdir = NULL)
cat("Wrote", out, "\n")
