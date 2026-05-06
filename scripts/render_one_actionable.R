#!/usr/bin/env Rscript
# scripts/render_one_actionable.R
#
# Render the Foundation One–style actionable report for ONE sample.
# Reuses the cached `.rds` artefacts from a previous pipeline run; if those
# don't exist yet, falls back to running run_one_sample_tso500.R first.
#
# Usage:
#   Rscript --vanilla scripts/render_one_actionable.R <SAMPLE_ID> <TUMOR_CODE>
# Example:
#   Rscript --vanilla scripts/render_one_actionable.R M26-0165R UCEC

.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library",
            .libPaths()))
setwd("/Users/htlin/ngs-tertiary-analysis-skills")
Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  R_LIBS_USER = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: render_one_actionable.R <SAMPLE_ID> <TUMOR_CODE>")
sample_id  <- args[[1]]
tumor_type <- args[[2]]
Sys.setenv(SAMPLE_ID = sample_id, TUMOR_TYPE = tumor_type)

t0 <- Sys.time()
cat(sprintf("\n==== %s (%s) actionable render ====\n", sample_id, tumor_type))

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tibble); library(stringr); library(glue)
  library(htmltools); library(jsonlite); library(yaml); library(fs); library(logger)
  library(tidyr); library(purrr); library(httr2)
})

source(here("R/utils.R"))
source(here("R/tso500_parser.R"))
source(here("R/oncokb_helpers.R"))
source(here("R/foundation_helpers.R"))
source(here("R/foundation_sections.R"))
source(here("R/foundation_report.R"))
source(here("R/foundation_plots.R"))

load_env()
config <- load_config()

data_dir <- here("reports", sample_id, "08-report", "data")
batch_tsv <- here("inputs", "TSO500-HRD", "Batch_7_20260505",
                   paste0(sample_id, "_CombinedVariantOutput.tsv"))

if (!dir_exists(data_dir) || !file.exists(file.path(data_dir, "oncokb.rds"))) {
  cat("No cached .rds found — running full pipeline first…\n")
  driver <- here("scripts/run_one_sample_tso500.R")
  status <- system2("Rscript", c("--vanilla", driver, sample_id,
                                 batch_tsv, tumor_type),
                    stdout = "", stderr = "")
  if (status != 0) stop("run_one_sample_tso500.R exited with status ", status)
}

# ── Load cached artefacts ───────────────────────────────────────────────────
parsed <- parse_tso500_combined_output(batch_tsv)
oncokb     <- readRDS(file.path(data_dir, "oncokb.rds"))
civic      <- readRDS(file.path(data_dir, "civic.rds"))
amp        <- readRDS(file.path(data_dir, "amp.rds"))
escat_raw  <- readRDS(file.path(data_dir, "escat.rds"))
literature <- readRDS(file.path(data_dir, "literature.rds"))

escat <- if (is.data.frame(escat_raw)) {
  escat_raw
} else if (is.list(escat_raw)) {
  escat_raw$escat_results %||% tibble()
} else {
  tibble()
}

# ── Static figures ──────────────────────────────────────────────────────────
fig_dir <- here("reports", sample_id, "08-report", "figures")
dir_create(fig_dir, recurse = TRUE)
vaf_png <- file.path(fig_dir, "top_vaf.png")
cnv_png <- file.path(fig_dir, "cnv_log2.png")
plot_vaf_static(parsed$variants, min_vaf = config$variant_calling$min_vaf %||% 0.05,
                output_png = vaf_png, top_n = 30)
plot_cnv_static(parsed$cnv, output_png = cnv_png, top_n = 60)

# ── Render HTML ─────────────────────────────────────────────────────────────
html_path <- render_actionable_report(
  sample_id = sample_id, parsed = parsed,
  oncokb = oncokb, civic = civic, amp = amp, escat = escat,
  literature = literature, config = config,
  vaf_png = vaf_png, cnv_png = cnv_png
)
cat(sprintf("HTML : %s (%.1f KB)\n", html_path, file.size(html_path)/1024))

# ── HTML → PDF via headless Google Chrome ───────────────────────────────────
pdf_path <- sub("\\.html$", ".pdf", html_path)
chrome <- "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
if (!file.exists(chrome)) {
  cat("Chrome not found at expected path, skipping PDF.\n")
} else {
  args_pdf <- c("--headless=new", "--disable-gpu", "--no-sandbox",
                "--no-pdf-header-footer",
                paste0("--print-to-pdf=", pdf_path),
                paste0("file://", normalizePath(html_path)))
  status <- system2(chrome, args_pdf, stdout = FALSE, stderr = FALSE)
  if (status == 0 && file.exists(pdf_path)) {
    cat(sprintf("PDF  : %s (%.1f KB)\n", pdf_path, file.size(pdf_path)/1024))
  } else {
    cat("PDF render failed (exit code ", status, ")\n")
  }
}

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("==== %s done in %.1fs ====\n", sample_id, elapsed))
