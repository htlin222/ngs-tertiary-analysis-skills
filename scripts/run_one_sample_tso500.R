#!/usr/bin/env Rscript
# scripts/run_one_sample_tso500.R
#
# Per-sample driver: parses the Illumina TSO500 CombinedVariantOutput.tsv,
# calls OncoKB + CiVIC + ESCAT + AMP + literature, renders the Quarto report,
# and prints to PDF via headless chromium. Bypasses the targets DAG (and VEP)
# because TSO500 already provides full annotation.
#
# Usage:
#   Rscript --vanilla scripts/run_one_sample_tso500.R \
#     <SAMPLE_ID> <CombinedVariantOutput.tsv> <OncoKB-tumor-code>
#
# Example:
#   Rscript --vanilla scripts/run_one_sample_tso500.R \
#     M26-0165R inputs/.../M26-0165R_CombinedVariantOutput.tsv UCEC

.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library",
            .libPaths()))
setwd("/Users/htlin/ngs-tertiary-analysis-skills")

# The project's renv is out-of-sync — make sure quarto's child R also skips it.
Sys.setenv(
  RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
  R_LIBS_USER = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_one_sample_tso500.R <SAMPLE_ID> <CombinedVariantOutput.tsv> <TUMOR_CODE>")
}
sample_id <- args[[1]]
tsv_path  <- args[[2]]
tumor_type <- args[[3]]
Sys.setenv(SAMPLE_ID = sample_id, TUMOR_TYPE = tumor_type)

cat(sprintf("\n==== %s | %s | %s ====\n", sample_id, tumor_type, tsv_path))
t0 <- Sys.time()

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr); library(glue)
  library(readr); library(purrr); library(jsonlite); library(logger); library(here)
  library(yaml); library(httr2); library(fs)
})

source(here("R/utils.R"))
source(here("R/tso500_parser.R"))
source(here("R/oncokb_helpers.R"))
source(here("R/api_clients.R"))
source(here("R/civic_client.R"))
source(here("R/amp_classification.R"))
source(here("R/esmo_helpers.R"))
source(here("06-clinical-annotation/query_oncokb.R"))
source(here("06-clinical-annotation/query_civic.R"))
source(here("06-clinical-annotation/classify_escat.R"))
source(here("07-literature/generate_narrative.R"))
source(here("08-report/render_report.R"))

load_env()
config <- load_config()  # picks up TUMOR_TYPE / SAMPLE_ID env overrides

# ── 1. Parse TSO500 deliverable ─────────────────────────────────────────────
parsed <- parse_tso500_combined_output(tsv_path)
log_info("Parsed TSO500 — variants:{nrow(parsed$variants)} cnv:{nrow(parsed$cnv)} fusions:{nrow(parsed$fusions)}")
cat(sprintf("Variants: %d | CNV: %d | Fusions: %d | TMB: %s | MSI: %s | GIS: %s\n",
            nrow(parsed$variants), nrow(parsed$cnv), nrow(parsed$fusions),
            parsed$tmb$tmb_score, parsed$msi$msi_status, parsed$hrd$hrd_score))

# Persist parsed inputs into the per-sample reports tree (audit trail)
sample_dir <- here("reports", sample_id)
dir_create(file.path(sample_dir, c("00-input", "02-annotation",
                                    "03-cnv", "04-fusions", "05-biomarkers")),
           recurse = TRUE)
write_tsv(parsed$variants, file.path(sample_dir, "02-annotation/merged_annotations.tsv"))
write_tsv(parsed$cnv,      file.path(sample_dir, "03-cnv/cnv_calls.tsv"))
write_tsv(parsed$fusions,  file.path(sample_dir, "04-fusions/fusions.tsv"))
write_json(parsed$tmb, file.path(sample_dir, "05-biomarkers/tmb_result.json"),
           pretty = TRUE, auto_unbox = TRUE)
write_json(parsed$msi, file.path(sample_dir, "05-biomarkers/msi_result.json"),
           pretty = TRUE, auto_unbox = TRUE)
write_json(parsed$hrd, file.path(sample_dir, "05-biomarkers/hrd_result.json"),
           pretty = TRUE, auto_unbox = TRUE)

# ── 2. Clinical annotation: OncoKB + CiVIC ─────────────────────────────────
oncokb <- tryCatch(
  query_oncokb(parsed$variants, parsed$cnv, parsed$fusions, config, sample_id),
  error = function(e) {
    log_warn("OncoKB query failed: {e$message}")
    list(mutations = list(), cnas = list(), fusions = list())
  }
)
civic <- tryCatch(
  query_civic(parsed$variants, parsed$cnv, parsed$fusions, config, sample_id),
  error = function(e) {
    log_warn("CiVIC query failed: {e$message}")
    list(variant_evidence = tibble(), gene_summaries = list(), assertions = tibble())
  }
)

# ── 3. ESCAT + AMP classification ──────────────────────────────────────────
escat_raw <- tryCatch(
  classify_escat(oncokb, config, sample_id),
  error = function(e) {
    log_warn("ESCAT failed: {e$message}")
    NULL
  }
)
# classify_escat returns a list of tibbles (escat_results, tier_summary, ...);
# the report template expects the per-alteration tibble.
escat <- if (is.data.frame(escat_raw)) {
  escat_raw
} else if (is.list(escat_raw) && !is.null(escat_raw$escat_results)) {
  escat_raw$escat_results
} else {
  tibble(gene = character(), alteration = character(), type = character(),
         oncogenic = character(), sensitive_level = character(),
         escat_tier = character(), escat_description = character())
}

variants_for_amp <- parsed$variants
if (length(oncokb$mutations) > 0) {
  oncokb_lookup <- bind_rows(map(oncokb$mutations, function(m) {
    if (is.list(m) && !is.null(m$gene)) {
      tibble(gene = m$gene,
             alteration_short = m$alteration %||% NA_character_,
             oncogenic = m$oncogenic %||% NA_character_,
             sensitive_level = m$highest_sensitive_level %||% NA_character_)
    } else NULL
  }))
  # Variants carry HGVS.p in 3-letter form (e.g. "p.Gly12Ser") but OncoKB
  # normalises to 1-letter ("G12S"); join on the short form.
  variants_for_amp <- variants_for_amp |>
    mutate(alteration_short = hgvsp_to_short(as.character(hgvsp))) |>
    left_join(oncokb_lookup, by = c("gene", "alteration_short"))
} else {
  variants_for_amp <- variants_for_amp |>
    mutate(oncogenic = NA_character_, sensitive_level = NA_character_)
}

amp <- tryCatch(
  classify_all_amp(variants_for_amp,
                   if (!is.null(civic$assertions) && nrow(civic$assertions) > 0)
                     civic$assertions else NULL),
  error = function(e) {
    log_warn("AMP classification failed: {e$message}")
    tibble(gene = character(), hgvsp = character(),
           amp_tier = character(), amp_level = character(), amp_evidence = character())
  }
)

# ── 4. Literature narrative ────────────────────────────────────────────────
literature <- tryCatch(
  generate_narrative(parsed$variants, parsed$cnv, parsed$fusions,
                     config, sample_id, oncokb = oncokb),
  error = function(e) {
    log_warn("Literature generation failed: {e$message}")
    list(pubmed = tibble(), scopus = tibble(), narratives = list())
  }
)

# ── 5. Render report ───────────────────────────────────────────────────────
html_path <- render_report(
  sample_id = sample_id, config = config,
  qc = parsed$qc, variants = parsed$variants,
  cnv = parsed$cnv, fusions = parsed$fusions,
  tmb = parsed$tmb, msi = parsed$msi, hrd = parsed$hrd,
  oncokb = oncokb, escat = escat, literature = literature,
  civic = civic, amp = amp
)
cat("HTML : ", html_path, " (", round(file.size(html_path) / 1024),
    " KB, self-contained — opens in any browser)\n", sep = "")

# ── 6. Also produce the actionable Foundation One–style HTML+PDF ──────────
source(here("R/foundation_helpers.R"))
source(here("R/foundation_sections.R"))
source(here("R/foundation_report.R"))
source(here("R/foundation_plots.R"))

fig_dir <- here("reports", sample_id, "08-report", "figures")
fs::dir_create(fig_dir, recurse = TRUE)
vaf_png <- file.path(fig_dir, "top_vaf.png")
cnv_png <- file.path(fig_dir, "cnv_log2.png")
plot_vaf_static(parsed$variants, min_vaf = config$variant_calling$min_vaf %||% 0.05,
                output_png = vaf_png, top_n = 30)
plot_cnv_static(parsed$cnv, output_png = cnv_png, top_n = 60)

action_html <- render_actionable_report(
  sample_id = sample_id, parsed = parsed,
  oncokb = oncokb, civic = civic, amp = amp, escat = escat,
  literature = literature, config = config,
  vaf_png = vaf_png, cnv_png = cnv_png
)
cat("HTML+: ", action_html, " (",
    round(file.size(action_html) / 1024),
    " KB, actionable Foundation One layout)\n", sep = "")

chrome <- "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
if (file.exists(chrome)) {
  pdf_path <- sub("\\.html$", ".pdf", action_html)
  rc <- system2(chrome, c("--headless=new", "--disable-gpu", "--no-sandbox",
                          "--no-pdf-header-footer",
                          paste0("--print-to-pdf=", pdf_path),
                          paste0("file://", normalizePath(action_html))),
                stdout = FALSE, stderr = FALSE)
  if (rc == 0 && file.exists(pdf_path)) {
    cat("PDF  : ", pdf_path, " (", round(file.size(pdf_path) / 1024),
        " KB)\n", sep = "")
  }
}

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("==== %s done in %.1f min ====\n", sample_id, as.numeric(elapsed)))
