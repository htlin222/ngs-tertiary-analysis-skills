#!/usr/bin/env Rscript
# scripts/run_one_sample.R — Run pipeline for a single TSO500 VCF and emit HTML+PDF.
#
# Reads VCF_PATH, SAMPLE_ID, TUMOR_TYPE from env (or args).
# Skips renv library (which is out-of-sync) and uses the system R library.

.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library",
            .libPaths()))

setwd("/Users/htlin/ngs-tertiary-analysis-skills")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) Sys.setenv(SAMPLE_ID = args[[1]])
if (length(args) >= 2) Sys.setenv(VCF_PATH = args[[2]])
if (length(args) >= 3) Sys.setenv(TUMOR_TYPE = args[[3]])

sample_id  <- Sys.getenv("SAMPLE_ID")
vcf_path   <- Sys.getenv("VCF_PATH")
tumor_type <- Sys.getenv("TUMOR_TYPE")

stopifnot(nzchar(sample_id), nzchar(vcf_path), nzchar(tumor_type))
if (!file.exists(vcf_path)) stop("VCF not found: ", vcf_path)

cat(sprintf("\n==== %s | %s | %s ====\n", sample_id, tumor_type, vcf_path))

t0 <- Sys.time()
suppressPackageStartupMessages({
  library(targets)
})

# Sequential controller — avoids crew worker churn for single-sample runs
options(targets.error = "null")

targets::tar_make(reporter = "summary", as_job = FALSE,
                  callr_function = NULL)  # in-process, share env

elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\n==== %s pipeline finished in %.1f min ====\n",
            sample_id, as.numeric(elapsed)))

html_path <- file.path("reports", sample_id, "08-report", "clinical_report.html")
if (!file.exists(html_path)) stop("Report HTML not produced: ", html_path)
cat("HTML: ", html_path, "\n", sep = "")

# ── Render PDF via headless chromium ────────────────────────────────────────
pdf_path <- sub("\\.html$", ".pdf", html_path)

# pagedown::chrome_print requires Chrome; use chromote directly (chromium ok)
suppressPackageStartupMessages(library(chromote))
chromote::set_chrome_args(c(chromote::default_chrome_args(),
                            "--no-sandbox", "--disable-gpu"))
Sys.setenv(CHROMOTE_CHROME = "/opt/homebrew/bin/chromium")

b <- chromote::ChromoteSession$new()
on.exit(try(b$close(), silent = TRUE), add = TRUE)
b$Page$navigate(paste0("file://", normalizePath(html_path)))
b$Page$loadEventFired(timeout_ = 60)
Sys.sleep(2)  # let plotly/JS settle
pdf_b64 <- b$Page$printToPDF(printBackground = TRUE,
                              preferCSSPageSize = TRUE,
                              marginTop = 0.4, marginBottom = 0.4,
                              marginLeft = 0.4, marginRight = 0.4)$data
writeBin(jsonlite::base64_dec(pdf_b64), pdf_path)
cat("PDF : ", pdf_path, " (", file.size(pdf_path), " B)\n", sep = "")
