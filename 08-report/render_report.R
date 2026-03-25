# 08-report/render_report.R — Render the clinical genomic profiling report

suppressPackageStartupMessages({
  library(fs)
  library(glue)
  library(logger)
})

source(here::here("R/utils.R"))

#' Render the clinical genomic profiling report
#'
#' Saves all pipeline stage results as .rds files for the Quarto document to
#' consume, renders the parameterized Quarto report, and copies the output
#' HTML to the final reports directory.
#'
#' @param sample_id Character. Sample identifier.
#' @param config List. Pipeline configuration (from load_config).
#' @param qc List or tibble. QC metrics from stage 00.
#' @param variants Tibble. Annotated somatic variants from stages 01-02.
#' @param cnv Tibble. Copy number alterations from stage 03.
#' @param fusions Tibble. Gene fusions from stage 04.
#' @param tmb List. TMB results from stage 05 (score, classification, details).
#' @param msi List. MSI results from stage 05 (score, status, details).
#' @param hrd List. HRD results from stage 05 (score, status, details).
#' @param oncokb List. OncoKB annotation results from stage 06.
#' @param escat List. ESCAT classification results from stage 06.
#' @param literature List or tibble. Literature narratives from stage 07.
#' @param civic List. CiVIC community evidence results (default NULL).
#' @param amp Tibble. AMP/ASCO/CAP classification results (default NULL).
#'
#' @return Character. Path to the rendered HTML report.
#'
#' @details
#' The function:
#' 1. Creates a data directory under `reports/{sample_id}/08-report/data/`
#' 2. Saves each input as an .rds file for the Quarto document to read
#' 3. Renders `08-report/clinical_report.qmd` with quarto::quarto_render()
#' 4. Copies the output HTML to `reports/{sample_id}/08-report/clinical_report.html`
#' 5. Returns the path to the final HTML file
#'
#' @export
render_report <- function(sample_id,
                          config,
                          qc,
                          variants,
                          cnv,
                          fusions,
                          tmb,
                          msi,
                          hrd,
                          oncokb,
                          escat,
                          literature,
                          civic = NULL,
                          amp = NULL) {
  log_info("Rendering clinical report for sample: {sample_id}")


  # ── 1. Prepare data directory ───────────────────────────────────────────────
  reports_dir <- here::here("reports")
  data_dir <- path(reports_dir, sample_id, "08-report", "data")
  dir_create(data_dir, recurse = TRUE)
  log_info("Report data directory: {data_dir}")

  # ── 2. Save all input data as .rds files ────────────────────────────────────
  datasets <- list(
    config     = config,
    qc         = qc,
    variants   = variants,
    cnv        = cnv,
    fusions    = fusions,
    tmb        = tmb,
    msi        = msi,
    hrd        = hrd,
    oncokb     = oncokb,
    escat      = escat,
    literature = literature,
    civic      = civic,
    amp        = amp
  )

  for (name in names(datasets)) {
    rds_path <- path(data_dir, glue("{name}.rds"))
    saveRDS(datasets[[name]], rds_path)
    log_debug("Saved {name}.rds ({object.size(datasets[[name]])} bytes)")
  }
  log_info("Saved {length(datasets)} data files to {data_dir}")

  # ── 3. Render the Quarto document ──────────────────────────────────────────
  qmd_path <- here::here("08-report", "clinical_report.qmd")
  if (!file_exists(qmd_path)) {
    log_error("Quarto template not found: {qmd_path}")
    stop(glue("Quarto template not found: {qmd_path}"))
  }

  output_dir <- path(reports_dir, sample_id, "08-report")
  dir_create(output_dir, recurse = TRUE)

  render_params <- list(
    sample_id   = sample_id,
    config_path = here::here("config/default.yaml"),
    reports_dir = as.character(reports_dir)
  )

  log_info("Rendering Quarto report with params: sample_id={sample_id}")

  tryCatch(
    {
      quarto::quarto_render(
        input       = as.character(qmd_path),
        output_file = "clinical_report.html",
        execute_params = render_params,
        quiet       = TRUE
      )
    },
    error = function(e) {
      log_error("Quarto render failed: {e$message}")
      stop(glue("Report rendering failed for {sample_id}: {e$message}"))
    }
  )

  # ── 4. Move rendered HTML to final location ────────────────────────────────
  # quarto_render outputs next to the .qmd by default; move to reports dir
  rendered_path <- path(here::here("08-report"), "clinical_report.html")
  final_path <- path(output_dir, "clinical_report.html")

  if (file_exists(rendered_path)) {
    file_copy(rendered_path, final_path, overwrite = TRUE)
    file_delete(rendered_path)
    log_info("Report moved to: {final_path}")
  } else if (file_exists(final_path)) {
    log_info("Report already at final location: {final_path}")
  } else {
    log_warn("Rendered HTML not found at expected location; checking alternatives")
    # quarto may output to the working directory
    alt_path <- path(getwd(), "clinical_report.html")
    if (file_exists(alt_path)) {
      file_copy(alt_path, final_path, overwrite = TRUE)
      file_delete(alt_path)
      log_info("Report moved from working directory to: {final_path}")
    } else {
      log_error("Could not locate rendered report HTML")
      stop("Rendered report HTML not found after quarto_render")
    }
  }

  # ── 5. Return final path ───────────────────────────────────────────────────
  log_info("Clinical report rendered successfully: {final_path}")
  as.character(final_path)
}
