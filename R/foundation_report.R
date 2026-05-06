# R/foundation_report.R — Foundation One–style actionable HTML report.
#
# Pure-R, no Quarto, no plotly/ggiraph. Self-contained HTML, then headless
# Chrome converts it to PDF (see scripts/render_one_actionable.R).

suppressPackageStartupMessages({
  library(htmltools); library(glue); library(fs); library(jsonlite)
})

source(here::here("R/foundation_helpers.R"))
source(here::here("R/foundation_sections.R"))
source(here::here("R/tso500_parser.R"))     # tumor_type_display()

#' Render the actionable HTML report for a single sample.
#'
#' All inputs are objects already on disk (the `.rds` artefacts produced by
#' the existing pipeline). No network calls inside this function.
#'
#' @param sample_id   Sample identifier (string).
#' @param parsed      Output of `parse_tso500_combined_output()`.
#' @param oncokb      OncoKB results list (mutations / cnas / fusions).
#' @param civic       CiVIC results list (variant_evidence / assertions).
#' @param amp         AMP/ASCO/CAP per-variant tibble (with `amp_evidence`).
#' @param escat       ESCAT classifications tibble (unused at the moment but
#'                    accepted for future per-tier badging).
#' @param literature  Literature narratives list (pubmed / scopus / narratives).
#' @param config      Loaded pipeline config (uses `sample$tumor_type`,
#'                    `report$institution`).
#' @param vaf_png,cnv_png  Optional PNG paths from `R/foundation_plots.R`.
#' @param output_path Override output HTML path; defaults to
#'                    `reports/<sample_id>/08-report/clinical_actionable.html`.
#'
#' @return Path to the rendered HTML file (always self-contained).
render_actionable_report <- function(sample_id, parsed, oncokb, civic, amp,
                                     escat, literature, config,
                                     vaf_png = NULL, cnv_png = NULL,
                                     output_path = NULL) {
  tumor_label <- tumor_type_display(config$sample$tumor_type %||% "Unknown")
  therapies <- extract_therapies(oncokb)

  css_path <- here::here("R/foundation_styles.css")
  css <- if (file.exists(css_path)) {
    paste(readLines(css_path, warn = FALSE), collapse = "\n")
  } else ""

  data_version <- tryCatch({
    v <- oncokb$mutations[[1]]$raw$dataVersion
    if (length(v) > 0) v[[1]] else "—"
  }, error = function(e) "—")

  doc <- tags$html(lang = "en",
    tags$head(
      tags$meta(charset = "utf-8"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tags$title(glue("Clinical Profile — {sample_id}")),
      tags$style(HTML(css))
    ),
    tags$body(
      section_cover(sample_id, tumor_label, parsed, config),
      section_biomarkers(parsed$tmb, parsed$msi, parsed$hrd),
      section_therapies_hero(therapies),
      section_resistance(therapies),
      section_variant_cards(amp, oncokb, civic, parsed = parsed),
      section_lit_review(sample_id),
      section_cnv(parsed$cnv),
      section_figures(vaf_png, cnv_png),
      section_gene_narratives(literature, oncokb),
      section_vus(amp),
      tags$footer(class = "report-footer",
        tags$p(class = "fine",
          glue("Generated {format(Sys.time(), '%Y-%m-%d %H:%M %Z')} · ",
               "OncoKB v{data_version} · CiVIC GraphQL · NCBI PubMed/Scopus")),
        tags$p(class = "fine",
          "TruSight Oncology 500 — Research Use Only. Tumor fraction, ploidy, ",
          "absolute copy number and LoH are beta features per Illumina guidance.")
      )
    )
  )

  if (is.null(output_path)) {
    output_path <- here::here("reports", sample_id, "08-report",
                              "clinical_actionable.html")
  }
  fs::dir_create(dirname(output_path), recurse = TRUE)
  htmltools::save_html(doc, output_path, libdir = NULL)
  invisible(output_path)
}
