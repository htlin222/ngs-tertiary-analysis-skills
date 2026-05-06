# R/foundation_plots.R — Static (no plotly / no ggiraph) figure helpers
# for the Foundation One–style actionable report.
#
# Each helper renders a PNG via ggsave() and returns the file path; the
# report module embeds the PNG as a base64 data URI so the final HTML stays
# self-contained.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(stringr); library(tibble)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

PATHO_COLORS <- c(
  "pathogenic"             = "#c62828",
  "likely_pathogenic"      = "#e53935",
  "uncertain_significance" = "#f9a825",
  "VUS"                    = "#f9a825",
  "likely_benign"          = "#66bb6a",
  "benign"                 = "#2e7d32"
)

CNV_COLORS <- c(
  "AMPLIFICATION" = "#c62828",
  "GAIN"          = "#ef6c00",
  "DELETION"      = "#1565c0",
  "LOSS"          = "#42a5f5",
  "NEUTRAL"       = "#bdbdbd"
)

#' Static lollipop of top-N variants by VAF.
#'
#' @param variants Tibble — needs gene, hgvsp, vaf, consequence, pathogenicity.
#' @param min_vaf  Numeric — drawn as a dashed reference line.
#' @param output_png Character — destination file path.
#' @param top_n  Integer — keep the top N rows by VAF (after deduplication).
plot_vaf_static <- function(variants, min_vaf = 0.05, output_png, top_n = 30) {
  if (!is.data.frame(variants) || nrow(variants) == 0) {
    df <- tibble(label = "(no variants)", vaf = 0, pathogenicity = "VUS")
  } else {
    df <- variants |>
      filter(!is.na(vaf), consequence != "amplification") |>
      arrange(desc(vaf)) |>
      slice_head(n = top_n) |>
      mutate(label = paste0(gene, " ", coalesce(hgvsp, consequence)),
             label = make.unique(label, sep = " #"),
             pathogenicity = ifelse(is.na(pathogenicity),
                                     "uncertain_significance", pathogenicity)) |>
      mutate(label = factor(label, levels = rev(label)))
  }

  p <- ggplot(df, aes(x = vaf, y = label, color = pathogenicity)) +
    geom_segment(aes(x = 0, xend = vaf, y = label, yend = label),
                 color = "grey75", linewidth = 0.4) +
    geom_point(size = 3) +
    geom_vline(xintercept = min_vaf, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    scale_x_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.05))) +
    scale_color_manual(values = PATHO_COLORS, na.value = "grey60", drop = FALSE) +
    labs(x = "Variant Allele Frequency", y = NULL, color = "Classification") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))

  fs::dir_create(dirname(output_png), recurse = TRUE)
  ggsave(output_png, p, width = 8, height = max(4, min(0.22 * nrow(df) + 1.5, 9)),
         dpi = 150, bg = "white")
  invisible(output_png)
}

#' Static bar chart: gene-level CNV log2 fold change, ordered by magnitude.
#'
#' @param cnv_data Tibble from `parsed$cnv` — needs gene, log2_ratio, type.
#' @param output_png Destination path.
plot_cnv_static <- function(cnv_data, output_png, top_n = 60) {
  if (!is.data.frame(cnv_data) || nrow(cnv_data) == 0) {
    df <- tibble(gene = "(none)", log2_ratio = 0, type = "NEUTRAL")
  } else {
    df <- cnv_data |>
      filter(!is.na(log2_ratio)) |>
      arrange(desc(abs(log2_ratio))) |>
      slice_head(n = top_n) |>
      mutate(type = ifelse(is.na(type), "NEUTRAL", type)) |>
      mutate(gene = factor(gene, levels = gene[order(log2_ratio)]))
  }

  p <- ggplot(df, aes(x = gene, y = log2_ratio, fill = type)) +
    geom_col() +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = c(0.58, -1.0), linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    scale_fill_manual(values = CNV_COLORS, na.value = "grey60") +
    coord_flip() +
    labs(x = NULL, y = "Log2 Fold Change", fill = "Type") +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          axis.text.y = element_text(size = 7))

  fs::dir_create(dirname(output_png), recurse = TRUE)
  ggsave(output_png, p, width = 8, height = max(4, min(0.18 * nrow(df) + 1.5, 11)),
         dpi = 150, bg = "white")
  invisible(output_png)
}
