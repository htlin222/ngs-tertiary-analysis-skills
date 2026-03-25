# R/plot_helpers.R — Interactive visualization functions for clinical report

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggiraph)
  library(dplyr)
  library(glue)
  library(tidyr)
  library(stringr)
})

# ── Color palettes ──────────────────────────────────────────────────────────

pathogenicity_colors <- c(
  "pathogenic"        = "#c62828",
  "likely_pathogenic" = "#e53935",
  "VUS"               = "#f9a825",
  "likely_benign"     = "#66bb6a",
  "benign"            = "#2e7d32"
)

cnv_type_colors <- c(
  "AMPLIFICATION" = "#c62828",
  "GAIN"          = "#ef6c00",
  "DELETION"      = "#1565c0",
  "LOSS"          = "#42a5f5"
)

# ── GRCh38 chromosome sizes ────────────────────────────────────────────────

grch38_sizes <- c(
  "chr1" = 248956422, "chr2" = 242193529, "chr3" = 198295559,
  "chr4" = 190214555, "chr5" = 181538259, "chr6" = 170805979,
  "chr7" = 159345973, "chr8" = 145138636, "chr9" = 138394717,
  "chr10" = 133797422, "chr11" = 135086622, "chr12" = 133275309,
  "chr13" = 114364328, "chr14" = 107043718, "chr15" = 101991189,
  "chr16" = 90338345, "chr17" = 83257441, "chr18" = 80373285,
  "chr19" = 58617616, "chr20" = 64444167, "chr21" = 46709983,
  "chr22" = 50818468, "chrX" = 156040895, "chrY" = 57227415
)

# ── 1. VAF Distribution (ggiraph lollipop) ─────────────────────────────────

#' Plot variant allele frequency distribution as interactive lollipop chart
#' @param variants Tibble with gene, hgvsp, vaf, consequence, pathogenicity
#' @param min_vaf Minimum VAF threshold (dashed line)
#' @return girafe object
plot_vaf_distribution <- function(variants, min_vaf = 0.05) {
  if (nrow(variants) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No variants to display",
               size = 5, color = "grey50") +
      theme_void()
    return(girafe(ggobj = p))
  }

  df <- variants |>
    mutate(
      label = paste0(gene, " ", hgvsp),
      tooltip = glue(
        "<b>{gene} {hgvsp}</b><br>",
        "VAF: {round(vaf * 100, 1)}%<br>",
        "Consequence: {consequence}<br>",
        "Classification: {pathogenicity}"
      ),
      pathogenicity = factor(pathogenicity,
        levels = names(pathogenicity_colors))
    ) |>
    arrange(desc(vaf)) |>
    mutate(label = factor(label, levels = rev(label)))

  p <- ggplot(df, aes(x = vaf, y = label, color = pathogenicity)) +
    geom_segment(aes(x = 0, xend = vaf, y = label, yend = label),
                 color = "grey70", linewidth = 0.4) +
    geom_point_interactive(
      aes(tooltip = tooltip, data_id = label),
      size = 4
    ) +
    geom_vline(xintercept = min_vaf, linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
    scale_color_manual(values = pathogenicity_colors, drop = FALSE, na.value = "grey50") +
    labs(x = "Variant Allele Frequency", y = NULL, color = "Classification") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )

  girafe(ggobj = p,
    options = list(
      opts_hover(css = "fill:orange;stroke:orange;r:6pt;"),
      opts_tooltip(css = "background:white;padding:8px;border-radius:4px;border:1px solid #ccc;font-size:12px;")
    ))
}

# ── 2. Per-Gene Coverage (ggiraph bar) ─────────────────────────────────────

#' Plot per-gene coverage as interactive horizontal bar chart
#' @param coverage_df Tibble with gene, mean_coverage, min_coverage, pct_above_200x
#' @param min_threshold Coverage threshold (dashed line)
#' @return girafe object
plot_gene_coverage <- function(coverage_df, min_threshold = 200) {
  if (nrow(coverage_df) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No coverage data",
               size = 5, color = "grey50") +
      theme_void()
    return(girafe(ggobj = p))
  }

  # Keep top 30 by coverage + all failing genes
  failing <- coverage_df |> filter(mean_coverage < min_threshold)
  top30 <- coverage_df |> arrange(desc(mean_coverage)) |> head(30)
  df <- bind_rows(top30, failing) |>
    distinct(gene, .keep_all = TRUE) |>
    mutate(
      status = ifelse(mean_coverage >= min_threshold, "PASS", "FAIL"),
      tooltip = glue(
        "<b>{gene}</b><br>",
        "Mean: {round(mean_coverage, 0)}x<br>",
        "Min: {round(min_coverage, 0)}x<br>",
        "% >= {min_threshold}x: {round(pct_above_200x, 1)}%"
      )
    ) |>
    arrange(mean_coverage) |>
    mutate(gene = factor(gene, levels = gene))

  p <- ggplot(df, aes(x = mean_coverage, y = gene, fill = status)) +
    geom_col_interactive(
      aes(tooltip = tooltip, data_id = gene),
      width = 0.7
    ) +
    geom_vline(xintercept = min_threshold, linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    scale_fill_manual(values = c("PASS" = "#66bb6a", "FAIL" = "#ef5350")) +
    labs(x = "Mean Coverage (x)", y = NULL, fill = "Status") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  girafe(ggobj = p,
    options = list(
      opts_hover(css = "fill:orange;"),
      opts_tooltip(css = "background:white;padding:8px;border-radius:4px;border:1px solid #ccc;font-size:12px;")
    ))
}

# ── 3. CNV Genome View (plotly scatter) ─────────────────────────────────────

#' Plot genome-wide CNV profile as interactive plotly scatter
#' @param cnv_data Tibble with chromosome, start, end, gene, log2_ratio, copy_number, length, type
#' @return plotly object
plot_cnv_genome_interactive <- function(cnv_data) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly package required for plot_cnv_genome_interactive()")
  }

  if (nrow(cnv_data) == 0) {
    p <- plotly::plot_ly() |>
      plotly::layout(
        title = "No CNV data available",
        xaxis = list(visible = FALSE),
        yaxis = list(visible = FALSE)
      )
    return(p)
  }

  # Compute chromosome offsets
  chr_order <- paste0("chr", c(1:22, "X", "Y"))
  offsets <- cumsum(c(0, grch38_sizes[chr_order][-length(chr_order)]))
  names(offsets) <- chr_order

  df <- cnv_data |>
    mutate(
      chr_label = ifelse(grepl("^chr", chromosome), chromosome, paste0("chr", chromosome)),
      midpoint = (start + end) / 2,
      genome_pos = offsets[chr_label] + midpoint,
      hover_text = glue(
        "<b>{gene}</b><br>",
        "Type: {type}<br>",
        "{chr_label}:{format(start, big.mark=',')}—{format(end, big.mark=',')}<br>",
        "Log2 ratio: {round(log2_ratio, 2)}<br>",
        "Copy number: {round(copy_number, 1)}"
      ),
      color = cnv_type_colors[type]
    )

  # Chromosome midpoints for axis labels
  chr_mids <- offsets + grch38_sizes[chr_order] / 2
  chr_labels <- gsub("chr", "", chr_order)

  p <- plotly::plot_ly(df,
    x = ~genome_pos,
    y = ~log2_ratio,
    color = ~type,
    colors = cnv_type_colors,
    text = ~hover_text,
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    marker = list(size = 8, opacity = 0.8)
  ) |>
    plotly::layout(
      title = "Genome-Wide Copy Number Profile",
      xaxis = list(
        title = "",
        tickvals = as.numeric(chr_mids),
        ticktext = chr_labels,
        tickangle = 0
      ),
      yaxis = list(title = "Log2 Ratio"),
      shapes = list(
        list(type = "line", x0 = 0, x1 = max(df$genome_pos, na.rm = TRUE),
             y0 = 1.5, y1 = 1.5, line = list(dash = "dash", color = "#c62828", width = 1)),
        list(type = "line", x0 = 0, x1 = max(df$genome_pos, na.rm = TRUE),
             y0 = -1.0, y1 = -1.0, line = list(dash = "dash", color = "#1565c0", width = 1))
      ),
      legend = list(orientation = "h", y = -0.15)
    )

  p
}

# ── 4. Circos Plot (circlize static PNG) ────────────────────────────────────

#' Generate circos plot with CNV and fusion tracks
#' @param cnv_data Tibble with chromosome, start, end, gene, log2_ratio, type
#' @param fusions Tibble with gene_a, gene_b, chr_a, pos_a, chr_b, pos_b, supporting_reads
#' @param output_file Path for output PNG
#' @return output_file path (invisibly)
plot_circos <- function(cnv_data, fusions, output_file) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("circlize package required for plot_circos()")
  }

  # Ensure chr prefix
  ensure_chr <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))

  png(output_file, width = 2400, height = 2400, res = 300)
  on.exit(dev.off(), add = TRUE)

  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = 2)
  circlize::circos.initializeWithIdeogram(
    species = "hg38",
    chromosome.index = paste0("chr", c(1:22, "X", "Y")),
    plotType = c("ideogram", "labels")
  )

  # Track 1: CNV log2 ratio points
  if (nrow(cnv_data) > 0) {
    cnv_bed <- data.frame(
      chr   = ensure_chr(cnv_data$chromosome),
      start = as.numeric(cnv_data$start),
      end   = as.numeric(cnv_data$end),
      value = cnv_data$log2_ratio,
      type  = cnv_data$type,
      stringsAsFactors = FALSE
    )

    cnv_colors <- ifelse(cnv_bed$type == "AMPLIFICATION" | cnv_bed$type == "GAIN",
                         "#c62828", "#1565c0")

    circlize::circos.genomicTrack(
      cnv_bed[, 1:4],
      panel.fun = function(region, value, ...) {
        idx <- which(cnv_bed$chr == circlize::get.cell.meta.data("sector.index") &
                     cnv_bed$start == region[[1]] &
                     cnv_bed$end == region[[2]])
        col <- cnv_colors[idx]
        if (length(col) == 0) col <- "grey50"
        circlize::circos.genomicPoints(region, value, pch = 16, cex = 0.6, col = col)
      },
      track.height = 0.15,
      ylim = c(-3, 3)
    )
  }

  # Track 2: Fusion links
  if (nrow(fusions) > 0) {
    for (i in seq_len(nrow(fusions))) {
      f <- fusions[i, ]
      chr_a <- ensure_chr(f$chr_a)
      chr_b <- ensure_chr(f$chr_b)
      circlize::circos.link(
        chr_a, as.numeric(f$pos_a),
        chr_b, as.numeric(f$pos_b),
        col = adjustcolor("#7b1fa2", alpha.f = 0.5),
        lwd = 2
      )
    }
  }

  circlize::circos.clear()

  invisible(output_file)
}

# ── 5. Fusion Arcs (ggiraph) ───────────────────────────────────────────────

#' Plot fusion events as interactive points
#' @param fusions Tibble with gene_a, gene_b, supporting_reads, fusion_type, known_fusion
#' @return girafe object
plot_fusion_arcs <- function(fusions) {
  if (nrow(fusions) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No fusions detected",
               size = 5, color = "grey50") +
      theme_void()
    return(girafe(ggobj = p))
  }

  df <- fusions |>
    mutate(
      fusion_name = paste0(gene_a, "::", gene_b),
      known_label = ifelse(known_fusion, "Known", "Novel"),
      tooltip = glue(
        "<b>{fusion_name}</b><br>",
        "Supporting reads: {supporting_reads}<br>",
        "Type: {fusion_type}<br>",
        "Status: {known_label}"
      ),
      y = seq_len(n())
    )

  p <- ggplot(df, aes(x = supporting_reads, y = reorder(fusion_name, supporting_reads),
                       color = known_label)) +
    geom_point_interactive(
      aes(tooltip = tooltip, data_id = fusion_name),
      size = 5
    ) +
    geom_segment(aes(x = 0, xend = supporting_reads,
                     y = reorder(fusion_name, supporting_reads),
                     yend = reorder(fusion_name, supporting_reads)),
                 color = "grey70", linewidth = 0.4) +
    scale_color_manual(values = c("Known" = "#2e7d32", "Novel" = "#ef6c00")) +
    labs(x = "Supporting Reads", y = NULL, color = "Fusion Status") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )

  girafe(ggobj = p,
    options = list(
      opts_hover(css = "fill:orange;stroke:orange;r:7pt;"),
      opts_tooltip(css = "background:white;padding:8px;border-radius:4px;border:1px solid #ccc;font-size:12px;")
    ))
}

# ── 6. Biomarker Gauges (ggiraph faceted bars) ─────────────────────────────

#' Plot biomarker values as interactive gauge-style bar charts
#' @param tmb TMB result list with tmb_score, tmb_class, variant_count
#' @param msi MSI result list with msi_score, msi_status, unstable_sites, total_sites
#' @param hrd HRD result list with hrd_score, hrd_status, loh_score, tai_score, lst_score
#' @return girafe object
plot_biomarker_gauges <- function(tmb, msi, hrd) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

  # Extract scores with safe defaults

  tmb_score    <- tmb$tmb_score %||% 0
  tmb_class    <- tmb$tmb_class %||% "N/A"
  tmb_variants <- tmb$variant_count %||% 0
  msi_score    <- msi$msi_score %||% 0
  msi_status   <- msi$msi_status %||% "N/A"
  msi_unstable <- msi$unstable_sites %||% 0
  msi_total    <- msi$total_sites %||% 0
  hrd_score    <- hrd$hrd_score %||% 0
  hrd_status   <- hrd$hrd_status %||% "N/A"
  loh_score    <- hrd$loh_score %||% 0
  tai_score    <- hrd$tai_score %||% 0
  lst_score    <- hrd$lst_score %||% 0

  # If all scores are effectively empty, return a message plot
  if (tmb_score == 0 && msi_score == 0 && hrd_score == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No biomarker data",
               size = 5, color = "grey50") +
      theme_void()
    return(girafe(ggobj = p))
  }

  df <- tibble(
    biomarker = c("TMB", "MSI", "HRD"),
    value = c(tmb_score, msi_score, hrd_score),
    threshold = c(10, 20, 42),
    status = c(tmb_class, msi_status, hrd_status),
    tooltip = c(
      glue("<b>TMB: {tmb_score} mut/Mb</b><br>",
           "Class: {tmb_class}<br>",
           "Variants: {tmb_variants}<br>",
           "Threshold: 10 mut/Mb"),
      glue("<b>MSI: {msi_score}%</b><br>",
           "Status: {msi_status}<br>",
           "Unstable sites: {msi_unstable}/{msi_total}<br>",
           "Threshold: 20%"),
      glue("<b>HRD Score: {hrd_score}</b><br>",
           "Status: {hrd_status}<br>",
           "LOH: {loh_score} | TAI: {tai_score} | LST: {lst_score}<br>",
           "Threshold: 42")
    ),
    above = c(tmb_score >= 10, msi_score >= 20, hrd_score >= 42)
  ) |>
    mutate(
      fill_color = ifelse(above, "#c62828", "#2e7d32"),
      biomarker = factor(biomarker, levels = c("TMB", "MSI", "HRD"))
    )

  p <- ggplot(df, aes(x = biomarker, y = value)) +
    geom_col_interactive(
      aes(tooltip = tooltip, data_id = biomarker, fill = fill_color),
      width = 0.5
    ) +
    geom_hline(aes(yintercept = threshold), linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    geom_text(aes(label = paste0("threshold: ", threshold)),
              y = df$threshold, x = 0.5, hjust = 0, vjust = -0.5,
              size = 3, color = "grey40") +
    scale_fill_identity() +
    facet_wrap(~biomarker, scales = "free", nrow = 1) +
    labs(x = NULL, y = "Value") +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(size = 13, face = "bold"),
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank()
    )

  girafe(ggobj = p,
    options = list(
      opts_hover(css = "fill:orange;"),
      opts_tooltip(css = "background:white;padding:8px;border-radius:4px;border:1px solid #ccc;font-size:12px;")
    ))
}
