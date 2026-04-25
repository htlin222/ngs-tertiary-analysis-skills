#!/usr/bin/env Rscript
# benchmarks/manuscript/figures/make_figures.R
# Generates Figure 1 (stratified safety) and Figure 2 (tier-transition matrix)
# from agent_combined_summary.csv. Outputs PDF + PNG.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

results_dir <- here::here("benchmarks/results")
fig_dir     <- here::here("benchmarks/manuscript/figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(file.path(results_dir, "agent_combined_summary.csv"),
              stringsAsFactors = FALSE) |>
  mutate(
    valid = as.logical(valid),
    agrees_with_baseline = as.logical(agrees_with_baseline),
    civic_evidence_in_different_tumor = as.logical(civic_evidence_in_different_tumor)
  ) |>
  filter(valid)

# ── Figure 1: Stratified safety analysis ─────────────────────────────────────

# Define strata for the safety analysis on Tier I A baseline cases
fig1_data <- d |>
  filter(stratum == "cohort") |>
  mutate(
    safety_stratum = case_when(
      is.na(civic_amp_level)                    ~ "OncoKB sole source",
      civic_evidence_in_different_tumor == TRUE ~ "Cross-tumor CiVIC",
      civic_evidence_in_different_tumor == FALSE ~ "Same-tumor matched",
      TRUE ~ "Other"
    ),
    safety_stratum = factor(
      safety_stratum,
      levels = c("OncoKB sole source", "Same-tumor matched", "Cross-tumor CiVIC"))
  ) |>
  group_by(safety_stratum) |>
  summarise(
    n_total      = n(),
    n_overridden = sum(!agrees_with_baseline),
    n_agree      = sum(agrees_with_baseline),
    pct_override = 100 * n_overridden / n_total,
    .groups = "drop"
  )

# Add the stress cohort as a fourth stratum
stress_safety <- d |>
  filter(stratum == "stress") |>
  summarise(
    safety_stratum = factor("Stress (out-of-cohort)",
                            levels = c("OncoKB sole source",
                                       "Same-tumor matched",
                                       "Cross-tumor CiVIC",
                                       "Stress (out-of-cohort)")),
    n_total = n(),
    n_overridden = sum(!agrees_with_baseline),
    n_agree = sum(agrees_with_baseline),
    pct_override = 100 * n_overridden / n_total
  )

fig1_data <- bind_rows(
  fig1_data |> mutate(safety_stratum = factor(
    as.character(safety_stratum),
    levels = c("OncoKB sole source", "Same-tumor matched", "Cross-tumor CiVIC",
               "Stress (out-of-cohort)"))),
  stress_safety
)

fig1 <- ggplot(fig1_data, aes(x = safety_stratum, y = pct_override,
                              fill = safety_stratum)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d / %d", n_overridden, n_total)),
            vjust = -0.5, size = 4) +
  scale_y_continuous(limits = c(0, 50),
                     labels = label_percent(scale = 1)) +
  scale_fill_manual(values = c(
    "OncoKB sole source"     = "#2c7a51",  # green = safe
    "Same-tumor matched"     = "#5b8aaf",
    "Cross-tumor CiVIC"      = "#d97706",  # amber = the value-add zone
    "Stress (out-of-cohort)" = "#7e22ce"
  )) +
  labs(
    title = NULL,
    x = "Evidence stratum",
    y = "% agent overrode baseline (with N counts)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(fig_dir, "fig1_stratified_safety.pdf"), fig1,
       width = 6.5, height = 4.0, device = cairo_pdf)
ggsave(file.path(fig_dir, "fig1_stratified_safety.png"), fig1,
       width = 6.5, height = 4.0, dpi = 300)
cat("Wrote fig1_stratified_safety.{pdf,png}\n")

# ── Figure 2: Tier transition heatmap ────────────────────────────────────────

tier_levels <- c("Tier I Level A", "Tier I Level B",
                 "Tier II Level C", "Tier II Level D",
                 "Tier III VUS", "Tier IV Benign")

fig2_data <- d |>
  mutate(
    baseline_tier = factor(baseline_tier, levels = tier_levels),
    agent_tier    = factor(agent_tier, levels = tier_levels)
  ) |>
  count(baseline_tier, agent_tier, .drop = FALSE) |>
  complete(baseline_tier, agent_tier, fill = list(n = 0L))

fig2 <- ggplot(fig2_data,
               aes(x = agent_tier, y = baseline_tier, fill = n)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            color = "black", size = 4.5) +
  scale_fill_gradient(low = "#f7f4ef", high = "#1f4e5f",
                      guide = guide_colorbar(title = "Count")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE, limits = rev(tier_levels)) +
  labs(
    title = NULL,
    x = "Agent unified AMP tier",
    y = "Deterministic baseline AMP tier"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(fig_dir, "fig2_tier_transitions.pdf"), fig2,
       width = 6.5, height = 5.0, device = cairo_pdf)
ggsave(file.path(fig_dir, "fig2_tier_transitions.png"), fig2,
       width = 6.5, height = 5.0, dpi = 300)
cat("Wrote fig2_tier_transitions.{pdf,png}\n")
