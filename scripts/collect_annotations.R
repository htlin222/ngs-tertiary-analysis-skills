#!/usr/bin/env Rscript
# scripts/collect_annotations.R — Stage every produced merged_annotations.tsv
# (across all batches) into annotation/batch{N}/{sample}_{tumor}_annotation.tsv
# and write a single annotation/_manifest.tsv with sha256 + n_variants.
#
# Idempotent: re-runs overwrite the staged copies and refresh the manifest.
# annotation/ is gitignored — patient data, never commit.

.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library",
            .libPaths()))
setwd("/Users/htlin/ngs-tertiary-analysis-skills")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(glue)
  library(here); library(fs); library(purrr); library(tibble)
  library(digest)
})

manifest <- read_tsv(here("inputs/TSO500-HRD/_manifest.tsv"),
                     show_col_types = FALSE)

# Map "Batch_1_20250908" → "batch1"
manifest <- manifest |>
  mutate(
    batch_num = as.integer(str_extract(batch, "(?<=Batch_)\\d+")),
    batch_short = sprintf("batch%d", batch_num),
    src_path = here("reports", sample_id, "02-annotation",
                    "merged_annotations.tsv"),
    dest_rel = file.path("annotation", batch_short,
                         glue("{sample_id}_{oncokb_code}_annotation.tsv")),
    dest_path = here(dest_rel)
  )

# Stage destination batch dirs
dest_root <- here("annotation")
dir_create(dest_root, recurse = TRUE)
walk(unique(manifest$batch_short),
     ~ dir_create(file.path(dest_root, .x), recurse = TRUE))

cat(sprintf("[collect] %d samples in manifest, %d batches\n",
            nrow(manifest), length(unique(manifest$batch_short))))

# Copy + measure each file
results <- pmap_dfr(manifest, function(sample_id, batch, tumor_zh, oncokb_code,
                                      tsv_path, batch_num, batch_short,
                                      src_path, dest_rel, dest_path) {
  if (!file_exists(src_path)) {
    warning(sprintf("MISSING: %s", src_path), call. = FALSE)
    return(tibble(sample_id = sample_id, batch = batch_short,
                  oncokb_code = oncokb_code, tumor_zh = tumor_zh,
                  src_path = as.character(src_path),
                  dest_path = dest_rel,
                  n_variants = NA_integer_,
                  bytes = NA_integer_,
                  sha256 = NA_character_,
                  status = "missing"))
  }
  file_copy(src_path, dest_path, overwrite = TRUE)

  # n_variants = lines - 1 (header)
  n_lines <- length(readLines(dest_path, warn = FALSE))
  n_var <- max(0L, n_lines - 1L)

  tibble(
    sample_id   = sample_id,
    batch       = batch_short,
    oncokb_code = oncokb_code,
    tumor_zh    = tumor_zh,
    src_path    = as.character(path_rel(src_path, here())),
    dest_path   = dest_rel,
    n_variants  = n_var,
    bytes       = as.integer(file_info(dest_path)$size),
    sha256      = digest(file = dest_path, algo = "sha256"),
    status      = "ok"
  )
})

# Write manifest at annotation/_manifest.tsv
manifest_out <- here("annotation", "_manifest.tsv")
results <- results |>
  mutate(collected_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")) |>
  arrange(batch, sample_id)
write_tsv(results, manifest_out)

# Per-batch summary
summary_tbl <- results |>
  group_by(batch) |>
  summarise(n_samples  = n(),
            n_ok       = sum(status == "ok"),
            n_missing  = sum(status == "missing"),
            n_variants = sum(n_variants, na.rm = TRUE),
            .groups = "drop")
cat("\n[collect] per-batch summary:\n")
print(summary_tbl)

cat(sprintf("\n[collect] manifest: %s\n", manifest_out))
cat(sprintf("[collect] %d ok, %d missing\n",
            sum(results$status == "ok"),
            sum(results$status == "missing")))
