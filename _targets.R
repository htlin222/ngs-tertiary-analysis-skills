# _targets.R — Pipeline DAG for NGS Tertiary Analysis
#
# Run with: targets::tar_make()
# Visualize: targets::tar_visnetwork()

library(targets)
library(tarchetypes)

# Source all R files
tar_source("R/")
tar_source("00-qc/")
tar_source("01-variant-calling/")
tar_source("02-annotation/")
tar_source("03-cnv/")
tar_source("04-fusions/")
tar_source("05-biomarkers/")
tar_source("06-clinical-annotation/")
tar_source("07-literature/")
tar_source("08-report/")

# ── Pipeline configuration ───────────────────────────────────────────────────

tar_option_set(
  packages = c(
    "dplyr", "tidyr", "purrr", "stringr", "glue", "fs",
    "yaml", "logger", "jsonlite", "httr2",
    "VariantAnnotation", "GenomicRanges", "Rsamtools"
  ),
  error = "stop"
)

# ── Target definitions ───────────────────────────────────────────────────────

list(
  # ── Setup ──────────────────────────────────────────────────────────────────
  tar_target(config, {
    load_env()
    load_config()
  }),

  tar_target(bam_path, {
    bam <- Sys.getenv("BAM_PATH", unset = "")
    if (nchar(bam) == 0) stop("BAM_PATH environment variable not set")
    if (!file.exists(bam)) stop(glue("BAM file not found: {bam}"))
    normalizePath(bam)
  }, format = "file"),

  tar_target(sample_id, {
    sid <- Sys.getenv("SAMPLE_ID", unset = config$sample$id)
    if (is.null(sid) || nchar(sid) == 0) sid <- "SAMPLE_001"
    sid
  }),

  # ── Stage 0: Quality Control ───────────────────────────────────────────────
  tar_target(qc_results, {
    run_qc(bam_path, config, sample_id)
  }),

  # ── Stage 1: Variant Calling ───────────────────────────────────────────────
  tar_target(mutect2_vcf, {
    run_mutect2(bam_path, config, sample_id)
  }, format = "file"),

  tar_target(filtered_vcf, {
    filter_variants(mutect2_vcf, config, sample_id)
  }, format = "file"),

  # ── Stage 2: Annotation ────────────────────────────────────────────────────
  tar_target(annotated_vcf, {
    run_annotation(filtered_vcf, config, sample_id)
  }, format = "file"),

  tar_target(merged_annotations, {
    merge_annotations(annotated_vcf, config, sample_id)
  }),

  # ── Stage 3: Copy Number Variants ──────────────────────────────────────────
  tar_target(cnv_results, {
    run_cnvkit(bam_path, config, sample_id)
  }),

  tar_target(parsed_cnv, {
    parse_cnv(cnv_results, config, sample_id)
  }),

  # ── Stage 4: Fusions ──────────────────────────────────────────────────────
  tar_target(fusion_results, {
    run_fusions(bam_path, config, sample_id)
  }),

  tar_target(parsed_fusions, {
    parse_fusions(fusion_results, config, sample_id)
  }),

  # ── Stage 5: Biomarkers ───────────────────────────────────────────────────
  tar_target(tmb_result, {
    calc_tmb(filtered_vcf, config, sample_id)
  }),

  tar_target(msi_result, {
    calc_msi(bam_path, config, sample_id)
  }),

  tar_target(hrd_result, {
    calc_hrd(cnv_results, config, sample_id)
  }),

  # ── Stage 6: Clinical Annotation ──────────────────────────────────────────
  tar_target(oncokb_results, {
    query_oncokb(
      variants = merged_annotations,
      cnv = parsed_cnv,
      fusions = parsed_fusions,
      config = config,
      sample_id = sample_id
    )
  }),

  tar_target(escat_tiers, {
    classify_escat(oncokb_results, config, sample_id)
  }),

  # ── Stage 7: Literature ───────────────────────────────────────────────────
  tar_target(literature_results, {
    generate_narrative(
      variants = merged_annotations,
      cnv = parsed_cnv,
      fusions = parsed_fusions,
      config = config,
      sample_id = sample_id
    )
  }),

  # ── Stage 8: Report Generation ────────────────────────────────────────────
  tar_target(clinical_report, {
    render_report(
      sample_id = sample_id,
      config = config,
      qc = qc_results,
      variants = merged_annotations,
      cnv = parsed_cnv,
      fusions = parsed_fusions,
      tmb = tmb_result,
      msi = msi_result,
      hrd = hrd_result,
      oncokb = oncokb_results,
      escat = escat_tiers,
      literature = literature_results
    )
  }, format = "file")
)
