# _targets.R — Pipeline DAG for NGS Tertiary Analysis
#
# Supports two entry points:
#   BAM input: runs full pipeline (stages 0-8)
#   VCF input: skips QC + variant calling, enters at annotation (stages 2-8)
#
# Usage:
#   BAM_PATH=tumor.bam  SAMPLE_ID=P001  targets::tar_make()
#   VCF_PATH=somatic.vcf.gz SAMPLE_ID=P001  targets::tar_make()
#   INPUT_PATH=file.bam  SAMPLE_ID=P001  targets::tar_make()  # auto-detect
#
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
    "yaml", "logger", "jsonlite", "httr2", "here",
    "VariantAnnotation", "GenomicRanges", "Rsamtools"
  ),
  error = "null"
)

# ── Target definitions ───────────────────────────────────────────────────────

list(
  # ── Setup & Input Detection ────────────────────────────────────────────────
  tar_target(config, {
    load_env()
    load_config()
  }),

  tar_target(sample_id, {
    sid <- Sys.getenv("SAMPLE_ID", unset = config$sample$id)
    if (is.null(sid) || nchar(sid) == 0) sid <- "SAMPLE_001"
    sid
  }),

  # Auto-detect input: BAM or VCF
  tar_target(pipeline_input, {
    resolve_input()
  }),

  tar_target(input_type, pipeline_input$type),
  tar_target(input_path, pipeline_input$path, format = "file"),

  # ── BAM-only stages (skipped when input is VCF) ───────────────────────────

  # Stage 0: QC — BAM only
  tar_target(qc_results, {
    if (input_type != "bam") {
      log_info("Skipping QC (input is VCF, not BAM)")
      return(list(
        summary = list(total_reads = NA, mapped_reads = NA, mapping_rate = NA,
                       mean_coverage = NA, on_target_rate = NA,
                       duplicate_rate = NA, tumor_purity = NA),
        per_gene_coverage = tibble(gene = character(), mean_coverage = numeric(),
                                   min_coverage = numeric(), pct_above_200x = numeric()),
        pass = NA,
        skipped = TRUE,
        skip_reason = "VCF input — BAM-level QC not available"
      ))
    }
    run_qc(input_path, config, sample_id)
  }),

  # Stage 1: Variant Calling — BAM only
  tar_target(mutect2_vcf, {
    if (input_type != "bam") {
      log_info("Skipping Mutect2 (input is VCF)")
      return(NULL)
    }
    run_mutect2(input_path, config, sample_id)
  }, format = "file"),

  tar_target(filtered_vcf, {
    if (input_type == "vcf") {
      # VCF input goes directly to annotation — use the input VCF as "filtered"
      log_info("Using input VCF as filtered variants: {input_path}")
      return(input_path)
    }
    filter_variants(mutect2_vcf, config, sample_id)
  }, format = "file"),

  # ── Annotation (runs for both BAM and VCF) ────────────────────────────────

  # Stage 2: Annotation
  tar_target(annotated_vcf, {
    run_annotation(filtered_vcf, config, sample_id)
  }, format = "file"),

  tar_target(merged_annotations, {
    merge_annotations(annotated_vcf, config, sample_id)
  }),

  # Stage 3: CNV — BAM only (CNV from VCF not supported)
  tar_target(cnv_results, {
    if (input_type != "bam") {
      log_info("Skipping CNVkit (requires BAM input)")
      return(NULL)
    }
    run_cnvkit(input_path, config, sample_id)
  }),

  tar_target(parsed_cnv, {
    if (is.null(cnv_results)) {
      log_info("No CNV results (VCF input) — returning empty table")
      return(tibble(chromosome = character(), start = integer(), end = integer(),
                    gene = character(), log2_ratio = numeric(), copy_number = numeric(),
                    length = integer(), type = character()))
    }
    parse_cnv(cnv_results, config, sample_id)
  }),

  # Stage 4: Fusions — BAM only
  tar_target(fusion_results, {
    if (input_type != "bam") {
      log_info("Skipping fusion detection (requires BAM input)")
      return(NULL)
    }
    run_fusions(input_path, config, sample_id)
  }),

  tar_target(parsed_fusions, {
    if (is.null(fusion_results)) {
      log_info("No fusion results (VCF input) — returning empty table")
      return(tibble(gene_a = character(), gene_b = character(),
                    fusion_name = character(), exon_a = integer(), exon_b = integer(),
                    chr_a = character(), pos_a = integer(),
                    chr_b = character(), pos_b = integer(),
                    supporting_reads = integer(), fusion_type = character(),
                    known_fusion = logical()))
    }
    parse_fusions(fusion_results, config, sample_id)
  }),

  # ── Biomarkers ─────────────────────────────────────────────────────────────

  # Stage 5: TMB works from VCF, MSI needs BAM, HRD needs CNV
  tar_target(tmb_result, {
    calc_tmb(filtered_vcf, config, sample_id)
  }),

  tar_target(msi_result, {
    if (input_type != "bam") {
      log_info("Skipping MSI (requires BAM input for microsatellite analysis)")
      return(list(msi_score = NA, msi_status = "Unknown (VCF input)",
                  unstable_sites = NA, total_sites = NA))
    }
    calc_msi(input_path, config, sample_id)
  }),

  tar_target(hrd_result, {
    if (is.null(cnv_results)) {
      log_info("Skipping HRD (requires CNV data from BAM)")
      return(list(hrd_score = NA, hrd_status = "Unknown (no CNV data)",
                  loh_score = NA, tai_score = NA, lst_score = NA, is_reliable = FALSE))
    }
    calc_hrd(cnv_results, config, sample_id)
  }),

  # ── Clinical Annotation (runs for both) ───────────────────────────────────

  # Stage 6: OncoKB + ESCAT
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

  tar_target(civic_results, {
    query_civic(
      variants = merged_annotations,
      cnv = parsed_cnv,
      fusions = parsed_fusions,
      config = config,
      sample_id = sample_id
    )
  }),

  # AMP/ASCO/CAP Oncogenicity Classification
  tar_target(amp_results, {
    source(here::here("R/amp_classification.R"))
    civic_assertions <- if (!is.null(civic_results$assertions) &&
                            nrow(civic_results$assertions) > 0) {
      civic_results$assertions
    } else {
      NULL
    }
    # Build variants with OncoKB data for classification
    variants_enriched <- merged_annotations |>
      left_join(
        bind_rows(
          map_dfr(oncokb_results$mutations %||% list(), function(m) {
            tibble(gene = m$gene, alteration = m$alteration,
                   oncogenic = m$oncogenic,
                   sensitive_level = m$highest_sensitive_level)
          })
        ),
        by = c("gene", "hgvsp" = "alteration")
      )
    classify_all_amp(variants_enriched, civic_assertions)
  }),

  # Stage 7: Literature
  tar_target(literature_results, {
    generate_narrative(
      variants = merged_annotations,
      cnv = parsed_cnv,
      fusions = parsed_fusions,
      config = config,
      sample_id = sample_id
    )
  }),

  # ── Report Generation ─────────────────────────────────────────────────────

  # Stage 8: HTML Report
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
      literature = literature_results,
      civic = civic_results,
      amp = amp_results
    )
  }, format = "file")
)
