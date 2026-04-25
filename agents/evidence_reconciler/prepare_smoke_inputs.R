# agents/evidence_reconciler/prepare_smoke_inputs.R
# Phase 1 of the evidence_reconciler smoke test.
#
# Writes 5 synthetic discordant cases as input JSONs under
# agents/inputs/evidence_reconciler/. The Claude Code orchestrator then invokes
# the evidence_reconciler subagent once per case and writes responses under
# agents/outputs/evidence_reconciler/. Phase 3 (analyze_smoke_outputs.R) reads
# them back and compares to the deterministic baseline.
#
# Usage (from repo root):
#   Rscript agents/evidence_reconciler/prepare_smoke_inputs.R

suppressPackageStartupMessages({
  library(here)
  library(tibble)
  library(purrr)
  library(logger)
})

source(here::here("R/utils.R"))
source(here::here("R/amp_classification.R"))           # deterministic baseline
source(here::here("agents/io.R"))
source(here::here("agents/evidence_reconciler/format_input.R"))

setup_logging("INFO")

# ── 5 synthetic discordant cases (one per discordance type + 1 concordant) ───

cases <- tribble(
  ~case_id,               ~gene,   ~variant, ~tumor_type,
  "01_level_mismatch",    "EGFR",  "L858R",  "NSCLC",
  "02_tumor_specificity", "BRAF",  "V600E",  "COAD",
  "03_oncogenicity",      "TP53",  "P72R",   "BRCA",
  "04_coverage_only",     "POLE",  "P286R",  "COAD",
  "05_concordant",        "KRAS",  "G12C",   "NSCLC"
)

evidence <- list(
  "01_level_mismatch" = list(
    oncokb = list(level = "LEVEL_1", oncogenic = "Oncogenic",
                  mutation_effect = "Gain-of-function",
                  evidence_text = "EGFR L858R in NSCLC: osimertinib FDA-approved 1L (FLAURA trial)."),
    civic  = list(amp_level = "TIER_II_LEVEL_C", evidence_count = 8L,
                  tumor_context = "NSCLC",
                  evidence_text = "8 supporting assertions. Most cite single-arm phase II data. No TIER_I_LEVEL_A assertion found in CiVIC."),
    vep_impact = "MODERATE", clinvar = "pathogenic"),

  "02_tumor_specificity" = list(
    oncokb = list(level = "LEVEL_3B", oncogenic = "Oncogenic",
                  mutation_effect = "Gain-of-function",
                  evidence_text = "BRAF V600E in COAD: encorafenib + cetuximab FDA-approved (BEACON). Note: LEVEL_1 context is melanoma."),
    civic  = list(amp_level = "TIER_I_LEVEL_A", evidence_count = 192L,
                  tumor_context = "colorectal",
                  evidence_text = "192 assertions across colorectal; TIER_I_LEVEL_A from BEACON CRC regulatory submission."),
    vep_impact = "MODERATE", clinvar = "pathogenic"),

  "03_oncogenicity" = list(
    oncokb = list(level = NA, oncogenic = "Neutral",
                  mutation_effect = "Neutral",
                  evidence_text = "TP53 P72R: common polymorphism (rs1042522); OncoKB calls Neutral."),
    civic  = list(amp_level = "TIER_I_LEVEL_B", evidence_count = 14L,
                  tumor_context = "breast",
                  evidence_text = "14 assertions; disputed predictive association with chemotherapy response."),
    vep_impact = "LOW", clinvar = "benign"),

  "04_coverage_only" = list(
    oncokb = list(level = "LEVEL_2", oncogenic = "Oncogenic",
                  mutation_effect = "Loss-of-function",
                  evidence_text = "POLE P286R in COAD: MSI-H / TMB-high surrogate. Pembrolizumab LEVEL_2."),
    civic  = list(amp_level = NA, evidence_count = 1L,
                  tumor_context = NA,
                  evidence_text = "One supporting assertion but no AMP-level assignment in CiVIC."),
    vep_impact = "MODERATE", clinvar = "pathogenic"),

  "05_concordant" = list(
    oncokb = list(level = "LEVEL_1", oncogenic = "Oncogenic",
                  mutation_effect = "Gain-of-function",
                  evidence_text = "KRAS G12C in NSCLC: sotorasib FDA-approved 2L+ (CodeBreaK 100)."),
    civic  = list(amp_level = "TIER_I_LEVEL_A", evidence_count = 57L,
                  tumor_context = "NSCLC",
                  evidence_text = "57 assertions; TIER_I_LEVEL_A predictive for sotorasib and adagrasib."),
    vep_impact = "MODERATE", clinvar = "pathogenic")
)

# ── Compute baseline + write inputs ──────────────────────────────────────────

walk(transpose(cases), function(row) {
  ev <- evidence[[row$case_id]]

  baseline <- classify_amp_tier(
    oncokb_oncogenic = ev$oncokb$oncogenic,
    oncokb_level     = ev$oncokb$level,
    civic_amp_level  = ev$civic$amp_level,
    vep_impact       = ev$vep_impact,
    clinvar          = ev$clinvar
  )

  user_turn <- format_reconciler_input(
    gene       = row$gene,
    variant    = row$variant,
    tumor_type = row$tumor_type,
    oncokb     = ev$oncokb,
    civic      = ev$civic,
    vep_impact = ev$vep_impact,
    clinvar    = ev$clinvar
  )

  metadata <- list(
    gene            = row$gene,
    variant         = row$variant,
    tumor_type      = row$tumor_type,
    baseline_tier   = baseline$amp_tier,
    baseline_level  = baseline$amp_level,
    baseline_evidence = baseline$amp_evidence
  )

  write_agent_input(
    agent_name = "evidence_reconciler",
    case_id    = row$case_id,
    user_turn  = user_turn,
    metadata   = metadata
  )
})

log_info("All 5 smoke-test inputs written under agents/inputs/evidence_reconciler/")
log_info("Next: operator invokes the evidence_reconciler subagent once per case.")
