# agents/evidence_reconciler/format_input.R
# Agent-specific helper that formats one variant's evidence into the user turn
# for the evidence_reconciler agent. Kept separate from agents/runner.R so the
# runner stays agent-agnostic.

suppressPackageStartupMessages({
  library(glue)
})

#' Format a single variant's OncoKB + CiVIC evidence into a prompt string.
#'
#' Deliberately terse and structured; the model does better with bounded fields
#' than free-form paragraphs. Missing fields become "not available" so the model
#' can distinguish "no evidence" from "evidence says no."
#'
#' @param gene Gene symbol (e.g. "BRAF")
#' @param variant Protein change (e.g. "V600E")
#' @param tumor_type Tumor-type code used by OncoKB (e.g. "NSCLC")
#' @param oncokb List with fields: level, oncogenic, mutation_effect, evidence_text
#'   Any field may be NULL/NA.
#' @param civic List with fields: amp_level, evidence_count, evidence_text,
#'   tumor_context.  Any field may be NULL/NA.
#' @param vep_impact VEP HIGH/MODERATE/LOW/MODIFIER (optional)
#' @param clinvar ClinVar clinical significance (optional)
#' @return Character scalar user prompt
format_reconciler_input <- function(gene, variant, tumor_type,
                                    oncokb = list(), civic = list(),
                                    vep_impact = NA_character_,
                                    clinvar    = NA_character_) {
  na_or <- function(x, fallback = "not available") {
    if (is.null(x) || (length(x) == 1 && (is.na(x) || nchar(as.character(x)) == 0))) {
      fallback
    } else {
      as.character(x)
    }
  }

  glue::glue("
Variant under review:
  gene            : {gene}
  protein_change  : {variant}
  tumor_type      : {tumor_type}

Functional / ClinVar context:
  vep_impact              : {na_or(vep_impact)}
  clinvar_significance    : {na_or(clinvar)}

OncoKB evidence:
  therapeutic_level       : {na_or(oncokb$level)}
  oncogenic               : {na_or(oncokb$oncogenic)}
  mutation_effect         : {na_or(oncokb$mutation_effect)}
  evidence_text           : {na_or(oncokb$evidence_text)}

CiVIC evidence:
  highest_amp_level       : {na_or(civic$amp_level)}
  supporting_evidence_ct  : {na_or(civic$evidence_count)}
  tumor_context           : {na_or(civic$tumor_context)}
  evidence_text           : {na_or(civic$evidence_text)}

Task: Produce your reconciled AMP/ASCO/CAP tier call as JSON per the schema.
")
}
