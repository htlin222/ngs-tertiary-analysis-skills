# ESMO Congress 2026 abstract — agentic NGS tertiary analysis

Working draft. Numbers below are pulled from `benchmarks/results/agent_combined_summary.csv`,
`agent_consistency_summary.csv`, and `agent_cross_model_summary.csv`. Re-run
`Rscript scripts/extract_abstract_metrics.R` before submission to confirm.

Source: condensed from `benchmarks/manuscript/main.qmd` abstract.

**Submission deadline: 12 May 2026, 21:00 CEST.**

---

## ESMO compliance constraints (applied to draft below)

- **Hard limit: 2,000 characters excluding spaces** for title + body + table combined (regulation 6.8). Current draft body verified under limit; re-check after any edit with `make abstract-charcount` (added below).
- **Section labels NOT included in body** (saves chars; regulation 6.3 explicitly allows omission).
- **Title must not refer to results or conclusions** (6.4) — current title is descriptive.
- **No commercial names in title** (6.5) — "TruSight Oncology 500" / "TSO500" appears in Methods only; allowed.
- **Abbreviations defined on first use** (6.7) — applied for KB, LLM, VUS, MTB, TSO500.
- **Encore not accepted** (1) — confirm this work has not been presented elsewhere before submitting.
- **Generative-AI writing disclosure** (3.23) — not required since we describe the LLM agent in Methods (research-process use, which is the correct disclosure form per the regulation).
- **Max 20 authors** (6.10).
- **Suggested category:** primary `AI for diagnostics and profiling (pathology, radiology, molecular biology)`; backup `AI for clinical workflows and decision-making`.
- **Preferred presentation type:** Rapid Oral (data-rich, novel architecture); fall back to Poster.

---

## Title

A model-agnostic LLM agent for cross-knowledge-base evidence reconciliation in cancer NGS tertiary analysis

## Body (no section headers — for paste into ESMO submission tool)

Deterministic AMP/ASCO/CAP tiering (Li 2017) takes the higher tier between OncoKB and CiVIC by fixed rule, conflating genuine knowledge-base (KB) agreement with one-sided coverage and with high-tier CiVIC evidence in a different tumor type than queried, risking inappropriate Tier I assignment. Recent large language model (LLM) benchmarks on variant tiering report over-classification of weak-evidence cases, motivating an architecture restricting LLM use to KB-discordant variants only.

We built a targets-orchestrated R pipeline (DRAGEN VCF to printable report) integrating VEP, OncoKB, CiVIC, ESCAT, and AMP/ASCO/CAP tiering. Concordant variants are resolved deterministically; KB-discordant variants escalate to a markdown-specified, model-agnostic agent emitting a unified tier, confidence, and rationale (JSON-schema validated). We benchmarked it against the deterministic baseline on 112 real variant-tumor cases (99 COSMIC plus 13 cross-tumor stress), measured within-model (3 runs) and cross-model reproducibility (Opus 4.7, Sonnet 4.6, Haiku 4.5), and deployed the pipeline on 8 pan-cancer patients profiled with TruSight Oncology 500 (TSO500).

112 of 112 agent calls produced schema-valid JSON. Agent-baseline tier agreement was 80.4% (90 of 112), concentrated at documented rule blind spots: in the 25-case OncoKB-only subset the agent never overrode baseline (0 of 25); when CiVIC carried different-tumor evidence, the agent correctly down-tiered RET M918T NSCLC from Tier IA to IIC. The dominant disagreement (20 of 112) was the Tier III variant of uncertain significance to Tier IID boundary for oncogenic variants lacking therapeutic level. Within-model reproducibility was perfect (Fleiss kappa 1.000); cross-model was substantial (0.717), with the RET catch and concordant Tier IA calls preserved in all three. All 8 clinical cases completed end-to-end with auditable per-case reconciliation logs.

Confining LLM invocation to KB-discordant variants matches deterministic reproducibility on safety-critical cases while detecting tumor-specificity errors the rule cannot encode, supporting deployment as second-opinion decision support in molecular tumor board workflows.

---

## How to verify before submission

```
Rscript scripts/extract_abstract_metrics.R   # confirm cited numbers
make abstract-charcount                       # verify <=2000 chars no-space
```

Numbers to confirm:

- 112 of 112 schema-valid → `agent_combined_summary.csv` rows where `valid == TRUE`
- 80.4% agent-baseline agreement → combined cohort+stress in `agent_combined_per_stratum.csv`
- 0 of 25 OncoKB-only override → cohort cases where `civic_amp_level` is NA
- Tier III VUS to Tier IID boundary count → `baseline_tier == "Tier III VUS" & agent_tier == "Tier II Level D"`
- Fleiss kappa within-model → `agent_consistency_summary.csv` `fleiss_kappa_tier`
- Fleiss kappa cross-model → `agent_cross_model_summary.csv`

## Pre-submission TODOs

1. Confirm work has not been presented elsewhere (encore policy, regulation 1).
2. Add author block (max 20, no degrees/titles, ordered as published).
3. Confirm at least one independent author can be presenter (regulation 3.11 — sponsor/AI-company employees cannot present clinical/translational data).
4. Complete Declaration of Interest for every author (3.17).
5. Pick category: primary `AI for diagnostics and profiling (pathology, radiology, molecular biology)`, backup `AI for clinical workflows and decision-making`.
6. Pick preferred presentation type: Rapid Oral (recommended) / Poster.
7. Decide on ESMO Merit Award application (first author under 40, ESMO member; due at submission).
8. Run `make abstract-charcount` and confirm under 2000.
9. Spell-check / proof-read in good English (6.14 — submitted abstracts cannot be edited).
