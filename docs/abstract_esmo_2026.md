# ESMO 2026 abstract — agentic NGS tertiary analysis

Working draft. Numbers below are pulled from `benchmarks/results/agent_combined_summary.csv`,
`agent_consistency_summary.csv`, and `agent_cross_model_summary.csv`. Re-run
`Rscript scripts/extract_abstract_metrics.R` before submission to confirm.

Source: condensed from `benchmarks/manuscript/main.qmd` abstract.

---

## Title

**A model-agnostic LLM agent for evidence reconciliation in cancer NGS tertiary analysis: cross-model benchmark and pan-cancer clinical deployment**

## Background

Deterministic AMP/ASCO/CAP classification (Li et al., 2017) takes the higher
tier between OncoKB and CiVIC by fixed rule, but cannot distinguish (i) genuine
inter-knowledge-base agreement from (ii) one-sided coverage gaps or (iii)
high-tier CiVIC evidence carried in a *different* tumor type than queried.
All three collapse to the same tier under the rule, but only the first
warrants high-confidence reporting; the third risks an inappropriate Tier I
call. Recent benchmarks of large language models (LLMs) on variant tiering
report a tendency to over-classify weak-evidence variants
(npj Precis Oncol 2025), motivating an architecture in which the LLM is
restricted to the discordance cases the rule cannot resolve.

## Methods

We implemented a `targets`-orchestrated R pipeline (DRAGEN VCF →
comprehensive genomic profiling–style printable report) that integrates
VEP, OncoKB, CiVIC, ESCAT v1.1, and AMP/ASCO/CAP tiering. Concordant variants
are resolved by a deterministic rule engine; cross-knowledge-base–discordant
variants are escalated to an evidence-reconciliation agent specified via the
**SKILL protocol** — a declarative, markdown-based agent specification
**decoupled from any specific LLM backend**, enabling portability across
commercial and local models. The agent emits a unified tier, confidence,
discordance category, and human-readable rationale, validated against a JSON
schema. We benchmarked the agent against the deterministic baseline on
**112 real variant–tumor cases** (99 well-characterized COSMIC variants
across 6 tumor types + 13 cross-tumor stress cases), measured within-model
reproducibility (3 independent runs, n=10) and cross-model reproducibility
(Opus 4.7, Sonnet 4.6, Haiku 4.5), and deployed the pipeline on
**8 pan-cancer patients** profiled with TruSight Oncology 500
(4 endometrial, 2 breast, 2 colorectal).

## Results

Output structure was reliable: **112/112 (100%)** of agent calls produced
schema-valid JSON. Agent–baseline tier agreement was **80.4% (90/112)**, with
disagreements structurally concentrated at documented rule blind spots: in
the **safety-critical 25-case subset where OncoKB was the sole evidence
source, the agent never overrode the baseline (0/25)**. When CiVIC carried
evidence in a different tumor type, the agent correctly down-tiered
RET M918T NSCLC (Tier I A → Tier II Level C, flagged
`tumor_specificity_mismatch`, citing CiVIC's medullary-thyroid context). The
dominant disagreement (20/112) was a documented Tier III VUS ↔ Tier II
Level D boundary for oncogenic variants lacking explicit therapeutic levels.
Within-model reproducibility was perfect (**Fleiss κ = 1.000**, 3 runs);
cross-model reproducibility was substantial (**Fleiss κ = 0.717**) with
the RET tumor-specificity catch and all concordant Tier I A calls preserved
across all three models. All eight clinical cases completed end-to-end with
auditable per-case agent reconciliation logs (TMB 3.9–10.5 mut/Mb,
MSI 0.0–3.3%).

## Conclusions

A SKILL-protocol, model-agnostic agentic architecture that confines LLM
invocation to knowledge-base–discordant variants delivers reproducibility
matching deterministic systems on safety-critical cases while detecting
tumor-specificity errors the rule cannot encode. Source code, prompts,
schemas, and benchmark scripts are openly released. Prospective validation
against expert molecular tumor board consensus is ongoing.

---

## How to fill / verify numbers

```
Rscript scripts/extract_abstract_metrics.R
```

Numbers to confirm before submission:

- 112/112 schema-valid → `agent_combined_summary.csv` rows where `valid == TRUE`
- 80.4% agent–baseline agreement → `pct_agree` in `agent_consistency_summary.csv` strata aggregate, currently 78.8% (cohort) + 92.3% (stress) = 80.4% combined
- 0/25 OncoKB-only override → cohort cases where `civic_amp_level` is NA
- Tier III VUS ↔ II D boundary count → `baseline_tier == "Tier III VUS" & agent_tier == "Tier II Level D"`
- Fleiss κ within-model → `agent_consistency_summary.csv` `fleiss_kappa_tier`
- Fleiss κ cross-model → `agent_cross_model_summary.csv`
- 8-patient TSO500 ranges → `reports/batch7/` index

## What still needs doing before submit

1. Run `extract_abstract_metrics.R` and confirm any rounding deltas.
2. Pick ESMO track (RWD & Digital Oncology recommended; main Congress riskier).
3. Add author block + affiliation.
4. Trim to track word limit (ESMO Congress is 350 words for body; this draft is ~340).
5. Decide whether to keep "SKILL protocol" branding or rephrase as
   "declarative markdown-based agent specification" if reviewers may not
   recognize the term — currently kept with an inline definition.
