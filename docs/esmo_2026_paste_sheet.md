# ESMO Congress 2026 — submission-page paste sheet (with table)

Final version. Body split into the four ESMO portal sections (regulation 6.3)
plus a standalone HTML table for the reproducibility + safety numbers
(regulation 6.8 allows one table, counted into the 2000-char budget).
Total title + body + table visible text = ~1940 / 2000 chars excl. spaces.

Deadline: **12 May 2026, 21:00 CEST** (= 03:00 Wed 13 May Taipei).

---

## Title

A model-agnostic LLM agent for cross-knowledge-base evidence reconciliation in cancer NGS tertiary analysis

---

## Background

Deterministic AMP/ASCO/CAP tiering (Li 2017) takes the higher tier between OncoKB and CiVIC, conflating knowledge-base (KB) agreement with one-sided coverage or with high-tier CiVIC evidence in a different tumor type, risking false Tier I calls. Recent large language model (LLM) variant-tiering benchmarks report over-classification, motivating an architecture restricting LLM use to KB-discordant variants only.

## Methods

We built a targets-orchestrated R pipeline (DRAGEN VCF to printable report) integrating VEP, OncoKB, CiVIC, ESCAT, and AMP/ASCO/CAP tiering. Concordant variants resolve deterministically; KB-discordant variants escalate to a markdown-specified, model-agnostic agent emitting a unified tier, confidence, and rationale (JSON-schema validated). We benchmarked it against the deterministic baseline on 112 real variant-tumor cases (99 COSMIC plus 13 cross-tumor stress), measured within-model (3 runs) and cross-model reproducibility (Opus 4.7, Sonnet 4.6, Haiku 4.5), and prospectively deployed it on 62 consecutive late-stage cancer patients (8 monthly TruSight Oncology 500 batches, 11 cancer types).

## Results

Output structure and reproducibility are summarised in the Table. Disagreements concentrated at documented rule blind spots: when CiVIC carried different-tumor evidence the agent correctly down-tiered RET M918T NSCLC from Tier IA to IIC; the dominant disagreement (20 of 112) was the Tier III VUS to Tier IID boundary for oncogenic variants lacking therapeutic level. The RET catch and Tier IA concordances were preserved across all three models. All 62 deployed cases completed end-to-end with auditable per-case reconciliation logs.

## Conclusions

Confining LLM invocation to KB-discordant variants matches deterministic reproducibility on safety-critical cases while detecting tumor-specificity errors the rule cannot encode, supporting deployment as second-opinion decision support in molecular tumor boards. Code, prompts, and schemas are open source.

---

## Table (paste into the portal's table-HTML field — standalone, copy as-is)

```html
<table border="1" cellpadding="3" cellspacing="0" style="border-collapse:collapse">
 <tbody>
  <tr>
   <td><b>Metric</b></td>
   <td><b>Value</b></td>
  </tr>
  <tr>
   <td>Schema-valid JSON output</td>
   <td>112 / 112</td>
  </tr>
  <tr>
   <td>Agent-baseline tier agreement</td>
   <td>90 / 112 (80.4%)</td>
  </tr>
  <tr>
   <td>OncoKB-only Tier I/II override</td>
   <td>0 / 70</td>
  </tr>
  <tr>
   <td>Within-model reproducibility (Fleiss kappa)</td>
   <td>1.000</td>
  </tr>
  <tr>
   <td>Cross-model reproducibility (Opus / Sonnet / Haiku)</td>
   <td>0.717</td>
  </tr>
 </tbody>
</table>
```

---

## Submission portal checklist (do these in the same session)

- [ ] Category / track: **AI for diagnostics and profiling (pathology, radiology, molecular biology)** (backup: AI for clinical workflows and decision-making)
- [ ] Presentation type preference: **Mini Oral** primary, **Proffered Paper (Oral)** fallback
- [ ] Author list (≤20, no degrees, no titles, ordered as you want it published)
- [ ] Affiliations (numbered, linked to authors)
- [ ] Declaration of Interest from every author (no DOI → author dropped)
- [ ] Encore policy: not previously presented anywhere — confirm
- [ ] Generative-AI disclosure: already covered in Methods body (LLM agent described as part of methodology, regulation 3.22 satisfied)
- [ ] ESMO Merit Award checkbox (if first author < 40 + ESMO member)
- [ ] **Click "Finish Submission"** — without this, abstract is NOT submitted (regulation 6.2)
- [ ] Save the confirmation email screenshot

---

## Numbers cited in the body — re-verify with `make abstract-metrics` morning of 12 May

| Number | Source |
|---|---|
| 112 of 112 schema-valid JSON | `agent_combined_summary.csv` rows where `valid == TRUE` |
| 80.4 % agreement (90 / 112) | combined cohort+stress |
| 0 of 70 OncoKB-only Tier I/II override | `civic_amp_level` is NA AND `baseline_tier` starts with "Tier I" / "Tier II" |
| 20 of 112 Tier III VUS ↔ Tier II Level D | `baseline_tier == "Tier III VUS" & agent_tier == "Tier II Level D"` |
| Fleiss κ 1.000 within-model | `agent_consistency_summary.csv` |
| Fleiss κ 0.717 cross-model | `agent_cross_model_summary.csv` |
| RET M918T NSCLC, IA → IIC | `agent_combined_summary.csv` filter `tumor_specificity_mismatch & !agrees_with_baseline` |
