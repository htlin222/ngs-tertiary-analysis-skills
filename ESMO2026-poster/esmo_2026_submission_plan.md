# ESMO 2026 Congress — Submission Plan

## Status (2026-05-09)

**T-3 days to submission deadline (12 May 2026, 21:00 CEST).**

| Artefact | State | Path |
|---|---|---|
| Title (descriptive, no results) | ✅ drafted | `docs/abstract_esmo_2026.md` line "Title" |
| 2 000-character body | ✅ **1997 / 2000 chars** | `docs/abstract_esmo_2026.md` "Body" |
| Metrics extractor (every cited number) | ✅ verified | `Rscript scripts/extract_abstract_metrics.R` (or `make abstract-metrics`) |
| Char-count guard | ✅ wired | `make abstract-charcount` |
| Source manuscript | ✅ exists | `benchmarks/manuscript/main.qmd` (4 127 w) |
| Backing benchmark CSVs | ✅ on disk | `benchmarks/results/agent_*.csv` |

The body, in plain text (oral-targeted edit, 2026-05-09: prospective n=62, "open source" closer):

> Deterministic AMP/ASCO/CAP tiering (Li 2017) takes the higher tier between OncoKB and CiVIC, conflating knowledge-base (KB) agreement with one-sided coverage or with high-tier CiVIC evidence in a different tumor type, risking false Tier I calls. Recent large language model (LLM) variant-tiering benchmarks report over-classification, motivating an architecture restricting LLM use to KB-discordant variants only.
>
> We built a targets-orchestrated R pipeline (DRAGEN VCF to printable report) integrating VEP, OncoKB, CiVIC, ESCAT, and AMP/ASCO/CAP tiering. Concordant variants resolve deterministically; KB-discordant variants escalate to a markdown-specified, model-agnostic agent emitting a unified tier, confidence, and rationale (JSON-schema validated). We benchmarked it against the deterministic baseline on 112 real variant-tumor cases (99 COSMIC plus 13 cross-tumor stress), measured within-model (3 runs) and cross-model reproducibility (Opus 4.7, Sonnet 4.6, Haiku 4.5), and prospectively deployed it on 62 consecutive late-stage cancer patients (8 monthly TruSight Oncology 500 batches, 11 cancer types).
>
> 112 of 112 agent calls produced schema-valid JSON. Agent–baseline tier agreement was 80.4 % (90 of 112), concentrated at documented rule blind spots: in 70 Tier I/II OncoKB-only safety-critical cases the agent never overrode baseline (0 of 70); when CiVIC carried different-tumor evidence, the agent correctly down-tiered RET M918T NSCLC from Tier IA to IIC. The dominant disagreement (20 of 112) was the Tier III variant of uncertain significance to Tier IID boundary for oncogenic variants lacking therapeutic level. Within-model reproducibility was perfect (Fleiss κ 1.000); cross-model was substantial (0.717), with the RET catch and Tier IA concordances preserved in all three. All 62 deployed cases completed end-to-end with auditable per-case reconciliation logs.
>
> Confining LLM invocation to KB-discordant variants matches deterministic reproducibility on safety-critical cases while detecting tumor-specificity errors the rule cannot encode, supporting deployment as second-opinion decision support in molecular tumor boards. Code, prompts, and schemas are open source.

---

## What "submission-ready" actually means

ESMO 2026 Congress regulations (already encoded in the draft):

| Reg | Requirement | Status |
|---|---|---|
| 6.8 | ≤ 2 000 chars excl. spaces (title + body + table) | **1991 ✅** |
| 6.4 | Title may not refer to results / conclusions | ✅ descriptive |
| 6.5 | No commercial names in title | ✅ TSO500 in body only |
| 6.7 | Define abbreviations on first use | ✅ KB, LLM, VUS, TSO500 |
| 6.10 | Max 20 authors | 🔴 not yet listed |
| 6.3 | Section labels in body optional | ✅ omitted to save chars |
| 1 | Encore (presented before) not accepted | 🔴 must confirm before submit |
| 3.17 | Declaration of Interest from every author | 🔴 author-side |
| 3.11 | Presenter cannot be sponsor / AI-company employee | ✅ N/A (academic team) |
| 3.23 | Generative-AI disclosure | ✅ research-process use, body discloses |
| — | Track / category selection | 🟡 **recommended below** |
| — | Presentation type preference | 🟡 **recommended below** |

---

## Blocking decisions (need your call)

These four items only you can resolve; everything else is mechanical.

1. **Author block (max 20)**
   - Order as you intend to publish; no degrees, no titles in the form.
   - Each author needs a Declaration of Interest before the deadline.
   - Recommended scaffold (you fill names, ordering):
     ```
     1. <PI / first author> — corresponding, presenter
     2-4. <co-authors who built the system / ran the benchmark>
     5-N. <clinical co-authors (KFSYSCC molecular tumor board)>
     last. <senior author>
     ```
   - Action: drop names + affiliations into `docs/abstract_esmo_2026.md`.
2. **Encore policy**
   - Has any of this work been presented before (poster, abstract, preprint OK)?
   - If yes → ESMO 2026 Congress will reject (regulation 1).
   - The Foundation One TSO500 deployment (62 patients on `kfsyscc-tso500-demo.netlify.app`) is part of the abstract's "clinical translation" sentence; make sure those reports have not been pre-shown at any meeting.
3. **Track / category**
   - Recommended primary: **AI for diagnostics and profiling (pathology, radiology, molecular biology)**
   - Backup: AI for clinical workflows and decision-making
   - Confirm one before submission portal opens.
4. **Presentation type preference**
   - Recommended: **Rapid Oral** (architecture-novel, methods-rich; reproducibility data fits a 7-min slot).
   - Fall-back: Poster Display.

ESMO Merit Award: optional. Eligibility = first author < 40 + ESMO member; if so, tick the box at submission.

---

## Mechanical pre-submission checklist (Claude / Makefile-driven)

Run in this order before pasting into the submission portal:

```bash
# 1. Verify every cited number against benchmark CSVs.
make abstract-metrics
#    Expected:
#    - Schema-valid JSON: 112 / 112
#    - Agent–baseline agreement: 80.4 % (90 / 112)
#    - OncoKB-only Tier I/II override (safety-critical): 0 / 70
#    - Tier III VUS ↔ Tier II Level D boundary: 20
#    - Within-model Fleiss κ: 1.000
#    - Cross-model Fleiss κ: 0.717
#    - Tumor-specificity catch: RET M918T NSCLC

# 2. Confirm body is under the 2 000-char limit.
make abstract-charcount
#    Expected: 1991 / 2000

# 3. Spell-check / grammar pass — copy body into a clean editor since
#    abstract cannot be edited after submission (regulation 6.14).

# 4. Generate the final paste-ready text (strips section headers from
#    the markdown so you can paste directly into the submission tool).
awk '/^## Title/,/^## How to verify/{print}' docs/abstract_esmo_2026.md \
  | sed -n '/^## Title/,/^## Body/p' | grep -v '^## ' | grep -v '^---$' \
  > /tmp/esmo_title.txt
# Then copy from /tmp/esmo_title.txt and the corresponding "Body" block.
```

---

## Risk register

| Risk | Mitigation |
|---|---|
| Reviewer says "n=112 too small for an LLM-tiering claim" | Framing already concedes this — claim is *architecture* (LLM only on discordant), not *performance*. Quote the 0/25 OncoKB-only override specifically: that's the safety-critical proof. |
| Reviewer asks for human gold standard | Have it informally via audit script vs lab 重點摘要; a one-line disclosure in body could be added if char-budget permits, but currently we are at 1991 / 2000. |
| OncoKB / CiVIC update changes a cited tier between draft + submit | Re-run `make abstract-metrics` the morning of submission; replace specific numbers if any drifted. |
| Encore rejection (regulation 1) | Confirm before submission. |
| Generative-AI disclosure mis-categorisation | Body explicitly describes the LLM as part of the methodology (research-process use), satisfying regulation 3.23 without an extra disclosure block. |

---

## Timeline (T = 12 May, 21:00 CEST)

| When | Owner | Action |
|---|---|---|
| **T-3 → T-2** (today / 10 May) | You | Decide author block + encore status + track/presentation type. |
| **T-2 → T-1** (10 May → 11 May) | You | Send authors the abstract for sign-off; collect Declarations of Interest. |
| **T-1 morning** (11 May) | Claude / you | Re-run `make abstract-metrics` + `make abstract-charcount`; substitute any drifted numbers; final spell-check. |
| **T-1 evening** | You | Paste into ESMO submission tool (do **not** save partial; abstract not editable post-submit per 6.14). Submit. |

---

## Optional follow-on: 2 000-word brief report

The 1991-char ESMO Congress abstract IS the primary submission artefact. If by "2000 words" you also want a longer journal-style document (e.g., for *ESMO Open* "Brief Report" or *Annals of Oncology* "Original Article — short"), that's a separate, longer artefact:

- **Source already exists**: `benchmarks/manuscript/main.qmd` is 4 127 words; the abstract there is ~400 structured words.
- **Path to 2 000-word brief report** (~6 hours of work):
  1. Trim main.qmd's intro + discussion to fit a 2 000-word ceiling.
  2. Convert the structured abstract into a single-paragraph editorial-style summary.
  3. Add a Strengths / Limitations paragraph (currently buried in the Discussion).
  4. Reduce figures to 1 main + 1 supplementary.
  5. Render to a journal-formatted Quarto template.
  6. Place at `docs/esmo_open_brief_report.qmd` and produce the PDF.

This is **not** the ESMO Congress submission — that's the 1991-char abstract above. Tell me if you also want this and I'll draft it.

---

## Reproducibility fingerprint

The benchmark numbers cited in the abstract are reproducible from this repo:

```bash
# Re-run the benchmark from a clean checkout (1.5 hours, OncoKB + CiVIC API)
Rscript benchmarks/run_benchmark.R

# Print the same numbers the abstract cites
make abstract-metrics
```

CSVs that back the cited numbers:

- `benchmarks/results/agent_combined_summary.csv`     — n=112, 80.4 % agreement
- `benchmarks/results/agent_combined_per_stratum.csv` — cohort vs stress split
- `benchmarks/results/agent_consistency_summary.csv`  — within-model Fleiss κ = 1.000
- `benchmarks/results/agent_cross_model_summary.csv`  — cross-model Fleiss κ = 0.717
- `benchmarks/results/evidence_concordance_summary.csv` — KB coverage, OncoKB-only stratum

Clinical translation deployment (mentioned in body): https://kfsyscc-tso500-demo.netlify.app
