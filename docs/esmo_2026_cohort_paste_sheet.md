# ESMO Congress 2026 — submission #2 (cohort poster)

Second abstract for the same submission session. Focus: pure single-center
real-world cohort, no AI claims. Targets **Poster Display** primary,
**Mini Oral** fallback.

Same deadline: **12 May 2026, 21:00 CEST** (= 03:00 Wed 13 May Taipei).

---

## Title

Real-world clinical actionability of comprehensive genomic profiling in advanced solid tumors: a prospective single-center pan-cancer cohort

---

## Background

Real-world actionability of comprehensive genomic profiling (CGP) in advanced solid tumors ranges 20-40% across published cohorts, mostly large US/European consortia; Asian single-center first-pass implementation data are scarce. We characterised the biomarker landscape and therapy-match rate at a Taiwanese tertiary cancer center.

## Methods

We prospectively enrolled 62 consecutive late-stage solid-tumor patients profiled with a 523-gene panel (TruSight Oncology 500, DRAGEN v2.6) over 8 monthly batches. Tumor mutational burden, microsatellite instability, genomic instability score (HRD proxy), small variants, copy-number alterations and large rearrangements were annotated by OncoKB and CiVIC and tiered by AMP/ASCO/CAP (Li 2017). Therapy-match rates were stratified into OncoKB on-label (Level 1/2/3A) versus off-label (3B/4), with resistance flags (R1/R2) tracked.

## Results

Eleven cancer types were represented (endometrial 16, breast 14, colorectal 9, NSCLC 8, ovarian 7, others 8). Biomarker and therapy-match rates appear in the Table. Most recurrent oncogenic driver genes were TP53 (26 patients), KRAS (16), PIK3CA (16), PTEN (11), ARID1A (9), APC (7), ESR1 (6). HRD-positive cases clustered in breast and ovarian; a marquee finding was a BRCA1 large rearrangement (exon 7 LOSS, fold change 0.15) in an HRD-positive breast cancer, supporting on-label PARP inhibition. OncoKB resistance variants (KRAS G12D, NRAS G13D/Q61L) flagged anti-EGFR ineligibility in 6 colorectal cases.

## Conclusions

Routine pan-cancer CGP at a Taiwanese tertiary cancer center yielded an on-label therapy match in roughly 1 in 3 advanced solid tumors, consistent with international cohorts. Biomarker-defined subsets (TMB-High, HRD-positive, MSI-High) identified additional immunotherapy- or PARPi-eligible patients. Patient reports openly browseable.

---

## Table (standalone HTML — paste into portal table-HTML field)

```html
<table border="1" cellpadding="3" cellspacing="0" style="border-collapse:collapse">
 <tbody>
  <tr>
   <td><b>Cohort metric</b></td>
   <td><b>n / 62 (%)</b></td>
  </tr>
  <tr>
   <td>TMB-High (≥10 mut/Mb)</td>
   <td>8 (12.9%)</td>
  </tr>
  <tr>
   <td>MSI-High</td>
   <td>2 (3.2%)</td>
  </tr>
  <tr>
   <td>HRD-positive (GIS ≥42)</td>
   <td>6 (9.7%)</td>
  </tr>
  <tr>
   <td>≥1 OncoKB on-label match (Level 1/2/3A)</td>
   <td>21 (33.9%)</td>
  </tr>
  <tr>
   <td>≥1 OncoKB off-label match (Level 3B/4)</td>
   <td>38 (61.3%)</td>
  </tr>
  <tr>
   <td>≥1 OncoKB resistance flag (R1/R2)</td>
   <td>6 (9.7%)</td>
  </tr>
 </tbody>
</table>
```

---

## Submission portal checklist (same session as primary abstract)

- [ ] Category / track: **Real-world data and patient-reported outcomes** (backup: Molecular oncology / Translational research)
- [ ] Presentation type preference: **Poster Display** primary, **Mini Oral** fallback (don't compete with AI abstract's mini-oral slot)
- [ ] Author list (≤20, same as primary abstract or rotated first/last)
- [ ] Affiliations
- [ ] Declaration of Interest from every author
- [ ] Encore policy: not previously presented — confirm
- [ ] Generative-AI disclosure: AI was used in pipeline annotation (mentioned in Methods of primary abstract); cohort abstract Methods describes only deterministic CGP + OncoKB/CiVIC/AMP tiering, no LLM mention here
- [ ] **Click "Finish Submission"** for THIS abstract specifically (separate from primary)
- [ ] **First-author note**: ESMO presenter cap is 2 abstracts/person — if you present both, you're at the limit; if a different first author for this one, both authors get a presenter slot

---

## Why this could land Poster (very likely) → Mini Oral (possible)

- Real-world cohort descriptive abstracts have ~75-85% accept rate at ESMO; n=62 is on the smaller side but acceptable for an implementation pilot
- The **33.9% on-label match rate** sits squarely within the 20-40% literature band — easy to defend
- The **BRCA1 LGR case** + the resistance-flag count (6/62) give the poster two specific clinical narratives during Q&A
- No outcome data and single-center are the main weaknesses — keep expectations at Poster, not Oral

---

## Numbers cited — re-verify with `Rscript scripts/extract_cohort_metrics.R` morning of 12 May

| Number | Source |
|---|---|
| n = 62 patients, 8 monthly batches, 11 cancer types | `inputs/TSO500-HRD/_manifest.tsv` |
| TMB-High 12.9% (8/62) | `_cohort_summary.json` → rates.tmb_high |
| MSI-High 3.2% (2/62) | `_cohort_summary.json` → rates.msi_high |
| HRD-positive 9.7% (6/62) | `_cohort_summary.json` → rates.hrd_positive |
| On-label match 33.9% (21/62) | `_cohort_summary.json` → rates.any_on_label |
| Off-label match 61.3% (38/62) | `_cohort_summary.json` → rates.any_off_label |
| Resistance flag 9.7% (6/62) | `_cohort_summary.json` → rates.any_resistance |
| Top genes: TP53 26, KRAS 16, PIK3CA 16, PTEN 11, ARID1A 9, APC 7, ESR1 6 | `_cohort_summary.json` → top_genes_ge2 (excludes KDR/AURKA known germline polymorphisms) |
| BRCA1 LGR FC 0.15 in HRD-positive breast | sample M26-0234R `[Large Rearrangements]` |
