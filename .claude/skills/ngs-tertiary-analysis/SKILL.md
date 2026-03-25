---
name: ngs-tertiary-analysis
description: >
  Run NGS tertiary analysis pipeline from BAM/VCF to ESMO-compliant clinical HTML report.
  Use when user has a cancer panel BAM/VCF file (TSO500) and needs somatic variant calling,
  annotation, OncoKB + CiVIC clinical actionability, ESCAT + AMP/ASCO/CAP classification,
  interactive visualizations, literature evidence, and a publication-quality genomics report.
  Two modes: /ngs-quick (one-shot) or /ngs-interactive (step-by-step with choices).
---

## Two Modes

This skill has two invocation modes. Choose based on context:

### `/ngs-quick` — One-Shot Pipeline Run
Use when the user already knows what they want. No questions asked — just run the pipeline and deliver the report.

**Trigger**: User provides a BAM/VCF path + tumor type, or says "run the pipeline"

**Behavior**:
1. Validate inputs (BAM/VCF exists, config valid, API keys present)
2. Run `targets::tar_make()` end-to-end
3. Report results summary + path to HTML report
4. No questions, no choices — just execute

### `/ngs-interactive` — Step-by-Step with Choices
Use when exploring data, learning the pipeline, or when the user wants control over each stage.

**Trigger**: User says "analyze this sample", "help me interpret", or asks about specific genes/variants

**Behavior**:
1. Ask about input (BAM or VCF?)
2. Ask about tumor type (with suggestions)
3. Run QC and show results — ask if acceptable
4. Show variants and let user select which to investigate
5. Present OncoKB + CiVIC findings per gene
6. Show AMP/ASCO/CAP classification and discuss
7. Generate visualizations on demand
8. Render final report with user's chosen options

---

## Pipeline Overview (v0.2.0)

| Stage | Module | What it Does |
|-------|--------|-------------|
| 0 | QC | Coverage stats, purity, on-target rate |
| 1 | Variant Calling | GATK Mutect2 somatic SNV/indel calling |
| 2 | Annotation | VEP + ClinVar + COSMIC + gnomAD |
| 3 | CNV | CNVkit copy number detection |
| 4 | Fusions | Manta structural variant / fusion calling |
| 5 | Biomarkers | TMB, MSI, HRD calculation |
| 6a | OncoKB | Clinical actionability + ESCAT tiers |
| 6b | CiVIC | Community evidence + AMP/ASCO/CAP assertions |
| 6c | AMP | Automated oncogenicity four-tier classification |
| 7 | Literature | PubMed + Scopus + treatment narratives |
| 8 | Report | ESMO 2024 HTML report with interactive plots |

## Knowledge Sources

| Source | Type | API Key Required |
|--------|------|-----------------|
| OncoKB | Clinical actionability | Yes (ONCOKB_API_KEY) |
| CiVIC | Community evidence + AMP assertions | No (free) |
| ClinVar | Pathogenicity | Via VEP |
| COSMIC | Somatic mutation catalog | Via VEP |
| gnomAD | Population frequency | Via VEP |
| PubMed | Literature evidence | Optional (PUBMED_API_KEY) |
| Scopus | Citation-ranked literature | Yes (SCOPUS_API_KEY) |
| Unpaywall | Open access PDFs | Yes (UNPAYWALL_EMAIL) |

## Interactive Visualizations (v0.2.0)

| Plot | Library | Interactive |
|------|---------|------------|
| VAF distribution | ggiraph | Hover: gene, variant, VAF%, classification |
| Gene coverage | ggiraph | Hover: mean/min coverage, % above threshold |
| CNV genome view | plotly | Zoom, pan, hover across all chromosomes |
| Circos plot | circlize | Static PNG at 300 DPI (publication quality) |
| Fusion arcs | ggiraph | Known vs novel, supporting reads |
| Biomarker gauges | ggiraph | TMB/MSI/HRD with threshold indicators |

## Classification Systems

### ESCAT (from OncoKB)
- Tier I: Standard of care (LEVEL_1/2)
- Tier II: Investigational with strong evidence (LEVEL_3A)
- Tier III: Other tumor types (LEVEL_3B)
- Tier IV: Preclinical (LEVEL_4)
- Tier X: No actionability

### AMP/ASCO/CAP (from OncoKB + CiVIC + VEP + ClinVar)
- Tier I Level A: FDA-approved / guideline-concordant
- Tier I Level B: Well-powered studies
- Tier II Level C: Different tumor type evidence
- Tier II Level D: Preclinical / oncogenic + HIGH impact
- Tier III: VUS (insufficient evidence)
- Tier IV: Benign / likely benign

## Report Security

Set `report.password` in config or `REPORT_PASSWORD` env var to encrypt the HTML report with a password gate (browser-native Web Crypto API).

## Quick Reference

```r
# Full pipeline from BAM
BAM_PATH=tumor.bam SAMPLE_ID=P001 Rscript -e 'targets::tar_make()'

# From VCF (skips QC, calling, CNV, fusions, MSI)
VCF_PATH=somatic.vcf.gz SAMPLE_ID=P001 Rscript -e 'targets::tar_make()'

# E2E test with mock ovarian cancer data
Rscript run_e2e_test.R

# Run specific targets
targets::tar_make(names = c("civic_results"))
targets::tar_make(names = c("amp_results"))
```

## Supporting Files

- `reference.md` — API docs, ESMO requirements, config schema
- `examples.md` — Example runs and tumor-specific configurations
- `R/api_clients.R` — OncoKB, PubMed, Scopus API wrappers
- `R/civic_client.R` — CiVIC GraphQL API client
- `R/amp_classification.R` — AMP/ASCO/CAP four-tier logic
- `R/plot_helpers.R` — Interactive visualization functions
- `R/esmo_helpers.R` — ESCAT classification and formatting
- `R/report_security.R` — Password-protected HTML encryption
