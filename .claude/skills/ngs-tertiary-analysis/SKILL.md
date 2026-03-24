---
name: ngs-tertiary-analysis
description: >
  Run NGS tertiary analysis pipeline from BAM to ESMO-compliant clinical HTML report.
  Use when user has a cancer panel BAM file (TSO500) and needs somatic variant calling,
  annotation, OncoKB clinical actionability, ESCAT classification, literature evidence,
  and a publication-quality genomics report.
---

## Overview

End-to-end NGS tertiary analysis pipeline for cancer genomics. Takes a BAM file from
TSO500 (or similar comprehensive genomic profiling panel) and produces an HTML clinical
report following ESMO 2024 guidelines.

## Quick Start

```r
# Set BAM path and sample config
Sys.setenv(BAM_PATH = "path/to/tumor.bam")
Sys.setenv(SAMPLE_ID = "PATIENT_001")

# Run full pipeline
targets::tar_make()

# Or run specific stages
targets::tar_make(names = c("qc_results"))
targets::tar_make(names = c("oncokb_results"))
```

## Pipeline Stages

| Stage | Dir | Input | Output |
|-------|-----|-------|--------|
| QC | `00-qc/` | BAM | coverage stats, purity estimate |
| Variant Calling | `01-variant-calling/` | BAM | filtered somatic VCF |
| Annotation | `02-annotation/` | VCF | annotated VCF (VEP, ClinVar, COSMIC) |
| CNV | `03-cnv/` | BAM | copy number segments |
| Fusions | `04-fusions/` | BAM | fusion calls |
| Biomarkers | `05-biomarkers/` | VCF, BAM | TMB, MSI, HRD |
| Clinical | `06-clinical-annotation/` | all variants | OncoKB + ESCAT tiers |
| Literature | `07-literature/` | key variants | PubMed/Scopus narratives |
| Report | `08-report/` | all above | ESMO HTML report |

## Configuration

Edit `config/default.yaml` to customize:
- `sample.tumor_type`: OncoKB tumor type (e.g., "NSCLC", "COAD", "BRCA")
- `qc.min_mean_coverage`: minimum acceptable coverage (default: 200x)
- `variant_calling.min_vaf`: VAF threshold (default: 0.05)
- `biomarkers.tmb.panel_coding_size_mb`: panel coding size for TMB denominator

## API Keys Required

Set in `.env` file:
- `ONCOKB_API_KEY` - OncoKB API token (https://www.oncokb.org/apiAccess)
- `PUBMED_API_KEY` - NCBI E-utilities API key
- `SCOPUS_API_KEY` - Elsevier Scopus API key
- `UNPAYWALL_EMAIL` - Email for Unpaywall API
- `EMBASE_API_KEY` - Embase API key (optional)

## Output Structure

All outputs go to `reports/{sample_id}/` (gitignored for patient privacy):

```
reports/PATIENT_001/
├── 00-qc/           # Coverage stats, QC plots
├── 01-variant-calling/  # VCF files
├── ...              # One dir per stage
└── 08-report/
    └── clinical_report.html  # Final ESMO report
```

## System Dependencies

Run `make setup-all` or use Docker:
- GATK 4.x, samtools, bcftools
- Ensembl VEP with GRCh38 cache
- CNVkit
- Quarto CLI

## Troubleshooting

- **OncoKB 401**: Check ONCOKB_API_KEY in .env, ensure approved account
- **VEP not found**: Run `make setup-vep` or use Docker
- **Low coverage warning**: Check BAM quality, may need to adjust `qc.min_mean_coverage`
- **Rate limiting**: Pipeline uses httr2 retry; increase wait in `R/api_clients.R`

## Reference Files

- `reference.md` - Full API docs, ESMO report requirements, config schema
- `examples.md` - Example runs, sample outputs, common configurations
- `R/api_clients.R` - OncoKB/PubMed/Scopus API wrapper functions
- `R/esmo_helpers.R` - ESCAT classification and ESMO formatting logic
