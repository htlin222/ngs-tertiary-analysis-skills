---
name: ngs-quick
description: >
  One-shot NGS pipeline run. No questions — validate, execute, deliver.
  Use when user provides BAM/VCF + tumor type and wants the report fast.
  Runs full targets pipeline, generates ESMO HTML report with all features.
---

## Behavior

You are a pipeline executor. Run the NGS tertiary analysis pipeline end-to-end with zero interaction. Do not ask questions — infer from context or use defaults.

## Steps

### 1. Detect Input

Check environment variables and user message for:
- `BAM_PATH` or `VCF_PATH` or `INPUT_PATH` — the input file
- `SAMPLE_ID` — sample identifier
- Tumor type — map to OncoKB code

If not provided, check if `run_e2e_test.R` mock data should be used.

### 2. Validate

Run these checks silently:
```r
source("R/validate_inputs.R")
# Check: input file exists, config valid, API keys present
```

If validation fails, report the specific issue and suggest fix. Do not proceed.

### 3. Configure

Update `config/default.yaml` if tumor type differs from default:
```r
config <- load_config()
config$sample$tumor_type <- "{user_tumor_type}"
config$sample$id <- "{sample_id}"
```

### 4. Execute

```r
Sys.setenv(BAM_PATH = "{path}", SAMPLE_ID = "{id}")
targets::tar_make()
```

Or for VCF:
```r
Sys.setenv(VCF_PATH = "{path}", SAMPLE_ID = "{id}")
targets::tar_make()
```

### 5. Report

After completion, provide:

```
Pipeline complete for {SAMPLE_ID}

Variants: {N} somatic | {N} CNV | {N} fusions
Biomarkers: TMB {score} ({class}) | MSI {status} | HRD {status}
Actionable: {N} ESCAT Tier I-III | {N} AMP Tier I-II
Knowledge: OncoKB + CiVIC ({N} evidence items)
Literature: {N} PubMed + {N} Scopus articles

Report: reports/{sample_id}/08-report/clinical_report.html
```

Then offer: "Open the report? Deploy to GitHub Pages? Run on another sample?"

## If No Input File

If the user doesn't have a BAM/VCF, offer:
1. Run `run_e2e_test.R` with mock ovarian cancer data (real API calls)
2. Download test data: `bash testdata/download_testdata.sh --small`

## Defaults

| Parameter | Default |
|-----------|---------|
| Tumor type | From config (NSCLC) |
| Min coverage | 200x |
| Min VAF | 5% |
| CiVIC | Enabled |
| Report password | None (set REPORT_PASSWORD to enable) |
