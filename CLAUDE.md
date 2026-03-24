# NGS Tertiary Analysis Pipeline

## Project Overview
R-based pipeline: BAM → ESMO 2024-compliant HTML clinical report for cancer genomics.
Uses `targets` for orchestration, Quarto for reporting.

## Commands
- `make setup-all` — install all system + R dependencies
- `targets::tar_make()` — run full pipeline
- `make test` — run testthat suite
- `BAM_PATH=x SAMPLE_ID=y Rscript -e 'targets::tar_make()'` — run on specific sample

## Structure
- `R/` — shared utilities (utils.R, api_clients.R, esmo_helpers.R)
- `00-qc/` through `08-report/` — pipeline stages
- `config/default.yaml` — all pipeline parameters
- `reports/{sample_id}/` — per-sample outputs (gitignored, patient data)
- `.claude/skills/ngs-tertiary-analysis/` — Claude Code skill files

## Key Conventions
- All outputs go to `reports/{sample_id}/{stage}/`
- External tools (GATK, VEP, CNVkit) called via `run_tool()` in R/utils.R
- API keys loaded from `.env` via `load_env()` / `get_api_key()`
- Config loaded from `config/default.yaml` via `load_config()`
- Use `httr2` for all HTTP requests (not httr)
- Use `gt` for publication tables, `ggplot2` for figures
- Bioconductor packages: `Rsamtools`, `VariantAnnotation`, `GenomicRanges`

## Testing
- `tests/` uses testthat
- `tests/fixtures/` has small BAM/VCF for unit tests
- `testdata/` has download script for full-size test BAM

## Do NOT
- Commit anything in `reports/` (patient data)
- Hardcode file paths — use config or here::here()
- Use httr (deprecated) — use httr2
- Skip error checking on system2() calls — use run_tool()
