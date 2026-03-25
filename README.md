# NGS Tertiary Analysis Pipeline

R-based pipeline for cancer genomics: **BAM &rarr; ESMO 2024-compliant HTML clinical report**.

Annotates somatic variants with [OncoKB](https://www.oncokb.org/) + [CiVIC](https://civicdb.org/), classifies clinical actionability using ESCAT and [AMP/ASCO/CAP](https://doi.org/10.1016/j.jmoldx.2016.10.002) four-tier systems, generates interactive visualizations, retrieves treatment-focused literature from PubMed/Scopus, and renders a publication-quality Quarto report.

**[Live demo report](https://htlin222.github.io/ngs-tertiary-analysis-skills/)** (mock ovarian cancer case)

## Pipeline Overview

```
BAM file (TSO500 panel)          VCF file (pre-called variants)
  |                                |
  +-- 00-qc/        (BAM only)    |
  +-- 01-calling/   (BAM only)    |
  |                                |
  +------------ OR ----------------+
  |                                |
  +-- 02-annotation/    Ensembl VEP (HGVS, ClinVar, COSMIC, gnomAD)
  +-- 03-cnv/           CNVkit (BAM only, skipped for VCF)
  +-- 04-fusions/       Manta (BAM only, skipped for VCF)
  +-- 05-biomarkers/    TMB (both), MSI (BAM only), HRD (BAM only)
  +-- 06-clinical/      OncoKB + CiVIC + ESCAT + AMP/ASCO/CAP classification
  +-- 07-literature/    PubMed & Scopus, treatment narratives, AMA citations
  +-- 08-report/        Quarto HTML report with interactive plots (ESMO 2024)
  |
  v
reports/{sample_id}/08-report/clinical_report.html
```

Orchestrated by the [`targets`](https://docs.ropensci.org/targets/) R package. Each stage reads from the previous and writes to `reports/{sample_id}/{stage}/`.

## Quick Start

### 1. Install dependencies

```bash
# System tools (macOS)
make setup-system      # samtools, bcftools, htslib

# R packages
Rscript -e 'renv::restore()'

# Optional: GATK, VEP, CNVkit (for stages 0-4)
make setup-all
```

### 2. Set up API keys

Create a `.env` file in the project root:

```env
ONCOKB_API_KEY=your-oncokb-token
PUBMED_API_KEY=your-ncbi-api-key
SCOPUS_API_KEY=your-elsevier-key
UNPAYWALL_EMAIL=your-email@example.com
EMBASE_API_KEY=your-embase-key
```

- **OncoKB**: Register at [oncokb.org/apiAccess](https://www.oncokb.org/apiAccess)
- **PubMed**: Get key at [ncbi.nlm.nih.gov/account/settings](https://www.ncbi.nlm.nih.gov/account/settings/)
- **Scopus**: Apply at [dev.elsevier.com](https://dev.elsevier.com/)

### 3. Run the pipeline

**From BAM** (full pipeline — stages 0-8):

```bash
BAM_PATH=path/to/tumor.bam SAMPLE_ID=PATIENT_001 make run
```

**From VCF** (skips QC + calling — stages 2-8):

```bash
VCF_PATH=path/to/somatic.vcf.gz SAMPLE_ID=PATIENT_001 make run
```

**Auto-detect** (BAM or VCF):

```bash
INPUT_PATH=path/to/file SAMPLE_ID=PATIENT_001 make run
```

Or in R:

```r
Sys.setenv(BAM_PATH = "path/to/tumor.bam", SAMPLE_ID = "PATIENT_001")
# or: Sys.setenv(VCF_PATH = "path/to/somatic.vcf.gz", SAMPLE_ID = "PATIENT_001")
targets::tar_make()
```

When VCF is provided, the pipeline auto-skips QC, variant calling, CNV, fusion, and MSI stages (which require BAM), and enters at annotation. TMB is still calculated from the VCF.

### 4. Run the demo (no BAM needed)

The demo uses mock ovarian cancer variants with **real API calls** to OncoKB, PubMed, and Scopus:

```bash
Rscript run_e2e_test.R
```

This produces a complete report at `reports/EGA_OV_TSO500_001/08-report/clinical_report.html`.

## Output

All outputs go to `reports/{sample_id}/` (gitignored for patient privacy):

```
reports/PATIENT_001/
  00-qc/                  coverage_stats.tsv, qc_summary.tsv
  01-variant-calling/     mutect2_raw.vcf.gz, filtered.vcf.gz
  02-annotation/          vep_annotated.vcf.gz, merged_annotations.tsv
  03-cnv/                 cnvkit_segments.tsv, cnv_plot.png
  04-fusions/             fusions.tsv
  05-biomarkers/          tmb_result.json, msi_result.json, hrd_result.json
  06-clinical-annotation/ oncokb_results.json, escat_tiers.csv, civic_results.rds
  07-literature/          pubmed_hits.json, scopus_hits.json, narratives.json
  08-report/              clinical_report.html   <-- final report
```

## Report Features

The HTML report follows [ESMO 2024 NGS reporting guidelines](https://www.annalsofoncology.org/article/S0923-7534(24)00111-X/fulltext):

- **Patient & sample info** with sequencing platform and specimen details
- **QC dashboard** with pass/fail thresholds for coverage, mapping, purity
- **Annotation methodology** table (VEP, SIFT, PolyPhen-2, gnomAD, ClinVar, CiVIC, AMP)
- **Somatic variants** with HGVS nomenclature, Ensembl transcript IDs, VAF, read depth
- **Interactive VAF distribution** (hover for gene, variant, classification)
- **Mutation confidence ranking** based on VAF vs detection threshold
- **Copy number alterations** with interactive genome-wide CNV plot (plotly zoom/pan)
- **Circos plot** showing CNV tracks + fusion arcs (publication quality, 300 DPI)
- **Gene fusions** with interactive visualization (known vs novel)
- **Biomarker gauges** for TMB, MSI, HRD with threshold indicators
- **ESCAT clinical actionability** with color-coded tiers, OncoKB + CiVIC sources
- **AMP/ASCO/CAP oncogenicity classification** (four-tier system per Li et al. 2017)
- **CiVIC community evidence** with variant evidence and AMP assertions
- **Treatment-focused literature** with AMA-style references
- **Per-gene coverage** interactive heatmap (ESMO Level A requirement)
- **Password protection** (optional, browser-native access control)

Rendered with Quarto (cosmo theme, left-side TOC, `gt` tables, ggiraph/plotly/circlize plots, self-contained HTML).

## Knowledge Sources

| Source | Type | API Key |
|--------|------|---------|
| [OncoKB](https://www.oncokb.org/) | Clinical actionability, ESCAT tiers | Required |
| [CiVIC](https://civicdb.org/) | Community evidence, AMP/ASCO/CAP assertions | Free (no key) |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Pathogenicity classification | Via VEP |
| [COSMIC](https://cancer.sanger.ac.uk/cosmic) | Somatic mutation catalog | Via VEP |
| [gnomAD](https://gnomad.broadinstitute.org/) | Population allele frequency | Via VEP |
| [PubMed](https://pubmed.ncbi.nlm.nih.gov/) | Literature evidence | Optional |
| [Scopus](https://www.scopus.com/) | Citation-ranked literature | Required |
| [Unpaywall](https://unpaywall.org/) | Open access PDF links | Required |

## Configuration

Edit `config/default.yaml` to customize:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sample.tumor_type` | `NSCLC` | OncoKB tumor type code (e.g., HGSOC, COAD, BRCA) |
| `qc.min_mean_coverage` | `200` | Minimum acceptable mean coverage |
| `variant_calling.min_vaf` | `0.05` | VAF threshold for variant filtering |
| `biomarkers.tmb.panel_coding_size_mb` | `1.94` | TSO500 coding region size |
| `biomarkers.msi.msi_h_threshold` | `20.0` | % unstable sites for MSI-H |
| `biomarkers.hrd.score_threshold` | `42` | HRD-positive threshold |
| `literature.pubmed.max_results` | `10` | Max PubMed articles per variant |
| `clinical_annotation.civic.enabled` | `true` | Enable CiVIC community evidence |
| `report.password` | `null` | Password-protect the HTML report |

## Project Structure

```
R/                      Shared utilities
  utils.R               Config loading, logging, .env, external tool wrappers
  api_clients.R         OncoKB, PubMed, Scopus, Unpaywall API wrappers
  civic_client.R        CiVIC GraphQL API client (free, no key)
  esmo_helpers.R        ESCAT tier classification, ESMO formatting
  amp_classification.R  AMP/ASCO/CAP four-tier oncogenicity classification
  plot_helpers.R        Interactive visualizations (ggiraph, plotly, circlize)
  report_security.R     Password-protected HTML report

00-qc/ .. 08-report/    Pipeline stage scripts (one dir per stage)
config/                 Pipeline config (default.yaml, panel BED)
tests/                  testthat unit tests (169 pass, 0 fail)
docs/                   SETUP.md, ESMO_GUIDELINES.md, live demo report
.claude/skills/         Claude Code skills (ngs-quick, ngs-interactive)

_targets.R              Pipeline DAG definition
Makefile                System dependency setup + pipeline commands
Dockerfile              Reproducible container environment
run_e2e_test.R          End-to-end demo with real API calls
run_standalone_report.R Fallback HTML report generator (no Quarto)
```

## Running Specific Stages

```r
# Visualize the pipeline DAG
targets::tar_visnetwork()

# Run only QC
targets::tar_make(names = c("qc_results"))

# Run from annotation onward
targets::tar_make(names = c("merged_annotations", "oncokb_results",
                             "escat_tiers", "literature_results",
                             "clinical_report"))
```

## Tests

```bash
make test
# or
Rscript -e 'testthat::test_dir("tests")'
```

Tests cover: utility functions, ESCAT mapping, OncoKB response parsing, CiVIC API, AMP/ASCO/CAP classification, interactive visualizations, report encryption. API tests are skipped unless keys are set.

## Docker

```bash
docker build -t ngs-tertiary .
docker run -v $(pwd)/reports:/app/reports \
           --env-file .env \
           -e BAM_PATH=/data/tumor.bam \
           -e SAMPLE_ID=PATIENT_001 \
           ngs-tertiary
```

## Key APIs

| API | Purpose | Endpoint |
|-----|---------|----------|
| [OncoKB](https://www.oncokb.org/) | Mutation/CNA/fusion clinical actionability | `www.oncokb.org/api/v1/annotate/*` |
| [CiVIC](https://civicdb.org/) | Community variant evidence + AMP assertions | `civicdb.org/api/graphql` |
| [PubMed](https://www.ncbi.nlm.nih.gov/home/develop/api/) | Literature search (E-utilities) | `eutils.ncbi.nlm.nih.gov` |
| [Scopus](https://dev.elsevier.com/) | Citation-ranked literature | `api.elsevier.com/content/search/scopus` |
| [Unpaywall](https://unpaywall.org/products/api) | Open access PDF links | `api.unpaywall.org/v2/{doi}` |

## References

- [ESMO 2024 NGS Recommendations](https://www.annalsofoncology.org/article/S0923-7534(24)00111-X/fulltext) — Recommendations for NGS use in advanced cancer
- [ESMO Clinical Reporting Guidelines](https://www.annalsofoncology.org/article/S0923-7534(24)01011-1/fulltext) — How to report NGS results
- [ESCAT Classification](https://www.annalsofoncology.org/article/S0923-7534(19)34179-1/fulltext) — Scale for Clinical Actionability of molecular Targets
- [AMP/ASCO/CAP Guidelines](https://doi.org/10.1016/j.jmoldx.2016.10.002) — Li et al., J Mol Diagn 2017; somatic variant classification
- [CiVIC Knowledgebase](https://civicdb.org/) — Open-access community clinical interpretation
- [HGVS Nomenclature](https://varnomen.hgvs.org/) — Variant naming standards
- [OncoKB API Docs](https://api.oncokb.org/) — Precision oncology knowledge base

## License

Research use only. OncoKB annotations require an [academic license](https://www.oncokb.org/apiAccess).
