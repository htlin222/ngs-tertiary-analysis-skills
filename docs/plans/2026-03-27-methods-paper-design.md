# Methods Paper: Open-Source NGS Tertiary Analysis Pipeline — Design Document

## Paper Target

**Journal**: Bioinformatics (Oxford) or BMC Bioinformatics
**Type**: Application Note / Software Paper
**Title (working)**: "ngs-tertiary-analysis: An open-source pipeline for automated clinical interpretation of cancer panel sequencing with multi-source evidence integration"

## Key Selling Points (vs. Illumina Connected Insights)

1. **Open-source** — fully reproducible, no vendor lock-in
2. **Multi-source evidence** — OncoKB + CiVIC + VEP/ClinVar/COSMIC (vs. proprietary CKB)
3. **Dual classification** — ESCAT + AMP/ASCO/CAP (not just one system)
4. **Interactive visualizations** — publication-ready plots embedded in HTML
5. **Local data processing** — patient data never leaves infrastructure
6. **Automated end-to-end** — BAM/VCF → clinical report in one command

## Validation Strategy

### Dataset 1: SEQC2 HCC1395 (Variant Detection Benchmark)

**Source**: BioProject PRJNA489865 (SRR8861483), well-characterized breast cancer cell line
**Truth set**: SEQC2 consortium published known somatic variants (Fang et al., Nature Biotechnology 2021)
**What it validates**: Variant calling sensitivity/specificity at different VAF thresholds

Metrics:
- Sensitivity (TP / TP+FN) for SNVs and indels separately
- Precision (TP / TP+FP)
- F1 score
- VAF correlation (called vs truth)
- Detection limit analysis (sensitivity by VAF bucket: >20%, 10-20%, 5-10%)

### Dataset 2: ClinVar Pathogenicity Concordance

**Source**: ClinVar VCF (public), filtered to TSO500 panel genes
**What it validates**: AMP/ASCO/CAP classification accuracy
**Approach**:
- Extract ClinVar variants in TSO500 panel genes with definitive classifications
- Run through pipeline's AMP classification
- Compare our Tier assignments against ClinVar expert-reviewed pathogenicity
- Expected: pathogenic/likely_pathogenic → Tier I-II, benign → Tier IV, VUS → Tier III

Metrics:
- Classification concordance rate
- Cohen's kappa for tier agreement
- Discordance analysis (which tiers disagree and why)

### Dataset 3: OncoKB vs CiVIC Evidence Concordance

**Source**: Run pipeline on a set of well-known oncogenic variants across multiple tumor types
**What it validates**: Multi-source evidence integration adds value beyond single-source
**Approach**:
- Select 100 well-characterized variants from COSMIC (top 100 most frequent)
- Query OncoKB and CiVIC for each
- Measure: overlap, unique findings per source, tier agreement

Metrics:
- Venn diagram: variants with evidence in OncoKB-only, CiVIC-only, both
- Cases where CiVIC evidence changes the AMP tier (added value)
- ESCAT vs AMP tier concordance

### Dataset 4: Runtime Benchmarks

**Source**: SEQC2 BAM + synthetic VCFs of varying sizes (10, 50, 100, 500 variants)
**What it validates**: Pipeline scalability and optimization impact

Metrics:
- Wall-clock time per stage
- Before/after optimization (crew, batch API, VEP fork)
- API call count reduction
- Memory usage peak

## Manuscript Structure (Bioinformatics Application Note, 2 pages)

### Abstract (150 words)
- Problem: clinical interpretation bottleneck, commercial tools expensive
- Solution: open-source pipeline with multi-source evidence
- Results: X% sensitivity on SEQC2, Y% AMP concordance, Z min runtime
- Availability: GitHub URL + Docker

### Introduction (300 words)
- NGS panels becoming standard in precision oncology
- ESMO 2024 guidelines mandate structured reporting
- Gap: free tools lack multi-source evidence + automated classification
- We present: open-source pipeline, TSO500, ESMO 2024 compliant

### Methods (500 words)
- Pipeline architecture (8 stages, targets DAG)
- Knowledge sources (OncoKB, CiVIC, VEP, ClinVar, COSMIC, PubMed, Scopus)
- Classification systems (ESCAT + AMP/ASCO/CAP)
- Benchmarking: SEQC2 truth set, ClinVar concordance, multi-source comparison
- Technology: R, Bioconductor, Quarto, Docker

### Results (400 words)
- Table 1: Variant detection performance (SEQC2)
- Table 2: AMP classification concordance (ClinVar)
- Figure 1: Pipeline overview diagram
- Figure 2: OncoKB vs CiVIC evidence Venn + added value
- Figure 3: Runtime benchmark (stages, optimization impact)
- Supplementary: interactive report screenshot, full concordance tables

### Discussion (200 words)
- Comparable to commercial solutions for tertiary analysis
- Limitations: panel-based HRD, no germline, no RNA fusions without RNA-seq
- CiVIC adds value in X% of cases beyond OncoKB alone
- Future: liquid biopsy, RNA integration, ClinGen allele registry

## Implementation: Benchmark Framework

### Directory Structure

```
benchmarks/
  datasets/
    seqc2/           # SEQC2 truth set VCF + pipeline results
    clinvar/          # ClinVar panel-filtered variants
    cosmic_top100/    # Top 100 COSMIC variants for concordance
  scripts/
    01_download_seqc2.R        # Download + prepare SEQC2 data
    02_run_seqc2_benchmark.R   # Run pipeline on SEQC2
    03_variant_accuracy.R      # Calculate sensitivity/specificity
    04_clinvar_concordance.R   # AMP vs ClinVar comparison
    05_evidence_concordance.R  # OncoKB vs CiVIC overlap
    06_runtime_benchmark.R     # Timing across stages
    07_generate_figures.R      # Publication figures
    08_generate_tables.R       # Publication tables
  results/
    figures/          # Publication-ready figures (300 DPI PNG + PDF)
    tables/           # LaTeX/CSV tables for manuscript
  manuscript/
    main.qmd          # Quarto manuscript (renders to PDF/DOCX)
    references.bib     # BibTeX references
    bioinformatics.csl # Journal citation style
```

### Key Scripts

**01_download_seqc2.R**: Download SEQC2 truth set from supplementary data, filter to TSO500 regions, prepare comparison VCF.

**03_variant_accuracy.R**: Compare pipeline VCF vs truth VCF using `vcfR` or `VariantAnnotation`. Calculate TP/FP/FN per variant type, generate ROC-like curve by VAF threshold.

**04_clinvar_concordance.R**: Download latest ClinVar VCF, intersect with TSO500 BED, run through AMP classification, compare tiers. Cohen's kappa.

**05_evidence_concordance.R**: For top 100 COSMIC variants across 5 tumor types, query OncoKB + CiVIC, tabulate evidence overlap. Generate Venn diagram with `ggVennDiagram`.

**07_generate_figures.R**: Use `R/plot_helpers.R` functions + `patchwork` for multi-panel figures. Pipeline diagram with `DiagrammeR`. All at 300 DPI.

### Quarto Manuscript

Use Quarto for reproducible manuscript rendering:
- `manuscript/main.qmd` → PDF (journal submission) + DOCX (review)
- Inline R code chunks that read from `results/` directory
- Automatic figure/table numbering
- BibTeX bibliography with `bioinformatics.csl`

## Timeline

| Week | Task |
|------|------|
| 1 | Set up benchmark framework, download SEQC2 data |
| 2 | Run SEQC2 benchmark, calculate variant accuracy |
| 3 | ClinVar concordance + OncoKB/CiVIC comparison |
| 4 | Runtime benchmarks, generate all figures/tables |
| 5 | Write manuscript draft in Quarto |
| 6 | Internal review, revise, submit |
