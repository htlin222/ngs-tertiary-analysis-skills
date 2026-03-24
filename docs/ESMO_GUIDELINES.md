# ESMO 2024 NGS Clinical Reporting Guidelines

This document summarizes the 2024 ESMO guidelines for next-generation sequencing (NGS) clinical reports. These are the current best practice standards for variant interpretation and clinical actionability assessment.

**Reference**: Mateo et al. "Recommendations for NGS-based diagnostics in oncology: standardization of variant classification and reporting" Annals of Oncology (2024)
https://www.annalsofoncology.org/article/S0923-7534(24)00111-X/fulltext

## Overview

ESMO 2024 guidelines establish standardized requirements for NGS clinical reports to ensure consistency, accuracy, and clinical utility across laboratories. These guidelines focus on:

- Mandatory reporting of pathogenic variants with high confidence
- Structured assessment of clinical actionability
- Transparent documentation of assay performance
- Explicit reporting of genes with insufficient coverage
- Clear communication of genomic findings to clinicians

## Level A Requirements (Mandatory)

All clinical NGS reports **must** include these elements:

### 1. Patient and Sample Identification

- Full patient name (or secure identifier)
- Date of birth (age at testing)
- Specimen type (e.g., formalin-fixed paraffin-embedded tumor, fresh frozen, cfDNA)
- Specimen collection date
- Sample processing date
- Anatomical site of origin (primary tumor location)
- Histological or cytological diagnosis
- Unique sample identifier (accession number)

### 2. Assay Description and Performance Characteristics

**Assay Method**:
- Sequencing platform used (Illumina, Ion Torrent, PacBio, Oxford Nanopore)
- Library preparation method
- Capture method (if applicable): whole exome, targeted panel, whole genome
- Targeted genes and genomic regions (if panel-based)
- Depth of sequencing strategy

**Performance Metrics**:
- Mean sequencing depth across targeted region
- Minimum sequencing depth threshold applied
- Base quality cutoff
- Mapping quality cutoff
- Variant calling algorithm and parameters
- Bioinformatics pipeline version

**Validation Status**:
- Clinical Laboratory Improvement Amendments (CLIA) certification or equivalent
- Clinical validity evidence
- Analytic validity (sensitivity, specificity)
- Turnaround time

### 3. Quality Control Metrics and Sample-Specific Performance

**Sample-Level QC**:
- Percentage of target region with ≥20x coverage
- Percentage of target region with ≥50x coverage
- Percentage of target region with ≥100x coverage (if applicable)
- Mean coverage for genes of clinical interest
- On-target alignment percentage
- Duplicate read percentage
- Insert size distribution (for paired-end sequencing)

**Pass/Fail Criteria**:
- For whole genome: minimum average coverage (typically ≥30x)
- For exomes: minimum average coverage (typically ≥100x)
- For panels: minimum average coverage (typically ≥500x)
- Genes failing QC thresholds must be explicitly listed

### 4. Genomic Alterations with Functional Annotation

**Variant Reporting**:
- All detected pathogenic and likely pathogenic variants (ACMG classification)
- Variants of uncertain significance (VUS) if clinically relevant
- Do NOT report benign/likely benign variants (except if explicitly requested)

**For Each Variant, Report**:
- **Gene name** (HUGO Gene Nomenclature Committee approved)
- **Protein change** (HGVS p. notation, e.g., p.Val600Glu)
- **Nucleotide change** (HGVS c. notation, e.g., c.1799T>A)
- **Variant type** (missense, frameshift, nonsense, splice site, indel, etc.)
- **Allele frequency** in sample (variant allele fraction, VAF; e.g., 45.3%)
- **Genomic coordinates** (chromosome, start position, reference/alternate alleles)
- **dbSNP ID** or ClinVar accession (if known polymorphism)
- **Functional impact prediction** (in-silico tools: SIFT, PolyPhen-2, REVEL)
- **Population frequency** (gnomAD, 1000 Genomes if allele frequency <0.01)
- **Zygosity** (heterozygous, homozygous, compound heterozygous context)

**Functional Annotation**:
- Effect on protein structure/function
- Known pathogenicity from ClinVar, COSMIC, OncoKB
- Association with specific cancer phenotypes
- Transcript affected (NM accession number)
- Exon/intron location

### 5. Clinical Actionability Assessment Using ESCAT

All variants must be assessed for clinical actionability using the **ESCAT (Esmo Scale for Clinical Actionability of molecular Targets)** framework:

- **ESCAT Tier I** variants require action and have strong clinical evidence
- **ESCAT Tier II** variants should be considered with moderate evidence
- **ESCAT Tier III** variants are emerging with limited evidence
- **ESCAT Tier IV** variants are in preclinical research stage
- **ESCAT Tier X** variants have no clinical action

Full tier descriptions in [ESCAT Tiers](#escat-tiers) below.

**For Each Actionable Variant**:
- ESCAT tier assignment
- Associated therapy or clinical trial
- Evidence level (strong, moderate, limited, preclinical)
- Supporting publications (at minimum, key references)
- Whether therapy is approved, in clinical trial, or emerging

### 6. Genes Not Meeting Minimum Coverage (Explicit Requirement)

**Critical Addition**: ESMO 2024 mandates explicit reporting of any:

- Genes targeted by the assay but failing QC coverage thresholds
- Specific exons or regions with suboptimal coverage
- Explanation of why coverage was insufficient (e.g., difficult-to-sequence regions, sample quality)
- Recommended follow-up testing if high clinical suspicion

**Example Table**:
| Gene | Region | Target Coverage | Achieved Coverage | Status |
|------|--------|-----------------|-------------------|--------|
| BRCA1 | Exon 5 | ≥500x | 287x | FAIL |
| TP53 | Full gene | ≥500x | 521x | PASS |

## Level B Requirements (Recommended)

These elements enhance clinical utility and should be included when possible:

### 1. Germline Variant Flagging

- Flag variants with known germline pathogenicity (even in somatic analysis)
- Recommend genetic counseling for pathogenic germline variants
- Note: ESMO guidance recommends returning germline BRCA1/2 findings to patients
- Distinguish variants by origin (somatic vs. germline) where determinable

### 2. Tumor Mutational Burden (TMB)

- Total number of somatic mutations per megabase (mut/Mb)
- TMB categories: Low (<5), Intermediate (5-20), High (>20) mut/Mb
- Threshold for immunotherapy eligibility per regulatory guidelines
- Note: Established primarily for atezolizumab, pembrolizumab, nivolumab

### 3. Microsatellite Instability (MSI) Status

- MSI status (microsatellite stable [MSS] vs. microsatellite instability [MSI-H/dMMR])
- If MSI detected: mismatch repair (MMR) gene mutations identified
- Clinical implication: MSI-H is FDA-approved biomarker for checkpoint inhibitors
- Mononucleotide repeat analysis recommended for robust MSI detection

### 4. Homologous Recombination Deficiency (HRD)

- HRD score (genomic instability measure)
- BRCA1/2 status and functional impact
- Implication for platinum-based chemotherapy
- Recommendation: Consider HRD testing in platinum-sensitive tumors

### 5. Literature References

- Cite peer-reviewed publications supporting clinical interpretation
- Include PMID (PubMed ID) or DOI for each reference
- Prioritize recent guidelines (ESMO, NCCN, CAP) over older literature
- Reference databases: ClinVar, COSMIC, OncoKB

## ESCAT Tiers

The ESCAT scale provides structured assessment of clinical actionability. All variants reported should include their assigned tier.

### Tier I (Strong Clinical Evidence)

**Definition**: Variants with strong evidence of clinical actionability; therapy is FDA/EMA-approved with strong supporting evidence.

**Criteria**:
- Drug/therapy is approved by regulatory authority (FDA, EMA)
- Published evidence of clinical benefit in multiple trials (Phase III minimum)
- Variant-specific therapies with clear benefit

**Examples**:
- EGFR exon 19 deletions or L858R in lung cancer → erlotinib, gefitinib, afatinib (TKI therapy approved)
- BRAF V600E mutations in melanoma → vemurafenib, dabrafenib (approved with strong evidence)
- HER2 amplifications in breast cancer → trastuzumab (Herceptin), approved standard therapy
- BRCA1/2 pathogenic variants in ovarian cancer → platinum-based chemotherapy, PARP inhibitors

**Clinical Action**: Recommend approved targeted therapy or standard care modification

---

### Tier II (Moderate Clinical Evidence)

**Definition**: Variants with moderate clinical evidence; therapy is approved or in advanced clinical trials with reasonable supporting evidence.

**Criteria**:
- Therapy FDA/EMA-approved or in late-stage clinical trial (Phase II-III)
- Published evidence from controlled trials, though possibly single-arm or smaller cohorts
- Clear molecular and biochemical rationale

**Examples**:
- KRAS G12C mutations in non-small cell lung cancer → sotorasib, adagrasib (approved 2021-2023)
- MET exon 14 skipping mutations in lung cancer → capmatinib, tepotinib
- ROS1 fusions in lung cancer → crizotinib, entrectinib
- ALK translocations in lung cancer → alectinib, crizotinib
- NRAS Q61 mutations in melanoma → MEK inhibitors (binimetinib, trametinib)

**Clinical Action**: Discuss targeted therapy options; recommend enrollment in clinical trials when available

---

### Tier III (Emerging Evidence)

**Definition**: Variants with limited but promising clinical evidence; potential therapeutic options exist but evidence is incomplete.

**Criteria**:
- Therapy in early-stage clinical trials (Phase I-II) or preclinical development
- Limited published data, often from case reports or small series
- Mechanistic rationale supports therapeutic targeting
- Ongoing research expected to clarify utility

**Examples**:
- IDH1/2 mutations in glioma, AML → ivosidenib, enasidenib (some approvals, limited tumor types)
- FGFR fusions/mutations in various cancers → erdafitinib, rogaratinib (limited approvals)
- PTEN loss in multiple tumor types → various targeted approaches (preclinical to Phase II)
- Mutation burden (TMB-high) in microsatellite stable tumors → immune checkpoint inhibitors (emerging)
- PIK3CA mutations in breast cancer → alpelisib + endocrine therapy (recent Tier II upgrades)

**Clinical Action**: Consider for clinical trial enrollment; discuss emerging options with patient

---

### Tier IV (Preclinical Evidence)

**Definition**: Variants with preclinical or very early clinical evidence; mechanistic rationale exists but minimal human data.

**Criteria**:
- No FDA/EMA approvals
- Only preclinical evidence (cell lines, animal models) or very early Phase I data
- Potential targets identified in laboratory studies

**Examples**:
- Rare fusion partners without established therapies
- Novel mutations in genes with limited clinical history
- Variants with in-vitro drug sensitivity data but no clinical translation yet
- Cancer predisposition genes in somatic analysis without clear therapeutic implications

**Clinical Action**: Note for research interest; not recommended for current clinical decision-making; consider enrollment in investigational studies

---

### Tier X (No Clinical Action)

**Definition**: Variants with no established or foreseeable clinical actionability.

**Criteria**:
- No therapeutic options available
- Limited evidence of pathogenicity
- Incidental finding without clinical significance in this disease context

**Examples**:
- VUS in genes without strong disease association
- Known polymorphisms classified as benign
- Variants of unclear functional impact in non-actionable genes
- Frameshift mutations in tumor suppressor genes without approved therapies

**Clinical Action**: No specific action; may mention for completeness or research relevance only

---

## Variant Nomenclature

### HGVS Notation (Required)

All variants must use Human Genome Variation Society (HGVS) nomenclature:

#### Protein Nomenclature (p. notation)

**Format**: p.[protein change]

**Examples**:
- `p.Val600Glu` – Missense mutation (valine to glutamic acid at position 600)
- `p.Asp235fs` – Frameshift mutation (deletion/insertion causing frameshift at position 235)
- `p.Ter123Lys` – Nonsense mutation converted by suppressor tRNA
- `p.Ile1234del` – In-frame deletion
- `p.Asp235_Ile237del` – Multi-amino acid deletion

#### Nucleotide Nomenclature (c. notation)

**Format**: c.[DNA change]

**Examples**:
- `c.1799T>A` – Missense mutation (thymine to adenine at position 1799)
- `c.235_237delTAC` – In-frame deletion (three nucleotides deleted)
- `c.234+2T>G` – Splice site mutation (position +2 of intron)
- `c.1234_1235insAT` – Insertion of AT nucleotides
- `c.1234dupA` – Duplication of adenine nucleotide

### Additional Information

- **Transcript reference**: Always specify NCBI NM accession (e.g., NM_007294.4 for BRCA1)
- **GRCh38 coordinates**: Provide genomic location (chr17:g.41219625T>A for BRCA1 c.1799T>A)
- **Exon/intron**: Note affected exon or intron for clarity

## Fusion Nomenclature

### Standard Format

Fusions are reported as: `geneA::geneB` with exon numbers included:

**Format**: `geneA(exonX)::geneB(exonY)`

**Examples**:
- `EML4::ALK(exon1)::ALK(exon20)` – EML4-ALK fusion (common in lung cancer)
- `BCR::ABL(exon1)::ABL(exon2)` – BCR-ABL fusion (Philadelphia chromosome in CML)
- `TMPRSS2::ERG(exon2)::ERG(exon4)` – TMPRSS2-ERG fusion in prostate cancer
- `EWSR1::FLI1(exon5)::FLI1(exon9)` – EWSR1-FLI1 fusion in Ewing sarcoma

### Additional Details

- Report breakpoint at nucleotide level when possible
- Specify if in-frame or out-of-frame (critical for functional assessment)
- Include transcript accessions for both genes
- Note whether fusion is detectable by FISH, RT-PCR, or NGS only

## Report Structure

### Recommended Section Order

For clear clinical utility and standardized format, organize reports as follows:

1. **Header/Title Page**
   - Laboratory name and accreditation (CLIA, ISO 15189, etc.)
   - Report date and result date
   - Report ID (unique identifier)

2. **Patient and Specimen Information**
   - Patient demographics
   - Specimen details (type, site, collection date)
   - Clinical history/diagnosis

3. **Test Description**
   - Assay name and method
   - Genes analyzed
   - Performance characteristics

4. **Quality Control Metrics**
   - Sample-level QC results (coverage, contamination, etc.)
   - Pass/Fail status
   - Genes with insufficient coverage (explicitly listed)

5. **Findings**
   - Pathogenic/likely pathogenic variants
   - Variants of uncertain significance (if reported)
   - Copy number variations
   - Fusions/translocations
   - TMB, MSI, HRD status (if applicable)

6. **Clinical Interpretation**
   - ESCAT tier assignment for each actionable variant
   - Associated therapies and evidence
   - Clinical recommendations
   - References to supporting literature

7. **Limitations**
   - Assay limitations (coverage gaps, regions not analyzed)
   - Interpretation limitations (VUS, emerging variants)
   - Recommendation for additional testing if needed

8. **References**
   - Cited literature (PubMed IDs, DOIs)
   - ESMO guidelines reference
   - Links to key databases (ClinVar, OncoKB, etc.)

9. **Pathologist/Bioinformatician Sign-off**
   - Signature block
   - Date
   - Credentials

## Quality Metrics and Thresholds

### Minimum Coverage Requirements (by assay type)

| Assay Type | Minimum Mean Depth | Minimum % Target ≥30x |
|------------|-------------------|----------------------|
| Whole Genome | 30x | 95% |
| Whole Exome | 100x | 95% |
| Targeted Panel | 500x | 98% |

### Variant Calling Quality Thresholds

- **Variant Quality Score**: ≥20 (Phred scale minimum)
- **Variant Allele Fraction (VAF)**: ≥5% somatic, ≥20% germline (assay-dependent)
- **Read Depth**: Minimum 10 supporting reads for variants <5% VAF
- **Base Quality**: Minimum QUAL score 20 (1% error rate)

### Annotation Database Requirements

Variants must be cross-referenced with:
- **ClinVar**: Clinical significance classification
- **gnomAD v4.0+**: Population allele frequency
- **COSMIC**: Cancer-specific variant frequency
- **OncoKB**: Cancer-specific actionability
- **HPD/OpenTargets**: Therapeutic targets

## Additional Considerations

### Tumor Purity Assessment

- Estimate tumor cell percentage in sample
- Adjust VAF interpretation based on purity
- Note if sample contains significant normal cell contamination

### Clonal vs. Subclonal Variants

- Variants at VAF >30%: likely clonal (present in most tumor cells)
- Variants at VAF 5-30%: subclonal (present in subset of tumor cells)
- Note evolutionary implications for Tier III+ variants

### Incidental Findings

- Report pathogenic germline variants (per ESMO guidance on BRCA1/2)
- Offer genetic counseling referral
- Distinguish clearly from somatic findings in report

## Reference and Further Reading

**Primary Reference**:
- Mateo et al. "Recommendations for NGS-based diagnostics in oncology: standardization of variant classification and reporting" *Annals of Oncology*. 2024.
  https://www.annalsofoncology.org/article/S0923-7534(24)00111-X/fulltext

**ESCAT Scale**:
- Mateo et al. "ESCAT: A standardised framework for clinical actionability of molecular targets in cancer." *Nature Medicine*. 2023.

**ACMG Variant Classification**:
- Richards et al. "Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology." *Genetics in Medicine*. 2015.

**Other Key Guidelines**:
- NCCN Precision Medicine & Personalized Treatment Panels
- CAP Guidelines for Molecular Testing
- ISO 15189 Medical Laboratory Accreditation

**Clinical Databases**:
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
- OncoKB: https://www.oncokb.org/
- COSMIC: https://cancer.sanger.ac.uk/cosmic/
- gnomAD: https://gnomad.broadinstitute.org/
