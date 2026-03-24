# NGS Tertiary Analysis - Reference Guide

## OncoKB API Reference

### Authentication
- **Base URL**: https://www.oncokb.org/api/v1
- **Auth Method**: Bearer token in `Authorization` header
- **Header Format**: `Authorization: Bearer {ONCOKB_API_KEY}`
- **Rate Limits**: 10 requests per second per API key
- **Response Format**: JSON

### Core Endpoints

#### Mutation Annotation
```
GET /annotate/mutations/byProteinChange
  ?hugoSymbol={gene}
  &alteration={change}
  &tumorType={type}
  &consequence={consequence}
```

**Parameters**:
- `hugoSymbol` (required): HUGO gene symbol (e.g., "BRAF", "EGFR", "TP53")
- `alteration` (required): Protein change in HGVS format (e.g., "V600E", "L858R", "deletion")
- `tumorType` (required): OncoKB tumor type (e.g., "NSCLC", "COAD", "Breast Cancer")
- `consequence` (optional): VEP consequence term (e.g., "missense_variant", "frameshift_variant")

**Response Fields**:
```json
{
  "oncogenic": "Likely Oncogenic|Likely Loss-of-function|Inconclusive",
  "mutationEffect": "Activating|Loss-of-function|Switch-of-function|Neutral",
  "highestSensitiveLevel": "LEVEL_1|LEVEL_2|LEVEL_3A|LEVEL_3B|LEVEL_4|LEVEL_R1",
  "highestResistanceLevel": "LEVEL_R1|LEVEL_R2|null",
  "treatments": [
    {
      "drugs": [{"drugName": "String"}],
      "level": "LEVEL_1|LEVEL_2|...",
      "indication": "Indication for the drug",
      "isStandard": true|false
    }
  ],
  "citations": {
    "abstracts": ["String"],
    "pmids": [12345]
  }
}
```

#### Copy Number Alteration Annotation
```
GET /annotate/copyNumberAlterations
  ?hugoSymbol={gene}
  &copyNameAlterationType={type}
  &tumorType={tumorType}
```

**Parameters**:
- `hugoSymbol` (required): Gene symbol
- `copyNameAlterationType` (required): "AMPLIFICATION" or "DELETION"
- `tumorType` (required): OncoKB tumor type

**Response**: Same structure as mutation annotation above

#### Structural Variant Annotation
```
GET /annotate/structuralVariants
  ?hugoSymbolA={geneA}
  &hugoSymbolB={geneB}
  &structuralVariantType={type}
  &tumorType={tumorType}
```

**Parameters**:
- `hugoSymbolA`, `hugoSymbolB`: Genes involved in fusion
- `structuralVariantType`: "FUSION", "TRANSLOCATION", "INVERSION", "DELETION"
- `tumorType`: OncoKB tumor type

**Response**: Same structure as mutation annotation

### OncoKB Tumor Types (Selection)
Common tumor types for CGP panels:
- `NSCLC` - Non-small cell lung cancer
- `Small Cell Lung Cancer` - Small cell lung cancer
- `COAD` - Colorectal adenocarcinoma
- `Breast Cancer`
- `Melanoma`
- `Ovarian Cancer`
- `Pancreatic Cancer`
- `Prostate Cancer`
- `Renal Cell Carcinoma`
- `Thyroid Cancer`
- `Endometrial Cancer`

Full list: https://www.oncokb.org/api/v1/tumorTypes

### OncoKB Evidence Levels
**Sensitive Levels** (treatment response):
- **LEVEL_1**: FDA-approved drug for the variant in this cancer type
- **LEVEL_2**: FDA-approved drug for the variant in other cancer type(s) OR well-powered clinical trial with strong results
- **LEVEL_3A**: Compelling clinical evidence from case reports or small clinical trials
- **LEVEL_3B**: Compelling biological evidence, including well-characterized DBs or case reports in other tumor types
- **LEVEL_4**: Biological evidence only

**Resistance Levels**:
- **LEVEL_R1**: FDA-approved drug resistance mutation in this cancer type
- **LEVEL_R2**: Compelling clinical evidence of resistance from case reports or small trials

---

## ESMO 2024 Report Requirements

The clinical report must include the following 10 sections per ESMO precision oncology guidelines:

### Section 1: Patient & Sample Information (Level A - Required)
- Patient identifier (de-identified)
- Age at diagnosis
- Sex
- Tumor type and primary site
- Sample type (tissue, blood, other)
- Specimen collection date
- Sample receipt date
- DNA/RNA quantity and quality metrics

### Section 2: Test Summary (Level A - Required)
- Test name and version (e.g., "TSO500 v1.2")
- Genomic regions covered (coding size in Mb)
- Depth of sequencing (target and achieved coverage)
- Tumor purity/cellularity estimate
- QC pass/fail status
- Data analysis pipeline version
- Reference genome build (GRCh38)

### Section 3: QC Metrics (Level A - Required)
- Mean coverage across sequencing
- Coverage uniformity (% bases at target coverage)
- Aligned bases
- On-target rate
- Contamination assessment if applicable
- Purity estimate method and value
- Interpretation of QC results (pass/fail thresholds)

### Section 4: Somatic Variants (Level A - Required)
For each reportable variant:
- Gene name and HGNC ID
- Genomic coordinate (GRCh38) in HGVS notation
- cDNA change (e.g., c.1799T>A)
- Protein change (e.g., p.V600E)
- Variant allele frequency (VAF) with confidence interval
- Type (SNV, indel, etc.)
- Zygosity
- Functional impact (missense, frameshift, stop-gain, etc.)
- ClinVar classification and ID
- COSMIC ID and frequency
- dbSNP ID if applicable
- Population frequency (gnomAD)

### Section 5: Copy Number Alterations (Level A - Required)
- Gene symbol
- Alteration type (amplification, deletion, LOH)
- Log2 ratio or copy number
- Genomic coordinates
- Gene dosage effect interpretation

### Section 6: Structural Variants & Fusions (Level A - Required)
- Gene partners (5' and 3')
- Fusion type and breakpoints
- Predicted fusion protein impact
- Frequency in databases (COSMIC, ChimerDB)
- Functional consequences

### Section 7: Biomarkers (Level A - Required)
- **TMB** (tumor mutational burden):
  - Count per Mb coding sequence
  - Interpretation (high/intermediate/low per tumor type)
  - Denominator used (panel size in Mb)
  - Relevance to immunotherapy

- **MSI Status**:
  - Predicted status (MSS/MSI-L/MSI-H) and method
  - Relevant for Lynch syndrome and immunotherapy

- **HRD Signature** (if applicable):
  - Score and status
  - Relevance to PARP inhibitor response

### Section 8: Clinical Actionability (Level A - Required)
For each tier, list:

**Tier I - Guideline-concordant therapy**:
- Variant details
- Linked drug(s) - FDA-approved for this indication
- Level of evidence (OncoKB LEVEL_1 or LEVEL_2)
- Clinical trial status
- Recommended action: strongly recommend inclusion in treatment plan

**Tier II - Investigational with strong evidence**:
- Variant details
- Linked drug(s) and clinical trials
- Evidence base (LEVEL_3A)
- Recommended action: consider for clinical trial enrollment

**Tier III - Other tumor types**:
- Variant details
- Linked drugs with efficacy in other malignancies
- Evidence base (LEVEL_3B)
- Recommended action: may have potential relevance outside indicated tumor type

**Tier IV - Preclinical/Early stage**:
- Variant details
- Biological relevance
- Status: research findings only, not actionable

**Tier X - No actionability**:
- Variants with uncertain pathogenicity or unknown clinical significance

### Section 9: Literature & Evidence (Level B - Recommended)
For each Tier I and Tier II variant:
- Up to 3 most recent publications (PubMed)
- Citation format: Authors, Year, Journal, PMID
- 1-2 sentence summary of key finding
- Open access indicator (with DOI link if available)
- Optional: Scopus citation count for high-impact studies

### Section 10: Clinical Recommendations & Interpretation (Level A - Required)
- Summary of key pathogenic variants
- Summary of key treatment opportunities (Tier I/II)
- Assessment of treatment resistance mutations
- Recommended next steps:
  - Molecular tumor boards review?
  - Refer to specific oncology specialist?
  - Enrollment in basket trials?
  - Repeat testing recommendations (timeframe)
- Limitations of the test
- Confirmatory testing recommendations if appropriate
- Disclaimer regarding interpretation validity

---

## ESCAT Tier Definitions

Classification system for precision oncology recommendations:

### Tier I: Standard of Care
**Criteria**:
- FDA-approved drug specifically for the variant in this cancer type
- OR sufficient clinical evidence in this indication with standard of care status
- **OncoKB Level**: LEVEL_1, LEVEL_2

**Reporting**:
- Strongly recommend inclusion in treatment plan
- Link to approved drugs and dosing
- Note standard of care status

### Tier II: Investigational, Strong Evidence
**Criteria**:
- Compelling clinical evidence (LEVEL_3A)
- Active clinical trials available
- Drug mechanism well-established for variant

**Reporting**:
- Recommend consideration for clinical trial enrollment
- List actively recruiting trials (NCT numbers)
- Include evidence summary from key publications

### Tier III: Actionability in Other Tumor Types
**Criteria**:
- Drug approved for variant in different tumor type
- Strong biological relevance (LEVEL_3B)
- Limited direct evidence in index tumor type

**Reporting**:
- Note potential cross-tumor applicability
- Recommend discussion with tumor board
- Mark as "may have relevance outside indicated use"

### Tier IV: Preclinical or Early-Stage Evidence
**Criteria**:
- Biological/mechanistic evidence only
- No clinical trial data yet
- OncoKB LEVEL_4

**Reporting**:
- Label as research findings
- Not suitable for clinical decision-making
- Mention as potential future opportunity

### Tier X: No Actionability Evidence
**Criteria**:
- Variant is pathogenic but no actionable therapies exist
- OR variant has uncertain pathogenicity

**Reporting**:
- Note pathogenic status if known
- Recommend monitoring for future therapies
- Do not recommend any specific treatment based on this variant alone

---

## Configuration Schema

Full `config/default.yaml` reference:

### Sample Configuration
```yaml
sample:
  # Sample identifier (no spaces, patient-safe)
  sample_id: "PATIENT_001"

  # Tumor type for OncoKB (must match OncoKB controlled vocabulary)
  # Examples: "NSCLC", "COAD", "Breast Cancer", "Melanoma"
  tumor_type: "NSCLC"

  # Specimen information
  specimen:
    collection_date: "2024-01-15"  # ISO format
    tissue_type: "FFPE"  # FFPE, fresh, blood, etc.
    purity: 0.5  # Tumor purity fraction (0-1)
```

### QC Configuration
```yaml
qc:
  # Minimum acceptable mean coverage
  min_mean_coverage: 200

  # Minimum percentage of bases at target coverage
  min_on_target_pct: 90

  # Minimum alignment percentage
  min_aligned_pct: 95

  # Contamination threshold (flag if > threshold)
  contamination_threshold: 0.05
```

### Variant Calling Configuration
```yaml
variant_calling:
  # Minimum variant allele frequency to report
  min_vaf: 0.05

  # Minimum number of alt reads
  min_alt_depth: 3

  # Minimum total depth at variant site
  min_total_depth: 50

  # Use matched normal for filtering (if available)
  use_normal: false
```

### Annotation Configuration
```yaml
annotation:
  # VEP version (e.g., "110")
  vep_version: "110"

  # Restrict to canonical transcripts only
  canonical_only: true

  # Include VEP plugins (list)
  plugins:
    - "Consequence"
    - "SIFT"
    - "PolyPhen"
    - "FATHMM"
```

### Biomarkers Configuration
```yaml
biomarkers:
  # TMB calculation
  tmb:
    # Panel coding size in Mb (for denominator)
    panel_coding_size_mb: 1.2

    # High TMB threshold (mutations per Mb)
    high_threshold: 10

    # Intermediate TMB threshold
    intermediate_threshold: 5

  # MSI prediction
  msi:
    # Prediction method: "manta" or "msi-sensor"
    method: "manta"

  # HRD signature (if applicable)
  hrd:
    enabled: true
    method: "signature"
```

### Clinical Annotation Configuration
```yaml
clinical:
  # OncoKB API key (set via environment variable)
  oncokb_api_key: ${ONCOKB_API_KEY}

  # PubMed E-utilities API key
  pubmed_api_key: ${PUBMED_API_KEY}

  # Scopus API key (Elsevier)
  scopus_api_key: ${SCOPUS_API_KEY}

  # Unpaywall email for open access lookup
  unpaywall_email: ${UNPAYWALL_EMAIL}

  # Embase API key (optional)
  embase_api_key: ${EMBASE_API_KEY}
```

### Report Configuration
```yaml
report:
  # Quarto output format
  format: "html"

  # Include literature evidence section
  include_literature: true

  # Maximum number of citations per variant
  max_citations_per_variant: 3

  # Report footer/institution info
  institution: "Clinical Genomics Laboratory"

  # Include disclaimer and limitations
  include_disclaimers: true
```

---

## API Response Examples

### OncoKB: BRAF V600E in Melanoma
```json
{
  "oncogenic": "Likely Oncogenic",
  "mutationEffect": "Activating",
  "highestSensitiveLevel": "LEVEL_1",
  "highestResistanceLevel": null,
  "treatments": [
    {
      "drugs": [
        {"drugName": "Vemurafenib"},
        {"drugName": "Dabrafenib"}
      ],
      "level": "LEVEL_1",
      "indication": "Melanoma",
      "isStandard": true
    },
    {
      "drugs": [
        {"drugName": "Trametinib"}
      ],
      "level": "LEVEL_1",
      "indication": "BRAF-mutant melanoma (combination therapy)",
      "isStandard": true
    }
  ],
  "citations": {
    "abstracts": ["..."],
    "pmids": [18776081, 19321294]
  }
}
```

### OncoKB: EGFR L858R in NSCLC
```json
{
  "oncogenic": "Likely Oncogenic",
  "mutationEffect": "Activating",
  "highestSensitiveLevel": "LEVEL_1",
  "highestResistanceLevel": null,
  "treatments": [
    {
      "drugs": [
        {"drugName": "Erlotinib"},
        {"drugName": "Gefitinib"},
        {"drugName": "Afatinib"}
      ],
      "level": "LEVEL_1",
      "indication": "EGFR-mutant NSCLC",
      "isStandard": true
    },
    {
      "drugs": [
        {"drugName": "Osimertinib"}
      ],
      "level": "LEVEL_1",
      "indication": "EGFR-mutant NSCLC (T790M or first-line)",
      "isStandard": true
    }
  ],
  "citations": {
    "pmids": [15235091, 19721078]
  }
}
```

### OncoKB: TP53 R248Q in Colorectal Cancer
```json
{
  "oncogenic": "Likely Oncogenic",
  "mutationEffect": "Loss-of-function",
  "highestSensitiveLevel": null,
  "highestResistanceLevel": null,
  "treatments": [],
  "citations": {
    "pmids": []
  }
}
```
Note: TP53 mutations are pathogenic but generally not directly actionable with approved drugs; associated with poor prognosis.
