# NGS Tertiary Analysis - Examples & Use Cases

## Example 1: Basic TSO500 Analysis Run

### Setup
```bash
# Clone repository and install dependencies
git clone <repo-url>
cd ngs-tertiary-analysis
make setup-all  # Install GATK, VEP, CNVkit, etc.

# Set up environment variables
cat > .env << 'EOF'
ONCOKB_API_KEY=your_api_key_here
PUBMED_API_KEY=your_ncbi_api_key
SCOPUS_API_KEY=your_scopus_key
UNPAYWALL_EMAIL=your_email@institution.edu
EOF

chmod 600 .env  # Protect API keys
```

### Minimal Configuration
Edit `config/default.yaml`:
```yaml
sample:
  sample_id: "TSO500_001"
  tumor_type: "NSCLC"
  specimen:
    collection_date: "2024-03-15"
    tissue_type: "FFPE"
    purity: 0.6
```

### Run Full Pipeline
```r
# In R console
Sys.setenv(BAM_PATH = "/path/to/tumor.bam")
Sys.setenv(SAMPLE_ID = "TSO500_001")

library(targets)
tar_make()
```

### Check Results
```bash
# View QC metrics
less reports/TSO500_001/00-qc/qc_summary.txt

# Open final report
open reports/TSO500_001/08-report/clinical_report.html
```

---

## Example 2: Customize for Lung Cancer (NSCLC)

### Configuration for Actionable Lung Cancer Genes
```yaml
sample:
  sample_id: "NSCLC_024"
  tumor_type: "NSCLC"
  specimen:
    collection_date: "2024-03-10"
    tissue_type: "Fresh tissue"
    purity: 0.75

qc:
  min_mean_coverage: 500  # Higher coverage for sensitive detection
  min_on_target_pct: 92

variant_calling:
  min_vaf: 0.02  # Detect lower VAF variants in NSCLC

biomarkers:
  tmb:
    panel_coding_size_mb: 1.2
    high_threshold: 10
    intermediate_threshold: 5
  msi:
    method: "manta"
```

### Expected Actionable Variants
Common NSCLC drivers detected on TSO500:
- **EGFR**: L858R, exon 19 deletions, T790M (resistance)
- **ALK**: Fusions (EML4-ALK, TFG-ALK)
- **ROS1**: Fusions
- **BRAF**: V600E
- **KRAS**: G12C (newly actionable)
- **MET**: Exon 14 skipping, amplifications
- **PD-L1**: Expression level (if available)

### Output Interpretation
```
Tier I (Standard of care):
- EGFR L858R → Gefitinib, Erlotinib, Afatinib (Level 1)
- ALK fusion → Crizotinib, Alectinib, Brigatinib (Level 1)

Tier II (Clinical trials):
- KRAS G12C → Sotorasib, Adagrasib (Level 1)
- MET exon 14 → Tepotinib (Level 1)

Tier III (Other cancers):
- BRAF V600E → Vemurafenib (approved in melanoma)
```

---

## Example 3: Colorectal Cancer (CRC) Pipeline

### CRC-Specific Configuration
```yaml
sample:
  sample_id: "CRC_2024_045"
  tumor_type: "COAD"  # Colorectal adenocarcinoma
  specimen:
    collection_date: "2024-02-28"
    tissue_type: "FFPE"
    purity: 0.5

variant_calling:
  min_vaf: 0.03

biomarkers:
  tmb:
    panel_coding_size_mb: 1.2
    high_threshold: 10  # TMB >10 → immunotherapy candidate
  msi:
    method: "manta"  # Important for CRC!
```

### CRC Key Biomarkers & Interpretation

**MSI Status** (critical for CRC):
```
If MSI-H (high instability):
  → Suitable for checkpoint inhibitors (Pembrolizumab, Nivolumab)
  → Suggest Tier I recommendation
  → Consider Lynch syndrome screening

If MSS (microsatellite stable):
  → Tailor therapy to specific drivers (KRAS, BRAF, etc.)
  → MSI-H testing recommended if patient is young
```

**TMB in CRC**:
- TMB-H (>10 mutations/Mb) → Immunotherapy candidate (Tier II)
- Correlates with MSI status

### Common CRC Actionable Variants
- **KRAS**: G12C, G12V, Q61H (no standard drug yet for most)
- **BRAF**: V600E (EGFR inhibitors contraindicated, MEK inhibitor may help)
- **MMR genes**: MLH1, MSH2, MSH6, PMS2 (Lynch syndrome)
- **SMAD4**: Loss-of-function (prognostic)
- **TP53**: Mutations (common but not directly actionable)

### Sample Output Section
```markdown
## Biomarkers & Tumor Characteristics

**Microsatellite Instability (MSI)**: Predicted MSI-H (128 indels in microsatellite regions)
- Clinical relevance: High immunogenicity; suitable for checkpoint inhibitors
- Actionability: Tier I consideration for Pembrolizumab or Nivolumab

**Tumor Mutational Burden (TMB)**: 8.3 mutations/Mb (intermediate)
- Threshold for immunotherapy consideration: >10 mutations/Mb
- Status: Borderline; review with tumor board

**Lynch Syndrome Screening**:
- MLH1: No pathogenic mutations detected
- MSH2: No pathogenic mutations detected
- Interpretation: No inherited mismatch repair deficiency identified
```

---

## Example 4: Breast Cancer with HRD Assessment

### Configuration
```yaml
sample:
  sample_id: "BREAST_2024_103"
  tumor_type: "Breast Cancer"
  specimen:
    collection_date: "2024-03-05"
    tissue_type: "FFPE"
    purity: 0.65

biomarkers:
  hrd:
    enabled: true
    method: "signature"  # Genomic scar analysis
```

### Key Breast Cancer Variants
- **BRCA1/BRCA2**: Germline and somatic (pathognomonic for HRD)
- **TP53**: Frequent in triple-negative breast cancer (TNBC)
- **PIK3CA**: Hormone receptor-positive breast cancer
- **ESR1**: Ligand-binding domain mutations (endocrine resistance)
- **ERBB2**: Amplifications (HER2-positive breast cancer)
- **PTEN**: Loss (PI3K pathway activation)

### HRD-Positive Output
```markdown
## Homologous Recombination Deficiency (HRD)

**HRD Status**: Positive (genomic scar score 68)
- Method: Large-scale genomic rearrangement detection
- Supporting variants: BRCA1 frameshift (8del), PTEN loss

**Clinical Significance**:
- Patient is PARP inhibitor-sensitive
- Consider Olaparib or Rucaparib as treatment option
- Strong consideration for platinum-based chemotherapy

**Tier I Recommendation**:
- PARP inhibitor (Olaparib, Rucaparib, Talazoparib) for HRD-positive breast cancer
- Level: LEVEL_1 evidence for BRCA-mutant or HRD-positive status
```

---

## Example 5: Running Specific Pipeline Stages

### Run Only Variant Calling
```r
# Skip QC and annotation, focus on calling
targets::tar_make(names = c("variant_calling_results"))

# Output: reports/{sample_id}/01-variant-calling/
#   ├── somatic.vcf.gz
#   ├── somatic.vcf.gz.tbi
#   └── calling_metrics.txt
```

### Run Only OncoKB Annotation
```r
# Assumes VCF already exists
targets::tar_make(names = c("oncokb_results"))

# Output: reports/{sample_id}/06-clinical-annotation/
#   ├── variants_oncokb.tsv
#   └── oncokb_summary.txt
```

### Run Only Literature Review
```r
# Assumes OncoKB results exist
targets::tar_make(names = c("literature_results"))

# Output: reports/{sample_id}/07-literature/
#   ├── pubmed_evidence.tsv
#   ├── scopus_citations.tsv
#   └── open_access_links.txt
```

---

## Example 6: OncoKB API Response for BRAF V600E

### Request
```bash
curl -H "Authorization: Bearer ${ONCOKB_API_KEY}" \
  "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol=BRAF&alteration=V600E&tumorType=Melanoma"
```

### Response (JSON)
```json
{
  "query": {
    "hugoSymbol": "BRAF",
    "alteration": "V600E",
    "tumorType": "Melanoma"
  },
  "oncogenic": "Likely Oncogenic",
  "mutationEffect": "Activating",
  "mutationEffectPmids": [12829827, 15235091],
  "highestSensitiveLevel": "LEVEL_1",
  "highestResistanceLevel": "LEVEL_R1",
  "highestDiagnosticImplicationLevel": null,
  "highestPrognosticImplicationLevel": "LEVEL_2",
  "treatments": [
    {
      "drugs": [
        {
          "drugName": "Vemurafenib",
          "ncitId": "C100467",
          "synonyms": ["PLX4032"]
        }
      ],
      "approvalStatus": "FDA_APPROVED",
      "level": "LEVEL_1",
      "indication": "Melanoma",
      "pmids": [18776081, 19321294],
      "abstracts": [
        "Randomized phase III trial demonstrating vemurafenib efficacy in BRAF-mutant melanoma..."
      ],
      "isStandard": true
    },
    {
      "drugs": [
        {
          "drugName": "Dabrafenib",
          "ncitId": "C68809"
        }
      ],
      "approvalStatus": "FDA_APPROVED",
      "level": "LEVEL_1",
      "indication": "Melanoma",
      "pmids": [21639808],
      "isStandard": true
    },
    {
      "drugs": [
        {
          "drugName": "Trametinib",
          "ncitId": "C77529"
        }
      ],
      "approvalStatus": "FDA_APPROVED",
      "level": "LEVEL_1",
      "indication": "BRAF-mutant melanoma (combination with BRAF inhibitor)",
      "pmids": [20534541],
      "isStandard": true
    },
    {
      "drugs": [
        {
          "drugName": "Cobimetinib",
          "ncitId": "C85550"
        }
      ],
      "approvalStatus": "FDA_APPROVED",
      "level": "LEVEL_1",
      "indication": "BRAF-mutant melanoma (combination with BRAF inhibitor)",
      "pmids": [20534541],
      "isStandard": true
    }
  ],
  "resistanceMutations": [
    {
      "drugs": [
        {
          "drugName": "Vemurafenib"
        }
      ],
      "level": "LEVEL_R1",
      "pmids": [22472876],
      "alterations": ["T1799I", "M1849V", "A1846T"]
    }
  ],
  "citations": {
    "pmids": [18776081, 19321294, 21639808],
    "abstracts": []
  }
}
```

### How to Interpret in Report
```markdown
### BRAF V600E - Melanoma

**Classification**: Likely Oncogenic | Activating mutation

**Tier**: I - Standard of Care

**FDA-Approved Therapies**:
1. Vemurafenib (Zelboraf)
   - Single-agent or in combination with trametinib
   - Evidence: LEVEL_1 (randomized phase III trial)
   - PMID: 18776081, 19321294

2. Dabrafenib (Tafinlar)
   - Often combined with trametinib for enhanced response
   - Evidence: LEVEL_1 (randomized phase III trial)
   - PMID: 21639808

3. Combination Therapy (BRAF + MEK inhibitor)
   - Vemurafenib + Trametinib
   - Dabrafenib + Trametinib
   - Dabrafenib + Cobimetinib
   - Evidence: LEVEL_1 (improved survival)

**Resistance Mutations to Monitor**:
- T1799I, M1849V, A1846T (secondary BRAF mutations)
- Monitor for emerging resistance; consider repeat testing if disease progresses

**Recommendation**: Strongly recommend BRAF inhibitor-based therapy as first-line treatment.
Consider combination with MEK inhibitor for improved response durability.
```

---

## Example 7: Full Report HTML Structure

### Generated Report Outline
```html
<!DOCTYPE html>
<html>
<head>
  <title>Clinical Genomics Report - PATIENT_001</title>
  <style>/* ESMO-compliant styling */</style>
</head>
<body>

<!-- Section 1: Header & Patient Info -->
<div class="report-header">
  <h1>Clinical Genomics Report</h1>
  <p>Patient ID: PATIENT_001 | DOB: [redacted] | Age: 62</p>
  <p>Tumor Type: NSCLC | Specimen: FFPE tumor tissue | Collection: 2024-03-15</p>
</div>

<!-- Section 2: Test Summary -->
<div class="test-summary">
  <h2>Test Overview</h2>
  <p>Test: TSO500 v1.2 | Coverage: 1247x mean (target: 500x)</p>
  <p>On-target: 94.2% | Tumor Purity: 60% | QC Status: PASS</p>
</div>

<!-- Section 3: QC Metrics -->
<div class="qc-metrics">
  <h2>Quality Control Metrics</h2>
  <table>
    <tr><th>Metric</th><th>Value</th><th>Threshold</th><th>Status</th></tr>
    <tr><td>Mean Coverage</td><td>1247x</td><td>>200x</td><td>PASS</td></tr>
    <tr><td>On-target %</td><td>94.2%</td><td>>90%</td><td>PASS</td></tr>
    <tr><td>Alignment %</td><td>98.5%</td><td>>95%</td><td>PASS</td></tr>
  </table>
</div>

<!-- Section 4-7: Variants, CNV, Fusions, Biomarkers -->
<div class="somatic-variants">
  <h2>Somatic Variants</h2>
  <!-- Table with variant details -->
</div>

<div class="cnv">
  <h2>Copy Number Alterations</h2>
  <!-- CNV details -->
</div>

<div class="biomarkers">
  <h2>Biomarkers</h2>
  <p><strong>TMB</strong>: 7.4 mutations/Mb (intermediate)</p>
  <p><strong>MSI</strong>: Predicted MSS</p>
</div>

<!-- Section 8: Clinical Actionability (Core Section) -->
<div class="clinical-actionability">
  <h2>Clinical Actionability & Treatment Recommendations</h2>

  <h3>Tier I - Guideline-Concordant Therapy</h3>
  <div class="tier-i">
    <h4>EGFR L858R (c.2573T>A)</h4>
    <p><strong>OncoKB Classification</strong>: Likely Oncogenic</p>
    <p><strong>VAF</strong>: 28.5% (95 alt reads / 333 total)</p>
    <p><strong>Evidence Level</strong>: LEVEL_1</p>
    <p><strong>Recommended Drugs</strong>:
      <ul>
        <li>Gefitinib (Iressa) - FDA-approved for EGFR-mutant NSCLC</li>
        <li>Erlotinib (Tarceva) - Alternative first-line EGFR TKI</li>
        <li>Afatinib (Gilotrif) - Irreversible EGFR inhibitor</li>
      </ul>
    </p>
    <p><strong>Clinical Action</strong>: Strongly recommend EGFR-targeted therapy as first-line treatment.</p>
  </div>

  <h3>Tier II - Investigational, Strong Evidence</h3>
  <div class="tier-ii">
    <h4>ALK Fusion (EML4-ALK)</h4>
    <p>Active clinical trials available (NCT02968147, NCT03181620)</p>
  </div>

  <h3>Tier X - No Actionability</h3>
  <div class="tier-x">
    <h4>TP53 R248Q (somatic)</h4>
    <p>Pathogenic but no approved targeted therapy available.</p>
  </div>
</div>

<!-- Section 9: Literature Evidence -->
<div class="literature">
  <h2>Literature Evidence</h2>
  <h3>EGFR L858R in NSCLC</h3>
  <ol>
    <li>Lynch et al. (2004). "Activating mutations in the epidermal growth factor receptor
        underlying responsiveness of non-small-cell lung cancer to gefitinib."
        <em>N Engl J Med</em>. PMID: 15235091
        [<a href="https://doi.org/10.1056/NEJMoa040938">DOI</a>]
    </li>
  </ol>
</div>

<!-- Section 10: Interpretation & Recommendations -->
<div class="interpretation">
  <h2>Clinical Interpretation & Recommendations</h2>
  <p><strong>Summary</strong>: This 62-year-old patient with NSCLC carries activating EGFR
     L858R mutation with strong evidence for targeted therapy response.</p>
  <p><strong>Next Steps</strong>:
    <ul>
      <li>Initiate EGFR-targeted tyrosine kinase inhibitor (TKI) therapy</li>
      <li>Monitor for acquired T790M resistance mutation at progression</li>
      <li>Consider osimertinib for second-line if resistance develops</li>
      <li>Repeat genomic testing 12-18 months or at disease progression</li>
    </ul>
  </p>
  <p><strong>Limitations</strong>: This report is based on somatic analysis only.
     Germline counseling recommended for any identified hereditary cancer variants.</p>
</div>

</body>
</html>
```

---

## Example 8: Troubleshooting Common Issues

### Issue: Low Coverage Warning
```
Error: Mean coverage 120x < threshold 200x
```

**Solution**:
```yaml
# Option A: Lower threshold if coverage is adequate for variant calls
qc:
  min_mean_coverage: 100

# Option B: Accept warning and proceed (inspect QC plot manually)
# Output: reports/{sample_id}/00-qc/coverage_plot.pdf
```

### Issue: OncoKB 401 Unauthorized
```
Error: HTTP 401 - OncoKB API authentication failed
```

**Solution**:
```bash
# Check API key is set correctly
echo $ONCOKB_API_KEY  # Should print your key (not empty)

# Verify account access at https://www.oncokb.org/apiAccess
# Ensure you have approved account status (may take 1-2 business days)

# Refresh .env file
source .env
```

### Issue: VEP Not Found
```
Error: vep command not found
```

**Solution**:
```bash
# Install VEP with cache
make setup-vep

# Or use Docker (preferred)
docker pull ensemblorg/ensembl-vep:latest
```

### Issue: Rate Limiting from PubMed
```
Warning: PubMed API rate limit exceeded, retrying...
```

**Solution**:
```r
# Increase wait time in R/api_clients.R
# Default: 0.3 sec between requests
# Increase to: 0.5-1.0 sec

# Or use batching with longer interval
Sys.sleep(2)
```

---

## Example 9: Multi-Sample Batch Analysis

### Run Multiple Samples
```r
# Create sample manifest
samples <- data.frame(
  sample_id = c("TSO500_001", "TSO500_002", "TSO500_003"),
  bam_path = c("path1.bam", "path2.bam", "path3.bam"),
  tumor_type = c("NSCLC", "COAD", "Breast Cancer")
)

# Run all samples
for (i in 1:nrow(samples)) {
  Sys.setenv(
    SAMPLE_ID = samples$sample_id[i],
    BAM_PATH = samples$bam_path[i]
  )
  targets::tar_make()
  cat("Sample", samples$sample_id[i], "complete\n")
}

# Aggregate results
results <- lapply(
  samples$sample_id,
  function(id) {
    read.csv(file.path("reports", id, "06-clinical-annotation", "variants_oncokb.tsv"))
  }
)
combined <- do.call(rbind, results)
```

### Generate Summary Report
```r
# Table: All Tier I variants across cohort
tier_i <- combined %>%
  filter(escat_tier == "I") %>%
  group_by(gene, hgvs_c) %>%
  summarize(
    n_samples = n_distinct(sample_id),
    drugs = paste(unique(na.omit(drug_name)), collapse = ", ")
  )
```

