# NGS Tertiary Analysis Pipeline Setup Guide

This guide covers all steps to set up the NGS tertiary analysis pipeline for variant annotation, clinical interpretation, and report generation.

## Prerequisites

Before starting, ensure you have:

- **R >= 4.4**: [Download from CRAN](https://cran.r-project.org)
- **Operating System**: macOS or Linux
- **Package Manager**: Homebrew (macOS) or equivalent Linux package manager
- **Disk Space**: ~50 GB (for reference genomes, tool installations, and cache files)
- **Internet Connection**: For downloading dependencies and API access

## Quick Setup

For the fastest path forward, run:

```bash
make setup-all
```

This single command handles system dependencies, reference data, R packages, and API configuration. Skip to [Running the Pipeline](#running-the-pipeline) if this succeeds.

If you prefer more control or encounter issues, follow the [Manual Setup](#manual-setup) section below.

## Manual Setup

### Step 1: Install System Tools

#### macOS (via Homebrew)

```bash
brew install samtools bcftools htslib
```

#### Linux (Ubuntu/Debian)

```bash
sudo apt-get install -y samtools bcftools libhts-dev
```

#### Linux (CentOS/RHEL)

```bash
sudo yum install -y samtools bcftools htslib-devel
```

Verify installation:

```bash
samtools --version
bcftools --version
```

### Step 2: Install GATK 4.6.1.0

GATK must be installed separately for variant calling and filtration:

```bash
mkdir -p ~/tools
cd ~/tools

# Download from GitHub releases
wget https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip
unzip gatk-4.6.1.0.zip
rm gatk-4.6.1.0.zip

# Add to PATH
export PATH="$HOME/tools/gatk-4.6.1.0:$PATH"

# Verify installation
gatk --version
```

Add the PATH export to your shell configuration (~/.bashrc, ~/.zshrc, or ~/.profile):

```bash
echo 'export PATH="$HOME/tools/gatk-4.6.1.0:$PATH"' >> ~/.zshrc
```

### Step 3: Install Ensembl VEP and Download Cache

VEP provides functional variant annotation:

```bash
# Clone VEP repository
mkdir -p ~/tools
cd ~/tools
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep

# Install Perl dependencies and VEP
perl INSTALL.pl --NO_TEST --NO_BIOPERL --CACHEDIR ~/vep_cache

# This takes 5-10 minutes. Follow any prompts to install required Perl modules.
```

Download GRCh38 cache (required for annotations):

```bash
mkdir -p ~/vep_cache
cd ~/vep_cache

# Download the cache file (~10 GB)
wget https://ftp.ensembl.org/pub/release-111/variation/vep/homo_sapiens_vep_111_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_111_GRCh38.tar.gz
rm homo_sapiens_vep_111_GRCh38.tar.gz
```

Add VEP to PATH:

```bash
export PATH="$HOME/tools/ensembl-vep:$PATH"
echo 'export PATH="$HOME/tools/ensembl-vep:$PATH"' >> ~/.zshrc

vep --version  # Verify
```

### Step 4: Install CNVkit

CNVkit is used for copy number variation analysis:

```bash
pip install cnvkit
```

Or if using a virtual environment:

```bash
uv pip install cnvkit
```

Verify:

```bash
cnvkit.py --version
```

### Step 5: Install/Verify Quarto CLI

Quarto is required for report generation:

#### macOS

```bash
brew install quarto
```

#### Linux

```bash
# Download from https://quarto.org/docs/get-started/
# Or use your package manager if available

# For Ubuntu:
sudo apt-get install -y quarto
```

Verify:

```bash
quarto --version
```

### Step 6: Initialize R Environment

Install R package dependencies using renv:

```bash
cd /Users/htlin/ngs-tertiary-analysis-skills

# Restore packages from lockfile
R --vanilla << 'EOF'
renv::restore()
EOF
```

This installs all required packages listed in `renv.lock`. This may take 10-15 minutes for first-time installation.

### Step 7: Download Reference Genome

Download GRCh38 reference sequence:

```bash
mkdir -p ~/references
cd ~/references

# Download GRCh38 from NCBI (full reference ~3.2 GB compressed)
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38.p14/seqs_for_alignment_pipelines.ucsc_ids/GRCh38_full_analysis_set.refasta.gz

# Decompress
gunzip GRCh38_full_analysis_set.refasta.gz
mv GRCh38_full_analysis_set.refasta GRCh38.fa

# Create index files (required by samtools and GATK)
samtools faidx GRCh38.fa
```

Create BWA index (if performing read alignment):

```bash
bwa index GRCh38.fa
```

Note: This takes ~2 hours. You can also use pre-built indices from Broad Institute.

## Docker Setup

For containerized deployment:

### Build Docker Image

```bash
cd /Users/htlin/ngs-tertiary-analysis-skills
docker build -t ngs-tertiary-analysis:latest .
```

### Run Pipeline in Docker

```bash
docker run --rm -v /path/to/data:/data ngs-tertiary-analysis:latest \
  Rscript /pipeline/main.R --input /data/sample.bam --output /data/report.html
```

### Run Interactive Session

```bash
docker run --rm -it -v /path/to/data:/data ngs-tertiary-analysis:latest bash
```

## API Keys Configuration

Several optional features require API keys. Create a `.env` file in the project root:

```bash
cd /Users/htlin/ngs-tertiary-analysis-skills
touch .env
```

### OncoKB API Key

1. Register at https://www.oncokb.org/account/register
2. Go to Account Settings → API Access
3. Copy your API token
4. Add to `.env`:

```
ONCOKB_API_KEY=your_token_here
```

### PubMed API Key (NCBI)

1. Register for NCBI account at https://www.ncbi.nlm.nih.gov/account/
2. Go to Settings → API Keys
3. Create new API key
4. Add to `.env`:

```
NCBI_API_KEY=your_key_here
```

### Scopus API Key

1. Register at https://dev.elsevier.com/
2. Create new API key in Developer Portal
3. Add to `.env`:

```
SCOPUS_API_KEY=your_key_here
```

### Unpaywall API Key

1. Register at https://unpaywall.org/products/api
2. Receive API token via email
3. Add to `.env`:

```
UNPAYWALL_API_KEY=your_email@example.com
```

### Embase API Access

Contact Elsevier for institutional access. Add credentials:

```
EMBASE_USERNAME=your_username
EMBASE_PASSWORD=your_password
```

**Note**: API keys are optional. The pipeline will skip annotation steps if keys are unavailable.

## Test Data

Download test BAM files for pipeline validation:

```bash
cd /Users/htlin/ngs-tertiary-analysis-skills
bash testdata/download_testdata.sh
```

This script downloads:
- Sample BAMs (exome capture sequencing)
- Expected variant calls for comparison
- Reference files (~2 GB total)

Files are stored in `testdata/bams/` directory.

## Running the Pipeline

### Basic Usage

Analyze a single sample:

```bash
Rscript main.R \
  --input path/to/sample.bam \
  --sample-id SAMPLE001 \
  --output results/SAMPLE001_report.html
```

### Advanced Usage

With custom parameters:

```bash
Rscript main.R \
  --input path/to/sample.bam \
  --sample-id SAMPLE001 \
  --output results/SAMPLE001_report.html \
  --min-depth 100 \
  --min-vaf 0.05 \
  --oncokb-key $ONCOKB_API_KEY \
  --threads 8
```

### Batch Processing

Process multiple samples:

```bash
for bam in data/samples/*.bam; do
  sample_id=$(basename "$bam" .bam)
  Rscript main.R --input "$bam" --sample-id "$sample_id" --output "results/${sample_id}_report.html"
done
```

### Pipeline Outputs

For each sample, the pipeline generates:

- `{SAMPLE_ID}_variants.vcf` - Raw variant calls
- `{SAMPLE_ID}_annotated.vcf` - Annotated variants
- `{SAMPLE_ID}_report.html` - Clinical report with ESMO guidelines compliance
- `{SAMPLE_ID}_qc_metrics.json` - Quality control metrics
- `{SAMPLE_ID}_cnv_calls.tsv` - Copy number variations

## Troubleshooting

### Issue: GATK not found

**Solution**: Verify GATK is in PATH and properly installed:

```bash
which gatk
# Should return /home/user/tools/gatk-4.6.1.0/gatk

# If not found, add to PATH manually
export PATH="$HOME/tools/gatk-4.6.1.0:$PATH"
```

### Issue: VEP cache not found

**Solution**: Ensure cache directory is accessible and set VEP_CACHEDIR:

```bash
export VEP_CACHEDIR="$HOME/vep_cache"
vep --cache --help | grep cache  # Verify VEP finds cache
```

### Issue: R packages fail to install

**Solution**: Update R and renv, then restore:

```bash
R --vanilla << 'EOF'
update.packages()
renv::restore(prompt=FALSE)
EOF
```

For persistent issues, check system dependencies:

```bash
# macOS
brew list | grep -E "samtools|bcftools|htslib"

# Linux
apt list --installed | grep -E "samtools|bcftools|libhts"
```

### Issue: Insufficient disk space

**Solution**: Check available space and relocate reference data:

```bash
df -h ~/references  # Check free space

# If low on space, use symbolic links to external drive:
mv ~/references /mnt/external/references
ln -s /mnt/external/references ~/references
```

### Issue: BAM file is corrupt or unindexed

**Solution**: Reindex the BAM file:

```bash
samtools index -b sample.bam sample.bam.bai
```

If BAM is corrupt, try to recover:

```bash
samtools quickcheck sample.bam  # Test file integrity
```

### Issue: Out of memory during variant calling

**Solution**: Reduce thread count or increase memory:

```bash
# Reduce threads
Rscript main.R --threads 2

# Or increase available memory in Docker/environment
export MALLOC_TRIM_THRESHOLD_=131072
```

## Reference Genome Notes

### Which Reference to Use?

- **GRCh38** (current, recommended): Human reference genome build 38
  - Used for modern NGS pipelines
  - Required for clinical reporting with ESMO guidelines
  - VEP cache available at release 111+

- **GRCh37** (deprecated): Older build, not recommended
  - Use only if legacy data requires it

### Pre-built Indices

To save time, download pre-built indices:

```bash
# From Broad Institute (recommended)
cd ~/references
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.gz
gunzip Homo_sapiens_assembly38.fasta.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
```

### Updating Reference Data

To update to a newer GRCh38 patch:

```bash
# Download new version
cd ~/references
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38.pXX/...

# Reindex
samtools faidx new_reference.fa
bwa index new_reference.fa
```

Update pipeline config to point to new reference path.

## Next Steps

After setup completes:

1. Run `bash testdata/download_testdata.sh` to get test data
2. Execute a test analysis with test data
3. Review the generated report at `results/test_sample_report.html`
4. Check [ESMO_GUIDELINES.md](ESMO_GUIDELINES.md) for reporting standards
5. Configure API keys in `.env` for enhanced annotations

For detailed information about clinical report standards, see [ESMO_GUIDELINES.md](ESMO_GUIDELINES.md).
