.PHONY: setup-all setup-system setup-gatk setup-vep setup-cnvkit setup-r setup-quarto run test clean \
        help batch-pipeline batch-render batch-index batch-audit batch-deploy batch-all \
        oe-queue oe-citations clean-cache \
        manifest all-batches all-index all-deploy

# Setup all dependencies
setup-all: setup-system setup-gatk setup-vep setup-cnvkit setup-r setup-quarto
	@echo "All dependencies installed successfully"

# Install system dependencies using Homebrew
setup-system:
	@echo "Installing system dependencies via Homebrew..."
	brew install samtools bcftools htslib
	@echo "System dependencies installed"

# Download and setup GATK 4.6.1.0
setup-gatk:
	@echo "Setting up GATK 4.6.1.0..."
	mkdir -p tools
	cd tools && \
	if [ ! -f gatk-4.6.1.0.zip ]; then \
		curl -L -o gatk-4.6.1.0.zip https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip && \
		unzip -o gatk-4.6.1.0.zip && \
		rm gatk-4.6.1.0.zip; \
	fi
	@echo "GATK 4.6.1.0 installed to tools/"

# Setup Ensembl VEP
setup-vep:
	@echo "Installing Ensembl VEP..."
	brew install ensembl-vep || git clone https://github.com/Ensembl/ensembl-vep.git tools/vep
	@echo "VEP installed"

# Setup CNVkit via pip
setup-cnvkit:
	@echo "Installing CNVkit..."
	pip install cnvkit
	@echo "CNVkit installed"

# Setup R environment with renv
setup-r:
	@echo "Restoring R environment with renv..."
	Rscript -e 'renv::restore()'
	@echo "R environment restored"

# Check and setup Quarto CLI
setup-quarto:
	@echo "Checking Quarto installation..."
	if ! command -v quarto &> /dev/null; then \
		echo "Quarto not found, installing via Homebrew..."; \
		brew install quarto; \
	else \
		echo "Quarto is already installed"; \
	fi
	@echo "Quarto ready"

# Run the pipeline
run:
	@echo "Running targets pipeline..."
	Rscript -e 'targets::tar_make()'

# Run tests
test:
	@echo "Running tests..."
	Rscript -e 'testthat::test_dir("tests")'

# Clean targets cache
clean:
	@echo "Cleaning _targets/ directory..."
	rip _targets/
	@echo "Clean complete"

# ─────────────────────────────────────────────────────────────────────────────
#  Foundation One–style actionable reports (TSO500 batches)
# ─────────────────────────────────────────────────────────────────────────────
#
# Default flow for a new batch directory layout:
#   inputs/TSO500-HRD/<BATCH>/
#     <SID>_CombinedVariantOutput.tsv
#     <SID>_dna.hard-filtered.vcf
#     <SID>_dna_TMB_Trace.tsv
#     20260428_癌末TSO 500 HRD <BATCH>_重點摘要.xlsx     (curated mapping)
#
# Per-batch sample/tumor mapping currently lives inline in
# scripts/render_actionable_batch7.sh. To handle a new batch, copy that
# script and edit the SAMPLES array, then override BATCH below.
#
# Usage:
#   make batch-all                            # pipeline + render + index
#   make batch-deploy SITE=kfsyscc-tso500-demo
#   make batch-audit
#
# OpenEvidence MCP queries are run *interactively* from a Claude session;
# `make oe-queue` produces the deduplicated query list, then run one
# oe_ask per sample, save each answer to reports/.openevidence_cache/
# sample_<SID>.md, then `make oe-citations` to extract the AMA refs JSON.

BATCH       ?= Batch_7_20260505
BATCH_DIR   ?= inputs/TSO500-HRD/$(BATCH)
HANDOFF_DIR ?= reports/_handoff_$(BATCH)
SITE_NAME   ?= kfsyscc-tso500-demo
SITE_ID     ?= dbb8cb56-44fa-4dcb-899f-3fe65959d3c9

help: ## list batch targets
	@printf '\n  Foundation One TSO500 — batch make targets\n'
	@printf '  Override BATCH=<dir>, SITE_NAME=<...>, SITE_ID=<...> as needed.\n\n'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
	  | awk -F: '{ split($$0,a,"## "); printf "  \033[36m%-22s\033[0m %s\n", a[1], a[2] }'
	@echo

batch-pipeline: ## run TSO500 → OncoKB/CiVIC/AMP/literature for every sample
	bash scripts/render_actionable_batch7.sh

batch-render: ## re-render Foundation One HTMLs + PDFs from cached .rds
	bash scripts/render_actionable_batch7.sh

batch-index: ## rewrite reports/_handoff_*/index.html
	Rscript --vanilla scripts/build_handoff_index.R

batch-audit: ## reconcile reports against curated 重點摘要
	Rscript --vanilla scripts/audit_batch7.R

batch-all: batch-pipeline batch-index ## end-to-end: pipeline + render + index
	@echo "Done. Open: $(HANDOFF_DIR)/index.html"

oe-queue: ## emit deduplicated OpenEvidence query list
	Rscript --vanilla scripts/openevidence_queue.R

oe-citations: ## extract AMA reference JSONs from persisted OE responses
	bash scripts/extract_oe_citations.sh

batch-deploy: ## netlify deploy --prod (uses SITE_ID / SITE_NAME)
	@if ! command -v netlify >/dev/null 2>&1; then \
	  echo "ERROR: netlify CLI not installed (brew install netlify-cli)"; exit 1; \
	fi
	netlify deploy \
	  --site $(SITE_ID) \
	  --dir $(HANDOFF_DIR) \
	  --prod \
	  --message "Batch $(BATCH) — TSO500 actionable reports"

clean-cache: ## wipe OncoKB / OpenEvidence / API caches (asks before deleting)
	@echo "Will delete: reports/.oncokb_cache/ reports/.openevidence_cache/ .api_cache/"
	@read -p "Confirm? (y/N) " ans && [ "$$ans" = "y" ] && \
	  rip reports/.oncokb_cache reports/.openevidence_cache .api_cache && \
	  echo "Caches cleared." || echo "Aborted."
# ─────────────────────────────────────────────────────────────────────────────
#  Multi-batch operation (all batches under inputs/TSO500-HRD/)
# ─────────────────────────────────────────────────────────────────────────────

ALL_HANDOFF_DIR ?= reports/_handoff_all_batches

manifest: ## (re)build the master sample-→-tumor manifest by parsing every batch xlsx
	uvx --with openpyxl python3 scripts/build_batch_manifest.py

all-batches: manifest ## run pipeline + Foundation One render for every sample, skipping completed
	bash scripts/run_all_batches.sh

all-index: ## stage every produced report into _handoff_all_batches/ + grouped index
	Rscript --vanilla scripts/build_unified_index.R

all-deploy: ## netlify deploy --prod (all-batches handoff)
	@if ! command -v netlify >/dev/null 2>&1; then \
	  echo "ERROR: netlify CLI not installed (brew install netlify-cli)"; exit 1; \
	fi
	netlify deploy \
	  --site $(SITE_ID) \
	  --dir $(ALL_HANDOFF_DIR) \
	  --prod \
	  --message "All batches — TSO500 actionable reports"

# ─────────────────────────────────────────────────────────────────────────────
#  ESMO abstract submission helpers
# ─────────────────────────────────────────────────────────────────────────────

abstract-charcount: ## count ESMO abstract chars (limit 2000, excl spaces)
	@count=$$(awk '/^## Title/,/^## How to verify/' docs/abstract_esmo_2026.md \
	  | sed '/^---$$/d;/^## /d;/^$$/d' \
	  | tr -d ' \n' | wc -c); \
	  echo "Abstract chars (excl spaces): $$count / 2000"; \
	  if [ $$count -gt 2000 ]; then echo "OVER LIMIT"; exit 1; fi

abstract-metrics: ## print every number cited in the ESMO abstract from the bench CSVs
	Rscript --vanilla scripts/extract_abstract_metrics.R
