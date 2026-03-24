.PHONY: setup-all setup-system setup-gatk setup-vep setup-cnvkit setup-r setup-quarto run test clean

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
