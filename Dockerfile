FROM bioconductor/bioconductor_docker:3.20

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    samtools \
    bcftools \
    tabix \
    curl \
    libxml2-dev \
    libcurl4-openssl-dev \
    git \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install GATK 4.6.1.0
RUN mkdir -p /app/tools && \
    cd /app/tools && \
    curl -L -o gatk-4.6.1.0.zip https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip && \
    unzip -o gatk-4.6.1.0.zip && \
    rm gatk-4.6.1.0.zip && \
    chmod +x gatk-4.6.1.0/gatk

# Install Ensembl VEP
RUN apt-get update && apt-get install -y \
    ensembl-vep \
    && rm -rf /var/lib/apt/lists/*

# Install CNVkit via pip
RUN pip3 install cnvkit

# Install Quarto CLI
RUN curl -L -o /tmp/quarto.deb https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.542/quarto-1.4.542-linux-amd64.deb && \
    dpkg -i /tmp/quarto.deb && \
    rm /tmp/quarto.deb

# Copy project files
COPY . /app/

# Restore R environment using renv
RUN Rscript -e 'renv::restore()'

# Default command runs the targets pipeline
CMD ["Rscript", "-e", "targets::tar_make()"]
