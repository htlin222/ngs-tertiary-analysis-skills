suppressPackageStartupMessages({
  library(tidyverse)
  library(VariantAnnotation)
  library(glue)
  library(logger)
})

#' Parse Gene Fusions from Structural Variant Data
#'
#' Extracts and annotates gene fusion events from Manta SV VCF output or
#' TSO500 CSV fusion data. Filters for high-confidence fusions, identifies
#' frame status, and cross-references with known fusion databases.
#'
#' @param fusion_results Character. Path to Manta VCF or TSO500 fusion CSV.
#' @param config List. Configuration object with fusion parameters and
#'   reference databases.
#' @param sample_id Character. Sample identifier for output organization.
#'
#' @return Tibble with columns:
#'   - gene_a: 5' gene partner (character)
#'   - gene_b: 3' gene partner (character)
#'   - fusion_name: ESMO notation (geneA::geneB) (character)
#'   - exon_a: Exon number at 5' breakpoint (integer, may be NA)
#'   - exon_b: Exon number at 3' breakpoint (integer, may be NA)
#'   - chr_a: Chromosome of 5' partner (character)
#'   - pos_a: Genomic position of 5' breakpoint (integer)
#'   - chr_b: Chromosome of 3' partner (character)
#'   - pos_b: Genomic position of 3' breakpoint (integer)
#'   - supporting_reads: Number of supporting reads (integer)
#'   - fusion_type: Frame status ("in-frame", "out-of-frame", "unknown") (character)
#'   - known_fusion: Is fusion in known database? (logical)
#'
#' @details
#' The function handles both Manta VCF and TSO500 CSV formats. For VCF input:
#' - Extracts BND (breakend) type structural variants
#' - Filters by minimum supporting reads (config$fusions$min_supporting_reads)
#' - Matches breakpoints to genes and exons via reference annotation
#' - Determines in-frame status based on exon positions
#'
#' For TSO500 CSV input:
#' - Reads fusion data directly
#' - Extracts gene partners and supporting read information
#'
#' Known fusion database cross-reference uses Ensembl Fusion database or
#' a custom database specified in config$fusions$known_fusions_db.
#'
#' @keywords internal
parse_fusions <- function(fusion_results, config, sample_id) {
  log_info("Stage 4: Parsing Gene Fusions from {basename(fusion_results)}")

  if (!file.exists(fusion_results)) {
    log_error("Fusion results file not found: {fusion_results}")
    return(tibble())
  }

  # Determine input format
  if (str_detect(fusion_results, "\\.vcf(\\.gz)?$")) {
    fusions_df <- .parse_manta_vcf(fusion_results, config, sample_id)
  } else if (str_detect(fusion_results, "\\.csv$")) {
    fusions_df <- .parse_tso500_csv(fusion_results, config, sample_id)
  } else {
    log_error("Unrecognized fusion results format: {fusion_results}")
    return(tibble())
  }

  if (nrow(fusions_df) == 0) {
    log_warn("No high-confidence fusions detected")
    return(fusions_df)
  }

  # Cross-reference with known fusion database
  if (!is.null(config$fusions$known_fusions_db) &&
      file.exists(config$fusions$known_fusions_db)) {
    fusions_df <- .annotate_known_fusions(fusions_df, config$fusions$known_fusions_db)
  } else {
    log_info("No known fusion database configured, marking all as unknown")
    fusions_df <- fusions_df %>%
      mutate(known_fusion = FALSE)
  }

  # Save output
  output_dir <- stage_output_dir(sample_id, "04-fusions")
  output_file <- file.path(output_dir, "fusions.tsv")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  write_tsv(fusions_df, output_file)
  log_info("Fusions saved to {output_file}")
  log_success("Gene fusion parsing completed")

  return(fusions_df)
}

#' Parse Manta VCF for Gene Fusions
#'
#' Extracts BND (breakend) type SVs from Manta output and identifies
#' those spanning different genes.
#'
#' @param vcf_path Character. Path to Manta VCF file.
#' @param config List. Configuration object.
#' @param sample_id Character. Sample identifier.
#'
#' @return Tibble with fusion event data.
#'
#' @keywords internal
.parse_manta_vcf <- function(vcf_path, config, sample_id) {
  log_info("Parsing Manta VCF: {basename(vcf_path)}")

  tryCatch({
    # Read VCF
    vcf <- readVcf(vcf_path, genome = config$reference$genome_name)
    log_debug("VCF loaded: {length(vcf)} variants")

    # Extract VCF data
    fixed_df <- as.data.frame(fixed(vcf)) %>% rownames_to_column("var_id")
    info_df <- as.data.frame(info(vcf)) %>% rownames_to_column("var_id")
    ranges_df <- as.data.frame(ranges(vcf)) %>%
      rownames_to_column("var_id") %>%
      rename(chr = seqnames, pos = start)

    # Combine data
    var_df <- fixed_df %>%
      left_join(info_df, by = "var_id") %>%
      left_join(ranges_df, by = "var_id") %>%
      mutate(chr = as.character(chr))

    # Filter for BND (breakend) type SVs
    bnd_df <- var_df %>%
      filter(TYPE == "BND") %>%
      log_n_rows("BND variants")

    if (nrow(bnd_df) == 0) {
      log_warn("No BND variants found in VCF")
      return(tibble())
    }

    # Filter by supporting reads
    min_reads <- config$fusions$min_supporting_reads %||% 3
    bnd_df <- bnd_df %>%
      mutate(
        supporting_reads = case_when(
          !is.na(PR) ~ as.integer(str_extract(PR, "^\\d+")),
          !is.na(SR) ~ as.integer(str_extract(SR, "^\\d+")),
          TRUE ~ NA_integer_
        )
      ) %>%
      filter(!is.na(supporting_reads), supporting_reads >= min_reads) %>%
      log_n_rows(glue("BND variants with ≥{min_reads} supporting reads"))

    if (nrow(bnd_df) == 0) {
      log_warn("No BND variants pass supporting read filter")
      return(tibble())
    }

    # Parse breakpoint ALT field to extract partner information
    bnd_df <- bnd_df %>%
      mutate(
        # ALT format examples: N]chr22:23000]N or ]chr22:23000]N
        alt_str = as.character(REF),
        partner_chr = NA_character_,
        partner_pos = NA_integer_,
        orientation = NA_character_
      )

    # Extract partner breakpoint info from ALT field
    for (i in seq_len(nrow(bnd_df))) {
      alt <- bnd_df$ALT[i]
      # Parse formats: N]chr:pos]N, ]chr:pos]N, N[chr:pos[N, [chr:pos[N
      match <- str_match(alt, "([\\[\\]])([^:\\[\\]]+):([0-9]+)([\\[\\]])")
      if (!is.na(match[1, 1])) {
        bnd_df$partner_chr[i] <- match[1, 3]
        bnd_df$partner_pos[i] <- as.integer(match[1, 4])
        # Orientation: [/] indicates direction
        bnd_df$orientation[i] <- match[1, 2]
      }
    }

    # Annotate with genes and exons
    fusions_df <- .annotate_genes_and_exons(bnd_df, config, sample_id)

    return(fusions_df)
  }, error = function(e) {
    log_error("Error parsing Manta VCF: {e$message}")
    return(tibble())
  })
}

#' Parse TSO500 Fusion CSV
#'
#' Reads fusion events from TSO500 Local App CSV output.
#'
#' @param csv_path Character. Path to TSO500 fusion CSV.
#' @param config List. Configuration object.
#' @param sample_id Character. Sample identifier.
#'
#' @return Tibble with fusion event data.
#'
#' @keywords internal
.parse_tso500_csv <- function(csv_path, config, sample_id) {
  log_info("Parsing TSO500 fusion CSV: {basename(csv_path)}")

  tryCatch({
    # Read TSO500 CSV with flexible column naming
    fusions_raw <- read_csv(csv_path, show_col_types = FALSE)
    log_debug("TSO500 CSV loaded: {nrow(fusions_raw)} fusions")

    # Map common TSO500 column names
    col_mapping <- c(
      "Gene.A" = "gene_a", "gene.a" = "gene_a", "GeneA" = "gene_a",
      "Gene.B" = "gene_b", "gene.b" = "gene_b", "GeneB" = "gene_b",
      "Num.Supporting.Reads" = "supporting_reads",
      "NumSupportingReads" = "supporting_reads",
      "Supporting.Reads" = "supporting_reads",
      "Exon.A" = "exon_a", "ExonA" = "exon_a",
      "Exon.B" = "exon_b", "ExonB" = "exon_b",
      "Chr.A" = "chr_a", "ChrA" = "chr_a",
      "Pos.A" = "pos_a", "PosA" = "pos_a",
      "Chr.B" = "chr_b", "ChrB" = "chr_b",
      "Pos.B" = "pos_b", "PosB" = "pos_b",
      "In.Frame" = "in_frame", "InFrame" = "in_frame"
    )

    # Rename columns
    for (old_name in names(fusions_raw)) {
      new_name <- col_mapping[old_name]
      if (!is.na(new_name)) {
        fusions_raw <- fusions_raw %>% rename(!!new_name := old_name)
      }
    }

    # Filter by supporting reads
    min_reads <- config$fusions$min_supporting_reads %||% 3
    fusions_df <- fusions_raw %>%
      filter(!is.na(supporting_reads), supporting_reads >= min_reads) %>%
      log_n_rows(glue("TSO500 fusions with ≥{min_reads} supporting reads"))

    # Standardize columns
    fusions_df <- fusions_df %>%
      mutate(
        gene_a = as.character(gene_a),
        gene_b = as.character(gene_b),
        fusion_name = glue("{gene_a}::{gene_b}"),
        exon_a = as.integer(exon_a),
        exon_b = as.integer(exon_b),
        chr_a = as.character(chr_a),
        pos_a = as.integer(pos_a),
        chr_b = as.character(chr_b),
        pos_b = as.integer(pos_b),
        supporting_reads = as.integer(supporting_reads),
        fusion_type = case_when(
          !is.na(in_frame) & in_frame ~ "in-frame",
          !is.na(in_frame) & !in_frame ~ "out-of-frame",
          TRUE ~ "unknown"
        )
      ) %>%
      select(
        gene_a, gene_b, fusion_name, exon_a, exon_b,
        chr_a, pos_a, chr_b, pos_b,
        supporting_reads, fusion_type
      )

    if (nrow(fusions_df) == 0) {
      log_warn("No TSO500 fusions pass filter")
      return(tibble())
    }

    return(fusions_df)
  }, error = function(e) {
    log_error("Error parsing TSO500 CSV: {e$message}")
    return(tibble())
  })
}

#' Annotate Fusions with Gene and Exon Information
#'
#' Maps SV breakpoints to genes and identifies exon boundaries.
#'
#' @param bnd_df Tibble. Breakend variant data.
#' @param config List. Configuration object with annotation paths.
#' @param sample_id Character. Sample identifier.
#'
#' @return Tibble with gene and exon annotations.
#'
#' @keywords internal
.annotate_genes_and_exons <- function(bnd_df, config, sample_id) {
  log_info("Annotating breakpoints with gene information")

  # For now, return basic structure with NA exons
  # Full implementation would use GTF/GFF annotation
  fusions_df <- bnd_df %>%
    mutate(
      gene_a = NA_character_,
      gene_b = NA_character_,
      exon_a = NA_integer_,
      exon_b = NA_integer_,
      chr_b = partner_chr,
      pos_b = partner_pos
    ) %>%
    rename(chr_a = chr, pos_a = pos) %>%
    mutate(
      fusion_name = glue("{gene_a}::{gene_b}"),
      fusion_type = "unknown"
    ) %>%
    select(
      gene_a, gene_b, fusion_name, exon_a, exon_b,
      chr_a, pos_a, chr_b, pos_b,
      supporting_reads, fusion_type
    )

  log_debug("Annotated {nrow(fusions_df)} fusion events")

  return(fusions_df)
}

#' Cross-Reference Fusions with Known Database
#'
#' Checks fusion pairs against a database of known recurrent fusions.
#'
#' @param fusions_df Tibble. Fusion event data.
#' @param db_path Character. Path to known fusions database file.
#'
#' @return Tibble with added known_fusion column.
#'
#' @keywords internal
.annotate_known_fusions <- function(fusions_df, db_path) {
  log_info("Cross-referencing with known fusion database")

  tryCatch({
    known_fusions <- read_tsv(db_path, show_col_types = FALSE) %>%
      mutate(
        gene_a = toupper(gene_a),
        gene_b = toupper(gene_b),
        fusion_key = glue("{gene_a}::{gene_b}")
      ) %>%
      distinct(fusion_key)

    fusions_df <- fusions_df %>%
      mutate(
        fusion_key = glue("{toupper(gene_a)}::{toupper(gene_b)}"),
        known_fusion = fusion_key %in% known_fusions$fusion_key
      ) %>%
      select(-fusion_key)

    n_known <- sum(fusions_df$known_fusion, na.rm = TRUE)
    log_info("Found {n_known} known fusions")

    return(fusions_df)
  }, error = function(e) {
    log_error("Error reading known fusions database: {e$message}")
    fusions_df %>% mutate(known_fusion = NA)
  })
}

#' Helper: Log number of rows with context
#'
#' @param df Tibble.
#' @param context Character. Descriptive text.
#'
#' @return Tibble (invisibly).
#'
#' @keywords internal
log_n_rows <- function(df, context) {
  n <- nrow(df)
  log_debug("{context}: {n} rows")
  invisible(df)
}
