# R/tso500_parser.R — Parse Illumina TSO500 CombinedVariantOutput.tsv
#
# TSO500 already provides full annotation (gene, HGVS, consequence) plus
# pre-computed TMB, MSI, GIS, CNV, LoH, fusion calls. This parser converts
# those sections into the tibbles the rest of the pipeline expects, so we
# can skip VEP and the targets DAG.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readr)
})

#' Read a TSO500 CombinedVariantOutput.tsv and split into named sections.
#' @return Named list: section header (without brackets) → character vector of body lines
read_tso500_sections <- function(tsv_path) {
  lines <- read_lines(tsv_path)
  is_header <- grepl("^\\[.+\\]\\s*$", str_trim(str_replace_all(lines, "\\t+$", "")))
  header_idx <- which(is_header)
  headers <- str_match(str_trim(str_replace_all(lines[header_idx], "\\t+$", "")),
                       "^\\[(.+)\\]")[, 2]
  ends <- c(header_idx[-1] - 1L, length(lines))
  out <- Map(function(s, e) lines[(s + 1L):e], header_idx, ends)
  names(out) <- headers
  lapply(out, function(x) x[!grepl("^\\s*$", x)])
}

# ── Section parsers ──────────────────────────────────────────────────────────

parse_kv_section <- function(section_lines) {
  if (length(section_lines) == 0) return(list())
  parts <- str_split_fixed(section_lines, "\\t", 2)
  vals <- str_trim(parts[, 2])
  names(vals) <- str_trim(parts[, 1])
  as.list(vals)
}

parse_table_section <- function(section_lines) {
  if (length(section_lines) < 2) return(tibble())
  if (length(section_lines) == 2 && trimws(section_lines[[2]]) %in% c("", "NA")) {
    return(tibble())
  }
  hdr <- str_split_fixed(section_lines[[1]], "\\t", Inf)[1, ]
  body <- section_lines[-1]
  body <- body[!grepl("^\\s*NA\\s*$", body) & !grepl("^\\s*$", body)]
  if (length(body) == 0) return(tibble())
  mat <- str_split_fixed(body, "\\t", length(hdr))
  colnames(mat) <- str_trim(hdr)
  as_tibble(mat) |> mutate(across(everything(), ~ na_if(str_trim(.x), "")))
}

# ── Specific section converters ──────────────────────────────────────────────

tso500_qc <- function(sections) {
  tmb <- parse_kv_section(sections[["TMB"]] %||% character())
  list(
    summary = list(
      total_reads = NA_integer_, mapped_reads = NA_integer_,
      mapping_rate = NA_real_, mean_coverage = NA_real_,
      on_target_rate = NA_real_, duplicate_rate = NA_real_,
      tumor_purity = suppressWarnings(as.numeric(
        parse_kv_section(sections[["GIS"]] %||% character())[["Tumor Fraction"]] %||% NA))
    ),
    per_gene_coverage = tibble(),
    pass = TRUE,
    skipped = TRUE,
    skip_reason = "TSO500 deliverable input — pre-computed metrics shown elsewhere"
  )
}

tso500_tmb <- function(sections) {
  kv <- parse_kv_section(sections[["TMB"]] %||% character())
  tmb_score <- suppressWarnings(as.numeric(kv[["Total TMB"]]))
  variant_count <- suppressWarnings(as.integer(kv[["Number of Passing Eligible Variants"]]))
  panel_mb <- suppressWarnings(as.numeric(kv[["Coding Region Size in Megabases"]]))
  tmb_class <- if (is.na(tmb_score)) "Unknown"
               else if (tmb_score >= 10) "TMB-High"
               else if (tmb_score >= 6) "TMB-Intermediate" else "TMB-Low"
  list(tmb_score = tmb_score, tmb_class = tmb_class,
       variant_count = variant_count %||% NA_integer_,
       panel_size_mb = panel_mb %||% NA_real_)
}

tso500_msi <- function(sections) {
  kv <- parse_kv_section(sections[["MSI"]] %||% character())
  pct <- suppressWarnings(as.numeric(kv[["Percent Unstable MSI Sites"]]))
  status <- if (is.na(pct)) "Unknown"
            else if (pct >= 20) "MSI-High"
            else if (pct >= 10) "MSI-Intermediate" else "MSS"
  list(msi_score = pct,
       msi_status = status,
       unstable_sites = suppressWarnings(as.integer(kv[["Total MSI Sites Unstable"]])),
       total_sites    = suppressWarnings(as.integer(kv[["Usable MSI Sites"]])))
}

tso500_hrd <- function(sections) {
  kv <- parse_kv_section(sections[["GIS"]] %||% character())
  gis <- suppressWarnings(as.numeric(kv[["Genomic Instability Score"]]))
  status <- if (is.na(gis)) "Unknown"
            else if (gis >= 42) "HRD-Positive (GIS≥42)"
            else if (gis >= 30) "HRD-Intermediate"
            else "HRD-Negative"
  list(hrd_score = gis, hrd_status = status,
       loh_score = NA_real_, tai_score = NA_real_, lst_score = NA_real_,
       is_reliable = !is.na(gis),
       tumor_fraction = suppressWarnings(as.numeric(kv[["Tumor Fraction"]])),
       ploidy         = suppressWarnings(as.numeric(kv[["Ploidy"]])))
}

tso500_cnv <- function(sections) {
  empty <- tibble(chromosome = character(), start = integer(), end = integer(),
                  gene = character(), log2_ratio = numeric(),
                  copy_number = numeric(), length = integer(), type = character())

  cnv <- parse_table_section(sections[["Copy Number Variants"]] %||% character())
  cnv_part <- if (nrow(cnv) == 0) empty else {
    fold <- suppressWarnings(as.numeric(cnv[["Fold Change"]]))
    cn   <- suppressWarnings(as.numeric(cnv[["Absolute Copy Number"]]))
    type_raw <- str_replace_all(cnv[["Copy Number Variant"]] %||% rep("", nrow(cnv)),
                                "[<>]", "")
    type <- case_when(
      type_raw == "DUP" & fold >= 2.5 ~ "AMPLIFICATION",
      type_raw == "DUP"               ~ "GAIN",
      type_raw == "DEL" & fold <= 0.5 ~ "DELETION",
      type_raw == "DEL"               ~ "LOSS",
      TRUE                            ~ "NEUTRAL"
    )
    tibble(
      chromosome = NA_character_,
      start = NA_integer_, end = NA_integer_,
      gene = cnv[["Gene"]],
      log2_ratio = log2(pmax(fold, 1e-3)),
      copy_number = cn,
      length = NA_integer_,
      type = type
    )
  }

  # [Large Rearrangements] section also describes gene-level CNVs (most often
  # BRCA1/BRCA2 LGRs). Merge them so OncoKB CNA querying picks them up.
  lr <- parse_table_section(sections[["Large Rearrangements"]] %||% character())
  lr_part <- if (nrow(lr) == 0) empty else {
    fold_lr <- suppressWarnings(as.numeric(lr[["Fold Change"]]))
    chrom_lr <- lr[["Chromosome"]]
    start_lr <- suppressWarnings(as.integer(lr[["Start"]]))
    end_lr   <- suppressWarnings(as.integer(lr[["Stop"]]))
    type_lr_raw <- toupper(str_trim(lr[["CNV Type"]] %||% ""))
    type_lr <- case_when(
      type_lr_raw == "LOSS" & fold_lr <= 0.5 ~ "DELETION",
      type_lr_raw == "LOSS"                  ~ "LOSS",
      type_lr_raw == "GAIN" & fold_lr >= 2.5 ~ "AMPLIFICATION",
      type_lr_raw == "GAIN"                  ~ "GAIN",
      TRUE                                    ~ "LOSS"
    )
    tibble(
      chromosome = chrom_lr,
      start = start_lr, end = end_lr,
      gene = lr[["Gene"]],
      log2_ratio = log2(pmax(fold_lr, 1e-3)),
      copy_number = NA_real_,
      length = ifelse(is.na(start_lr) | is.na(end_lr), NA_integer_, end_lr - start_lr),
      type = type_lr
    ) |> filter(!is.na(gene), gene != "")
  }

  bind_rows(cnv_part, lr_part) |>
    distinct(gene, type, .keep_all = TRUE)
}

tso500_large_rearr <- function(sections) {
  lr <- parse_table_section(sections[["Large Rearrangements"]] %||% character())
  if (nrow(lr) == 0) return(tibble())
  lr |>
    mutate(
      Start = suppressWarnings(as.integer(Start)),
      Stop  = suppressWarnings(as.integer(Stop)),
      `Fold Change` = suppressWarnings(as.numeric(`Fold Change`))
    )
}

tso500_fusions <- function(sections) {
  empty <- tibble(gene_a = character(), gene_b = character(),
                  fusion_name = character(), exon_a = integer(), exon_b = integer(),
                  chr_a = character(), pos_a = integer(),
                  chr_b = character(), pos_b = integer(),
                  supporting_reads = integer(), fusion_type = character(),
                  known_fusion = logical())
  fus <- parse_table_section(sections[["Fusions"]] %||% character())
  if (nrow(fus) == 0) return(empty)
  parts <- str_split_fixed(fus[["Gene Pair"]], "[-:]", 2)
  bp1 <- str_split_fixed(fus[["Breakpoint 1"]] %||% rep("", nrow(fus)), ":", 2)
  bp2 <- str_split_fixed(fus[["Breakpoint 2"]] %||% rep("", nrow(fus)), ":", 2)
  tibble(
    gene_a = parts[, 1], gene_b = parts[, 2],
    fusion_name = fus[["Gene Pair"]],
    exon_a = NA_integer_, exon_b = NA_integer_,
    chr_a = bp1[, 1], pos_a = suppressWarnings(as.integer(bp1[, 2])),
    chr_b = bp2[, 1], pos_b = suppressWarnings(as.integer(bp2[, 2])),
    supporting_reads = suppressWarnings(as.integer(fus[["Fusion Supporting Reads"]])),
    fusion_type = "in-frame",
    known_fusion = TRUE
  )
}

# Map TSO500 [Small Variants] into the merged_annotations schema.
tso500_variants <- function(sections) {
  sv <- parse_table_section(sections[["Small Variants"]] %||% character())
  empty <- tibble(
    chrom = character(), pos = integer(), ref = character(), alt = character(),
    qual = numeric(), filter = character(),
    gene = character(), transcript_id = character(),
    hgvsc = character(), hgvsp = character(),
    consequence = character(), impact = character(),
    sift_prediction = character(), polyphen_prediction = character(),
    gnomad_af = character(), gnomad_af_numeric = numeric(),
    clinvar_significance = character(), cosmic_id = character(),
    pathogenic_class = character(), pathogenicity = character(),
    vaf = numeric(), alt_depth = integer(), total_depth = integer(),
    confidence = character(), confidence_reason = character(),
    classification = character()
  )
  if (nrow(sv) == 0) return(empty)

  cdot <- sv[["C-Dot Notation"]] %||% rep(NA_character_, nrow(sv))
  pdot <- sv[["P-Dot Notation"]] %||% rep(NA_character_, nrow(sv))
  transcript_id <- str_extract(cdot, "^[A-Z]+_[0-9]+\\.[0-9]+")
  hgvsc <- str_extract(cdot, "c\\..+$")
  hgvsp <- str_extract(pdot, "p\\.\\(?[A-Za-z0-9_*=>]+\\)?")
  hgvsp <- str_replace_all(hgvsp, "[()]", "")

  consequence <- sv[["Consequence(s)"]] %||% rep(NA_character_, nrow(sv))
  impact <- case_when(
    is.na(consequence) ~ "MODIFIER",
    str_detect(consequence, "stop_gained|stop_lost|frameshift|start_lost|splice_acceptor|splice_donor") ~ "HIGH",
    str_detect(consequence, "missense|inframe|protein_altering") ~ "MODERATE",
    str_detect(consequence, "synonymous|splice_region|5_prime_UTR|3_prime_UTR") ~ "LOW",
    TRUE ~ "MODIFIER"
  )
  pathogenicity <- case_when(
    impact == "HIGH"     ~ "likely_pathogenic",
    impact == "MODERATE" ~ "uncertain_significance",
    impact == "LOW"      ~ "likely_benign",
    TRUE                 ~ "uncertain_significance"
  )

  vaf <- suppressWarnings(as.numeric(sv[["Allele Frequency"]]))
  total_depth <- suppressWarnings(as.integer(sv[["Depth"]]))
  alt_depth <- as.integer(round(vaf * total_depth))

  confidence <- case_when(
    is.na(vaf) | is.na(total_depth) ~ "Uncertain",
    total_depth >= 200 & vaf >= 0.1 ~ "High confidence",
    total_depth >= 100 & vaf >= 0.05 ~ "Moderate confidence",
    TRUE ~ "Questionable"
  )
  confidence_reason <- glue::glue("VAF {round(vaf*100,1)}% at {total_depth}x depth")

  tibble(
    chrom = sv[["Chromosome"]],
    pos = suppressWarnings(as.integer(sv[["Genomic Position"]])),
    ref = sv[["Reference Call"]],
    alt = sv[["Alternative Call"]],
    qual = NA_real_, filter = "PASS",
    gene = sv[["Gene"]],
    transcript_id = transcript_id,
    hgvsc = hgvsc, hgvsp = hgvsp,
    consequence = consequence, impact = impact,
    sift_prediction = NA_character_, polyphen_prediction = NA_character_,
    gnomad_af = NA_character_, gnomad_af_numeric = NA_real_,
    clinvar_significance = NA_character_, cosmic_id = NA_character_,
    pathogenic_class = pathogenicity,
    pathogenicity = pathogenicity,
    classification = pathogenicity,
    vaf = vaf, alt_depth = alt_depth, total_depth = total_depth,
    confidence = confidence,
    confidence_reason = as.character(confidence_reason)
  ) |>
    filter(!is.na(gene) & gene != "")
}

#' Top-level: parse a TSO500 CombinedVariantOutput.tsv into the data shapes the
#' downstream pipeline expects.
parse_tso500_combined_output <- function(tsv_path) {
  stopifnot(file.exists(tsv_path))
  sections <- read_tso500_sections(tsv_path)
  list(
    qc        = tso500_qc(sections),
    variants  = tso500_variants(sections),
    cnv       = tso500_cnv(sections),
    fusions   = tso500_fusions(sections),
    tmb       = tso500_tmb(sections),
    msi       = tso500_msi(sections),
    hrd       = tso500_hrd(sections),
    large_rearrangements = tso500_large_rearr(sections),
    analysis_details = parse_kv_section(sections[["Analysis Details"]] %||% character()),
    sequencing       = parse_kv_section(sections[["Sequencing Run Details"]] %||% character())
  )
}

# Helper: OncoKB tumor type code → display label (for the report header)
tumor_type_display <- function(code) {
  map <- c(
    UCEC = "Endometrial Carcinoma (UCEC)",
    BRCA = "Breast Cancer (BRCA)",
    COAD = "Colorectal Adenocarcinoma (COAD)",
    READ = "Rectal Adenocarcinoma (READ)",
    LUAD = "Lung Adenocarcinoma (LUAD)",
    NSCLC = "Non-Small Cell Lung Cancer (NSCLC)",
    OV   = "Ovarian Cancer (OV)",
    HGSOC = "High-Grade Serous Ovarian Carcinoma (HGSOC)",
    PAAD = "Pancreatic Adenocarcinoma (PAAD)",
    STAD = "Stomach Adenocarcinoma (STAD)"
  )
  unname(map[code]) %||% code
}
