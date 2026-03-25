# R/utils.R — Shared utilities for NGS tertiary analysis pipeline

# ── Package loading ──────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(yaml)
  library(fs)
  library(logger)
  library(glue)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

# ── Environment ──────────────────────────────────────────────────────────────

#' Load .env file into environment variables
#' @param env_file Path to .env file (default: project root .env)
load_env <- function(env_file = here::here(".env")) {
  if (!file.exists(env_file)) {
    log_warn("No .env file found at {env_file}")
    return(invisible(NULL))
  }
  lines <- readLines(env_file, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines) & nchar(trimws(lines)) > 0]
  for (line in lines) {
    eq_pos <- regexpr("=", line, fixed = TRUE)
    if (eq_pos < 0) next
    key <- trimws(substr(line, 1, eq_pos - 1))
    value <- trimws(substr(line, eq_pos + 1, nchar(line)))
    if (nchar(key) > 0) {
      # Strip surrounding quotes
      value <- str_remove_all(value, '^["\']|["\']$')
      do.call(Sys.setenv, setNames(list(value), key))
    }
  }
  log_info("Loaded {length(lines)} env vars from {env_file}")
}

#' Get a required API key from environment
#' @param key Environment variable name
#' @return The API key value
#' @examples get_api_key("ONCOKB_API_KEY")
get_api_key <- function(key) {
  value <- Sys.getenv(key, unset = "")
  if (nchar(value) == 0) {
    log_error("{key} not set. Add it to .env file.")
    stop(glue("{key} is required but not set"))
  }
  value
}

# ── Configuration ────────────────────────────────────────────────────────────

#' Load pipeline configuration from YAML
#' @param config_path Path to YAML config file
#' @param overrides Named list of config overrides (dot-separated keys)
#' @return Nested list of config values
load_config <- function(config_path = here::here("config/default.yaml"),
                        overrides = list()) {
  if (!file.exists(config_path)) {
    stop(glue("Config file not found: {config_path}"))
  }
  config <- read_yaml(config_path)

  # Apply overrides (dot-separated keys like "qc.min_mean_coverage")
  for (key in names(overrides)) {
    parts <- str_split_1(key, "\\.")
    config <- set_nested(config, parts, overrides[[key]])
  }

  log_info("Loaded config from {config_path}")
  config
}

#' Set a value in a nested list by key path
#' @param lst Nested list
#' @param keys Character vector of keys
#' @param value Value to set
#' @return Modified list
set_nested <- function(lst, keys, value) {
  if (length(keys) == 1) {
    lst[[keys[1]]] <- value
    return(lst)
  }
  if (is.null(lst[[keys[1]]])) lst[[keys[1]]] <- list()
  lst[[keys[1]]] <- set_nested(lst[[keys[1]]], keys[-1], value)
  lst
}

# ── Output directories ───────────────────────────────────────────────────────

#' Get or create the output directory for a pipeline stage
#' @param sample_id Sample identifier
#' @param stage Stage directory name (e.g., "00-qc")
#' @param base_dir Base reports directory
#' @return Path to the stage output directory (created if needed)
stage_output_dir <- function(sample_id, stage,
                             base_dir = here::here("reports")) {
  out_dir <- path(base_dir, sample_id, stage)
  dir_create(out_dir, recurse = TRUE)
  as.character(out_dir)
}

#' Write a stage result as RDS for targets caching
#' @param data R object to save
#' @param sample_id Sample identifier
#' @param stage Stage name
#' @param filename Output filename (without path)
save_stage_result <- function(data, sample_id, stage, filename) {
  out_dir <- stage_output_dir(sample_id, stage)
  out_path <- path(out_dir, filename)
  saveRDS(data, out_path)
  log_info("Saved {stage}/{filename} for {sample_id}")
  as.character(out_path)
}

# ── External tool wrappers ───────────────────────────────────────────────────

#' Run an external command with logging and error checking
#' @param cmd Command name
#' @param args Character vector of arguments
#' @param description Human-readable description for logging
#' @return stdout as character vector (invisible)
run_tool <- function(cmd, args, description = cmd) {
  log_info("Running: {description}")
  log_debug("Command: {cmd} {paste(args, collapse = ' ')}")

  result <- system2(cmd, args, stdout = TRUE, stderr = TRUE)
  exit_code <- attr(result, "status") %||% 0L

  if (exit_code != 0) {
    log_error("{description} failed with exit code {exit_code}")
    log_error("Output: {paste(result, collapse = '\\n')}")
    stop(glue("{description} failed (exit {exit_code})"))
  }

  log_info("{description} completed successfully")
  invisible(result)
}

#' Check if an external tool is available
#' @param tool Tool name
#' @return TRUE if found, FALSE otherwise
check_tool <- function(tool) {
  found <- nchar(Sys.which(tool)) > 0
  if (!found) log_warn("Tool not found: {tool}")
  found
}

# ── Logging setup ────────────────────────────────────────────────────────────

#' Initialize logging for the pipeline
#' @param level Log level (default: INFO)
#' @param log_file Optional log file path
setup_logging <- function(level = "INFO", log_file = NULL) {
  log_threshold(level)
  if (!is.null(log_file)) {
    dir_create(path_dir(log_file))
    log_appender(appender_tee(log_file))
  }
  log_info("Pipeline logging initialized at {level} level")
}

# ── Input detection ──────────────────────────────────────────────────────────

#' Detect input file type (BAM or VCF)
#' @param input_path Path to input file
#' @return "bam" or "vcf"
detect_input_type <- function(input_path) {
  ext <- tolower(path_ext(input_path))
  # Handle .vcf.gz, .vcf, .bam
  if (ext == "gz") {
    base_ext <- tolower(path_ext(path_ext_remove(input_path)))
    if (base_ext == "vcf") return("vcf")
    if (base_ext == "bam") return("bam")  # unlikely but handle
  }
  if (ext == "vcf") return("vcf")
  if (ext == "bam") return("bam")

  # Fallback: check magic bytes
  con <- file(input_path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", n = 4)

  # BAM magic: 42 41 4d 01
  if (length(magic) >= 4 && identical(magic[1:3], charToRaw("BAM"))) return("bam")

  # VCF magic: ##fileformat (text) or gzip magic 1f 8b
  if (length(magic) >= 2 && magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)) {
    return("vcf")  # gzipped, assume VCF if not BAM
  }
  if (length(magic) >= 2 && identical(magic[1:2], charToRaw("##"))) return("vcf")

  stop(glue("Cannot detect input type for: {input_path}. Expected .bam or .vcf/.vcf.gz"))
}

#' Resolve pipeline input — accepts BAM_PATH or VCF_PATH env vars
#' @return List with: path, type ("bam" or "vcf")
resolve_input <- function() {
  bam <- Sys.getenv("BAM_PATH", unset = "")
  vcf <- Sys.getenv("VCF_PATH", unset = "")
  input <- Sys.getenv("INPUT_PATH", unset = "")  # generic

  # Priority: INPUT_PATH > VCF_PATH > BAM_PATH
  if (nchar(input) > 0) {
    if (!file.exists(input)) stop(glue("INPUT_PATH not found: {input}"))
    path <- normalizePath(input)
    type <- detect_input_type(path)
  } else if (nchar(vcf) > 0) {
    if (!file.exists(vcf)) stop(glue("VCF_PATH not found: {vcf}"))
    path <- normalizePath(vcf)
    type <- "vcf"
  } else if (nchar(bam) > 0) {
    if (!file.exists(bam)) stop(glue("BAM_PATH not found: {bam}"))
    path <- normalizePath(bam)
    type <- "bam"
  } else {
    stop("No input file specified. Set BAM_PATH, VCF_PATH, or INPUT_PATH env var.")
  }

  log_info("Input detected: {type} — {path}")
  list(path = path, type = type)
}

# ── Variant helpers ──────────────────────────────────────────────────────────

#' Format HGVS protein notation
#' @param gene Gene symbol
#' @param protein_change Protein change (e.g., "V600E")
#' @return Formatted string (e.g., "BRAF p.Val600Glu")
format_hgvs_protein <- function(gene, protein_change) {
  glue("{gene} p.{protein_change}")
}

#' Classify variant by VAF and context
#' @param vaf Variant allele frequency
#' @param tumor_purity Estimated tumor purity
#' @return Adjusted VAF accounting for purity
adjust_vaf_for_purity <- function(vaf, tumor_purity) {
  if (tumor_purity <= 0 || tumor_purity > 1) {
    log_warn("Invalid tumor purity: {tumor_purity}, returning raw VAF")
    return(vaf)
  }
  # Cancer cell fraction estimation (assuming diploid, heterozygous)
  ccf <- (vaf * 2) / tumor_purity
  pmin(ccf, 1.0)  # Cap at 1.0
}
