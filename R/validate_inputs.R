# R/validate_inputs.R — Input validation for NGS tertiary analysis pipeline
# Validates BAM files, configuration, API keys, and reference files

suppressPackageStartupMessages({
  library(logger)
  library(glue)
  library(fs)
  library(httr2)
})

#' Validate pipeline inputs
#'
#' @param bam_path Path to BAM file
#' @param config Configuration list
#'
#' @return List with components:
#'   - `valid` (logical): Overall validation status
#'   - `errors` (character): Vector of error messages
#'   - `warnings` (character): Vector of warning messages
#'
#' @details
#' Performs comprehensive validation:
#' - BAM file existence, size, indexing, and magic bytes
#' - Configuration required keys and values
#' - API key availability and connectivity (OncoKB)
#' - Reference files (optional, warnings only)
#'
#' @export
validate_pipeline_inputs <- function(bam_path, config) {
  errors <- character()
  warnings <- character()

  log_info("Starting pipeline input validation")

  # ── Validate BAM file ───────────────────────────────────────────────────────

  log_info("Validating BAM file: {bam_path}")

  # Check file exists
  if (!file.exists(bam_path)) {
    errors <- c(errors, glue("BAM file not found: {bam_path}"))
  } else {
    log_debug("BAM file exists")

    # Check file is readable
    if (file.access(bam_path, mode = 4) != 0) {
      errors <- c(errors, glue("BAM file not readable: {bam_path}"))
    } else {
      log_debug("BAM file is readable")
    }

    # Check file size > 0
    bam_size <- file.size(bam_path)
    if (is.na(bam_size) || bam_size == 0) {
      errors <- c(errors, glue("BAM file is empty or unreadable size: {bam_path}"))
    } else {
      log_debug("BAM file size: {format(as.numeric(bam_size), big.mark=',')} bytes")
    }

    # Check BAM magic bytes (BAM\1)
    tryCatch({
      con <- file(bam_path, "rb")
      magic <- readBin(con, "raw", 4)
      close(con)

      expected_magic <- as.raw(c(0x42, 0x41, 0x4d, 0x01)) # BAM\1
      if (!identical(magic, expected_magic)) {
        errors <- c(errors,
          glue("BAM magic bytes invalid. Expected BAM\\1, got {paste(magic, collapse=',')}"))
      } else {
        log_debug("BAM magic bytes valid")
      }
    }, error = function(e) {
      errors <<- c(errors, glue("Failed to read BAM magic bytes: {e$message}"))
    })

    # Check BAM index
    bai_path_1 <- paste0(bam_path, ".bai")
    bai_path_2 <- file.path(dirname(bam_path), paste0(basename(bam_path), ".bai"))

    if (!file.exists(bai_path_1) && !file.exists(bai_path_2)) {
      log_warn("BAM index not found, attempting to create with samtools")

      tryCatch({
        result <- system2("samtools", c("index", bam_path),
                         stdout = TRUE, stderr = TRUE)
        exit_code <- attr(result, "status") %||% 0L

        if (exit_code != 0) {
          errors <- c(errors,
            glue("Failed to create BAM index. samtools returned: {exit_code}"))
        } else {
          log_info("BAM index created successfully")
        }
      }, error = function(e) {
        errors <<- c(errors, glue("samtools not available or failed: {e$message}"))
      })
    } else {
      log_debug("BAM index found")
    }
  }

  # ── Validate configuration ──────────────────────────────────────────────────

  log_info("Validating configuration")

  if (!is.list(config)) {
    errors <- c(errors, "config must be a list")
  } else {
    # Required top-level keys
    required_keys <- c("reference", "variant_calling", "biomarkers")
    missing_keys <- setdiff(required_keys, names(config))

    if (length(missing_keys) > 0) {
      errors <- c(errors,
        glue("Missing required config sections: {paste(missing_keys, collapse=', ')}"))
    }

    # Validate reference section
    if ("reference" %in% names(config)) {
      ref <- config$reference

      # Check required reference keys
      if (is.null(ref$genome) || nchar(as.character(ref$genome)) == 0) {
        errors <- c(errors, "reference.genome is required and must be non-empty")
      } else {
        log_debug("reference.genome: {ref$genome}")
      }

      if (is.null(ref$panel_bed) || nchar(as.character(ref$panel_bed)) == 0) {
        errors <- c(errors, "reference.panel_bed is required and must be non-empty")
      } else {
        # Check panel BED file exists
        panel_bed <- ref$panel_bed
        if (!file.exists(panel_bed)) {
          errors <- c(errors, glue("Panel BED file not found: {panel_bed}"))
        } else {
          log_debug("Panel BED file exists: {panel_bed}")
        }
      }

      # Validate optional reference files (warnings only)
      if (!is.null(ref$fasta) && nchar(as.character(ref$fasta)) > 0) {
        if (!file.exists(ref$fasta)) {
          warnings <- c(warnings, glue("Reference FASTA not found: {ref$fasta}"))
        } else {
          log_debug("Reference FASTA found: {ref$fasta}")
        }
      }

      if (!is.null(ref$vep_cache_dir) && nchar(as.character(ref$vep_cache_dir)) > 0) {
        if (!dir.exists(ref$vep_cache_dir)) {
          warnings <- c(warnings, glue("VEP cache directory not found: {ref$vep_cache_dir}"))
        } else {
          log_debug("VEP cache directory found: {ref$vep_cache_dir}")
        }
      }
    }

    # Validate variant_calling section
    if ("variant_calling" %in% names(config)) {
      vc <- config$variant_calling

      if (is.null(vc$min_vaf) || !is.numeric(vc$min_vaf)) {
        errors <- c(errors, "variant_calling.min_vaf is required and must be numeric")
      } else if (vc$min_vaf < 0 || vc$min_vaf > 1) {
        errors <- c(errors,
          glue("variant_calling.min_vaf must be between 0 and 1, got: {vc$min_vaf}"))
      } else {
        log_debug("variant_calling.min_vaf: {vc$min_vaf}")
      }
    }

    # Validate biomarkers section
    if ("biomarkers" %in% names(config)) {
      bm <- config$biomarkers

      if (!is.null(bm$tmb) && !is.null(bm$tmb$panel_coding_size_mb)) {
        pcm <- bm$tmb$panel_coding_size_mb
        if (!is.numeric(pcm) || pcm <= 0) {
          errors <- c(errors,
            glue("biomarkers.tmb.panel_coding_size_mb must be positive, got: {pcm}"))
        } else {
          log_debug("biomarkers.tmb.panel_coding_size_mb: {pcm} Mb")
        }
      }
    }
  }

  # ── Validate API keys ───────────────────────────────────────────────────────

  log_info("Validating API keys")

  # OncoKB API key (required)
  oncokb_api_key <- Sys.getenv("ONCOKB_API_KEY", unset = "")
  if (nchar(oncokb_api_key) == 0) {
    errors <- c(errors, "ONCOKB_API_KEY environment variable not set")
  } else {
    log_debug("ONCOKB_API_KEY is set")

    # Test OncoKB connectivity
    tryCatch({
      log_info("Testing OncoKB connectivity")
      response <- httr2::request("https://www.oncokb.org/api/v1") %>%
        httr2::req_url_path_append("utils/allCuratedGenes") %>%
        httr2::req_headers("Authorization" = glue("Bearer {oncokb_api_key}")) %>%
        httr2::req_timeout(10) %>%
        httr2::req_perform()

      if (httr2::resp_status(response) == 200) {
        log_info("OncoKB connectivity test passed")
      } else {
        warnings <- c(warnings,
          glue("OncoKB API returned status {httr2::resp_status(response)}"))
      }
    }, error = function(e) {
      warnings <<- c(warnings, glue("OncoKB connectivity test failed: {e$message}"))
    })
  }

  # PubMed API key (optional, warning if missing)
  pubmed_api_key <- Sys.getenv("PUBMED_API_KEY", unset = "")
  if (nchar(pubmed_api_key) == 0) {
    warnings <- c(warnings, "PUBMED_API_KEY not set (optional)")
  } else {
    log_debug("PUBMED_API_KEY is set")
  }

  # Scopus API key (optional, warning if missing)
  scopus_api_key <- Sys.getenv("SCOPUS_API_KEY", unset = "")
  if (nchar(scopus_api_key) == 0) {
    warnings <- c(warnings, "SCOPUS_API_KEY not set (optional)")
  } else {
    log_debug("SCOPUS_API_KEY is set")
  }

  # ── Final validation status ─────────────────────────────────────────────────

  valid <- length(errors) == 0

  log_info("Validation complete: valid={valid}, errors={length(errors)}, warnings={length(warnings)}")

  if (!valid) {
    log_error("Validation failed with {length(errors)} error(s)")
    error_msg <- glue("Pipeline input validation failed:\n{paste('  - ', errors, collapse='\n')}")
    stop(error_msg)
  }

  if (length(warnings) > 0) {
    for (w in warnings) {
      log_warn(w)
    }
  }

  # Return validation results
  list(
    valid = valid,
    errors = errors,
    warnings = warnings
  )
}
