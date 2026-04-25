# agents/io.R
# Filesystem bridge between R benchmark drivers and the Claude Code orchestrator.
#
# R scripts write structured INPUTS (one JSON per case) under
#   agents/inputs/<agent_name>/<case_id>.json
# The Claude Code session (the operator) invokes the corresponding subagent
# for each input and writes the response to
#   agents/outputs/<agent_name>/<case_id>.json
# R analysis scripts then validate outputs against the agent's JSON schema and
# compare to the deterministic baseline.

suppressPackageStartupMessages({
  library(jsonlite)
  library(here)
  library(fs)
  library(logger)
  library(glue)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# ── Paths ────────────────────────────────────────────────────────────────────

agent_input_dir <- function(agent_name) {
  d <- here::here("agents", "inputs", agent_name)
  dir_create(d)
  d
}

agent_output_dir <- function(agent_name) {
  d <- here::here("agents", "outputs", agent_name)
  dir_create(d)
  d
}

agent_schema_path <- function(agent_name) {
  here::here("agents", agent_name, "schema.json")
}

# ── Writing inputs ────────────────────────────────────────────────────────────

#' Write one input JSON for an agent case.
#'
#' The JSON has two top-level fields:
#'   - case_id   : identifier the orchestrator will echo in the output filename
#'   - user_turn : the text the operator should pass to the subagent as the
#'                 user message (already formatted by the agent-specific
#'                 format_input() helper)
#'   - metadata  : free-form context not shown to the agent (e.g. gene,
#'                 variant, tumor_type, baseline result) — used by R analysis
write_agent_input <- function(agent_name, case_id, user_turn, metadata = list()) {
  payload <- list(
    case_id   = case_id,
    user_turn = user_turn,
    metadata  = metadata
  )
  path <- file.path(agent_input_dir(agent_name), paste0(case_id, ".json"))
  write_json(payload, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  log_info("Wrote input: {path}")
  invisible(path)
}

# ── Reading outputs ───────────────────────────────────────────────────────────

#' Read and validate one agent output.
#'
#' The output file is the raw JSON the subagent produced. This function parses
#' it and checks that all keys in the schema's `required` list are present.
#' Schema-level type/enum validation is opportunistic — we verify enums here
#' and leave deeper validation to downstream analysis.
#'
#' @return list with fields: case_id, output (parsed), valid (logical),
#'         errors (character vector)
read_agent_output <- function(agent_name, case_id) {
  out_path <- file.path(agent_output_dir(agent_name), paste0(case_id, ".json"))
  if (!file.exists(out_path)) {
    return(list(case_id = case_id, output = NULL, valid = FALSE,
                errors = glue("Output file missing: {out_path}")))
  }

  parsed <- tryCatch(
    fromJSON(out_path, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(parsed)) {
    return(list(case_id = case_id, output = NULL, valid = FALSE,
                errors = glue("Not valid JSON: {out_path}")))
  }

  schema <- fromJSON(agent_schema_path(agent_name), simplifyVector = FALSE)
  errors <- character()

  # Required keys
  required <- unlist(schema$required %||% list())
  missing_keys <- setdiff(required, names(parsed))
  if (length(missing_keys) > 0) {
    errors <- c(errors, glue("Missing required key(s): {paste(missing_keys, collapse=', ')}"))
  }

  # Enum check
  for (key in names(parsed)) {
    enum <- schema$properties[[key]]$enum
    if (!is.null(enum)) {
      allowed <- unlist(enum)
      if (!parsed[[key]] %in% allowed) {
        errors <- c(errors, glue("{key}='{parsed[[key]]}' not in enum [{paste(allowed, collapse=', ')}]"))
      }
    }
  }

  list(
    case_id = case_id,
    output  = parsed,
    valid   = length(errors) == 0,
    errors  = if (length(errors) == 0) character() else errors
  )
}
