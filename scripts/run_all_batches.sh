#!/usr/bin/env bash
# scripts/run_all_batches.sh — Iterate the master manifest and run the full
# pipeline + Foundation One render for any sample that doesn't already have
# a report. Per-sample errors do not abort the batch.
#
# Usage:
#   bash scripts/run_all_batches.sh                # all samples in manifest
#   bash scripts/run_all_batches.sh M25-1093R      # one sample (regex match)

set -uo pipefail
cd /Users/htlin/ngs-tertiary-analysis-skills

MANIFEST=inputs/TSO500-HRD/_manifest.tsv
LOG_DIR=reports/_run_logs
mkdir -p "$LOG_DIR"
ONLY="${1:-}"

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: manifest not found at $MANIFEST. Run scripts/build_batch_manifest.py first."
  exit 1
fi

T0=$(date +%s)
NSKIP=0
NRUN=0
NFAIL=0
declare -a FAILED=()

# Header: sample_id batch tumor_zh oncokb_code tsv_path
while IFS=$'\t' read -r sid batch zh code tsv; do
  [[ "$sid" == "sample_id" ]] && continue
  [[ -z "$sid" || -z "$code" ]] && continue
  [[ -n "$ONLY" && "$sid" != *"$ONLY"* ]] && continue

  out=reports/${sid}/08-report/clinical_actionable.html
  if [[ -f "$out" ]]; then
    NSKIP=$((NSKIP+1))
    continue
  fi

  if [[ ! -f "$tsv" ]]; then
    echo "  [MISS-TSV] $sid (expected $tsv)"
    NFAIL=$((NFAIL+1)); FAILED+=("$sid:no-tsv")
    continue
  fi

  log="$LOG_DIR/${sid}.log"
  echo "[$NRUN] $sid ($batch / $code)"
  if Rscript --vanilla scripts/run_one_sample_tso500.R "$sid" "$tsv" "$code" \
       > "$log" 2>&1; then
    NRUN=$((NRUN+1))
    if [[ ! -f "$out" ]]; then
      echo "    WARN no actionable.html produced — see $log"
      NFAIL=$((NFAIL+1)); FAILED+=("$sid:no-html")
    fi
  else
    echo "    FAIL ($?) — see $log"
    NFAIL=$((NFAIL+1)); FAILED+=("$sid:driver-exit")
  fi
done < "$MANIFEST"

ELAPSED=$(( $(date +%s) - T0 ))
echo
echo "========================================================================"
echo "  RUN: $NRUN   SKIP: $NSKIP   FAIL: $NFAIL   ELAPSED: ${ELAPSED}s"
echo "========================================================================"
if [[ ${#FAILED[@]} -gt 0 ]]; then
  printf '    failed: %s\n' "${FAILED[@]}"
fi
