#!/usr/bin/env bash
# scripts/render_actionable_batch7.sh
#
# End-to-end Foundation One–style render for Batch 7 (8 patients).
#
#   1. (re-)runs the full pipeline per sample if the per-sample .rds cache
#      is missing — picks up the OncoKB encoding fix automatically.
#   2. Renders the new actionable HTML + PDF.
#   3. Stages clinical_actionable.{html,pdf} into reports/_handoff_Batch_7…
#      and rewrites that folder's index.html.
#
# Usage:
#   bash scripts/render_actionable_batch7.sh                # all 8
#   bash scripts/render_actionable_batch7.sh M26-0165R      # just one

set -uo pipefail
cd /Users/htlin/ngs-tertiary-analysis-skills

BATCH=inputs/TSO500-HRD/Batch_7_20260505
HANDOFF=reports/_handoff_Batch_7_20260505
mkdir -p "$HANDOFF"

declare -a SAMPLES=(
  "M26-0165R UCEC"
  "M26-0234R BRCA"
  "M26-0284R COAD"
  "M26-0320R BRCA"
  "M26-0400R UCEC"
  "M26-0401R COAD"
  "M26-0410R UCEC"
  "M26-0414R UCEC"
)
ONLY="${1:-}"

declare -a STATUS=()
TOTAL_T0=$(date +%s)

for entry in "${SAMPLES[@]}"; do
  sid="${entry% *}"
  tt="${entry#* }"
  [[ -n "$ONLY" && "$ONLY" != "$sid" ]] && continue

  tsv="$BATCH/${sid}_CombinedVariantOutput.tsv"
  if [[ ! -f "$tsv" ]]; then
    STATUS+=("MISS  $sid — no input TSV at $tsv")
    continue
  fi

  pipeline_log="$HANDOFF/${sid}.pipeline.log"
  render_log="$HANDOFF/${sid}.render.log"

  echo
  echo "########################################################################"
  echo "## $sid ($tt)"
  echo "########################################################################"

  # Step 1: full pipeline (only if .rds artefacts are missing OR --force).
  rds=reports/${sid}/08-report/data/oncokb.rds
  if [[ ! -f "$rds" ]]; then
    echo "  [1/2] running full pipeline (OncoKB / CiVIC / AMP / literature)…"
    if ! Rscript --vanilla scripts/run_one_sample_tso500.R "$sid" "$tsv" "$tt" \
         > "$pipeline_log" 2>&1; then
      STATUS+=("PIPE-FAIL $sid — log $pipeline_log")
      continue
    fi
  else
    echo "  [1/2] cached .rds present — skipping pipeline rerun"
  fi

  # Step 2: actionable render
  echo "  [2/2] rendering actionable report…"
  if Rscript --vanilla scripts/render_one_actionable.R "$sid" "$tt" \
       > "$render_log" 2>&1; then
    src_html="reports/${sid}/08-report/clinical_actionable.html"
    src_pdf="reports/${sid}/08-report/clinical_actionable.pdf"
    cp "$src_html" "$HANDOFF/${sid}_${tt}_actionable.html"
    [[ -f "$src_pdf" ]] && cp "$src_pdf" "$HANDOFF/${sid}_${tt}_actionable.pdf"
    STATUS+=("OK   $sid ($tt)  HTML=$(stat -f %z "$src_html")B  PDF=$([[ -f "$src_pdf" ]] && stat -f %z "$src_pdf" || echo MISS)B")
    echo "  [OK]  $sid"
  else
    STATUS+=("RND-FAIL $sid — log $render_log")
    echo "  [FAIL] render failed — see $render_log"
  fi
done

echo
echo "########################################################################"
echo "## Summary"
echo "########################################################################"
for s in "${STATUS[@]}"; do echo "  $s"; done
echo "  Handoff: $HANDOFF/"
echo "  Total elapsed: $(( $(date +%s) - TOTAL_T0 ))s"
