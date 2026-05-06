#!/usr/bin/env bash
# scripts/run_batch7.sh — Run the TSO500 driver for all 8 patients in Batch 7
# and stage self-contained HTML reports under reports/_handoff_Batch_7/.

set -uo pipefail
cd /Users/htlin/ngs-tertiary-analysis-skills

BATCH=inputs/TSO500-HRD/Batch_7_20260505
HANDOFF=reports/_handoff_Batch_7_20260505
mkdir -p "$HANDOFF"

# Sample → OncoKB tumor type mapping (from 重點摘要.xlsx)
declare -a SAMPLES=(
  "M26-0165R UCEC"   # 子宮內膜癌
  "M26-0234R BRCA"   # 乳癌
  "M26-0284R COAD"   # 大腸癌
  "M26-0320R BRCA"   # 乳癌
  "M26-0400R UCEC"   # 子宮內膜癌
  "M26-0401R COAD"   # 乙狀結腸 (sigmoid colon)
  "M26-0410R UCEC"   # endometrial ca
  "M26-0414R UCEC"   # endometrial ca
)

declare -a STATUS=()
TOTAL_T0=$(date +%s)

for entry in "${SAMPLES[@]}"; do
  sid="${entry% *}"
  tt="${entry#* }"
  tsv="$BATCH/${sid}_CombinedVariantOutput.tsv"
  log="$HANDOFF/${sid}.log"

  echo
  echo "########################################################################"
  echo "## $sid ($tt)"
  echo "########################################################################"

  if Rscript --vanilla scripts/run_one_sample_tso500.R "$sid" "$tsv" "$tt" \
       > "$log" 2>&1; then
    html_src="reports/${sid}/08-report/clinical_report.html"
    if [[ -f "$html_src" ]]; then
      cp "$html_src" "$HANDOFF/${sid}_${tt}_clinical_report.html"
      STATUS+=("OK   $sid ($tt)")
      echo "  [OK]  $sid"
    else
      STATUS+=("MISS $sid ($tt) — log: $log")
      echo "  [MISS] HTML not produced — see $log"
    fi
  else
    STATUS+=("FAIL $sid ($tt) — log: $log")
    echo "  [FAIL] driver exited non-zero — see $log"
  fi
done

echo
echo "########################################################################"
echo "## Summary"
echo "########################################################################"
for s in "${STATUS[@]}"; do echo "  $s"; done
echo "  Handoff: $HANDOFF/"
echo "  Total elapsed: $(( $(date +%s) - TOTAL_T0 ))s"
