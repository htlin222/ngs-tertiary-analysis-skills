#!/usr/bin/env bash
# Extract unique citations + DOI metadata from each persisted OpenEvidence
# response file and write a slim JSON per sample at
# reports/.openevidence_cache/sample_<SID>.citations.json.
#
# Mapping of timestamp-named tool-results files → sample id is positional and
# was recorded when the queries were issued. If the OE responses are
# regenerated, regenerate the mapping and re-run this script.

set -euo pipefail
cd /Users/htlin/ngs-tertiary-analysis-skills

D=/Users/htlin/.claude/projects/-Users-htlin-ngs-tertiary-analysis-skills/6b73cd97-99ff-4248-8e5d-052f5d2704c3/tool-results
OUT=reports/.openevidence_cache
mkdir -p "$OUT"

declare -a MAP=(
  "M26-0165R 1778054603799"
  "M26-0234R 1778054678263"
  "M26-0284R 1778054756598"
  "M26-0320R 1778054874617"
  "M26-0400R 1778055006793"
  "M26-0401R 1778055115191"
  "M26-0410R 1778055164170"
  "M26-0414R 1778055233185"
)

for entry in "${MAP[@]}"; do
  sid="${entry% *}"
  ts="${entry#* }"
  src="$D/mcp-openevidence-oe_ask-${ts}.txt"
  if [[ ! -f "$src" ]]; then
    echo "  MISSING $sid → $src"; continue
  fi
  dst="$OUT/sample_${sid}.citations.json"

  # Pull unique citation_detail blocks; keep only fields we render.
  jq '[.. | objects | select(has("citation")) | .metadata.citation_detail]
        | map(select(. != null))
        | unique_by(.title // .href)
        | map({
            title: .title,
            authors: .authors_string,
            journal: (.journal_name // .repository),
            journal_short: (.journal_short_name // ""),
            citation_string: (.publication_info_string // ""),
            doi: (.doi // null),
            pmid: (.pmid // null),
            href: (.href // null),
            year: (.dt_published // null),
            publication_types: (.publication_types // [])
          })' "$src" > "$dst"

  count=$(jq 'length' "$dst")
  echo "  ${sid}: ${count} citations → ${dst}"
done
