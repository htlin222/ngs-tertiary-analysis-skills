# run_standalone_report.R — Standalone HTML report generator
# Generates ESMO-compliant clinical report directly in R (no Quarto dependency)

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(gt)
  library(logger)
  library(htmltools)
  library(purrr)
  library(stringr)
  library(fs)
})

source(here::here("R/esmo_helpers.R"))

generate_standalone_report <- function(sample_id, config, qc, variants, cnv,
                                       fusions, tmb, msi, hrd, escat,
                                       literature, output_path) {
  log_info("Generating standalone HTML report for {sample_id}")

  report_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  tumor_type <- config$sample$tumor_type %||% "Unknown"

  # ── CSS Styles ──────────────────────────────────────────────────────────
  css <- '
  <style>
    * { box-sizing: border-box; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
           max-width: 1100px; margin: 0 auto; padding: 20px; color: #333; line-height: 1.6; }
    h1 { color: #1a237e; border-bottom: 3px solid #1a237e; padding-bottom: 10px; }
    h2 { color: #283593; border-bottom: 2px solid #e8eaf6; padding-bottom: 8px; margin-top: 40px; }
    h3 { color: #3949ab; }
    .header-box { background: linear-gradient(135deg, #1a237e, #283593); color: white; padding: 30px;
                  border-radius: 8px; margin-bottom: 30px; }
    .header-box h1 { color: white; border: none; margin: 0; }
    .header-box p { margin: 5px 0; opacity: 0.9; }
    .badge { display: inline-block; padding: 3px 10px; border-radius: 4px; font-weight: bold;
             font-size: 0.85em; color: white; }
    .badge-pass { background: #2e7d32; }
    .badge-fail { background: #c62828; }
    .badge-high { background: #c62828; }
    .badge-intermediate { background: #f57c00; }
    .badge-low { background: #2e7d32; }
    .badge-positive { background: #c62828; }
    .badge-negative { background: #2e7d32; }
    .badge-mss { background: #2e7d32; }
    .badge-msi-h { background: #c62828; }
    .escat-I { background: #d32f2f; color: white; }
    .escat-II { background: #f57c00; color: white; }
    .escat-III { background: #fbc02d; color: #333; }
    .escat-IV { background: #7cb342; color: white; }
    .escat-X { background: #9e9e9e; color: white; }
    .card-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 20px 0; }
    .card { background: #f5f5f5; border-radius: 8px; padding: 20px; text-align: center;
            border-left: 4px solid #1a237e; }
    .card .value { font-size: 2em; font-weight: bold; color: #1a237e; }
    .card .label { font-size: 0.9em; color: #666; margin-top: 5px; }
    table { width: 100%; border-collapse: collapse; margin: 15px 0; }
    th { background: #e8eaf6; color: #1a237e; padding: 10px; text-align: left; font-weight: 600; }
    td { padding: 8px 10px; border-bottom: 1px solid #e0e0e0; }
    tr:nth-child(even) { background: #fafafa; }
    tr:hover { background: #f0f0ff; }
    .pathogenic { background-color: #ffebee !important; }
    .actionable-summary { background: #e3f2fd; border-radius: 8px; padding: 20px; margin: 20px 0;
                          border-left: 4px solid #1565c0; }
    .narrative { background: #fff8e1; border-radius: 8px; padding: 15px; margin: 15px 0;
                 border-left: 4px solid #f9a825; }
    .warning-box { background: #fff3e0; border: 1px solid #ff9800; border-radius: 8px; padding: 15px;
                   margin: 15px 0; }
    .disclaimer { background: #f5f5f5; border-radius: 8px; padding: 15px; margin-top: 40px;
                  font-size: 0.85em; color: #666; }
    .toc { background: #f5f5f5; border-radius: 8px; padding: 15px 25px; margin: 20px 0; }
    .toc a { color: #1a237e; text-decoration: none; }
    .toc a:hover { text-decoration: underline; }
    @media print { body { max-width: 100%; } .header-box { break-inside: avoid; } }
  </style>'

  # ── Helper functions ────────────────────────────────────────────────────
  escat_badge <- function(tier) {
    css_class <- glue("escat-{str_replace(tier, '-', '')}")
    glue('<span class="badge {css_class}">ESCAT {tier}</span>')
  }

  status_badge <- function(value, type = "pass") {
    css_class <- glue("badge-{tolower(value)}")
    glue('<span class="badge {css_class}">{value}</span>')
  }

  # ── Build HTML sections ─────────────────────────────────────────────────

  # 1. Header
  header <- glue('
  <div class="header-box">
    <h1>Comprehensive Genomic Profiling Report</h1>
    <p><strong>Panel:</strong> TruSight Oncology 500 (TSO500) | <strong>ESMO 2024 Compliant</strong></p>
    <p><strong>Report Date:</strong> {report_date}</p>
  </div>')

  # 2. Table of Contents
  toc <- '
  <div class="toc">
    <h3>Contents</h3>
    <ol>
      <li><a href="#patient">Patient &amp; Sample Information</a></li>
      <li><a href="#qc">Quality Control</a></li>
      <li><a href="#variants">Somatic Variants</a></li>
      <li><a href="#cnv">Copy Number Alterations</a></li>
      <li><a href="#fusions">Gene Fusions</a></li>
      <li><a href="#biomarkers">Genomic Biomarkers</a></li>
      <li><a href="#actionability">Clinical Actionability (ESCAT)</a></li>
      <li><a href="#literature">Literature Evidence</a></li>
      <li><a href="#coverage">Coverage Gaps</a></li>
    </ol>
  </div>'

  # 3. Patient Info
  patient_section <- glue('
  <h2 id="patient">1. Patient &amp; Sample Information</h2>
  <table>
    <tr><th>Field</th><th>Value</th></tr>
    <tr><td>Sample ID</td><td><strong>{sample_id}</strong></td></tr>
    <tr><td>Tumor Type</td><td>{tumor_type} (High-Grade Serous Ovarian Carcinoma)</td></tr>
    <tr><td>Panel</td><td>TruSight Oncology 500 (TSO500)</td></tr>
    <tr><td>Specimen Type</td><td>FFPE Tumor Tissue</td></tr>
    <tr><td>Sequencing Platform</td><td>Illumina NovaSeq X</td></tr>
    <tr><td>Institution</td><td>{config$report$institution %||% "Research Institution"}</td></tr>
    <tr><td>Report Date</td><td>{report_date}</td></tr>
    <tr><td>Data Source</td><td>EGA EGAD50000001451 (Ovarian Cancer TSO500 Cohort)</td></tr>
  </table>')

  # 4. QC Section
  qc_pass <- if (qc$pass) status_badge("PASS", "pass") else status_badge("FAIL", "fail")
  qc_section <- glue('
  <h2 id="qc">2. Quality Control</h2>
  <p>Overall QC Status: {qc_pass}</p>
  <table>
    <tr><th>Metric</th><th>Value</th><th>Threshold</th><th>Status</th></tr>
    <tr><td>Total Reads</td><td>{format(qc$summary$total_reads, big.mark=",")}</td><td>-</td><td>-</td></tr>
    <tr><td>Mapped Reads</td><td>{format(qc$summary$mapped_reads, big.mark=",")} ({round(qc$summary$mapping_rate*100,1)}%)</td><td>&gt;95%</td><td>{if(qc$summary$mapping_rate>0.95) "PASS" else "FAIL"}</td></tr>
    <tr><td>Mean Coverage</td><td>{qc$summary$mean_coverage}x</td><td>&gt;200x</td><td>{if(qc$summary$mean_coverage>200) "PASS" else "FAIL"}</td></tr>
    <tr><td>On-Target Rate</td><td>{round(qc$summary$on_target_rate*100,1)}%</td><td>&gt;80%</td><td>{if(qc$summary$on_target_rate>0.80) "PASS" else "FAIL"}</td></tr>
    <tr><td>Duplicate Rate</td><td>{round(qc$summary$duplicate_rate*100,1)}%</td><td>&lt;30%</td><td>{if(qc$summary$duplicate_rate<0.30) "PASS" else "FAIL"}</td></tr>
    <tr><td>Estimated Tumor Purity</td><td>{round(qc$summary$tumor_purity*100,0)}%</td><td>&gt;10%</td><td>{if(qc$summary$tumor_purity>0.10) "PASS" else "FAIL"}</td></tr>
  </table>')

  # 4b. Annotation Methodology
  annotation_section <- '
  <h3>2b. Annotation Methodology</h3>
  <table>
    <tr><th>Component</th><th>Tool / Database</th><th>Version / Details</th></tr>
    <tr><td>Variant Calling</td><td>GATK Mutect2</td><td>v4.6.1.0, tumor-only mode</td></tr>
    <tr><td>Variant Annotation</td><td>Ensembl VEP</td><td>GRCh38, canonical transcripts, --everything flag</td></tr>
    <tr><td>Transcript Database</td><td>Ensembl / GENCODE</td><td>Canonical transcripts per gene (one consequence per variant via --pick)</td></tr>
    <tr><td>Functional Prediction</td><td>SIFT + PolyPhen-2</td><td>Protein impact prediction for missense variants</td></tr>
    <tr><td>Population Frequency</td><td>gnomAD</td><td>v4.0, germline filtering at AF &gt;1%</td></tr>
    <tr><td>Clinical Databases</td><td>ClinVar + COSMIC</td><td>Pathogenicity and somatic mutation databases</td></tr>
    <tr><td>Clinical Actionability</td><td>OncoKB API</td><td>Therapeutic levels mapped to ESCAT tiers</td></tr>
    <tr><td>Nomenclature</td><td>HGVS</td><td>Protein (p.) and coding (c.) notation per <a href="https://varnomen.hgvs.org/">HGVS recommendations</a></td></tr>
    <tr><td>CNV Detection</td><td>CNVkit</td><td>Hybrid capture mode for panel data</td></tr>
    <tr><td>Fusion Detection</td><td>Manta</td><td>Structural variant caller for gene fusions</td></tr>
  </table>
  <p style="font-size:0.85em;color:#666"><em>Note: Multiple transcripts may arise from the same genomic locus
  (coding, noncoding, splice variants). This report uses Ensembl canonical transcripts.
  Where library fragments do not cover all exons, coverage gaps are reported in Section 9.</em></p>'

  # 5. Somatic Variants — now with transcript, depth, SIFT/PolyPhen, and confidence ranking
  variant_rows <- map_chr(seq_len(nrow(variants)), function(i) {
    v <- variants[i, ]
    if (v$consequence == "amplification") return("")
    row_class <- if (v$pathogenicity %in% c("pathogenic", "likely_pathogenic")) ' class="pathogenic"' else ""
    vaf_text <- if (!is.na(v$vaf)) glue("{round(v$vaf*100,1)}%") else "N/A"
    depth_text <- if (!is.na(v$alt_depth)) glue("{v$alt_depth}/{v$total_depth}") else "N/A"
    transcript <- if ("transcript_id" %in% names(v) && !is.na(v$transcript_id)) v$transcript_id else "-"
    sift <- if ("sift_prediction" %in% names(v) && !is.na(v$sift_prediction)) v$sift_prediction else "-"
    polyphen <- if ("polyphen_prediction" %in% names(v) && !is.na(v$polyphen_prediction)) v$polyphen_prediction else "-"
    glue('<tr{row_class}>',
         '<td><strong>{v$gene}</strong><br><span style="font-size:0.75em;color:#666">{transcript}</span></td>',
         '<td>{v$hgvsp}</td><td>{v$hgvsc %||% "-"}</td>',
         '<td>{vaf_text}<br><span style="font-size:0.75em;color:#666">{depth_text}</span></td>',
         '<td>{v$consequence}</td><td>{v$impact}</td>',
         '<td>{sift}<br>{polyphen}</td>',
         '<td>{v$clinvar_significance %||% "-"}</td>',
         '<td>{v$pathogenicity}</td></tr>')
  })
  variant_rows <- paste(variant_rows[nchar(variant_rows) > 0], collapse = "\n    ")

  # Mutation confidence ranking section (Table 3 per best practices)
  has_confidence <- "confidence" %in% names(variants)
  confidence_section <- ""
  if (has_confidence) {
    conf_variants <- variants |> filter(consequence != "amplification", !is.na(confidence))
    conf_rows <- map_chr(seq_len(nrow(conf_variants)), function(i) {
      cv <- conf_variants[i, ]
      badge_class <- switch(cv$confidence,
        "High confidence" = "badge-pass",
        "Questionable but potentially impactful" = "badge-high",
        "badge-intermediate"
      )
      glue('<tr><td><strong>{cv$gene}</strong> {cv$hgvsp}</td>',
           '<td><span class="badge {badge_class}">{cv$confidence}</span></td>',
           '<td style="font-size:0.85em">{cv$confidence_reason}</td></tr>')
    })
    confidence_section <- glue('
    <h3>3b. Mutation Confidence Ranking</h3>
    <p>Variants are ranked by call confidence based on VAF relative to detection threshold
    ({round(config$variant_calling$min_vaf * 100, 0)}%), read depth, and database concordance
    (following HGVS nomenclature, <a href="https://varnomen.hgvs.org/">varnomen.hgvs.org</a>).</p>
    <table>
      <tr><th>Variant</th><th>Confidence</th><th>Rationale</th></tr>
      {paste(conf_rows, collapse = "\\n      ")}
    </table>
    <p style="font-size:0.85em;color:#666"><em>Questionable variants with high clinical impact should be confirmed
    by orthogonal methods (e.g., Sanger sequencing, ddPCR) before clinical action.</em></p>')
  }

  variants_section <- glue('
  <h2 id="variants">3. Somatic Variants</h2>
  <p>{nrow(variants |> filter(consequence != "amplification"))} somatic variants detected.
  Annotated using Ensembl VEP (GRCh38) with canonical transcripts.
  Variants reported using <a href="https://varnomen.hgvs.org/">HGVS nomenclature</a>.</p>
  <table>
    <tr><th>Gene<br><span style="font-size:0.75em">Transcript</span></th><th>Protein</th><th>Coding</th><th>VAF<br><span style="font-size:0.75em">Alt/Total</span></th><th>Consequence</th><th>Impact</th><th>SIFT<br>PolyPhen</th><th>ClinVar</th><th>Classification</th></tr>
    {variant_rows}
  </table>
  {confidence_section}')

  # 6. CNV Section
  cnv_rows <- map_chr(seq_len(nrow(cnv)), function(i) {
    c <- cnv[i, ]
    color <- if (c$type == "AMPLIFICATION") "color:#c62828;font-weight:bold" else "color:#1565c0;font-weight:bold"
    glue('<tr><td><strong>{c$gene}</strong></td><td style="{color}">{c$type}</td>',
         '<td>{round(c$log2_ratio, 2)}</td><td>{c$copy_number}</td>',
         '<td>{c$chromosome}:{format(c$start, big.mark=",")}-{format(c$end, big.mark=",")}</td></tr>')
  })

  cnv_section <- glue('
  <h2 id="cnv">4. Copy Number Alterations</h2>
  <p>{nrow(cnv)} clinically significant copy number alterations detected.</p>
  <table>
    <tr><th>Gene</th><th>Type</th><th>Log2 Ratio</th><th>Copy Number</th><th>Coordinates</th></tr>
    {paste(cnv_rows, collapse = "\n    ")}
  </table>')

  # 7. Fusions
  fusion_rows <- map_chr(seq_len(nrow(fusions)), function(i) {
    f <- fusions[i, ]
    known_badge <- if (f$known_fusion) '<span class="badge badge-high">Known</span>' else "Novel"
    glue('<tr><td><strong>{f$fusion_name}</strong></td><td>exon {f$exon_a}::exon {f$exon_b}</td>',
         '<td>{f$chr_a}:{format(f$pos_a, big.mark=",")} — {f$chr_b}:{format(f$pos_b, big.mark=",")}</td>',
         '<td>{f$supporting_reads}</td><td>{f$fusion_type}</td><td>{known_badge}</td></tr>')
  })

  fusions_section <- glue('
  <h2 id="fusions">5. Gene Fusions</h2>
  <p>{nrow(fusions)} gene fusion(s) detected.</p>
  <table>
    <tr><th>Fusion</th><th>Exons</th><th>Breakpoints</th><th>Supporting Reads</th><th>Type</th><th>Status</th></tr>
    {paste(fusion_rows, collapse = "\n    ")}
  </table>')

  # 8. Biomarkers
  tmb_badge <- if (tmb$tmb_class == "TMB-High") "badge-high" else if (tmb$tmb_class == "TMB-Intermediate") "badge-intermediate" else "badge-low"
  msi_badge <- if (msi$msi_status == "MSI-H") "badge-msi-h" else "badge-mss"
  hrd_badge <- if (hrd$hrd_status == "HRD-positive") "badge-positive" else "badge-negative"

  biomarkers_section <- glue('
  <h2 id="biomarkers">6. Genomic Biomarkers</h2>
  <div class="card-grid">
    <div class="card">
      <div class="value">{tmb$tmb_score}</div>
      <div class="label">TMB (mutations/Mb)</div>
      <div><span class="badge {tmb_badge}">{tmb$tmb_class}</span></div>
      <div class="label" style="margin-top:8px">{tmb$variant_count} variants / {tmb$panel_size_mb} Mb panel</div>
    </div>
    <div class="card">
      <div class="value">{msi$msi_score}%</div>
      <div class="label">MSI Score</div>
      <div><span class="badge {msi_badge}">{msi$msi_status}</span></div>
      <div class="label" style="margin-top:8px">{msi$unstable_sites}/{msi$total_sites} unstable sites</div>
    </div>
    <div class="card">
      <div class="value">{hrd$hrd_score}</div>
      <div class="label">HRD Score</div>
      <div><span class="badge {hrd_badge}">{hrd$hrd_status}</span></div>
      <div class="label" style="margin-top:8px">LOH:{hrd$loh_score} TAI:{hrd$tai_score} LST:{hrd$lst_score}</div>
    </div>
  </div>')

  # 9. Clinical Actionability
  actionable_count <- nrow(escat |> filter(escat_tier %in% c("I", "II", "III")))
  escat_rows <- map_chr(seq_len(nrow(escat)), function(i) {
    e <- escat[i, ]
    badge <- escat_badge(e$escat_tier)
    glue('<tr><td>{badge}</td><td><strong>{e$gene}</strong></td><td>{e$alteration}</td>',
         '<td>{e$type}</td><td>{e$oncogenic}</td><td>{e$sensitive_level %||% "-"}</td>',
         '<td>{e$escat_description}</td></tr>')
  })

  actionability_section <- glue('
  <h2 id="actionability">7. Clinical Actionability (ESCAT Classification)</h2>
  <div class="actionable-summary">
    <strong>{actionable_count} actionable alteration(s) identified</strong> (ESCAT Tier I-III)
    across {nrow(escat)} total classified alterations.
  </div>
  <table>
    <tr><th>ESCAT</th><th>Gene</th><th>Alteration</th><th>Type</th><th>Oncogenic</th><th>Level</th><th>Description</th></tr>
    {paste(escat_rows, collapse = "\n    ")}
  </table>')

  # 10. Literature
  narrative_html <- map_chr(names(literature$narratives), function(key) {
    text <- literature$narratives[[key]]
    gene_name <- str_extract(key, "^[^_]+")
    glue('<div class="narrative"><h4>{gene_name}</h4><p>{text}</p></div>')
  })

  pubmed_count <- if (is.data.frame(literature$pubmed)) nrow(literature$pubmed) else 0
  scopus_count <- if (is.data.frame(literature$scopus)) nrow(literature$scopus) else 0

  literature_section <- glue('
  <h2 id="literature">8. Literature Evidence</h2>
  <p>Literature search yielded <strong>{pubmed_count}</strong> PubMed and <strong>{scopus_count}</strong> Scopus articles for key alterations.</p>
  {paste(narrative_html, collapse = "\n  ")}')

  # 11. Coverage Gaps
  gap_genes <- qc$per_gene_coverage |> filter(mean_coverage < 200)
  if (nrow(gap_genes) > 0) {
    gap_rows <- map_chr(seq_len(nrow(gap_genes)), function(i) {
      g <- gap_genes[i, ]
      glue('<tr><td><strong>{g$gene}</strong></td><td>{round(g$mean_coverage, 0)}x</td><td>200x</td><td><span class="badge badge-fail">BELOW THRESHOLD</span></td></tr>')
    })
    coverage_section <- glue('
    <h2 id="coverage">9. Coverage Gaps</h2>
    <div class="warning-box">
      <strong>Warning:</strong> {nrow(gap_genes)} gene(s) did not meet the minimum 200x coverage threshold.
      Results for these genes should be interpreted with caution.
    </div>
    <table>
      <tr><th>Gene</th><th>Mean Coverage</th><th>Required</th><th>Status</th></tr>
      {paste(gap_rows, collapse = "\n      ")}
    </table>')
  } else {
    coverage_section <- '
    <h2 id="coverage">9. Coverage Gaps</h2>
    <p>All genes met the minimum 200x coverage threshold. No gaps to report.</p>'
  }

  # 12. Disclaimer
  disclaimer <- '
  <div class="disclaimer">
    <h3>Disclaimer</h3>
    <p>This report is generated for research purposes following ESMO 2024 guidelines for NGS clinical reporting.
    Clinical decisions should not be based solely on this report without confirmation by a certified molecular pathology laboratory.</p>
    <p>OncoKB annotations are provided under academic license. ESCAT classifications follow the ESMO Scale for Clinical
    Actionability of molecular Targets (Mateo et al., Ann Oncol, 2018).</p>
    <p><strong>Pipeline:</strong> ngs-tertiary-analysis v0.1.0 | <strong>Reference:</strong> GRCh38 | <strong>Panel:</strong> TSO500 (523 genes, 1.94 Mb coding)</p>
    <p><em>Report generated: {report_date}</em></p>
  </div>'

  # ── Assemble full HTML ──────────────────────────────────────────────────
  html <- glue('<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>CGP Report — {sample_id}</title>
  {css}
</head>
<body>
  {header}
  {toc}
  {patient_section}
  {qc_section}
  {annotation_section}
  {variants_section}
  {cnv_section}
  {fusions_section}
  {biomarkers_section}
  {actionability_section}
  {literature_section}
  {coverage_section}
  {disclaimer}
</body>
</html>')

  writeLines(html, output_path)
  log_info("Standalone report written to {output_path}")
  output_path
}
