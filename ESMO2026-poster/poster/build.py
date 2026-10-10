#!/usr/bin/env python3
"""Build the ESMO 2026 posters: src/*.html -> out/*.html + true-size PDF + preview PNG.

    Rscript --vanilla ESMO2026-poster/poster/extract_cohort.R   # once, cohort data
    uv run --no-project ESMO2026-poster/poster/build.py

Tokens in src/*.html: {{css}}, {{fig:<name>}}, {{qr:<url>}}.
Figures are drawn from repo data; the numbers quoted in the poster text are
checked against the same data in check_numbers() so text and figures cannot drift.
The same single-page PDF serves as paper poster (print at 100 %) and e-Poster.
"""

import base64
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import fig_agent
import fig_cohort

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
MM_PT = 72 / 25.4
MAC_CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
CHROME = (
    MAC_CHROME
    if Path(MAC_CHROME).exists()
    else shutil.which("google-chrome") or shutil.which("chromium")
)

FIGS = {
    f.__name__: f
    for f in (
        fig_agent.evidence_states,
        fig_agent.architecture,
        fig_agent.transition_matrix,
        fig_agent.strata_bars,
        fig_agent.repro_grid,
        fig_cohort.tumor_bars,
        fig_cohort.landscape,
        fig_cohort.match_by_tumor,
        fig_cohort.onlabel_table,
    )
}

# Reports how much vertical room each [data-fit] block has left (mm); negative = overflow.
FIT_PROBE = """<script>
addEventListener('load', () => {
  const mm = px => Math.round(px / 96 * 25.4);
  document.documentElement.dataset.fit = [...document.querySelectorAll('[data-fit]')].map(el => {
    const used = [...el.children].reduce((s, c) => s + c.getBoundingClientRect().height, 0);
    const wide = el.scrollWidth > el.clientWidth + 2 ? ' WIDE' : '';
    return el.dataset.fit + '=' + mm(el.clientHeight - used) + wide;
  }).join(';');
});
</script>"""


def check_numbers():
    """Every number hand-typed in src/*.html, asserted against the data."""
    a = fig_agent.metrics()
    assert (a["n"], a["valid"], a["agree"]) == (112, 112, 90), a
    assert a["strata"] == {
        "oncokb": [70, 0],
        "same": [18, 1],
        "cross": [3, 1],
        "none": [21, 20],
    }
    assert a["transitions"]["Tier III VUS", "Tier II Level D"] == 20
    assert a["transitions"]["Tier I Level A", "Tier II Level C"] == 2
    assert (a["kappa_within"], a["kappa_cross"]) == (1.0, 0.717)

    c = fig_cohort.metrics()
    assert (c["n"], c["types"], c["on"], c["off_any"]) == (62, 11, 21, 38), c
    assert Counter(c["tmb_high"]) == {
        "BRCA": 3,
        "UCEC": 2,
        "PAAD": 1,
        "NSCLC": 1,
        "CESC": 1,
    }
    assert Counter(c["hrd_pos"]) == {"OV": 3, "BRCA": 2, "COAD": 1}
    assert Counter(c["msi_high"]) == {"UCEC": 1, "CHOL": 1}
    assert c["resistance_pts"] == 6
    assert Counter(t for t, _, _ in c["resistance"]) == {"COAD": 5, "NSCLC": 1}
    assert ("NSCLC", "EGFR", "T790M") in c["resistance"]
    assert {g for t, g, _ in c["resistance"] if t == "COAD"} == {"KRAS", "NRAS"}
    assert c["match_by_group"]["Ovarian"] == {"none": 7}
    assert dict(c["group_n"]) == {
        "Endometrial": 16,
        "Breast": 14,
        "Colon": 9,
        "NSCLC": 8,
        "Ovarian": 7,
        "Other": 8,
    }
    assert [c["gene_n"][g] for g in fig_cohort.GENES[:7]] == [26, 16, 16, 11, 9, 7, 6]
    brca1 = [
        p
        for p in fig_cohort.patients()
        if any(x["gene"] == "BRCA1" for x in p["on_label"])
    ]
    assert len(brca1) == 1 and brca1[0]["tumor"] == "BRCA" and brca1[0]["gis"] >= 42


def qr(url):
    out = subprocess.run(
        ["qrencode", "-t", "SVG", "-m", "0", "-l", "M", "-o", "-", url],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    out = out[out.index("<svg") :]
    return re.sub(r'<svg width="[^"]*" height="[^"]*"', "<svg", out, count=1)


def data_uri(path):
    mime = {"png": "image/png", "svg": "image/svg+xml", "jpg": "image/jpeg"}[
        path.suffix[1:]
    ]
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode()


def render(src):
    html = src.read_text()
    html = html.replace("{{css}}", (HERE / "poster.css").read_text())
    html = re.sub(r"\{\{fig:(\w+)\}\}", lambda m: FIGS[m.group(1)](), html)
    html = re.sub(r"\{\{qr:([^}]+)\}\}", lambda m: qr(m.group(1)), html)
    html = re.sub(r"\{\{img:([^}]+)\}\}", lambda m: data_uri(HERE / m.group(1)), html)
    return html.replace("</body>", FIT_PROBE + "\n</body>")


# Paper poster is 190 x 110 cm (ESMO Poster Presenters Instructions); the e-Poster
# must be a landscape 16:9 PDF, so the same layout is re-paged on a wider sheet.
EPOSTER_W = 1100 * 16 / 9
VARIANTS = [
    ("", 1900, ""),
    (
        "_ePoster_16x9",
        EPOSTER_W,
        f"@page {{ size: {EPOSTER_W:.2f}mm 1100mm; }} body, .poster {{ width: {EPOSTER_W:.2f}mm; }} "
        f".poster {{ padding-left: {45 + (EPOSTER_W - 1900) / 2:.2f}mm; "
        f"padding-right: {45 + (EPOSTER_W - 1900) / 2:.2f}mm; }}",
    ),
]


def export(html, html_path, pdf, width_mm):
    """Print one HTML variant to PDF and verify page count, size and column fit."""
    html_path.write_text(html)
    base = [CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer"]
    subprocess.run(
        base + [f"--print-to-pdf={pdf}", html_path.as_uri()],
        check=True,
        capture_output=True,
    )
    dom = subprocess.run(
        base + ["--dump-dom", html_path.as_uri()],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    fit = re.search(r'data-fit="([^"]*)"', dom[:2000])
    info = subprocess.run(
        ["pdfinfo", str(pdf)], check=True, capture_output=True, text=True
    ).stdout
    pages = int(re.search(r"Pages:\s+(\d+)", info).group(1))
    w, h = map(float, re.search(r"Page size:\s+([\d.]+) x ([\d.]+)", info).groups())
    subprocess.run(
        [
            "pdftoppm",
            "-r",
            "20",
            "-png",
            "-singlefile",
            str(pdf),
            str(html_path.with_name(html_path.stem + "_preview")),
        ],
        check=True,
    )
    problems = []
    if pages != 1:
        problems.append(f"{pages} pages (must be 1)")
    if abs(w - width_mm * MM_PT) > 2 or abs(h - 1100 * MM_PT) > 2:
        problems.append(
            f"page is {w / MM_PT:.0f} x {h / MM_PT:.0f} mm "
            f"(must be {width_mm:.0f} x 1100)"
        )
    room = dict(kv.split("=") for kv in fit.group(1).split(";")) if fit else {}
    problems += [
        f"'{k}' overflows ({v} mm)"
        for k, v in room.items()
        if v.startswith("-") or v.endswith("WIDE")
    ]
    print(
        f"{pdf.name}: {pages} page, {w / MM_PT:.0f} x {h / MM_PT:.0f} mm "
        f"({w / h:.3f}:1), {pdf.stat().st_size / 1e6:.1f} MB"
    )
    print("  free room (mm):", ", ".join(f"{k} {v}" for k, v in room.items()))
    for t in sorted(set(re.findall(r"\[\[[^\]]+\]\]", html))):
        print("  TO FILL before printing:", t)
    for prob in problems:
        print("  PROBLEM:", prob)
    return not problems


def build(src):
    html = render(src)
    # PDF name comes from <meta name="pdf-name">; the ESMO portal rejects names over 100 chars
    slug = re.search(r'<meta name="pdf-name" content="([^"]+)"', html).group(1)
    ok = True
    for suffix, width_mm, css in VARIANTS:
        pdf_name = f"{slug}{suffix}.pdf"
        assert len(pdf_name) <= 100 and re.fullmatch(r"[\w.-]+", pdf_name), pdf_name
        variant = html.replace("</style>", css + "</style>", 1) if css else html
        ok &= export(
            variant, OUT / f"{src.stem}{suffix}.html", OUT / pdf_name, width_mm
        )
    return ok


if __name__ == "__main__":
    check_numbers()
    OUT.mkdir(exist_ok=True)
    ok = [build(src) for src in sorted((HERE / "src").glob("*.html"))]
    sys.exit(0 if all(ok) else 1)
