"""Figures for poster #6004 (pan-cancer cohort), built from reports/_poster_cohort.json.

Run extract_cohort.R first; the JSON is de-identified and gitignored.
"""

import json
from collections import Counter, defaultdict

from svgkit import (
    ACCENT,
    AXIS,
    BLUE,
    CELL,
    COL,
    FS,
    FS_S,
    INK,
    INK2,
    MUTED,
    ROOT,
    WARM,
    WIDE,
    hbar,
    line,
    rect,
    svg,
    swatch,
    text,
)

COHORT = ROOT / "reports" / "_poster_cohort.json"
GROUP = {
    "UCEC": "Endometrial",
    "BRCA": "Breast",
    "COAD": "Colon",
    "NSCLC": "NSCLC",
    "OV": "Ovarian",
}
ORDER = ["Endometrial", "Breast", "Colon", "NSCLC", "Ovarian", "Other"]
TUMOR_NAME = {
    "UCEC": "endometrial",
    "BRCA": "breast",
    "COAD": "colon",
    "NSCLC": "NSCLC",
    "OV": "ovarian",
    "PAAD": "pancreatic",
    "READ": "rectal",
    "STAD": "gastric",
    "CHOL": "biliary",
    "CESC": "cervical",
    "SARC": "sarcoma",
    "CANCER": "unspecified",
}
# Rows of the landscape: the seven most recurrent drivers reported in the
# abstract, then therapy-defining genes.
GENES = [
    "TP53",
    "KRAS",
    "PIK3CA",
    "PTEN",
    "ARID1A",
    "APC",
    "ESR1",
    "EGFR",
    "NRAS",
    "BRCA1/2",
    "ERBB2",
]
TMB_HIGH, GIS_POS = 10, 42
MATCH_FILL = {"on": ACCENT, "off": BLUE[250], "none": CELL}


def patients():
    with open(COHORT) as f:
        pts = json.load(f)
    for p in pts:
        p["group"] = GROUP.get(p["tumor"], "Other")
        p["match"] = "on" if p["n_on_label"] else "off" if p["n_off_label"] else "none"
        genes = set(p["oncogenic_genes"])
        if genes & {"BRCA1", "BRCA2"}:
            genes.add("BRCA1/2")
        p["genes"] = genes
    return pts


def metrics():
    pts = patients()
    by = defaultdict(Counter)
    for p in pts:
        by[p["group"]][p["match"]] += 1
    return {
        "n": len(pts),
        "types": len({p["tumor"] for p in pts} - {"CANCER"}),
        "tmb_high": [p["tumor"] for p in pts if p["tmb"] >= TMB_HIGH],
        "msi_high": [p["tumor"] for p in pts if p["msi_high"]],
        "hrd_pos": [p["tumor"] for p in pts if p["gis"] >= GIS_POS],
        "on": sum(p["match"] == "on" for p in pts),
        "off_any": sum(p["n_off_label"] > 0 for p in pts),
        "resistance": [
            (p["tumor"], a["gene"], a["alteration"])
            for p in pts
            for a in p["resistance"]
        ],
        "resistance_pts": sum(bool(p["resistance"]) for p in pts),
        "match_by_group": by,
        "group_n": Counter(p["group"] for p in pts),
        "gene_n": {g: sum(g in p["genes"] for p in pts) for g in GENES},
    }


def tumor_bars():
    n = metrics()["group_n"]
    left, pitch = 100, 29
    scale = (COL - left - 28) / max(n.values())
    b = []
    for i, g in enumerate(ORDER):
        y = i * pitch
        b.append(text(left - 8, y + 9, g, FS, "end", 500))
        b.append(hbar(left, y, n[g] * scale, 18, AXIS if g == "Other" else ACCENT))
        b.append(text(left + n[g] * scale + 5, y + 9, n[g], FS, weight=700))
    return svg(COL, len(ORDER) * pitch - 11, b, "Patients by cancer type")


def landscape():
    """Oncoprint-style cohort overview: one column per patient."""
    pts = patients()
    m = metrics()
    left, right, gap = 124, 62, 7
    pitch = (WIDE - left - right - gap * (len(ORDER) - 1)) / len(pts)
    cw = pitch - 1.5
    pts.sort(
        key=lambda p: (
            ORDER.index(p["group"]),
            ["on", "off", "none"].index(p["match"]),
            [g not in p["genes"] for g in GENES],
        )
    )
    xs, spans, x = [], {}, left
    for i, p in enumerate(pts):
        if i and p["group"] != pts[i - 1]["group"]:
            x += gap
        spans.setdefault(p["group"], [x, x])[1] = x + cw
        xs.append(x)
        x += pitch
    x_end = x - (pitch - cw)
    rx = x_end + 7
    b = []
    for g, (a, z) in spans.items():
        b.append(text((a + z) / 2, 6, f"{g} ({m['group_n'][g]})", FS_S, "middle", 600))
        b.append(line(a, 13.5, z, 13.5, INK, 0.7))

    def track(y, h, key, vmax, thr, name, unit, count):
        s = [
            text(left - 9, y + h / 2 - 5.4, name, FS, "end", 600),
            text(left - 9, y + h / 2 + 5.6, unit, FS_S, "end", fill=MUTED),
            line(left, y + h, x_end, y + h, AXIS, 0.5),
        ]
        for px, p in zip(xs, pts):
            bh = max(p[key] / vmax * h, 0.6)
            s.append(
                rect(px, y + h - bh, cw, bh, ACCENT if p[key] >= thr else AXIS, rx=0.8)
            )
        ty = y + h - thr / vmax * h
        s.append(line(left, ty, x_end, ty, INK, 0.45))
        s.append(text(rx, ty, f"≥{thr}", FS_S, fill=INK2))
        s.append(text(rx, y + h - 5, f"n = {count}", FS_S, weight=600))
        return "".join(s)

    b.append(track(22, 58, "tmb", 26, TMB_HIGH, "TMB", "mut/Mb", len(m["tmb_high"])))
    b.append(track(90, 58, "gis", 80, GIS_POS, "GIS", "HRD proxy", len(m["hrd_pos"])))

    def row(y, name, fills, count, italic=False):
        s = [
            text(
                left - 9,
                y + 6.2,
                name,
                FS,
                "end",
                500,
                style="font-style:italic" if italic else "",
            )
        ]
        s += [rect(px, y, cw, 12.4, f, rx=0.8) for px, f in zip(xs, fills)]
        s.append(text(rx, y + 6.2, count, FS_S, weight=600))
        return "".join(s)

    y = 160
    b.append(
        row(
            y,
            "MSI-High",
            [ACCENT if p["msi_high"] else CELL for p in pts],
            len(m["msi_high"]),
        )
    )
    b.append(
        row(
            y + 15,
            "OncoKB match",
            [MATCH_FILL[p["match"]] for p in pts],
            f"{m['on']} on-label",
        )
    )
    b.append(
        row(
            y + 30,
            "Resistance",
            [WARM if p["resistance"] else CELL for p in pts],
            m["resistance_pts"],
        )
    )
    y0 = y + 52
    for i, g in enumerate(GENES):
        n = m["gene_n"][g]
        b.append(
            row(
                y0 + i * 15,
                g,
                [INK if g in p["genes"] else CELL for p in pts],
                f"{n} ({n / len(pts):.0%})",
                italic=True,
            )
        )
    ly = y0 + len(GENES) * 15 + 9
    lx = left
    for label, fill in [
        ("On-label match (Level 1/2/3A)", ACCENT),
        ("Off-label only (Level 3B/4)", BLUE[250]),
        ("Resistance flag (R1)", WARM),
        ("Oncogenic / likely oncogenic alteration", INK),
    ]:
        s, lx = swatch(lx, ly, label, fill)
        b.append(s)
    return svg(WIDE, ly + 6, b, "Genomic landscape of the 62-patient cohort")


def match_by_tumor():
    """Best OncoKB therapy level per patient, by cancer type (counts)."""
    m = metrics()
    left, pitch = 132, 23
    scale = (COL - left - 6) / max(m["group_n"].values())
    b, lx = [], 0
    for label, k in [
        ("On-label", "on"),
        ("Off-label only", "off"),
        ("No match", "none"),
    ]:
        s, lx = swatch(
            lx, 5, label, MATCH_FILL[k], stroke=AXIS if k == "none" else "none"
        )
        b.append(s)
    for i, g in enumerate(ORDER):
        y = 20 + i * pitch
        b.append(text(left - 8, y + 7.5, f"{g} ({m['group_n'][g]})", FS, "end", 500))
        x = left
        segs = [(k, m["match_by_group"][g][k]) for k in ("on", "off", "none")]
        segs = [(k, n) for k, n in segs if n]
        for j, (k, n) in enumerate(segs):
            w = n * scale - 1
            last = j == len(segs) - 1
            b.append(
                hbar(x, y, w, 15, MATCH_FILL[k])
                if last
                else rect(x, y, w, 15, MATCH_FILL[k])
            )
            if k == "none":
                b.append(line(x, y + 15, x + w, y + 15, AXIS, 0.5))
            if w >= 12:
                b.append(
                    text(
                        x + w / 2,
                        y + 7.5,
                        n,
                        FS_S,
                        "middle",
                        700,
                        "#ffffff" if k == "on" else INK,
                    )
                )
            x += n * scale
    return svg(COL, 20 + len(ORDER) * pitch - 6, b, "Therapy match by cancer type")


def onlabel_table():
    """HTML table: patients with an on-label (Level 1/2/3A) alteration, by gene."""
    by_gene = defaultdict(lambda: {"pts": 0, "tumors": Counter(), "levels": set()})
    for p in patients():
        for gene in {a["gene"] for a in p["on_label"]}:
            g = by_gene[gene]
            g["pts"] += 1
            g["tumors"][TUMOR_NAME[p["tumor"]]] += 1
            g["levels"] |= {
                a["level"].replace("LEVEL_", "")
                for a in p["on_label"]
                if a["gene"] == gene
            }
    ranked = sorted(by_gene.items(), key=lambda kv: (-kv[1]["pts"], kv[0]))
    rows = []
    for gene, g in ranked:
        tumors = ", ".join(
            t if n == g["pts"] else f"{t} {n}" for t, n in g["tumors"].most_common()
        )
        rows.append(
            f"<tr><td><i>{gene}</i></td><td class='num'>{g['pts']}</td>"
            f"<td>{tumors}</td><td>{' / '.join(sorted(g['levels']))}</td></tr>"
        )
    return (
        "<table class='data'><thead><tr><th>Gene</th><th class='num'>Patients</th>"
        "<th>Cancer type</th><th>Level</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
    )
