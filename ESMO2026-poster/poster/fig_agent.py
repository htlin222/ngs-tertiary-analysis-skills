"""Figures for poster #4317 (model-agnostic agent), built from benchmarks/results."""

import csv
from collections import Counter

from svgkit import (
    ACCENT,
    AXIS,
    BLUE,
    CELL,
    COL,
    FS,
    FS_S,
    HAIR,
    INK,
    INK2,
    MUTED,
    ROOT,
    TIER_FILL,
    TIER_LABEL,
    TIERS,
    WARM,
    WARM_SOFT,
    WIDE,
    arrow_right,
    chip,
    hbar,
    line,
    on,
    rect,
    svg,
    swatch,
    text,
)

RES = ROOT / "benchmarks" / "results"
TUMOR = {"BRCA": "breast", "OVT": "ovarian", "COAD": "colon", "MEL": "melanoma"}


def _rows(name):
    with open(RES / name, newline="") as f:
        return list(csv.DictReader(f))


def strata(rows):
    """Partition the 112 cases by the evidence the deterministic rule saw."""

    def key(r):
        if not r["civic_amp_level"]:
            return "none" if r["baseline_tier"] == "Tier III VUS" else "oncokb"
        return "cross" if r["civic_evidence_in_different_tumor"] == "TRUE" else "same"

    out = {k: [0, 0] for k in ("oncokb", "same", "cross", "none")}
    for r in rows:
        k = key(r)
        out[k][0] += 1
        out[k][1] += r["agrees_with_baseline"] != "TRUE"
    return out


def metrics():
    rows = _rows("agent_combined_summary.csv")
    cons = {
        r["metric"]: float(r["value"]) for r in _rows("agent_consistency_summary.csv")
    }
    cross = {
        r["metric"]: float(r["value"]) for r in _rows("agent_cross_model_summary.csv")
    }
    return {
        "n": len(rows),
        "valid": sum(r["valid"] == "TRUE" for r in rows),
        "agree": sum(r["agrees_with_baseline"] == "TRUE" for r in rows),
        "transitions": Counter((r["baseline_tier"], r["agent_tier"]) for r in rows),
        "strata": strata(rows),
        "kappa_within": cons["fleiss_kappa_tier"],
        "kappa_cross": cross["fleiss_kappa_5_raters"],
    }


def _verdict(x, y, ok):
    if ok:
        return (
            f'<path d="M{x - 5:.1f},{y:.1f} l3.6,4 l6.6,-8" fill="none" '
            f'stroke="{INK2}" stroke-width="1.5" stroke-linecap="round" '
            f'stroke-linejoin="round"/>'
        )
    return line(x - 4, y - 4, x + 4, y + 4, WARM, 1.7) + line(
        x - 4, y + 4, x + 4, y - 4, WARM, 1.7
    )


def evidence_states():
    """Three evidence states that the higher-tier rule maps to the same call."""
    X_OK, X_CV, X_RULE, X_V = 128, 216, 336, 414
    b = [
        text(0, 6, "Evidence state", FS_S, fill=MUTED),
        text(X_OK, 6, "OncoKB", FS_S, fill=MUTED),
        text(X_CV, 6, "CiVIC", FS_S, fill=MUTED),
        text(X_RULE, 6, "Rule output", FS_S, fill=MUTED),
        line(0, 13, COL, 13, HAIR),
    ]
    hi = BLUE[700]
    states = [
        (["Both KBs agree"], ("Level 1", hi, {}), ("Tier I A", hi, {}), hi, True),
        (
            ["One KB has", "no assertion"],
            ("Level 1", hi, {}),
            ("none", "#ffffff", dict(stroke=AXIS, dash="2 1.6", ink=MUTED)),
            hi,
            True,
        ),
        (
            ["Evidence from a", "different tumor"],
            ("Level 3B", BLUE[400], {}),
            ("Tier I A*", "#ffffff", dict(stroke=hi, ink=hi)),
            WARM,
            False,
        ),
    ]
    for i, (label, okb, civic, rule_fill, ok) in enumerate(states):
        y = 20 + i * 33
        cy = y + 12
        for j, ln in enumerate(label):
            b.append(
                text(0, cy + (j - (len(label) - 1) / 2) * 11.4, ln, FS, weight=500)
            )
        b.append(chip(X_OK, y, 80, 24, okb[0], okb[1], **okb[2]))
        b.append(chip(X_CV, y, 96, 24, civic[0], civic[1], **civic[2]))
        b.append(line(316, cy, 328, cy, INK, 0.7))
        b.append(f'<path d="M327,{cy - 2.2} l4,2.2 l-4,2.2 z" fill="{INK}"/>')
        b.append(
            chip(X_RULE, y, 64, 24, "Tier I A", rule_fill, ink=INK if not ok else None)
        )
        b.append(_verdict(X_V, cy, ok))
    b.append(
        text(
            COL,
            124,
            "* asserted in another tumor type: a false Tier I call",
            FS_S,
            anchor="end",
            fill=INK2,
        )
    )
    return svg(COL, 130, b, "Three evidence states collapse to the same rule output")


def architecture():
    """Routing: concordant variants stay deterministic, discordant ones go to the agent."""

    def box(x, y, w, h, lines, fill="#ffffff", stroke=INK, ink=INK, sub=INK2):
        s = rect(x, y, w, h, fill, rx=2, stroke=stroke, sw=0.7)
        y0 = y + h / 2 - (len(lines) - 1) * 5.7
        for i, ln in enumerate(lines):
            s += text(
                x + w / 2,
                y0 + i * 11.4,
                ln,
                FS if i == 0 else FS_S,
                "middle",
                600 if i == 0 else 400,
                ink if i == 0 else sub,
            )
        return s

    mid = 56
    b = [
        box(0, mid - 24, 150, 48, ["TSO500 variants", "VEP annotation"]),
        arrow_right(150, mid, 172),
        box(172, mid - 24, 158, 48, ["OncoKB + CiVIC", "evidence retrieval"]),
        arrow_right(330, mid, 352),
        rect(352, mid - 18, 150, 36, "#ffffff", rx=18, stroke=INK, sw=0.7),
        text(427, mid, "KBs concordant?", FS, "middle", 600),
        line(502, mid, 518, mid, INK, 0.7),
        line(518, 20, 518, 85, INK, 0.7),
        arrow_right(518, 20, 548),
        arrow_right(518, 85, 548),
        text(512, 24, "concordant", FS_S, "end", fill=INK2),
        text(512, 90, "discordant", FS_S, "end", fill=INK2),
        box(548, 0, 176, 40, ["Deterministic rule", "Li 2017 decision tree"]),
        box(
            548,
            58,
            176,
            54,
            ["LLM agent", "model-agnostic · no tools", "JSON-schema-validated"],
            fill=BLUE[100],
            stroke=ACCENT,
        ),
        line(724, 20, 738, 20, INK, 0.7),
        line(724, 85, 738, 85, INK, 0.7),
        line(738, 20, 738, 85, INK, 0.7),
        arrow_right(738, mid, 756),
        box(
            756,
            mid - 26,
            WIDE - 756,
            52,
            ["Unified tier", "MTB report", "+ audit log"],
            fill=INK,
            ink="#ffffff",
            sub="#ffffff",
        ),
    ]
    return svg(WIDE, 112, b, "Pipeline architecture and agent routing")


def _seq(count, peak):
    steps = [100, 200, 250, 350, 400, 500, 550, 600]
    return BLUE[steps[round((count / peak) ** 0.5 * (len(steps) - 1))]]


def transition_matrix():
    """Baseline (rows) versus agent (columns) tier for all benchmark cases."""
    m = metrics()
    left, top, ch = 80, 34, 50
    cw = (COL - left) / 5
    peak = max(m["transitions"].values())
    b = [
        text(
            left + cw * 2.5,
            7,
            "AGENT TIER",
            FS_S,
            "middle",
            600,
            MUTED,
            "letter-spacing:.08em",
        ),
        f'<text transform="translate(7,{top + ch * 2.5}) rotate(-90)" font-size="{FS_S}" '
        f'text-anchor="middle" font-weight="600" fill="{MUTED}" dominant-baseline="central" '
        f'style="letter-spacing:.08em">DETERMINISTIC BASELINE</text>',
    ]
    for j, t in enumerate(TIERS):
        b.append(text(left + cw * (j + 0.5), 24, TIER_LABEL[t], FS, "middle", 500))
    for i, base in enumerate(TIERS):
        y = top + i * ch
        b.append(text(left - 6, y + ch / 2, TIER_LABEL[base], FS, "end", 500))
        for j, agent in enumerate(TIERS):
            n = m["transitions"].get((base, agent), 0)
            x = left + j * cw
            if not n:
                fill = CELL
            elif i == j:
                fill = _seq(n, peak)
            else:
                fill = WARM if n > 5 else WARM_SOFT
            b.append(rect(x + 0.9, y + 0.9, cw - 1.8, ch - 1.8, fill, rx=1.6))
            if n:
                b.append(
                    text(
                        x + cw / 2,
                        y + ch / 2,
                        n,
                        17,
                        "middle",
                        700,
                        INK if i != j else on(fill),
                    )
                )
    ly = top + 5 * ch + 13
    s, x = swatch(left, ly, f"Agreement ({m['agree']})", ACCENT)
    b.append(s)
    b.append(swatch(x, ly, f"Agent re-tiered ({m['n'] - m['agree']})", WARM)[0])
    return svg(COL, ly + 7, b, "Tier transition matrix, baseline versus agent")


def strata_bars():
    """Agreement and re-tiering within each evidence stratum."""
    st = metrics()["strata"]
    labels = [
        ("oncokb", "OncoKB only · baseline Tier I/II"),
        ("same", "OncoKB + CiVIC · same tumor type"),
        ("cross", "CiVIC evidence from a different tumor"),
        ("none", "No therapeutic level · baseline Tier III"),
    ]
    scale = 300 / max(n for n, _ in st.values())
    b = []
    s, x = swatch(0, 5, "Agrees with baseline", ACCENT)
    b.append(s)
    b.append(swatch(x, 5, "Re-tiered by agent", WARM)[0])
    for i, (k, label) in enumerate(labels):
        n, over = st[k]
        y = 22 + i * 43
        b.append(text(0, y, label, FS, weight=500))
        w_ok, w_over = (n - over) * scale, over * scale
        by = y + 8.5
        if over:
            b.append(rect(0, by, w_ok, 14, ACCENT) if w_ok else "")
            b.append(hbar(w_ok + (1 if w_ok else 0), by, w_over, 14, WARM))
        else:
            b.append(hbar(0, by, w_ok, 14, ACCENT))
        end = w_ok + w_over + (1 if w_ok and over else 0)
        b.append(
            f'<text x="{end + 6:.2f}" y="{by + 7}" font-size="{FS}" font-weight="700" '
            f'fill="{INK}" dominant-baseline="central">{over} / {n}<tspan dx="3.5" '
            f'font-size="{FS_S}" font-weight="400" fill="{INK2}">re-tiered</tspan></text>'
        )
    return svg(COL, 22 + 4 * 43 - 8, b, "Re-tiering by evidence stratum")


def repro_grid():
    """Tier call per case across three Opus runs, Sonnet and Haiku."""
    cases = {r["case_id"]: r for r in _rows("agent_combined_summary.csv")}
    rows = _rows("agent_cross_model.csv")
    cols = ["opus_r1_tier", "opus_r2_tier", "opus_r3_tier", "sonnet_tier", "haiku_tier"]
    rows.sort(
        key=lambda r: (
            r["cross_model_unanimous"] != "TRUE",
            TIERS.index(r["opus_r1_tier"]),
        )
    )
    left, top, pitch = 140, 30, 20.5
    cw = (COL - left) / 5
    b = [
        text(left + cw * 1.5, 5, "Opus 4.7", FS_S, "middle", 600),
        line(left + 2, 12, left + cw * 3 - 2, 12, INK, 0.6),
        text(left + cw * 3.5, 5, "Sonnet 4.6", FS_S, "middle", 600),
        text(left + cw * 4.5, 5, "Haiku 4.5", FS_S, "middle", 600),
    ]
    for j, lab in enumerate(["run 1", "run 2", "run 3", "", ""]):
        b.append(text(left + cw * (j + 0.5), 20, lab, FS_S, "middle", fill=MUTED))
    for i, r in enumerate(rows):
        c = cases[r["case_id"]]
        y = top + i * pitch
        tumor = TUMOR.get(c["tumor_type"], c["tumor_type"])
        b.append(text(0, y + 9.5, f"{c['gene']} {c['variant']}", FS_S, weight=600))
        b.append(text(left - 8, y + 9.5, tumor, FS_S, "end", fill=INK2))
        for j, col in enumerate(cols):
            fill = TIER_FILL[r[col]]
            b.append(rect(left + j * cw + 0.8, y, cw - 1.6, 19, fill, rx=1.6))
            b.append(
                text(
                    left + cw * (j + 0.5),
                    y + 9.5,
                    TIER_LABEL[r[col]][5:],
                    FS_S,
                    "middle",
                    600,
                    on(fill),
                )
            )
    return svg(
        COL, top + len(rows) * pitch - 1.5, b, "Tier calls across runs and models"
    )
