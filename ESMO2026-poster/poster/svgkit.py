"""Tiny SVG helpers for the poster figures.

One user unit = 1 mm at true size: every figure's viewBox width equals the
width of the poster column it sits in, so font sizes below are real millimetres
(8.6 mm = 24 pt, the smallest label size used on the poster).
"""

from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# Palette (validated with the dataviz validator against a white surface)
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
HAIR, AXIS, CELL = "#e1e0d9", "#c3c2b7", "#f1f0ec"
BLUE = {
    100: "#cde2fb",
    200: "#9ec5f4",
    250: "#86b6ef",
    350: "#5598e7",
    400: "#3987e5",
    500: "#256abf",
    550: "#1c5cab",
    600: "#184f95",
    700: "#0d366b",
}
ACCENT = BLUE[600]
WARM, WARM_SOFT = "#eb6834", "#f6b99f"

FS, FS_S = 9.6, 8.6  # 27 pt / 24 pt
COL, WIDE = 425.5, 887.0  # single column / double column width (mm)

TIERS = [
    "Tier I Level A",
    "Tier I Level B",
    "Tier II Level C",
    "Tier II Level D",
    "Tier III VUS",
]
TIER_LABEL = {
    "Tier I Level A": "Tier I A",
    "Tier I Level B": "Tier I B",
    "Tier II Level C": "Tier II C",
    "Tier II Level D": "Tier II D",
    "Tier III VUS": "Tier III",
}
TIER_FILL = {
    "Tier I Level A": BLUE[700],
    "Tier I Level B": BLUE[550],
    "Tier II Level C": BLUE[400],
    "Tier II Level D": BLUE[250],
    "Tier III VUS": AXIS,
}
DARK_FILLS = {BLUE[700], BLUE[600], BLUE[550], BLUE[500], BLUE[400], INK}


def on(fill):
    """Text colour that clears contrast on a given fill."""
    return "#ffffff" if fill in DARK_FILLS else INK


def text(x, y, s, size=FS, anchor="start", weight=400, fill=INK, style=""):
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" '
        f'text-anchor="{anchor}" font-weight="{weight}" fill="{fill}" '
        f'dominant-baseline="central" style="{style}">{escape(str(s))}</text>'
    )


def rect(x, y, w, h, fill, rx=0.0, stroke="none", sw=0.0, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(w, 0):.2f}" '
        f'height="{h:.2f}" rx="{rx}" fill="{fill}" stroke="{stroke}" '
        f'stroke-width="{sw}"{d}/>'
    )


def line(x1, y1, x2, y2, stroke=AXIS, sw=0.5):
    return (
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
        f'stroke="{stroke}" stroke-width="{sw}" stroke-linecap="round"/>'
    )


def hbar(x, y, w, h, fill, r=2.0):
    """Horizontal bar: square at the baseline, rounded at the data end."""
    if w <= 0:
        return ""
    r = min(r, w, h / 2)
    return (
        f'<path d="M{x:.2f},{y:.2f} h{w - r:.2f} a{r},{r} 0 0 1 {r},{r} '
        f'v{h - 2 * r:.2f} a{r},{r} 0 0 1 {-r},{r} h{-(w - r):.2f} z" '
        f'fill="{fill}"/>'
    )


def chip(x, y, w, h, label, fill, size=FS, stroke="none", dash="", ink=None):
    return rect(
        x,
        y,
        w,
        h,
        fill,
        rx=2,
        stroke=stroke,
        sw=0.7 if stroke != "none" else 0,
        dash=dash,
    ) + text(
        x + w / 2,
        y + h / 2,
        label,
        size=size,
        anchor="middle",
        weight=600,
        fill=ink or on(fill),
    )


def arrow_right(x1, y, x2, stroke=INK, sw=0.7):
    return (
        line(x1, y, x2 - 2.6, y, stroke, sw)
        + f'<path d="M{x2 - 3.6:.2f},{y - 2.2:.2f} L{x2:.2f},{y:.2f} '
        f'L{x2 - 3.6:.2f},{y + 2.2:.2f} z" fill="{stroke}"/>'
    )


def swatch(x, y, label, fill, size=FS_S, stroke="none"):
    """Legend key: returns (svg, x after the label)."""
    s = rect(
        x,
        y - 3.6,
        7.2,
        7.2,
        fill,
        rx=1.2,
        stroke=stroke,
        sw=0.5 if stroke != "none" else 0,
    )
    s += text(x + 10.5, y, label, size=size, fill=INK2)
    return s, x + 10.5 + len(label) * size * 0.5 + 12


def svg(w, h, body, label):
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h:.1f}" '
        f'role="img" aria-label="{escape(label)}">{"".join(body)}</svg>'
    )
