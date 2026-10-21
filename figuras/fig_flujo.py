#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. 1 de A1 (flujo de evaluación), generada como SVG estático a partir del widget del borrador (pub 3, 19-10-2026)
y convertida a PDF vectorial y PNG con Chromium (Playwright). Sin título dentro de la figura: va en el pie.
"""
import os
from playwright.sync_api import sync_playwright

HERE = os.path.dirname(os.path.abspath(__file__))
EDGE, ACCENT, INK, QUIET = "#8a897f", "#2a78d6", "#0b0b0b", "#52514e"
Y0, H = 44, 524          # recorte: sin la línea de título del widget


def box(x, y, w, h, main=False):
    if main:
        return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{ACCENT}" fill-opacity="0.12" '
                f'stroke="{ACCENT}" stroke-width="2"/>')
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="none" stroke="{EDGE}" stroke-width="1.25"/>'


def txt(x, y, s, bold=False, size=13, color=INK, anchor="middle"):
    w = ' font-weight="600"' if bold else ""
    return f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}"{w} fill="{color}">{s}</text>'


def block(x, y, w, h, title, lines, main=False):
    cx = x + w / 2
    out = [box(x, y, w, h, main), txt(cx, y + 24, title, bold=True)]
    for i, line in enumerate(lines):
        out.append(txt(cx, y + 40 + 16 * i, line, size=11.5, color=INK if main else QUIET))
    return "\n".join(out)


def svg():
    arrows = ["M380 148V168", "M380 224V264", "M228 356V376", "M532 356V376", "M380 468V488", "M560 524H592",
              "M700 496V128"]
    plain = ["M132 128V148", "M380 128V148", "M628 128V148", "M132 148H628", "M380 336V356", "M228 356H532",
             "M228 448V468", "M532 448V468", "M228 468H532"]
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 {Y0} 760 {H}" width="760" height="{H}" '
             f'font-family="Inter, DejaVu Sans, Arial, sans-serif">',
             '<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" '
             f'orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="{EDGE}"/></marker></defs>',
             f'<g fill="none" stroke="{EDGE}" stroke-width="1.25">']
    parts += [f'<path d="{d}"/>' for d in plain]
    parts += [f'<path d="{d}" marker-end="url(#arr)"/>' for d in arrows]
    parts.append("</g>")
    parts.append(txt(692, 292, "resamples", size=11.5, color=QUIET, anchor="end"))
    parts.append(txt(692, 306, "coefficients", size=11.5, color=QUIET, anchor="end"))
    parts.append(block(24, 56, 216, 72, "ERA5 winds", ["Hourly, at 10 m and 100 m", "Along the route and at ports"]))
    parts.append(block(272, 56, 216, 72, "Ship and route", ["Power–speed curve and SFOC",
                                                            "Timetable and time at berth"]))
    parts.append(block(520, 56, 216, 72, "Coefficient database", ["Configurations C0 to C3",
                                                                  "Central values and ranges"]))
    parts.append(block(200, 168, 360, 56, "Hourly state", ["Ship position, true and apparent wind, spin ratio"]))
    parts.append(block(160, 264, 440, 72, "Benefit of each admissible mode",
                       ["Forces, spin power, turbine power and induced resistance",
                        "Every flow converted into fuel at its point of entry"]))
    parts.append(block(96, 376, 264, 72, "Optimal strategy", ["Best mode each hour; with a switching",
                                                              "cost, by dynamic programming"]))
    parts.append(block(400, 376, 264, 72, "Rule-based strategy", ["Favourable or unfavourable wind,",
                                                                  "as in the patent rules"]))
    parts.append(block(200, 488, 360, 72, "Annual fuel and CO₂ per topology",
                       ["At sea and at berth, relative to R1", "Ranking, cost of rules and switching, ε*"], main=True))
    parts.append(block(592, 496, 144, 56, "Monte Carlo", ["N = 1,000 samples"]))
    parts.append("</svg>")
    return "\n".join(parts)


def main():
    s = svg()
    open(os.path.join(HERE, "fig_flujo.svg"), "w", encoding="utf-8").write(s)
    html = ('<!doctype html><html><head><meta charset="utf-8"><style>@page{size:760px ' + str(H) + 'px;margin:0}'
            'html,body{margin:0;padding:0;background:#fff}svg{display:block}</style></head><body>' + s +
            '</body></html>')
    with sync_playwright() as p:
        b = p.chromium.launch()
        pg = b.new_page(viewport={"width": 760, "height": H}, device_scale_factor=3)
        pg.set_content(html)
        pg.pdf(path=os.path.join(HERE, "fig_flujo.pdf"), width="760px", height=f"{H}px", print_background=True,
               margin={"top": "0", "bottom": "0", "left": "0", "right": "0"})
        pg.screenshot(path=os.path.join(HERE, "fig_flujo.png"), full_page=True)
        b.close()


if __name__ == "__main__":
    main()
