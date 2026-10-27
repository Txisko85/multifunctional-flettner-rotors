#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Graphical abstract of A1 for Applied Energy (21-10-2026; version 4, 26-10-2026; version 5, 27-10-2026: April-October,
scenario A-VG).

13 cm x 5.2 cm: the 2.5:1 proportion of Elsevier (minimum 1328 x 531 pixels), legible at 13 cm wide.
Left: the case study (map of the study area and ship with four rotors) and the saving of the conventional rotors
relative to the ship without them (R1 against R0). Right: the fuel saved by each added function relative to the
conventional rotor, with four units, in the April-October service and coefficient scenario A: reference value
(columns) and 5-95 % interval of the Sobol analysis within that scenario (lines), over the drawings of the topologies
of Fig. 1 (fig_topologies.py), at the same scale as each other; next to the extension, what the comparators without a
new mechanism save (rotor 6 m taller and tiltable rotor).
All values are read from the results of the model (../route/route_results.json and ../route/uncertainty.json).
The map shows only the study area and draws no borders, like Fig. 4 (map policy of Elsevier).
Writes graphical_abstract.pdf, graphical_abstract.png (300 dpi) and graphical_abstract.tif (600 dpi).
"""
import json
import math
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROUTE = os.path.join(HERE, "..", "route")
sys.path.insert(0, HERE)
sys.path.insert(0, ROUTE)
import fig_topologies as ft          # noqa: E402  drawings of Fig. 1 (and the Inter typeface, if installed)
import fig_route_map as rm           # noqa: E402  study area and projection of Fig. 4

# Typeface: Arial, as the Elsevier guide asks for the graphical abstract, or Liberation Sans, with identical metrics,
# if Arial is not installed.
from matplotlib import font_manager  # noqa: E402
_FONTS = {f.name for f in font_manager.fontManager.ttflist}
FONT = next((f for f in ("Arial", "Liberation Sans") if f in _FONTS), "DejaVu Sans")

CM = 1 / 2.54
W, H = 13.0 * CM, 5.2 * CM           # inches
INK, MUTED, QUIET = "#1f1f1d", "#52514e", "#8a897f"
GRID, PANEL, RULE = "#e4e2db", "#f4f3ef", "#cfccc3"
SEA, LAND, COAST = "#eef3f7", "#e7e4dc", "#a9a59a"
HULL = "#4b5764"
LWS = 0.6                            # the lines of the drawings of Fig. 1, thinner at this size

_line = ft.line


def _thin_line(ax, pts, color=ft.OUTLINE, lw=ft.LW, zo=5, **kw):
    _line(ax, pts, color=color, lw=lw * LWS, zo=zo, **kw)


ft.line = _thin_line


def axes_in(fig, x, y, w, h):
    """Axes in inches from the lower left corner."""
    return fig.add_axes([x / W, y / H, w / W, h / H])


def text_in(fig, x, y, s, **kw):
    return fig.text(x / W, y / H, s, **kw)


def fmt(v):
    """+64, +9.8, −24: two significant figures and a typographic sign."""
    s = f"{abs(v):.0f}" if abs(v) >= 10 else f"{abs(v):.1f}"
    return ("+" if v > 0 else "−") + s


# ---------------------------------------------------------------- data
def load_results():
    res = json.load(open(os.path.join(ROUTE, "route_results.json"), encoding="utf-8"))
    mc = json.load(open(os.path.join(ROUTE, "uncertainty.json"), encoding="utf-8"))["scenarios"]["A"]["summary"]
    rows = {r[0]: r[3] for r in res["base"]["rows"]}          # total saving relative to R1, t/yr
    r1 = res["base"]["sea"]["R1"]
    return dict(rows=rows, mc=mc, r1=r1, r1_pct=100 * r1 / res["base"]["fuel_total"],
                r1_p5=mc["R1 %"]["p5"], r1_p95=mc["R1 %"]["p95"], eps=res["eps_star"],
                r130=res["base"]["sea"]["R1-30"], r1t=res["base"]["sea"]["R1t"],
                r1_vg_pct=res["scenarios"]["coefficients A-VG"]["R1_pct"])


# ---------------------------------------------------------------- left zone: case study
def mini_map(fig, x, y, w):
    route = json.load(open(os.path.join(ROUTE, "route.json"), encoding="utf-8"))
    land = json.load(open(os.path.join(ROUTE, "map_data", "ne_10m_land_clipped.geojson"), encoding="utf-8"))
    pts = list(rm.COORDS.values()) + [(p["lat"], p["lon"]) for p in route["points"] + route["ports"]]
    lon0, lon1 = min(p[1] for p in pts) - rm.PAD, max(p[1] for p in pts) + rm.PAD
    lat0, lat1 = min(p[0] for p in pts) - rm.PAD, max(p[0] for p in pts) + rm.PAD
    y0, y1 = float(rm.merc(lat0)), float(rm.merc(lat1))
    h = w * (y1 - y0) / (lon1 - lon0)
    ax = axes_in(fig, x, y, w, h)
    ax.set_facecolor(SEA)
    for feat in land["features"]:
        g = feat["geometry"]
        for poly in (g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]):
            ax.add_patch(PathPatch(rm.ring_path(poly), fc=LAND, ec=COAST, lw=0.3, zorder=1))
    lat = np.array([rm.COORDS[k][0] for k in rm.SEQ])
    lon = np.array([rm.COORDS[k][1] for k in rm.SEQ])
    yy = rm.merc(lat)
    ax.plot(lon[:3], yy[:3], color=INK, lw=0.8, ls=(0, (1, 1.4)), zorder=3)
    ax.plot(lon[5:], yy[5:], color=INK, lw=0.8, ls=(0, (1, 1.4)), zorder=3)
    ax.plot(lon[2:6], yy[2:6], color=INK, lw=1.0, zorder=3, solid_joinstyle="round")
    for p in route["ports"]:
        ax.scatter(p["lon"], rm.merc(p["lat"]), s=13, marker="s", c="white", edgecolors=INK, linewidths=0.6,
                   zorder=5)
    # names of the ports outside the frame, above and below each one
    lab = dict(color=INK, fontsize=6.3, clip_on=False)
    ax.text(rm.COORDS["Portsmouth"][1] + 0.25, y1 + 0.30, "Portsmouth", ha="right", va="bottom", **lab)
    ax.text(rm.COORDS["Bilbao"][1], y0 - 0.40, "Bilbao", ha="center", va="top", **lab)
    ax.set_xlim(lon0, lon1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.5)
        s.set_color(QUIET)
    return h


def ferry(fig, x, y, w):
    """Generic ro-pax in profile, bow to the right, with four rotors on the upper deck (base at 28 m)."""
    x0, x1, z0, z1 = -6.0, 222.0, -7.0, 56.0
    h = w * (z1 - z0) / (x1 - x0)
    ax = axes_in(fig, x, y, w, h)
    ax.set_xlim(x0, x1)
    ax.set_ylim(z0, z1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((x0, z0), x1 - x0, -z0, fc=SEA, ec="none", zorder=0))
    ax.plot([x0, x1], [0, 0], color="#b9c9d6", lw=0.5, zorder=0.5)
    ax.add_patch(Polygon([(4, 0), (0, 14), (196, 14), (215, 16.5), (207, 5), (199, 0)], closed=True, fc=HULL,
                         ec="none", zorder=1))
    ax.add_patch(Polygon([(24, 14), (24, 28), (166, 28), (180, 24), (187, 14)], closed=True, fc="#f3f1ec",
                         ec=QUIET, lw=0.3, zorder=2))
    for z in (18.5, 23.0):
        ax.plot([30, 176], [z, z], color="#8d9aa6", lw=0.9, zorder=2.5, solid_capstyle="butt")
    ax.add_patch(Polygon([(148, 28), (148, 32), (170, 32), (176, 28)], closed=True, fc="#f3f1ec", ec=QUIET, lw=0.3,
                         zorder=2))
    ax.add_patch(Polygon([(40, 28), (42, 40), (54, 40), (56, 28)], closed=True, fc=HULL, ec="none", zorder=2))
    for xc in (72.0, 92.0, 112.0, 132.0):
        ax.add_patch(Rectangle((xc - 2.0, 28), 4.0, 24.0, fc=ft.ROTOR, ec=ft.OUTLINE, lw=0.3, zorder=3))
        ax.add_patch(Rectangle((xc - 4.0, 52), 8.0, 1.4, fc="#dcdad3", ec=ft.OUTLINE, lw=0.3, zorder=3.1))
    return h


# ---------------------------------------------------------------- right zone: saving by added function
def column(ax, xc, v, width, color, rx, ry):
    """Column from 0 to v with the data end rounded (radii rx, ry in axes units)."""
    a, b = xc - width / 2, xc + width / 2
    s = 1.0 if v >= 0 else -1.0
    ry = min(ry, abs(v) / 2)
    k = 0.5523                                        # Bézier curve of a quarter circle
    top = v
    verts = [(a, 0), (a, top - s * ry), (a, top - s * ry * (1 - k)), (a + rx * (1 - k), top), (a + rx, top),
             (b - rx, top), (b - rx * (1 - k), top), (b, top - s * ry * (1 - k)), (b, top - s * ry), (b, 0), (a, 0)]
    codes = [Path.MOVETO, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4,
             Path.CURVE4, Path.LINETO, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), fc=color, ec="none", zorder=3))


def icons(fig, x, y, w, h, centers_in, scale):
    """Drawings of Fig. 1 at the same scale (inches per metre), centred below each column."""
    ax = axes_in(fig, x, y, w, h)
    ax.set_xlim(0, w / scale)
    ax.set_ylim(-2.4, -2.4 + h / scale)
    ax.set_aspect("equal")
    ax.axis("off")
    blue, orange, green = ft.CONF["C1"][0], ft.CONF["C2"][0], ft.CONF["C3"][0]
    c = [(ci - x) / scale for ci in centers_in]
    ft.extension(ax, c[0], orange)
    ft.contact_shadow(ax, c[1], 3.6)
    ft.savonius(ax, c[1], 0.0, ft.H_F, ft.R_F, green)
    ft.rotor(ax, c[2])
    ft.turbine(ax, c[2], ft.H_F + ft.T_E + 0.15, blue, ft.H_F + ft.T_E)
    ft.rotor(ax, c[3])
    ft.turbine(ax, c[3], ft.H_F + ft.T_E + 0.15, blue, ft.H_F + ft.T_E)
    for k in (0, 2):                                  # device stowed inside the rotor
        xa = c[k] + 3.6
        ax.annotate("", xy=(xa, ft.H_F + 6.4), xytext=(xa, ft.H_F + 0.8),
                    arrowprops=dict(arrowstyle="<|-|>", lw=0.5, color=MUTED, mutation_scale=3.2, shrinkA=0,
                                    shrinkB=0), zorder=9)
    for ci in c:
        ft.deck(ax, ci - 5.2, ci + 5.2)
    return ax


def main():
    d = load_results()
    plt.rcParams.update({"font.family": FONT, "font.size": 7.0})
    fig = plt.figure(figsize=(W, H))

    # --- case study
    mw = 0.62
    mh = mini_map(fig, 0.08, 0.0, mw)
    fig.axes[-1].set_position([0.08 / W, (H - mh) / 2 / H, mw / W, mh / H])
    tx = 0.83
    ferry(fig, tx, H - 0.05 - 0.25, 0.92)
    yt = H - 0.05 - 0.25 - 0.08
    text_in(fig, tx, yt, "Ro-pax ferry with\nfour Flettner rotors", ha="left", va="top", fontsize=7.6,
            fontweight="bold", color=INK, linespacing=1.12)
    text_in(fig, tx, yt - 0.30, "Portsmouth\u2013Bilbao,\nApril\u2013October,\nhourly ERA5 winds", ha="left", va="top",
            fontsize=7.2, color=MUTED, linespacing=1.12)
    text_in(fig, tx, 0.72, f"{d['r1_pct']:.1f}%", ha="left", va="bottom", fontsize=15, fontweight="bold",
            color=INK)
    text_in(fig, tx, 0.68, "less fuel with\nconventional rotors\n"
            f"({d['r1_p5']:.1f}\u2013{d['r1_p95']:.1f}%;\n{d['r1_vg_pct']:.1f}% if fitted\nto full scale)", ha="left", va="top", fontsize=7.2, color=MUTED,
            linespacing=1.12)
    fig.add_artist(plt.Line2D([1.90 / W, 1.90 / W], [0.10 / H, (H - 0.10) / H], color=RULE, lw=0.5))

    # --- saving by added function
    px = 2.30
    pw = W - px - 0.06
    py, ph = 0.88, 0.93
    ax = axes_in(fig, px, py, pw, ph)
    cols = [("T3x", "Telescopic\nextension", "C2"), ("T5", "Savonius\nconversion", "C3"),
            ("T4", "Retractable\nturbine", "C1"), ("T1", "Fixed\nturbine", "C1")]
    xl, xr = 0.45, 4.55
    ylo, yhi = -32.0, 124.0
    ax.set_xlim(xl, xr)
    ax.set_ylim(ylo, yhi)
    ax.axis("off")
    sx, sy = pw / (xr - xl), ph / (yhi - ylo)          # inches per axes unit
    for v in (50, 100):
        ax.plot([xl, xr], [v, v], color=GRID, lw=0.45, zorder=1)
    for v in (0, 50, 100):
        ax.text(xl - 0.04 / sx, v, f"{v}", ha="right", va="center", fontsize=6.6, color=MUTED)
    ax.plot([xl, xr], [0, 0], color=INK, lw=0.6, zorder=4)
    width = 0.17 / sx
    for i, (key, name, conf) in enumerate(cols, start=1):
        v = d["rows"][key]
        lo, hi = d["mc"][key]["p5"], d["mc"][key]["p95"]
        column(ax, i, v, width, ft.CONF[conf][0], 0.03 / sx, 0.03 / sy)
        ax.plot([i, i], [lo, hi], color=INK, lw=0.6, zorder=5, solid_capstyle="butt")
        for e in (lo, hi):
            ax.plot([i - 0.035 / sx, i + 0.035 / sx], [e, e], color=INK, lw=0.6, zorder=5)
        ax.text(i + width / 2 + 0.035 / sx, v + (0.015 / sy if v > 0 else 0), fmt(v), ha="left",
                va="bottom" if v > 0 else "center", fontsize=8.0, fontweight="bold", color=INK, zorder=6)
        if key == "T3x":                                # key of the lines, only once, and comparators
            ax.text(i, hi + 0.03 / sy, "5\u201395%", ha="center", va="bottom", fontsize=6.6, color=MUTED)
    centers = [px + (i - xl) * sx for i in range(1, 5)]

    # message: interaction needed for a turbine to match the extension
    bx0, bx1, by1 = 2.55, 4.55, 122.0
    pad = 0.05
    t = ax.text(bx0 + 0.07 / sx, by1 - pad / sy, "To match the extension,\na turbine on the rotor would\n"
                f"need {100 * d['eps']:.0f}% more rotor force", ha="left", va="top", fontsize=6.9, color=INK,
                linespacing=1.12, zorder=2)
    fig.canvas.draw()
    bb = t.get_window_extent().transformed(ax.transData.inverted())
    by0 = bb.y0 - pad / sy
    ax.add_patch(FancyBboxPatch((bx0, by0), bx1 - bx0, by1 - by0, boxstyle="round,pad=0,rounding_size=0.06",
                                fc=PANEL, ec="none", zorder=1.5, mutation_aspect=sx / sy))
    ax.text(1.32, -0.06 / sy, f"Rotor 6 m taller {fmt(d['r130'])} \u00b7 tiltable rotor {fmt(d['r1t'])}", ha="left", va="top", fontsize=6.3, color=MUTED, linespacing=1.12,
            zorder=6)

    text_in(fig, 2.02, H - 0.05, "Fuel saved vs conventional rotors (t/yr, four units)", ha="left", va="top",
            fontsize=7.6, fontweight="bold", color=INK)

    # drawings and names
    iy, ih = 0.31, 0.53
    scale = ih / 34.8                                  # from -2.4 m (deck) to 32.4 m
    icons(fig, px, iy, pw, ih, centers, scale)
    for c, (key, name, conf) in zip(centers, cols):
        text_in(fig, c, iy - 0.02, name, ha="center", va="top", fontsize=7.0, color=INK, linespacing=1.05)

    for ext, dpi in (("pdf", 300), ("png", 300)):
        fig.savefig(os.path.join(HERE, f"graphical_abstract.{ext}"), dpi=dpi)
    tif = os.path.join(HERE, "graphical_abstract.tif")
    fig.savefig(tif, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    Image.open(tif).convert("RGB").save(tif, compression="tiff_lzw", dpi=(600, 600))   # no alpha channel
    print({k: round(v, 1) for k, v in d["rows"].items()}, round(d["r1_pct"], 2), "font", FONT,
          "box", round(bx1 - bx0, 2), round(by0, 1), "text width (in)", round(bb.width * sx, 2))


if __name__ == "__main__":
    main()
