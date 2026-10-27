#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map of the Portsmouth–Bilbao route (A1, 21-10-2026): the reconstructed route (route.json, from route_geometry.py), its
two open-sea legs, the ten ERA5 sampling points and the two ports, coloured by the mean wind speed at 100 m in
2015–2024 (from the hourly ERA5 file).
The map shows only the study area, as the map policy of Elsevier asks: the rectangle that contains the route, its
waypoints and the two ports, enlarged by 0.25° on each side. It draws no borders.
Coastline: Natural Earth 1:10 m (public domain), clipped with clip_coastline.py (map_data/).
Mercator projection. Writes route_map.pdf and route_map.png.
Usage: python fig_route_map.py ../data/era5_hourly_route.csv.gz
"""
import csv
import gzip
import json
import math
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path

if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, QUIET = "#1f1f1d", "#52514e", "#8a897f"
SEA, LAND, COAST, GRID = "#eef3f7", "#e7e4dc", "#a9a59a", "#dfe6ec"
PAD = 0.25                      # margin around the study area, in degrees
# single-hue sequential ramp (blue), from light to dark
CMAP = LinearSegmentedColormap.from_list("blue", ["#d9e8f8", "#8fbaf0", "#2a78d6", "#123f7a"])
VMIN, VMAX = 4.0, 10.0
# waypoints of route_geometry.py (lat, lon)
COORDS = {"Portsmouth": (50.790, -1.110), "Nab": (50.670, -0.950), "St Catherine S": (50.500, -1.300),
          "Start Point S": (49.950, -3.500), "Ushant SW": (48.650, -5.750), "Bilbao approach": (43.450, -3.080),
          "Bilbao": (43.360, -3.030)}
SEQ = ["Portsmouth", "Nab", "St Catherine S", "Start Point S", "Ushant SW", "Bilbao approach", "Bilbao"]


def merc(lat):
    """Mercator ordinate in degrees (the abscissa is the longitude), for a 1:1 aspect."""
    return np.degrees(np.log(np.tan(np.pi / 4 + np.radians(np.asarray(lat, float)) / 2)))


def ring_path(rings):
    verts, codes = [], []
    for ring in rings:
        pts = [(x, float(merc(y))) for x, y in ring]
        verts += pts + [pts[0]]
        codes += [Path.MOVETO] + [Path.LINETO] * (len(pts) - 1) + [Path.CLOSEPOLY]
    return Path(verts, codes)


def mean_wind(csv_path):
    acc = {}
    with gzip.open(csv_path, "rt", encoding="utf-8") as f:
        for row in csv.reader(line for line in f if not line.startswith("#")):
            if row[0] == "point":
                continue
            s, n = acc.get(row[0], (0.0, 0))
            acc[row[0]] = (s + float(row[3]), n + 1)
    return {k: s / n for k, (s, n) in acc.items()}


def main():
    era5 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "data", "era5_hourly_route.csv.gz")
    u100 = mean_wind(era5)
    route = json.load(open(os.path.join(HERE, "route.json"), encoding="utf-8"))
    land = json.load(open(os.path.join(HERE, "map_data", "ne_10m_land_clipped.geojson"), encoding="utf-8"))

    # study area: route, waypoints, ERA5 points and ports, plus the margin
    pts = list(COORDS.values()) + [(p["lat"], p["lon"]) for p in route["points"] + route["ports"]]
    lon0, lon1 = min(p[1] for p in pts) - PAD, max(p[1] for p in pts) + PAD
    lat0, lat1 = min(p[0] for p in pts) - PAD, max(p[0] for p in pts) + PAD
    y0, y1 = float(merc(lat0)), float(merc(lat1))

    # layout in inches: map on the left, legend on the right
    mw = 2.62
    mh = mw * (y1 - y0) / (lon1 - lon0)
    ml, mb, mt, cw = 0.40, 0.27, 0.07, 1.50
    fw, fh = ml + mw + cw, mb + mh + mt
    plt.rcParams.update({"font.size": 7.0})
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([ml / fw, mb / fh, mw / fw, mh / fh])
    ax.set_facecolor(SEA)

    # graticule every degree
    glon = range(math.ceil(lon0), math.floor(lon1) + 1)
    glat = range(math.ceil(lat0), math.floor(lat1) + 1)
    for lon in glon:
        ax.plot([lon, lon], [y0, y1], color=GRID, lw=0.5, zorder=0.5)
    for lat in glat:
        ax.plot([lon0, lon1], [merc(lat)] * 2, color=GRID, lw=0.5, zorder=0.5)

    # land
    for feat in land["features"]:
        g = feat["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for poly in polys:
            ax.add_patch(PathPatch(ring_path(poly), fc=LAND, ec=COAST, lw=0.35, zorder=1))

    # route: port approaches (dotted) and open-sea legs
    lat = np.array([COORDS[k][0] for k in SEQ])
    lon = np.array([COORDS[k][1] for k in SEQ])
    y = merc(lat)
    ax.plot(lon[:3], y[:3], color=INK, lw=1.0, ls=(0, (1, 1.6)), zorder=3, solid_capstyle="round")
    ax.plot(lon[5:], y[5:], color=INK, lw=1.0, ls=(0, (1, 1.6)), zorder=3)
    ax.plot(lon[2:6], y[2:6], color=INK, lw=1.3, zorder=3, solid_joinstyle="round")
    # ends of the legs: start of the Channel, Ushant and end of the Bay of Biscay
    for k in ("St Catherine S", "Ushant SW", "Bilbao approach"):
        i = SEQ.index(k)
        ax.plot(lon[i], y[i], marker="o", ms=3.0, mfc="white", mec=INK, mew=0.8, zorder=4)

    # sampling points and ports, coloured by the mean wind at 100 m
    norm = Normalize(VMIN, VMAX)
    names = ["C1", "C2", "C3", "C4", "B1", "B2", "B3", "B4", "B5", "B6"]
    for name, p in zip(names, route["points"]):
        ax.scatter(p["lon"], merc(p["lat"]), s=34, c=[CMAP(norm(u100[name]))], edgecolors=INK, linewidths=0.55,
                   zorder=5)
    for p, key in ((route["ports"][0], "PORTSMOUTH"), (route["ports"][1], "BILBAO")):
        ax.scatter(p["lon"], merc(p["lat"]), s=46, marker="s", c=[CMAP(norm(u100[key]))], edgecolors=INK,
                   linewidths=0.7, zorder=6)

    # labels
    lab = dict(color=INK, fontsize=7.2, zorder=7)
    ax.text(-1.30, merc(50.90), "Portsmouth", ha="right", va="center", **lab)
    ax.text(-2.86, merc(43.31), "Bilbao", ha="left", va="center", **lab)
    ax.text(-5.06, merc(48.53), "Ushant", ha="center", va="bottom", color=MUTED, fontsize=6.6, zorder=7)
    sea = dict(color="#5d7891", fontsize=7.2, style="italic", ha="center", va="center", zorder=2)
    ax.text(-1.55, merc(50.06), "English\nChannel", linespacing=1.0, **sea)
    ax.text(-2.45, merc(45.55), "Bay of\nBiscay", linespacing=1.0, **sea)
    country = dict(color=QUIET, fontsize=6.4, ha="center", va="center", zorder=2)

    def spaced(t):
        return " ".join(" ".join(w) for w in t.split())

    ax.text(-3.55, merc(50.86), spaced("UNITED KINGDOM"), **country)
    ax.text(-1.35, merc(47.75), spaced("FRANCE"), **country)
    ax.text(-4.65, merc(43.24), spaced("SPAIN"), **country)

    def along(i0, i1, frac, off, text):
        """Label parallel to the leg i0-i1, shifted off degrees to its left (seen from i0 to i1)."""
        xa, ya, xb, yb = lon[i0], y[i0], lon[i1], y[i1]
        dx, dy = xb - xa, yb - ya
        ln = math.hypot(dx, dy)
        nx, ny = -dy / ln, dx / ln
        ang = math.degrees(math.atan2(dy, dx))
        if ang > 90 or ang < -90:
            ang -= 180 * (1 if ang > 0 else -1)
        ax.text(xa + frac * dx + off * nx, ya + frac * dy + off * ny, text, rotation=ang, rotation_mode="anchor",
                ha="center", va="center", color=INK, fontsize=6.8, zorder=7)

    along(3, 4, 0.42, -0.30, "Channel leg, 208 nm")
    along(4, 5, 0.47, -0.30, "Bay of Biscay leg, 331 nm")

    # scale bar of 100 nautical miles at 43.78°N, at sea west of the route
    lat_s, lon_s = 43.78, lon0 + 0.10
    dlon = (100 / 60) / math.cos(math.radians(lat_s))
    ys = float(merc(lat_s))
    for k in range(2):
        ax.add_patch(Rectangle((lon_s + k * dlon / 2, ys), dlon / 2, 0.09, fc=INK if k == 0 else "white", ec=INK,
                               lw=0.5, zorder=7))
    ax.text(lon_s, ys + 0.2, "0", ha="center", va="bottom", fontsize=6.4, color=INK, zorder=7)
    ax.text(lon_s + dlon, ys + 0.2, "100 nm", ha="right", va="bottom", fontsize=6.4, color=INK, zorder=7)

    # axes
    ax.set_xlim(lon0, lon1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks(list(glon))
    ax.set_xticklabels([f"{abs(v)}°W" if v < 0 else f"{v}°E" if v > 0 else "0°" for v in glon])
    ax.set_yticks([float(merc(v)) for v in glat])
    ax.set_yticklabels([f"{v}°N" for v in glat])
    ax.tick_params(labelsize=6.6, length=2.5, width=0.5, colors=MUTED)
    for s in ax.spines.values():
        s.set_linewidth(0.6)
        s.set_color(MUTED)

    # legend, to the right of the map (coordinates in inches)
    lx = ml + mw + 0.16
    top = mb + mh
    fig.text(lx / fw, (top - 0.02) / fh, "Mean wind speed at 100 m,\n2015–2024 (m/s)", ha="left", va="top",
             fontsize=6.6, color=INK, linespacing=1.15)
    cax = fig.add_axes([lx / fw, (top - 0.42) / fh, 1.10 / fw, 0.075 / fh])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_ticks([4, 6, 8, 10])
    cb.ax.tick_params(labelsize=6.4, length=2, width=0.5, colors=INK)
    cb.outline.set_linewidth(0.5)
    ov = fig.add_axes([0, 0, 1, 1])
    ov.set_xlim(0, fw)
    ov.set_ylim(0, fh)
    ov.axis("off")
    rows = [top - 0.86 - 0.2 * k for k in range(4)]
    tx = lx + 0.26
    ov.scatter([lx + 0.1], [rows[0]], s=30, c="white", edgecolors=INK, linewidths=0.55)
    ov.text(tx, rows[0], "ERA5 point on the route", va="center", fontsize=6.6, color=INK)
    ov.scatter([lx + 0.1], [rows[1]], s=40, marker="s", c="white", edgecolors=INK, linewidths=0.7)
    ov.text(tx, rows[1], "ERA5 point at port", va="center", fontsize=6.6, color=INK)
    ov.plot([lx, lx + 0.2], [rows[2]] * 2, color=INK, lw=1.0, ls=(0, (1, 1.6)))
    ov.text(tx, rows[2], "Port approach", va="center", fontsize=6.6, color=INK)
    ov.plot([lx + 0.1], [rows[3]], marker="o", ms=3.0, mfc="white", mec=INK, mew=0.8)
    ov.text(tx, rows[3], "End of an open-sea leg", va="center", fontsize=6.6, color=INK)

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"route_map.{ext}"), dpi=300)
    print("extent", round(lon0, 3), round(lon1, 3), round(lat0, 3), round(lat1, 3), "figure", round(fw, 3),
          round(fh, 3))
    print({k: round(v, 2) for k, v in u100.items()})


if __name__ == "__main__":
    main()
