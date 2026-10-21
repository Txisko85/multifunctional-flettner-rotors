#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapa de la ruta Portsmouth–Bilbao (A1, 21-10-2026): la ruta reconstruida (ruta.json, de ruta_geometria.py), sus dos
tramos de mar abierto, los diez puntos de muestreo de ERA5 y los dos puertos, coloreados por la velocidad media del
viento a 100 m en 2015–2024 (del fichero horario de ERA5).
Costa y fronteras: Natural Earth 1:10 m (dominio público), recortadas con ruta_recorte_costa.py (datos_mapa/).
Proyección de Mercator. Escribe fig_mapa.pdf y fig_mapa.png.
Uso: python ruta_figura_mapa.py ../datos/era5_horario_ruta.csv.gz
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
LON0, LON1, LAT0, LAT1 = -9.6, 1.0, 42.55, 51.3
# rampa secuencial de un solo tono (azul), de claro a oscuro
CMAP = LinearSegmentedColormap.from_list("azul", ["#d9e8f8", "#8fbaf0", "#2a78d6", "#123f7a"])
VMIN, VMAX = 4.0, 10.0


def merc(lat):
    """Ordenada de Mercator en grados (la abscisa es la longitud), para un aspecto 1:1."""
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
            if row[0] == "punto":
                continue
            s, n = acc.get(row[0], (0.0, 0))
            acc[row[0]] = (s + float(row[3]), n + 1)
    return {k: s / n for k, (s, n) in acc.items()}


def main():
    era5 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "datos", "era5_horario_ruta.csv.gz")
    u100 = mean_wind(era5)
    ruta = json.load(open(os.path.join(HERE, "ruta.json"), encoding="utf-8"))
    land = json.load(open(os.path.join(HERE, "datos_mapa", "ne_10m_tierra.geojson"), encoding="utf-8"))
    borders = json.load(open(os.path.join(HERE, "datos_mapa", "ne_10m_fronteras.geojson"), encoding="utf-8"))

    plt.rcParams.update({"font.size": 7.0})
    fig = plt.figure(figsize=(4.45, 5.25))
    ax = fig.add_axes([0.105, 0.075, 0.875, 0.905])
    ax.set_facecolor(SEA)

    # retícula
    for lon in range(-8, 1, 2):
        ax.plot([lon, lon], [merc(LAT0), merc(LAT1)], color=GRID, lw=0.5, zorder=0.5)
    for lat in range(44, 52, 2):
        ax.plot([LON0, LON1], [merc(lat)] * 2, color=GRID, lw=0.5, zorder=0.5)

    # tierra y fronteras
    for feat in land["features"]:
        g = feat["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for poly in polys:
            ax.add_patch(PathPatch(ring_path(poly), fc=LAND, ec=COAST, lw=0.35, zorder=1))
    for feat in borders["features"]:
        g = feat["geometry"]
        lines = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
        for line in lines:
            xy = np.array(line)
            ax.plot(xy[:, 0], merc(xy[:, 1]), color=COAST, lw=0.5, ls=(0, (3, 1.5)), zorder=1.2)

    # ruta: aproximaciones a puerto (punteadas) y tramos de mar abierto
    coords = {"Portsmouth": (50.790, -1.110), "Nab": (50.670, -0.950), "St Catherine S": (50.500, -1.300),
              "Start Point S": (49.950, -3.500), "Ushant SW": (48.650, -5.750), "Bilbao approach": (43.450, -3.080),
              "Bilbao": (43.360, -3.030)}       # los waypoints de ruta_geometria.py
    seq = ["Portsmouth", "Nab", "St Catherine S", "Start Point S", "Ushant SW", "Bilbao approach", "Bilbao"]
    lat = np.array([coords[k][0] for k in seq])
    lon = np.array([coords[k][1] for k in seq])
    y = merc(lat)
    ax.plot(lon[:3], y[:3], color=INK, lw=1.0, ls=(0, (1, 1.6)), zorder=3, solid_capstyle="round")
    ax.plot(lon[5:], y[5:], color=INK, lw=1.0, ls=(0, (1, 1.6)), zorder=3)
    ax.plot(lon[2:6], y[2:6], color=INK, lw=1.3, zorder=3, solid_joinstyle="round")
    # límites de los tramos: comienzo del Canal, Ushant y final del golfo
    for k in ("St Catherine S", "Ushant SW", "Bilbao approach"):
        i = seq.index(k)
        ax.plot(lon[i], y[i], marker="o", ms=3.0, mfc="white", mec=INK, mew=0.8, zorder=4)

    # puntos de muestreo y puertos, coloreados por el viento medio a 100 m
    norm = Normalize(VMIN, VMAX)
    names = ["C1", "C2", "C3", "C4", "B1", "B2", "B3", "B4", "B5", "B6"]
    for name, p in zip(names, ruta["puntos"]):
        ax.scatter(p["lon"], merc(p["lat"]), s=34, c=[CMAP(norm(u100[name]))], edgecolors=INK, linewidths=0.55,
                   zorder=5)
    for p, key in ((ruta["puertos"][0], "PORTSMOUTH"), (ruta["puertos"][1], "BILBAO")):
        ax.scatter(p["lon"], merc(p["lat"]), s=46, marker="s", c=[CMAP(norm(u100[key]))], edgecolors=INK,
                   linewidths=0.7, zorder=6)

    # rótulos
    lab = dict(color=INK, fontsize=7.2, zorder=7)
    ax.text(-0.83, merc(50.93), "Portsmouth", ha="left", va="bottom", **lab)
    ax.text(-2.72, merc(43.10), "Bilbao", ha="left", va="top", **lab)
    ax.text(-5.06, merc(48.53), "Ushant", ha="center", va="bottom", color=MUTED, fontsize=6.6, zorder=7)
    sea = dict(color="#5d7891", fontsize=7.2, style="italic", ha="center", va="center", zorder=2)
    ax.text(-1.15, merc(50.12), "English\nChannel", linespacing=1.0, **sea)
    ax.text(-2.55, merc(45.55), "Bay of\nBiscay", linespacing=1.0, **sea)
    ax.text(-7.6, merc(49.55), "Celtic Sea", **sea)
    country = dict(color=QUIET, fontsize=6.4, ha="center", va="center", zorder=2)

    def spaced(t):
        return "\u2002".join("\u2009".join(w) for w in t.split())

    ax.text(-3.0, merc(51.1), spaced("UNITED KINGDOM"), **country)
    ax.text(-0.3, merc(47.45), spaced("FRANCE"), **country)
    ax.text(-5.6, merc(42.92), spaced("SPAIN"), **country)

    def along(i0, i1, frac, off, text):
        """Rótulo paralelo al tramo i0-i1, desplazado off grados hacia su izquierda (vista de i0 a i1)."""
        x0, y0, x1, y1 = lon[i0], y[i0], lon[i1], y[i1]
        dx, dy = x1 - x0, y1 - y0
        ln = math.hypot(dx, dy)
        nx, ny = -dy / ln, dx / ln
        ang = math.degrees(math.atan2(dy, dx))
        if ang > 90 or ang < -90:
            ang -= 180 * (1 if ang > 0 else -1)
        ax.text(x0 + frac * dx + off * nx, y0 + frac * dy + off * ny, text, rotation=ang, rotation_mode="anchor",
                ha="center", va="center", color=INK, fontsize=6.8, zorder=7)

    along(3, 4, 0.42, -0.36, "Channel leg, 208 nm")
    along(4, 5, 0.47, -0.34, "Bay of Biscay leg, 331 nm")

    # barra de escala de 100 millas náuticas a 44,4°N
    lat_s, lon_s = 44.05, -9.25
    dlon = (100 / 60) / math.cos(math.radians(lat_s))
    ys = merc(lat_s)
    for k in range(2):
        ax.add_patch(Rectangle((lon_s + k * dlon / 2, ys), dlon / 2, 0.09, fc=INK if k == 0 else "white", ec=INK,
                               lw=0.5, zorder=7))
    ax.text(lon_s, ys + 0.2, "0", ha="center", va="bottom", fontsize=6.4, color=INK, zorder=7)
    ax.text(lon_s + dlon, ys + 0.2, "100 nm", ha="center", va="bottom", fontsize=6.4, color=INK, zorder=7)

    # barra de color
    cax = fig.add_axes([0.16, 0.30, 0.32, 0.018])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_ticks([4, 6, 8, 10])
    cb.ax.tick_params(labelsize=6.4, length=2, width=0.5, colors=INK)
    cb.outline.set_linewidth(0.5)
    cax.set_title("Mean wind speed at 100 m,\n2015–2024 (m/s)", fontsize=6.6, color=INK, pad=3, linespacing=1.1)
    # leyenda de símbolos
    lx, ly = -9.25, merc(47.45)
    ax.scatter([lx + 0.12], [ly], s=30, c="white", edgecolors=INK, linewidths=0.55, zorder=7)
    ax.text(lx + 0.4, ly, "ERA5 point on the route", va="center", fontsize=6.6, color=INK, zorder=7)
    ax.scatter([lx + 0.12], [ly - 0.42], s=40, marker="s", c="white", edgecolors=INK, linewidths=0.7, zorder=7)
    ax.text(lx + 0.4, ly - 0.42, "ERA5 point at port", va="center", fontsize=6.6, color=INK, zorder=7)
    ax.plot([lx - 0.05, lx + 0.3], [ly - 0.84] * 2, color=INK, lw=1.0, ls=(0, (1, 1.6)), zorder=7)
    ax.text(lx + 0.4, ly - 0.84, "Port approach", va="center", fontsize=6.6, color=INK, zorder=7)
    ax.plot([lx + 0.12], [ly - 1.26], marker="o", ms=3.0, mfc="white", mec=INK, mew=0.8, zorder=7)
    ax.text(lx + 0.4, ly - 1.26, "End of an open-sea leg", va="center", fontsize=6.6, color=INK, zorder=7)

    # ejes
    ax.set_xlim(LON0, LON1)
    ax.set_ylim(merc(LAT0), merc(LAT1))
    ax.set_aspect("equal")
    ax.set_xticks(range(-8, 1, 2))
    ax.set_xticklabels([f"{abs(v)}°W" if v < 0 else "0°" for v in range(-8, 1, 2)])
    ax.set_yticks([float(merc(v)) for v in range(44, 52, 2)])
    ax.set_yticklabels([f"{v}°N" for v in range(44, 52, 2)])
    ax.tick_params(labelsize=6.6, length=2.5, width=0.5, colors=MUTED)
    for s in ax.spines.values():
        s.set_linewidth(0.6)
        s.set_color(MUTED)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"fig_mapa.{ext}"), dpi=300)
    print({k: round(v, 2) for k, v in u100.items()})


if __name__ == "__main__":
    main()
