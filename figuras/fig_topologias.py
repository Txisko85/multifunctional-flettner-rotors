#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figura de las topologías de A1 (20-10-2026): esquemas en alzado, a la escala del caso de estudio, de R1 y de
T1-T5 en sus modos de operación, coloreados según la configuración aerodinámica de cada modo (C0-C3).
Inspirada en las realizaciones de la patente ES 2 579 555 B1 y simplificada: la turbina se dibuja como la
H-Darrieus de tres palas que modela el artículo y no se dibujan los mecanismos.
Geometría (apartado 2.4): rotor 4 m x 24 m con placa de 8 m; turbina y prolongación 3,4 m x 6 m; Savonius con
cubetas de cuerda d = D_F y solape e = 0,2 d (D_S = 7,2 m), cada mitad desplazada (d - e)/2 = 1,6 m.
Vista oblicua ligeramente desde arriba (elipses con aplastamiento K). Bajo los modos Savonius, la planta; en esos
modos no se dibujan las medias placas de extremo.
Escribe fig_topologias.pdf y fig_topologias.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle, Ellipse, Arc, Polygon, FancyBboxPatch
from matplotlib import font_manager

# Tipo de letra de las figuras del artículo (Inter); si no está instalado, el de Matplotlib por defecto.
if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, EDGE = "#1f1f1d", "#52514e", "#8a897f"
ROTOR, OUTLINE, PLATE = "#cfcdc4", "#3d3c38", "#f4f3ef"
CONF = {"C0": ("#6f6e68", "Rotor alone"), "C1": ("#2a78d6", "Rotor and turbine"),
        "C2": ("#e0702c", "Extended rotor"), "C3": ("#239271", "Savonius rotor")}
K = 0.30                      # aplastamiento de las elipses
R_F, H_F, R_E = 2.0, 24.0, 4.0
R_X, H_X = 1.7, 6.0
LW = 0.7


def mix(color, a):
    """Mezcla con blanco: a = 1 color puro, a = 0 blanco."""
    c = np.array(to_rgb(color))
    return tuple(1 - a * (1 - c))


def shaded_band(ax, x0, x1, z0, z1, color, convex=True, zo=2, n=28):
    """Banda vertical con sombreado de cilindro (luz desde delante a la izquierda)."""
    edges = np.linspace(x0, x1, n + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        u = ((a + b) / 2 - x0) / (x1 - x0) * 2 - 1           # -1 ... 1 a lo ancho
        th = np.arcsin(np.clip(u, -1, 1))
        light = np.cos(th + (0.6 if convex else -0.6))
        ax.add_patch(Rectangle((a, z0), b - a + 1e-3, z1 - z0, fc=mix(color, 0.55 + 0.45 * (1 - light)),
                               ec="none", zorder=zo))


def cylinder(ax, xc, z0, z1, r, color, edge=OUTLINE, zo=2):
    shaded_band(ax, xc - r, xc + r, z0, z1, color, zo=zo)
    ax.add_patch(Ellipse((xc, z0), 2 * r, 2 * r * K, fc=mix(color, 0.8), ec="none", zorder=zo - 0.05))
    ax.plot([xc - r, xc - r], [z0, z1], color=edge, lw=LW, zorder=zo + 0.1)
    ax.plot([xc + r, xc + r], [z0, z1], color=edge, lw=LW, zorder=zo + 0.1)
    ax.add_patch(Arc((xc, z0), 2 * r, 2 * r * K, theta1=180, theta2=360, color=edge, lw=LW, zorder=zo + 0.1))
    ax.add_patch(Ellipse((xc, z1), 2 * r, 2 * r * K, fc=mix(color, 0.45), ec=edge, lw=LW, zorder=zo + 0.2))


def plate(ax, xc, z, r, zo=4):
    t = 0.35
    ax.add_patch(Rectangle((xc - r, z), 2 * r, t, fc=PLATE, ec="none", zorder=zo))
    ax.add_patch(Arc((xc, z), 2 * r, 2 * r * K, theta1=180, theta2=360, color=OUTLINE, lw=LW, zorder=zo + 0.1))
    ax.plot([xc - r, xc - r], [z, z + t], color=OUTLINE, lw=LW, zorder=zo + 0.1)
    ax.plot([xc + r, xc + r], [z, z + t], color=OUTLINE, lw=LW, zorder=zo + 0.1)
    ax.add_patch(Ellipse((xc, z + t), 2 * r, 2 * r * K, fc=PLATE, ec=OUTLINE, lw=LW, zorder=zo + 0.2))


def half_plate(ax, xc, z, r, back, zo=4.5):
    th = np.radians(np.linspace(0, 180, 60) if back else np.linspace(180, 360, 60))
    pts = list(zip(xc + r * np.cos(th), z + 0.35 + r * K * np.sin(th)))
    ax.add_patch(Polygon(pts, closed=True, fc=PLATE, ec=OUTLINE, lw=LW, zorder=zo + (0 if back else 0.1)))


def turbine(ax, xc, z0, color, zo=5):
    """H-Darrieus de tres palas (3,4 m x 6 m) con su eje."""
    z1 = z0 + H_X
    ax.plot([xc, xc], [z0 - 0.3, z1 + 0.2], color=color, lw=1.0, zorder=zo)
    for zs in (z0 + 1.4, z1 - 1.4):
        ax.add_patch(Ellipse((xc, zs), 2 * R_X, 2 * R_X * K, fc="none", ec=color, lw=0.3, alpha=0.6,
                             zorder=zo - 0.5))
    for xb, front in ((-1.47, True), (1.47, True), (0.0, False)):
        for zs in (z0 + 1.4, z1 - 1.4):
            ax.plot([xc, xc + xb], [zs, zs + (-0.25 if front else 0.25)], color=color, lw=0.6, zorder=zo)
        ax.add_patch(Rectangle((xc + xb - 0.12, z0), 0.24, H_X, fc=color if front else "white", ec=color,
                               lw=0.5, zorder=zo + (0.5 if front else -0.6)))


def stowed(ax, xc, z0, z1, r, color):
    ax.add_patch(Rectangle((xc - r, z0), 2 * r, z1 - z0, fc="none", ec=color, lw=0.8, ls=(0, (2.2, 1.6)), zorder=6))


def savonius(ax, xc, z0, z1, r, color, plates=True, r_plate=R_E):
    """Rotor abierto: dos cubetas semicilíndricas desplazadas +-(d - e)/2 a lo largo de la cuerda; se ve la cara
    cóncava de la trasera (izquierda) y la convexa de la delantera (derecha)."""
    s = 0.8 * r
    shaded_band(ax, xc - s - r, xc - s + r, z0, z1, mix(color, 0.55), convex=False, zo=2)
    ax.plot([xc - s - r, xc - s - r], [z0, z1], color=OUTLINE, lw=LW, zorder=2.1)
    shaded_band(ax, xc + s - r, xc + s + r, z0, z1, color, convex=True, zo=3)
    for xx in (xc + s - r, xc + s + r):
        ax.plot([xx, xx], [z0, z1], color=OUTLINE, lw=LW, zorder=3.1)
    ax.add_patch(Arc((xc + s, z0), 2 * r, 2 * r * K, theta1=180, theta2=360, color=OUTLINE, lw=LW, zorder=3.1))
    # bordes superiores de las cubetas (la S en vista oblicua)
    th = np.radians(np.linspace(0, 180, 40))
    ax.plot(xc - s + r * np.cos(th), z1 + r * K * np.sin(th), color=OUTLINE, lw=LW, zorder=4)
    th = np.radians(np.linspace(180, 360, 40))
    ax.plot(xc + s + r * np.cos(th), z1 + r * K * np.sin(th), color=OUTLINE, lw=LW, zorder=4)
    if plates:
        half_plate(ax, xc - s, z1, r_plate, back=True)
        half_plate(ax, xc + s, z1, r_plate, back=False)


def savonius_plan(ax, xc, zc, r, color, scale=0.6):
    """Planta: dos cubetas semicirculares de cuerda 2r con solape 0,2 * 2r, dentro de D_S = 3,6 r."""
    sh, rr = 0.8 * r * scale, r * scale
    ax.add_patch(Ellipse((xc, zc), 2 * (rr + sh), 2 * (rr + sh), fc="none", ec=EDGE, lw=0.4, ls=(0, (1.5, 1.5)),
                         zorder=5))
    ax.add_patch(Arc((xc - sh, zc), 2 * rr, 2 * rr, theta1=0, theta2=180, color=color, lw=1.4, zorder=6))
    ax.add_patch(Arc((xc + sh, zc), 2 * rr, 2 * rr, theta1=180, theta2=360, color=color, lw=1.4, zorder=6))


def rotor(ax, xc, top_plate=True):
    cylinder(ax, xc, 0.0, H_F, R_F, ROTOR)
    if top_plate:
        plate(ax, xc, H_F, R_E)


def tag(ax, x, y, conf):
    col, _ = CONF[conf]
    ax.add_patch(FancyBboxPatch((x - 1.7, y - 1.0), 3.4, 2.0, boxstyle="round,pad=0.05,rounding_size=0.6",
                                fc=col, ec="none", zorder=7, clip_on=False))
    ax.text(x, y, conf, ha="center", va="center", fontsize=6.6, color="white", fontweight="bold", zorder=8,
            clip_on=False)


def main():
    plt.rcParams.update({"font.size": 6.5})
    fig = plt.figure(figsize=(7.48, 3.3))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_aspect("equal")
    ax.axis("off")
    blue, orange, green = CONF["C1"][0], CONF["C2"][0], CONF["C3"][0]

    GAP, SLOT, WIDE = 4.5, 10.0, 12.5
    x = 7.0
    groups = []

    def place(width):
        nonlocal x
        c = x + width / 2
        x += width
        return c

    # R1
    g0 = x; c = place(SLOT); rotor(ax, c)
    groups.append(("R1", "Conventional\nrotor", g0, x, [(c, "Rotor", "C0")])); x += GAP
    # T1
    g0 = x; c = place(SLOT); rotor(ax, c); turbine(ax, c, H_F + 0.65, blue)
    groups.append(("T1", "Fixed turbine", g0, x, [(c, "Rotor and\nturbine", "C1")])); x += GAP
    # T2
    g0 = x; c1 = place(SLOT); rotor(ax, c1); turbine(ax, c1, H_F + 0.65, blue)
    c2 = place(WIDE); savonius(ax, c2, 0.0, H_F, R_F, green, plates=False); turbine(ax, c2, H_F + 0.65, blue)
    groups.append(("T2", "Fixed turbine and Savonius", g0, x,
                   [(c1, "Rotor and\nturbine", "C1"), (c2, "Savonius and\nturbine", "C3")])); x += GAP
    # T3 (y T3x)
    g0 = x; c1 = place(SLOT); rotor(ax, c1); stowed(ax, c1, H_F - H_X - 0.4, H_F - 0.4, R_X, orange)
    c2 = place(SLOT); cylinder(ax, c2, 0.0, H_F, R_F, ROTOR); cylinder(ax, c2, H_F, H_F + H_X, R_X, orange, zo=3)
    plate(ax, c2, H_F + H_X, R_E)
    c3 = place(WIDE); savonius(ax, c3, 0.0, H_F, R_F, green, plates=False)
    savonius(ax, c3, H_F, H_F + H_X, R_X, green, plates=False)
    groups.append(("T3", "Extension and Savonius", g0, x,
                   [(c1, "Extension\nstowed", "C0"), (c2, "Extension\nraised", "C2"),
                    (c3, "Both cylinders\nas Savonius", "C3")]))
    t3x = (c1, c2); x += GAP
    # T4
    g0 = x; c1 = place(SLOT); rotor(ax, c1); stowed(ax, c1, H_F - H_X - 0.4, H_F - 0.4, R_X, blue)
    c2 = place(SLOT); rotor(ax, c2)
    ax.plot([c2, c2], [H_F + 0.35, H_F + 1.9], color=blue, lw=1.0, zorder=5)
    turbine(ax, c2, H_F + 1.9, blue)
    groups.append(("T4", "Retractable turbine", g0, x,
                   [(c1, "Turbine\nstowed", "C0"), (c2, "Turbine\nraised", "C1")])); x += GAP
    # T5
    g0 = x; c1 = place(SLOT); rotor(ax, c1)
    c2 = place(WIDE); savonius(ax, c2, 0.0, H_F, R_F, green, plates=False)
    groups.append(("T5", "Savonius conversion", g0, x, [(c1, "Rotor", "C0"), (c2, "Rotor opened\nas Savonius", "C3")]))
    x_end = x

    for name, title, a, b, modes in groups:
        ax.add_patch(Rectangle((a - 0.8, -1.0), b - a + 1.6, 1.0, fc="#efeee9", ec="none", zorder=0.5))
        ax.plot([a - 0.8, b + 0.8], [0, 0], color=EDGE, lw=1.1, zorder=1, solid_capstyle="butt")
        ax.text((a + b) / 2, 38.2, name, ha="center", va="bottom", fontsize=9.0, fontweight="bold", color=INK)
        ax.text((a + b) / 2, 37.4, title, ha="center", va="top", fontsize=7.0, color=MUTED, linespacing=1.05)
        for c, label, conf in modes:
            if conf == "C3":
                savonius_plan(ax, c, -3.7, R_F, green)
            ax.text(c, -6.9, label, ha="center", va="top", fontsize=6.7, color=INK, linespacing=1.1)
            tag(ax, c, -13.0, conf)
    a, b = t3x
    ax.plot([a - 3.5, a - 3.5, b + 3.5, b + 3.5], [-15.0, -15.5, -15.5, -15.0], color=MUTED, lw=0.6)
    ax.text((a + b) / 2, -15.9, "T3x: extension only", ha="center", va="top", fontsize=6.6, color=MUTED)

    # escala y cotas
    xs = 2.0
    ax.plot([xs, xs], [0, 30], color=INK, lw=0.6)
    for z in (0, 10, 20, 30):
        ax.plot([xs - 0.5, xs + 0.5], [z, z], color=INK, lw=0.6)
        ax.text(xs - 0.9, z, f"{z}", ha="right", va="center", fontsize=6.2, color=MUTED)
    ax.text(xs - 0.9, 32.4, "m", ha="right", va="center", fontsize=6.2, color=MUTED)
    c = groups[0][4][0][0]
    ax.annotate("", xy=(c - R_F, 6), xytext=(c + R_F, 6),
                arrowprops=dict(arrowstyle="<->", lw=0.5, color=INK, mutation_scale=4), zorder=9)
    ax.text(c, 6.5, "4 m", ha="center", va="bottom", fontsize=6.0, color=INK, zorder=9)
    ax.text(c, H_F + 1.6, "8 m endplate", ha="center", va="bottom", fontsize=6.0, color=MUTED)

    # leyenda
    names = [(conf, nm) for conf, (col, nm) in CONF.items()]
    total = x_end - 9.0
    for i, (conf, nm) in enumerate(names):
        xx = 9.0 + i * total / 4.0
        tag(ax, xx, -21.0, conf)
        ax.text(xx + 2.4, -21.0, nm, ha="left", va="center", fontsize=7.0, color=INK)

    ax.set_xlim(-2.5, x_end + 1.5)
    ax.set_ylim(-23.6, 42.6)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"fig_topologias.{ext}"), dpi=300)


if __name__ == "__main__":
    main()
