#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figura de las topologías de A1 (versión del 21-10-2026): esquemas en alzado, a la escala del caso de estudio, de R1
y de T1-T5 en sus modos de operación, coloreados según la configuración aerodinámica de cada modo (C0-C3).
Inspirada en las realizaciones de la patente ES 2 579 555 B1 y simplificada: la turbina se dibuja como la
H-Darrieus de tres palas que modela el artículo y no se dibujan los mecanismos.
Geometría (apartado 2.4): rotor 4 m x 24 m con placa de 8 m; turbina y prolongación 3,4 m x 6 m; Savonius con
cubetas de cuerda d = D_F y solape e = 0,2 d (D_S = 7,2 m), cada mitad desplazada (d - e)/2 = 1,6 m.
Vista oblicua ligeramente desde arriba (las circunferencias se ven como elipses de aplastamiento K) e iluminación
desde delante a la izquierda: cada superficie se sombrea con un modelo difuso más especular, en bandas verticales
recortadas a su contorno. Las palas, los brazos y el eje de la turbina se dibujan algo más gruesos que a escala
para que se vean. Bajo los modos Savonius, la planta; en esos modos no se dibujan las medias placas de extremo.
Escribe fig_topologias.pdf y fig_topologias.png.
"""
import math
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
ROTOR, OUTLINE, PLATE = "#c9c7be", "#3d3c38", "#e9e8e2"
DECK_TOP, DECK_FRONT = "#efeee9", "#d6d4cc"
CONF = {"C0": ("#6f6e68", "Rotor alone"), "C1": ("#2a78d6", "Rotor and turbine"),
        "C2": ("#e0702c", "Extended rotor"), "C3": ("#239271", "Savonius rotor")}
TAG = {"C0": "#6f6e68", "C1": "#2a78d6", "C2": "#c2581e", "C3": "#1d7d61"}   # algo más oscuros, por el texto blanco
K = 0.30                      # aplastamiento de las elipses (vista desde arriba)
R_F, H_F, R_E, T_E = 2.0, 24.0, 4.0, 0.35
R_X, H_X = 1.7, 6.0
LW = 0.55
NB = 90                       # bandas de sombreado por superficie
ALPHA_L, BETA_H = math.radians(40), math.radians(20)


# ---------------------------------------------------------------- color y sombreado
def mix(c1, c2, t):
    a, b = np.array(to_rgb(c1)), np.array(to_rgb(c2))
    return tuple(np.clip(a + (b - a) * t, 0, 1))


def tone(base, intensity):
    """Color de la superficie según la intensidad de luz (0 sombra, 0,75 tono base, 1,3 brillo)."""
    dark, light = mix(base, "#000000", 0.42), mix(base, "#ffffff", 0.78)
    if intensity <= 0.75:
        return mix(dark, base, intensity / 0.75)
    return mix(base, light, min(1.0, (intensity - 0.75) / 0.55))


def light(u, concave=False):
    """Intensidad en la posición u (-1 izquierda, 1 derecha) de un cilindro visto de frente; concave: cara interior."""
    u = float(np.clip(u, -0.999, 0.999))
    if not concave:
        phi = math.asin(u)
        d = max(0.0, math.cos(phi + ALPHA_L))
        s = max(0.0, math.cos(phi + BETA_H)) ** 40
        return 0.28 + 0.62 * d + 0.42 * s
    w = math.sqrt(1 - u * u)                     # cara interior de la mitad trasera: la normal mira al eje
    d = max(0.0, u * math.sin(ALPHA_L) + w * math.cos(ALPHA_L))
    return 0.85 * (0.24 + 0.62 * d)


def shade(ax, pts, x0, x1, base, zo, concave=False, n=NB, lift=1.0):
    """Rellena el polígono pts con bandas verticales sombreadas como un cilindro entre x0 y x1."""
    clip = Polygon(pts, closed=True, transform=ax.transData)
    zs = [p[1] for p in pts]
    z0, z1 = min(zs) - 0.05, max(zs) + 0.05
    edges = np.linspace(x0, x1, n + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        u = ((a + b) / 2 - x0) / (x1 - x0) * 2 - 1
        r = Rectangle((a, z0), (b - a) * 1.06, z1 - z0, fc=tone(base, light(u, concave) * lift), ec="none", zorder=zo)
        ax.add_patch(r)
        r.set_clip_path(clip)


def arc_pts(xc, z, r, t0, t1, n=72):
    th = np.radians(np.linspace(t0, t1, n))
    return list(zip(xc + r * np.cos(th), z + K * r * np.sin(th)))


def line(ax, pts, color=OUTLINE, lw=LW, zo=5, **kw):
    p = np.array(pts)
    ax.plot(p[:, 0], p[:, 1], color=color, lw=lw, zorder=zo, solid_capstyle="round", **kw)


# ---------------------------------------------------------------- piezas
def side(ax, xc, z0, z1, r, base, zo):
    """Cara lateral visible de un cilindro vertical (entre los arcos delanteros de abajo y de arriba)."""
    pts = arc_pts(xc, z0, r, 180, 360) + arc_pts(xc, z1, r, 360, 180)
    shade(ax, pts, xc - r, xc + r, base, zo)
    line(ax, [(xc - r, z0), (xc - r, z1)], zo=zo + 0.3)
    line(ax, [(xc + r, z0), (xc + r, z1)], zo=zo + 0.3)
    line(ax, arc_pts(xc, z0, r, 180, 360), zo=zo + 0.3)


def top_face(ax, xc, z, r, base, zo, edge=True):
    """Cara superior (elipse completa), algo más clara por delante a la izquierda."""
    pts = arc_pts(xc, z, r, 0, 360, 120)
    clip = Polygon(pts, closed=True, transform=ax.transData)
    n = 40
    for i in range(n):
        f = i / (n - 1)
        rr = r * (1 - 0.85 * f)
        e = Ellipse((xc - 0.25 * r * f, z - 0.25 * K * r * f), 2 * rr, 2 * rr * K,
                    fc=mix(base, "#ffffff", 0.15 + 0.45 * f), ec="none", zorder=zo)
        ax.add_patch(e)
        e.set_clip_path(clip)
    if edge:
        line(ax, pts, zo=zo + 0.2)


def seams(ax, xc, r, zs, zo):
    for z in zs:
        line(ax, arc_pts(xc, z, r, 182, 358), color="#000000", lw=0.35, zo=zo, alpha=0.30)
        line(ax, arc_pts(xc, z - 0.13, r, 184, 356), color="#ffffff", lw=0.35, zo=zo, alpha=0.45)


def plate_shadow(ax, xc, z, r, zo, depth=1.6):
    """Sombra de la placa de extremo sobre el cilindro que tiene debajo."""
    for i in range(10):
        f0, f1 = i / 10, (i + 1) / 10
        pts = arc_pts(xc, z - depth * f0, r, 180, 360) + arc_pts(xc, z - depth * f1, r, 360, 180)
        ax.add_patch(Polygon(pts, closed=True, fc="#000000", ec="none", alpha=0.16 * (1 - f0) ** 1.6, zorder=zo))


def endplate(ax, xc, z, r=R_E, zo=6):
    rim = arc_pts(xc, z, r, 180, 360) + arc_pts(xc, z + T_E, r, 360, 180)
    shade(ax, rim, xc - r, xc + r, PLATE, zo, n=60)
    line(ax, arc_pts(xc, z, r, 180, 360), zo=zo + 0.3)
    line(ax, [(xc - r, z), (xc - r, z + T_E)], zo=zo + 0.3)
    line(ax, [(xc + r, z), (xc + r, z + T_E)], zo=zo + 0.3)
    top_face(ax, xc, z + T_E, r, PLATE, zo + 0.4)


def contact_shadow(ax, xc, r, zo=1.4):
    for i in range(6):
        f = i / 5
        ax.add_patch(Ellipse((xc + 0.25 * r, 0.05), 2 * r * (1.55 - 0.45 * f), 2 * r * K * (1.55 - 0.45 * f),
                             fc="#000000", ec="none", alpha=0.035, zorder=zo))


def rotor(ax, xc, top_plate=True, zo=2):
    contact_shadow(ax, xc, R_F)
    side(ax, xc, 0.0, H_F, R_F, ROTOR, zo)
    seams(ax, xc, R_F, (6.0, 12.0, 18.0), zo + 0.2)
    if top_plate:
        plate_shadow(ax, xc, H_F, R_F, zo + 0.25)
        endplate(ax, xc, H_F)


def extension(ax, xc, color):
    """T3 con la prolongación subida: cara superior del rotor, prolongación y placa encima."""
    contact_shadow(ax, xc, R_F)
    side(ax, xc, 0.0, H_F, R_F, ROTOR, 2)
    seams(ax, xc, R_F, (6.0, 12.0, 18.0), 2.2)
    top_face(ax, xc, H_F, R_F, ROTOR, 2.5)
    side(ax, xc, H_F, H_F + H_X, R_X, color, 3)
    line(ax, arc_pts(xc, H_F + 0.12, R_X, 180, 360), color="#000000", lw=0.5, zo=3.4, alpha=0.35)
    plate_shadow(ax, xc, H_F + H_X, R_X, 3.25)
    endplate(ax, xc, H_F + H_X)


def stowed(ax, xc, z0, z1, r, color, blades=False):
    """Dispositivo recogido dentro del rotor, en transparencia (contorno discontinuo)."""
    pts = arc_pts(xc, z0, r, 180, 360) + arc_pts(xc, z1, r, 360, 180)
    ax.add_patch(Polygon(pts, closed=True, fc=color, ec="none", alpha=0.13, zorder=4))
    dash = dict(ls=(0, (2.2, 1.4)))
    line(ax, [(xc - r, z0), (xc - r, z1)], color=color, lw=0.8, zo=4.2, **dash)
    line(ax, [(xc + r, z0), (xc + r, z1)], color=color, lw=0.8, zo=4.2, **dash)
    line(ax, arc_pts(xc, z0, r, 180, 360), color=color, lw=0.8, zo=4.2, **dash)
    line(ax, arc_pts(xc, z1, r, 0, 360, 120), color=color, lw=0.8, zo=4.2, **dash)
    if blades:
        for th in (30.0, 150.0, 270.0):
            xb = xc + r * math.cos(math.radians(th))
            dz = K * r * math.sin(math.radians(th))
            line(ax, [(xb, z0 + 0.25 + dz), (xb, z1 - 0.25 + dz)], color=color, lw=1.2, zo=4.3, alpha=0.65)


def shaft(ax, xc, z0, z1, r=0.17, color="#7f8fa3", zo=5.0):
    side(ax, xc, z0, z1, r, color, zo)
    top_face(ax, xc, z1, r, color, zo + 0.35)


def turbine(ax, xc, z0, color, shaft_from, zb=3.0):
    """H-Darrieus de tres palas (3,4 m x 6 m): palas, dos niveles de brazos y eje, en vista oblicua. Se apoya
    directamente sobre la placa del rotor (o sobre las cubetas), como en las figuras 1 y 21 de la patente y como en el
    modelo (centro de la turbina a 55 m); zb eleva su orden de dibujo por encima de la placa."""
    z1 = z0 + H_X
    levels = (z0 + 1.25, z1 - 1.25)
    thetas = (30.0, 150.0, 270.0)                 # dos palas detrás, a los lados; una delante, sobre el eje
    c_draw = 0.52                                 # cuerda dibujada (m)
    for zs in levels:                             # recorrido de las puntas de los brazos
        line(ax, arc_pts(xc, zs, R_X, 0, 360, 120), color=color, lw=0.35, zo=4.6 + zb, alpha=0.45)

    def blade(th, zo, back):
        t = math.radians(th)
        xb, dz = xc + R_X * math.cos(t), K * R_X * math.sin(t)
        w = max(0.17, c_draw * abs(math.sin(t)))
        pts = [(xb - w / 2, z0 + dz), (xb + w / 2, z0 + dz), (xb + w / 2, z1 + dz), (xb - w / 2, z1 + dz)]
        shade(ax, pts, xb - w / 2, xb + w / 2, mix(color, "#ffffff", 0.32) if back else color, zo, n=14)
        ax.add_patch(Polygon(pts, closed=True, fc="none", ec=mix(color, "#000000", 0.45), lw=0.35, zorder=zo + 0.1))

    def struts(th, zo, back):
        t = math.radians(th)
        xb, dz = xc + R_X * math.cos(t), K * R_X * math.sin(t)
        for zs in levels:
            line(ax, [(xc, zs), (xb, zs + dz)], color=mix(color, "#000000", 0.25 if not back else 0.0), lw=1.0,
                 zo=zo, alpha=0.75 if back else 1.0)

    for th in thetas[:2]:
        struts(th, 4.7 + zb, True)
        blade(th, 4.75 + zb, True)
    shaft(ax, xc, shaft_from, z1 + 0.35, zo=5.0 + zb)
    for zs in levels:                              # cubos
        side(ax, xc, zs - 0.18, zs + 0.18, 0.3, mix(color, "#000000", 0.2), 5.5 + zb)
    struts(thetas[2], 5.8 + zb, False)
    blade(thetas[2], 5.9 + zb, False)


def savonius(ax, xc, z0, z1, r, color):
    """Rotor abierto: dos cubetas semicilíndricas desplazadas +-(d - e)/2 a lo largo de la cuerda. Se ve la cara
    cóncava de la de la izquierda (abierta hacia el observador) y la convexa de la de la derecha."""
    s = 0.8 * r
    xl, xr = xc - s, xc + s
    # cubeta izquierda: cara interior de su pared trasera
    inner = [(xl - r, z0), (xl - r, z1)] + arc_pts(xl, z1, r, 180, 0) + [(xl + r, z1), (xl + r, z0)] + \
        arc_pts(xl, z0, r, 0, 180)
    shade(ax, inner, xl - r, xl + r, color, 2.0, concave=True)
    line(ax, arc_pts(xl, z0, r, 0, 180), color=OUTLINE, lw=LW * 0.8, zo=2.2, alpha=0.6)
    line(ax, [(xl - r, z0), (xl - r, z1)], zo=2.3)
    rim = arc_pts(xl, z1, r, 180, 0)
    line(ax, rim, color=mix(color, "#ffffff", 0.7), lw=1.8, zo=2.35)
    line(ax, rim, zo=2.4)
    # cubeta derecha: cara exterior convexa, delante
    outer = arc_pts(xr, z0, r, 180, 360) + arc_pts(xr, z1, r, 360, 180)
    shade(ax, outer, xr - r, xr + r, color, 3.0)
    mouth = [(xr - r, z1), (xr + r, z1)] + arc_pts(xr, z1, r, 360, 180)
    ax.add_patch(Polygon(mouth, closed=True, fc=mix(color, "#000000", 0.55), ec="none", zorder=3.2))
    line(ax, [(xr - r, z0), (xr - r, z1)], zo=3.3)
    line(ax, [(xr + r, z0), (xr + r, z1)], zo=3.3)
    line(ax, arc_pts(xr, z0, r, 180, 360), zo=3.3)
    rim = arc_pts(xr, z1, r, 180, 360)
    line(ax, rim, color=mix(color, "#ffffff", 0.7), lw=1.8, zo=3.35)
    line(ax, rim, zo=3.4)
    line(ax, [(xr - r, z1), (xr + r, z1)], color=OUTLINE, lw=0.35, zo=3.4, alpha=0.5)


def savonius_plan(ax, xc, zc, r, color, scale=0.58):
    """Planta: dos cubetas semicirculares de cuerda 2r con solape 0,2 * 2r, dentro de D_S = 3,6 r."""
    sh, rr = 0.8 * r * scale, r * scale
    ax.add_patch(Ellipse((xc, zc), 2 * (rr + sh), 2 * (rr + sh), fc="white", ec=EDGE, lw=0.4, ls=(0, (1.5, 1.5)),
                         zorder=5))
    for cx, t0, t1 in ((xc - sh, 0, 180), (xc + sh, 180, 360)):
        ax.add_patch(Arc((cx, zc), 2 * rr, 2 * rr, theta1=t0, theta2=t1, color=mix(color, "#000000", 0.15),
                         lw=1.6, zorder=6, capstyle="round"))
    ax.add_patch(Ellipse((xc, zc), 0.28, 0.28, fc=OUTLINE, ec="none", zorder=6))


def tag(ax, x, y, conf):
    ax.add_patch(FancyBboxPatch((x - 1.7, y - 1.0), 3.4, 2.0, boxstyle="round,pad=0.05,rounding_size=0.6",
                                fc=TAG[conf], ec="none", zorder=7, clip_on=False))
    ax.text(x, y, conf, ha="center", va="center", fontsize=6.6, color="white", fontweight="bold", zorder=8,
            clip_on=False)


def deck(ax, a, b):
    """Cubierta: superficie vista desde arriba (de y = -4,5 a 4,5 m) y canto delantero."""
    zf, zb, t = -K * 4.5, K * 4.5, 0.9
    n = 12
    for i in range(n):
        f0, f1 = i / n, (i + 1) / n
        ax.add_patch(Rectangle((a, zf + (zb - zf) * f0), b - a, (zb - zf) / n * 1.05,
                               fc=mix(DECK_TOP, "#ffffff", 0.5 * (1 - f0)), ec="none", zorder=0.5))
    ax.add_patch(Rectangle((a, zf - t), b - a, t, fc=DECK_FRONT, ec="none", zorder=0.6))
    ax.plot([a, b], [zf, zf], color=EDGE, lw=0.6, zorder=0.7)
    ax.plot([a, b], [zf - t, zf - t], color=EDGE, lw=0.6, zorder=0.7)


def main():
    plt.rcParams.update({"font.size": 6.5})
    fig = plt.figure(figsize=(7.48, 3.55))
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
    g0 = x; c = place(SLOT); rotor(ax, c); turbine(ax, c, H_F + T_E + 0.15, blue, H_F + T_E)
    groups.append(("T1", "Fixed turbine", g0, x, [(c, "Rotor and\nturbine", "C1")])); x += GAP
    # T2
    g0 = x; c1 = place(SLOT); rotor(ax, c1); turbine(ax, c1, H_F + T_E + 0.15, blue, H_F + T_E)
    c2 = place(WIDE); contact_shadow(ax, c2, 3.6)
    savonius(ax, c2, 0.0, H_F, R_F, green); turbine(ax, c2, H_F + 0.15, blue, H_F)
    groups.append(("T2", "Fixed turbine and Savonius", g0, x,
                   [(c1, "Rotor and\nturbine", "C1"), (c2, "Savonius and\nturbine", "C3")])); x += GAP
    # T3 (y T3x)
    g0 = x; c1 = place(SLOT); rotor(ax, c1); stowed(ax, c1, H_F - H_X - 0.5, H_F - 0.5, R_X, orange)
    c2 = place(SLOT); extension(ax, c2, orange)
    c3 = place(WIDE); contact_shadow(ax, c3, 3.6)
    savonius(ax, c3, 0.0, H_F, R_F, green)
    savonius(ax, c3, H_F, H_F + H_X, R_X, green)
    groups.append(("T3", "Extension and Savonius", g0, x,
                   [(c1, "Extension\nstowed", "C0"), (c2, "Extension\nraised", "C2"),
                    (c3, "Both cylinders\nas Savonius", "C3")]))
    t3x = (c1, c2); x += GAP
    # T4
    g0 = x; c1 = place(SLOT); rotor(ax, c1); stowed(ax, c1, H_F - H_X - 0.5, H_F - 0.5, R_X, blue, blades=True)
    c2 = place(SLOT); rotor(ax, c2)
    turbine(ax, c2, H_F + T_E + 0.15, blue, H_F + T_E)
    groups.append(("T4", "Retractable turbine", g0, x,
                   [(c1, "Turbine\nstowed", "C0"), (c2, "Turbine\nraised", "C1")])); x += GAP
    # T5
    g0 = x; c1 = place(SLOT); rotor(ax, c1)
    c2 = place(WIDE); contact_shadow(ax, c2, 3.6); savonius(ax, c2, 0.0, H_F, R_F, green)
    groups.append(("T5", "Savonius conversion", g0, x, [(c1, "Rotor", "C0"), (c2, "Rotor opened\nas Savonius", "C3")]))
    x_end = x

    for name, title, a, b, modes in groups:
        deck(ax, a - 0.8, b + 0.8)
        ax.text((a + b) / 2, 40.2, name, ha="center", va="bottom", fontsize=9.0, fontweight="bold", color=INK)
        ax.text((a + b) / 2, 39.4, title, ha="center", va="top", fontsize=7.0, color=MUTED, linespacing=1.05)
        for c, label, conf in modes:
            if conf == "C3":
                savonius_plan(ax, c, -4.9, R_F, green)
            ax.text(c, -8.2, label, ha="center", va="top", fontsize=6.7, color=INK, linespacing=1.1)
            tag(ax, c, -14.3, conf)
    a, b = t3x
    ax.plot([a - 3.5, a - 3.5, b + 3.5, b + 3.5], [-16.3, -16.8, -16.8, -16.3], color=MUTED, lw=0.6)
    ax.text((a + b) / 2, -17.2, "T3x: extension only", ha="center", va="top", fontsize=6.6, color=MUTED)

    # escala y cotas
    xs = 2.0
    ax.plot([xs, xs], [0, 30], color=INK, lw=0.6)
    for z in (0, 10, 20, 30):
        ax.plot([xs - 0.5, xs + 0.5], [z, z], color=INK, lw=0.6)
        ax.text(xs - 0.9, z, f"{z}", ha="right", va="center", fontsize=6.2, color=MUTED)
    ax.text(xs - 0.9, 32.4, "m", ha="right", va="center", fontsize=6.2, color=MUTED)
    c = groups[0][4][0][0]
    ax.annotate("", xy=(c - R_F, 7.5), xytext=(c + R_F, 7.5),
                arrowprops=dict(arrowstyle="<->", lw=0.5, color=INK, mutation_scale=4), zorder=9)
    ax.text(c, 8.0, "4 m", ha="center", va="bottom", fontsize=6.0, color=INK, zorder=9)
    ax.text(c, H_F + T_E + 1.7, "8 m endplate", ha="center", va="bottom", fontsize=6.0, color=MUTED)

    # leyenda
    names = [(conf, nm) for conf, (col, nm) in CONF.items()]
    total = x_end - 9.0
    for i, (conf, nm) in enumerate(names):
        xx = 9.0 + i * total / 4.0
        tag(ax, xx, -22.3, conf)
        ax.text(xx + 2.4, -22.3, nm, ha="left", va="center", fontsize=7.0, color=INK)

    ax.set_xlim(-2.5, x_end + 1.5)
    ax.set_ylim(-24.9, 44.6)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"fig_topologias.{ext}"), dpi=300)


if __name__ == "__main__":
    main()
