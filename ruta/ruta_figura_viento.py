#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figura de las rosas de viento de A1 (versión 5, 27-10-2026): viento verdadero a 40 m respecto al rumbo del buque
(0 grados = de proa, 90 = de través, 180 = de popa; babor y estribor plegados) en las horas en que el buque navega
según el horario de 2026 (histogramas condicionados al horario, era5_horario_condicionado.py), ERA5 2015-2024: travesías
hacia Bilbao y hacia Portsmouth, de abril a octubre (caso base) y de diciembre a marzo. Hasta la versión 4, por tramo
(Canal, golfo) y estación (abril-noviembre, diciembre-marzo), con todas las horas y el rumbo hacia Bilbao.
Imprime también las estadísticas que cita el texto: velocidad media a 40 m, fracción de horas con viento verdadero y
aparente por delante del través en cada travesía, y fracción de horas con el rotor en marcha.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ruta_modelo as rm
from matplotlib import font_manager

# Tipo de letra de las figuras del artículo (Inter); si no está instalado, el de Matplotlib por defecto.
if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))
CLASSES = [(0, 4), (4, 8), (8, 12), (12, 16), (16, 30)]
LABELS = ["0–4", "4–8", "8–12", "12–16", "≥ 16"]
COLORS = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]   # rampa azul ordinal (250-650)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
DIR_NAMES = {"S": "To Bilbao", "N": "To Portsmouth"}
SEASON_NAMES = {0: "April–October", 1: "December–March"}


def rose(groups, direction, s):
    """Histograma medio (velocidad x ángulo respecto al rumbo propio) de las travesías de un sentido."""
    F = 0.0
    for name, dr, _ in rm.CROSSINGS:
        if dr == direction:
            Fi = groups["trav"][name][s]["F"]
            F = F + (Fi[:, ::-1] if direction == "N" else Fi)      # hacia Portsmouth: gamma -> 180 - gamma
    return F / F.sum()


def main():
    d, groups, nu, ng = rm.load_wind()
    u_mid = np.arange(nu) + 0.5
    theta = np.radians((np.arange(ng) + 0.5) * 15.0)
    width = np.radians(15.0)
    rmax = 0.0
    data = {}
    for g in ("S", "N"):
        for s in (0, 1):
            F = rose(groups, g, s)                                       # filas: velocidad; columnas: ángulo
            stack = np.array([F[(u_mid >= a) & (u_mid < b)].sum(axis=0) for a, b in CLASSES]) * 100
            data[g, s] = (F, stack)
            rmax = max(rmax, stack.sum(axis=0).max())
    rmax = np.ceil(rmax / 2.0) * 2.0
    rmax = max(rmax, 15.0)
    fig, axes = plt.subplots(1, 4, figsize=(7.48, 3.55), subplot_kw=dict(projection="polar"))
    panels = [("S", 0), ("N", 0), ("S", 1), ("N", 1)]
    for n, (g, s) in enumerate(panels):
        ax = axes[n]
        F, stack = data[g, s]
        bottom = np.zeros(ng)
        for c in range(len(CLASSES)):
            ax.bar(theta, stack[c], width=width, bottom=bottom, color=COLORS[c], edgecolor="white",
                   linewidth=0.5, label=f"{LABELS[c]} m/s" if n == 0 else None)
            bottom += stack[c]
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_thetamin(0)
        ax.set_thetamax(180)
        ax.set_ylim(0, rmax)
        ax.set_yticks([5, 10, 15])
        ax.set_yticklabels(["5%", "10%", "15%"] if n == 0 else [])
        ax.set_rlabel_position(90)
        ax.tick_params(axis="y", labelsize=5.5, labelcolor=MUTED, pad=0)
        ax.set_xticks(np.radians([0, 45, 90, 135, 180]))
        ax.set_xticklabels(["ahead", "", "beam", "", "astern"], fontsize=6.5, color=MUTED)
        ax.tick_params(axis="x", pad=0)
        ax.grid(color=GRID, linewidth=0.5)
        ax.spines["polar"].set_color("#c3c2b7")
        ax.spines["polar"].set_linewidth(0.6)
        u_mean = float((F.sum(axis=1) * u_mid).sum())
        fwd = float(F[:, :6].sum()) * 100
        ax.set_title(f"({'abcd'[n]}) {DIR_NAMES[g]}\n{SEASON_NAMES[s]}", fontsize=7.5, loc="left", color=INK,
                     pad=10, x=0.12)
        ax.text(0.5, -0.06, f"mean {u_mean:.1f} m/s\n{fwd:.0f}% forward of beam", transform=ax.transAxes,
                ha="left", va="top", fontsize=6.5, color=MUTED)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False, fontsize=7, handlelength=1.2,
               columnspacing=1.2, bbox_to_anchor=(0.5, -0.01), title="True wind speed at 40 m, hours at sea:", title_fontsize=7,
               alignment="left")
    fig.subplots_adjust(left=-0.04, right=0.99, top=0.84, bottom=0.25, wspace=-0.35)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_viento.{ext}"), dpi=300)

    # Estadísticas para el texto
    print("Velocidad media a 40 m y viento verdadero por delante del través (rumbo propio):")
    for (g, s), (F, _) in data.items():
        print(f"  {DIR_NAMES[g]:14s} {SEASON_NAMES[s]:15s} {float((F.sum(axis=1) * u_mid).sum()):.2f} m/s, "
              f"por delante del través {100 * float(F[:, :6].sum()):.0f} %")
    p = rm.P()
    print("Viento aparente por delante del través (cos beta > 0) y rotor en marcha, por travesía, abril-octubre:")
    for name, direction, hours in rm.crossings(p):
        v_kn, _ = rm.crossing_speed(hours, p)
        v = v_kn * rm.KN
        fwd, run, tot = 0.0, 0.0, 0.0
        for wg, F, ratio in rm.hist_entries(groups, 0, name):
            for iu in range(nu):
                for ig in range(ng):
                    f = F[iu, ig] * wg
                    if f == 0:
                        continue
                    gam = (ig + 0.5) * 15.0
                    if direction == "N":
                        gam = 180.0 - gam
                    t = rm.bin_terms(iu + 0.5, ratio[iu], gam, v, p)
                    tot += f
                    fwd += f * (t["cosb"] > 0)
                    run += f * t["running"]
        print(f"  {name:28s} {v_kn:5.2f} kn: aparente por delante del través {100 * fwd / tot:.0f} %, "
              f"rotor en marcha {100 * run / tot:.0f} %")


if __name__ == "__main__":
    main()
