#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. 3 de A1: rosas del viento verdadero respecto al rumbo hacia Bilbao (ERA5 2015-2024, 40 m), por tramo
(Canal de la Mancha, golfo de Bizkaia) y estación (abril-noviembre, diciembre-marzo). Las rosas son medias
porque babor y estribor se pliegan: 0 grados = viento de proa, 90 = de través, 180 = de popa. En las travesías
hacia Portsmouth los ángulos se reflejan respecto al través (gamma_N = 180 - gamma_S).
Imprime también las estadísticas que cita el texto: velocidad media a 40 m, fracción de horas con viento
verdadero y aparente por delante del través en cada travesía, y fracción de horas con el rotor en marcha.
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
GROUP_NAMES = {"C": "English Channel", "B": "Bay of Biscay"}
SEASON_NAMES = {0: "April–November", 1: "December–March"}


def main():
    d, groups, nu, ng = rm.load_wind()
    u_mid = np.arange(nu) + 0.5
    theta = np.radians((np.arange(ng) + 0.5) * 15.0)
    width = np.radians(15.0)
    rmax = 0.0
    data = {}
    for g in ("C", "B"):
        for s in (0, 1):
            F = groups[g]["F"][s] / groups[g]["F"][s].sum()            # filas: velocidad; columnas: ángulo
            stack = np.array([F[(u_mid >= a) & (u_mid < b)].sum(axis=0) for a, b in CLASSES]) * 100
            data[g, s] = (F, stack)
            rmax = max(rmax, stack.sum(axis=0).max())
    rmax = np.ceil(rmax / 2.0) * 2.0
    rmax = max(rmax, 15.0)
    fig, axes = plt.subplots(1, 4, figsize=(7.48, 3.55), subplot_kw=dict(projection="polar"))
    panels = [("C", 0), ("C", 1), ("B", 0), ("B", 1)]
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
        ax.set_title(f"({'abcd'[n]}) {GROUP_NAMES[g]}\n{SEASON_NAMES[s]}", fontsize=7.5, loc="left", color=INK,
                     pad=10, x=0.12)
        ax.text(0.5, -0.06, f"mean {u_mean:.1f} m/s\n{fwd:.0f}% forward of beam", transform=ax.transAxes,
                ha="left", va="top", fontsize=6.5, color=MUTED)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False, fontsize=7, handlelength=1.2,
               columnspacing=1.2, bbox_to_anchor=(0.5, -0.01), title="True wind speed at 40 m:", title_fontsize=7,
               alignment="left")
    fig.subplots_adjust(left=-0.04, right=0.99, top=0.84, bottom=0.25, wspace=-0.35)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_viento.{ext}"), dpi=300)

    # Estadísticas para el texto
    print("Velocidad media a 40 m y viento verdadero por delante del través (rumbo hacia Bilbao):")
    for (g, s), (F, _) in data.items():
        print(f"  {GROUP_NAMES[g]:16s} {SEASON_NAMES[s]:15s} {float((F.sum(axis=1) * u_mid).sum()):.2f} m/s, "
              f"por delante del través {100 * float(F[:, :6].sum()):.0f} %")
    p = rm.P()
    print("Viento aparente por delante del través (cos beta > 0) y rotor en marcha, por travesía, año completo:")
    w_tot = sum(rm.NM_OPEN.values())
    for name, direction, hours in rm.crossings(p):
        v_kn, _ = rm.crossing_speed(hours, p)
        v = v_kn * rm.KN
        fwd, run, tot = 0.0, 0.0, 0.0
        for s in (0, 1):
            ws = rm.DAYS[s] / sum(rm.DAYS.values())
            for g, gd in groups.items():
                F, ratio = gd["F"][s], gd["ratio"][s]
                wg = rm.NM_OPEN[g] / w_tot
                for iu in range(nu):
                    for ig in range(ng):
                        f = F[iu, ig] * wg * ws
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
