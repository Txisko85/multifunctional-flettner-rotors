#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figura del modelo horario (A1; versión 5, 27-10-2026): (a) ahorro de R1 y aumento de T3x por año, 2015-2024, con el
horario de 2026 de abril a octubre (caso base); (b) aumento sobre R1 de la prolongación (T3x), del rotor abatible
(R1t, comparador) y de la turbina retráctil (T4) frente al tiempo de cada cambio de modo, con control óptimo
(programación dinámica, conociendo el viento de la travesía), causal (persistencia L = 6 h), hora a hora y, a
tau = 0,25 h, control predictivo con previsión buena y mala (horizonte de 12 h).
Lee resultados_horarios.json (ruta_horaria.py). Escribe fig_horaria.png y fig_horaria.pdf.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Tipo de letra de las figuras del artículo (Inter); si no está instalado, el de Matplotlib por defecto.
if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
BLUE_D, BLUE, BLUE_L = "#104281", "#2a78d6", "#86b6ef"
ORANGE = "#c96a1b"
GREEN = "#2f9e6e"
L_SHOW = 6
H_SHOW = 12


def main():
    r = json.load(open(os.path.join(HERE, "resultados_horarios.json"), encoding="utf-8"))
    years = sorted(int(y) for y in r["R1"]["por_anio"])
    r1 = np.array([r["R1"]["por_anio"][str(y)] for y in years])
    t3x = np.array([r["incrementos_mar"]["T3x"]["por_anio inst"][str(y)] for y in years])
    taus = [0.0] + list(r["tau_h"])

    plt.rcParams.update({"font.size": 7, "axes.edgecolor": "#8a897f", "axes.linewidth": 0.6,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.48, 2.75), gridspec_kw=dict(width_ratios=[1.25, 1]))

    # (a) por año
    x = np.arange(len(years))
    ax1.bar(x, r1, width=0.68, color=BLUE, label="R1, four rotors, relative to R0")
    ax1.bar(x, t3x, width=0.68, bottom=r1, color=BLUE_L, label="T3x, extension, relative to R1")
    ax1.axhline(r1.mean(), color=BLUE_D, lw=0.8, ls="--", label="Mean of R1", zorder=0.9)
    ax1.axhline((r1 + t3x).mean(), color=MUTED, lw=0.8, ls=":", label="Mean of R1 + T3x", zorder=0.9)
    for xi, a, b in zip(x, r1, t3x):
        ax1.text(xi, a + 0.5 * b, f"{100 * b / a:.0f}%", ha="center", va="center", fontsize=5.6, color=BLUE_D)
    ax1.set_xticks(x)
    ax1.set_xticklabels([str(y) for y in years], fontsize=6.5)
    ax1.set_ylabel("Fuel saved at sea (t/yr)")
    ax1.set_ylim(0, 560)
    ax1.set_yticks(range(0, 451, 100))
    ax1.grid(axis="y", color=GRID, lw=0.5)
    ax1.set_axisbelow(True)
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.legend(loc="upper left", frameon=False, fontsize=6.3, handlelength=1.6, ncol=2, columnspacing=1.2)
    ax1.set_title("(a) Year by year, 2026 timetable", loc="left", fontsize=7.5, color=INK)

    # (b) coste de conmutación: etiquetas directas al final de cada curva
    series = (("T3x", r["incrementos_mar"]["T3x"], list(r["tau_h"]), BLUE_D),
              ("R1t", r["R1t"], list(r["tau_abatible"]), GREEN),
              ("T4", r["incrementos_mar"]["T4"], list(r["tau_h"]), ORANGE))
    ends = {}
    for topo, d, tl, col in series:
        inst = d["inst"]
        xs = [0.0] + tl
        opt = [inst] + [d[f"optima {t}"] for t in tl]
        cau = [inst] + [d[f"causal L{L_SHOW} {t}"] for t in tl]
        gre = [inst] + [d[f"voraz {t}"] for t in tl]
        ax2.plot(xs, opt, color=col, lw=1.2, marker="o", ms=3.2)
        ax2.plot(xs, cau, color=col, lw=1.0, ls="--", marker="s", ms=2.8, mfc="white")
        ax2.plot(xs, gre, color=col, lw=0.9, ls=":", marker="^", ms=3.0, mfc="white")
        for lev, mk in (("buena", "D"), ("mala", "d")):
            key = f"mpc {lev} H{H_SHOW} 0.25"
            if key in d:
                ax2.plot([0.25 + 0.012], [d[key]], mk, ms=3.6, color=col, mfc=col if lev == "buena" else "white",
                         mew=0.9, ls="none")
        ends[topo] = (xs[-1], opt[-1], cau[-1], gre[-1])
    ax2.axhline(0, color="#8a897f", lw=0.6)
    ax2.set_xlim(-0.02, 0.69)
    ax2.set_xticks([0.0, 0.1, 0.25, 0.5])
    ax2.set_xticklabels(["0", "0.1", "0.25", "0.5"])
    ax2.set_xlabel("Switching time per change of mode, τ (h)")
    ax2.set_ylabel("Fuel saved relative to R1 (t/yr)")
    ax2.set_ylim(-25, 90)
    ax2.grid(color=GRID, lw=0.5)
    ax2.set_axisbelow(True)
    ax2.spines[["top", "right"]].set_visible(False)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=BLUE_D, lw=1.4), Line2D([], [], color=GREEN, lw=1.4),
         Line2D([], [], color=ORANGE, lw=1.4),
         Line2D([], [], color=MUTED, lw=1.2, marker="o", ms=3.2),
         Line2D([], [], color=MUTED, lw=1.0, ls="--", marker="s", ms=2.8, mfc="white"),
         Line2D([], [], color=MUTED, lw=0.9, ls=":", marker="^", ms=3.0, mfc="white"),
         Line2D([], [], color=MUTED, marker="D", ms=3.6, ls="none"),
         Line2D([], [], color=MUTED, marker="d", ms=3.6, mfc="white", ls="none")]
    ax2.legend(h, ["T3x, extension", "R1t, tiltable rotor", "T4, retractable turbine",
                   "Optimal, perfect foresight", f"Causal, persistence {L_SHOW} h", "Hour by hour",
                   f"Predictive, good forecast ({H_SHOW} h)", f"Predictive, poor forecast ({H_SHOW} h)"],
               loc="lower left", bbox_to_anchor=(0.0, 0.22), frameon=False, fontsize=5.8, handlelength=1.8,
               ncol=2, columnspacing=0.9, labelspacing=0.35)
    ax2.set_title("(b) Cost of switching", loc="left", fontsize=7.5, color=INK)

    fig.subplots_adjust(left=0.07, right=0.99, top=0.91, bottom=0.17, wspace=0.25)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_horaria.{ext}"), dpi=300)
    print("T3x / R1 por año (%):", ", ".join(f"{y} {100 * b / a:.1f}" for y, a, b in zip(years, r1, t3x)))


if __name__ == "__main__":
    main()
