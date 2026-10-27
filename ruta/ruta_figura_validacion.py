#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figura del contraste a escala real (versión 4): (a) potencia propulsora equivalente al empuje de un rotor de 24 m x 4 m
en la ruta del Viking Grace a 21 kn, por escenario de coeficientes y con las sensibilidades del escenario A, frente a
los 207-282 kW verificados (Norsepower 2019); (b) desglose por rotor en el Galicia (abril-noviembre): empuje, giro,
rotor parado, resistencia inducida y neto, por escenario. Datos: validacion_viking_grace.json.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))
BLUE, ORANGE, INK, GREY, LIGHT = "#2a78d6", "#eb6834", "#0b0b0b", "#8a8984", "#e6e5e0"
GREEN, PURPLE = "#2f9e6e", "#7b5fc4"
SCEN = [("A+", "A+:"), ("A", "A:"), ("AB", "AB:"), ("B", "B:"), ("B′", "B':"), ("B″", "B'':")]


def pick(d, prefix):
    for k, v in d.items():
        if k.startswith(prefix):
            return v
    raise KeyError(prefix)


def main():
    d = json.load(open(os.path.join(HERE, "validacion_viking_grace.json"), encoding="utf-8"))
    lo, hi = d["ref_kw"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.48, 3.3), constrained_layout=True,
                                   gridspec_kw=dict(width_ratios=[1.15, 1.0]))
    for ax, lab in ((ax1, "(a)"), (ax2, "(b)")):
        ax.grid(True, axis="y", color=LIGHT, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_title(lab, loc="left", fontsize=9)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    # (a) Viking Grace
    names = [n for n, _ in SCEN]
    vals = [pick(d["escenarios"], pre)["empuje"] for _, pre in SCEN]
    x = np.arange(len(names))
    ax1.axhspan(lo, hi, color=GREEN, alpha=0.15, lw=0)
    ax1.text(len(names) - 0.55, hi + 4, "Verified, 207–282 kW", ha="right", va="bottom", fontsize=7, color=INK)
    cols = [BLUE, BLUE, GREY, ORANGE, ORANGE, ORANGE]
    ax1.bar(x, vals, width=0.62, color=cols, edgecolor="white", lw=2)
    for xi, v in zip(x, vals):
        ax1.text(xi, v + 4, f"{v:.0f}", ha="center", va="bottom", fontsize=7, color=INK)
    sens = d["sensibilidad"]
    marks = [("periodo del ensayo", "Measurement year", "o", 1.37), ("sin celdas con tierra", "Sea cells only", "s", 1.50),
             ("A+, sin celdas con tierra, viento", "A+, sea cells, wind ×1.1", "D", 0.42)]
    for key, lab, mk, xm in marks:
        v = pick(sens, key)["empuje"]
        ax1.plot(xm, v, mk, ms=5, mfc="white", mec=INK, mew=1.0, label=lab, ls="none")
    ax1.set_xticks(x)
    ax1.set_xticklabels(names)
    ax1.set_xlabel("Rotor coefficient scenario")
    ax1.set_ylabel("Thrust power / $\\eta_D$ per rotor (kW)")
    ax1.set_ylim(0, 300)
    ax1.legend(frameon=False, fontsize=6.5, loc="upper right", bbox_to_anchor=(1.0, 0.86))
    # (b) Galicia, desglose
    gal = d["galicia_abr_nov"]
    sc = [("A", "A:"), ("AB", "AB:"), ("B", "B:")]
    comp = [("_R1 empuje", "Thrust", BLUE, 1), ("_R1 giro", "Spin power", PURPLE, -1),
            ("_R1 parado", "Stopped-rotor drag", GREY, -1), ("_pen R1", "Induced resistance", ORANGE, -1)]
    xb = np.arange(len(sc))
    wdt = 0.16
    for j, (key, lab, col, sg) in enumerate(comp):
        v = [sg * pick(gal, pre)[key] for _, pre in sc]
        ax2.bar(xb + (j - 2) * wdt, v, width=wdt, color=col, label=lab, edgecolor="white", lw=0.8)
    net = [pick(gal, pre)["R1"] for _, pre in sc]
    ax2.bar(xb + 2 * wdt, net, width=wdt, color=INK, label="Net saving", edgecolor="white", lw=0.8)
    for xi, v in zip(xb, net):
        ax2.text(xi + 2 * wdt, v + (4 if v >= 0 else -12), f"{v:.0f}", ha="center", va="bottom", fontsize=7)
    ax2.axhline(0, color=INK, lw=0.6)
    ax2.set_xticks(xb)
    ax2.set_xticklabels([n for n, _ in sc])
    ax2.set_xlabel("Rotor coefficient scenario")
    ax2.set_ylabel("Galicia, April–November, kW per rotor")
    ax2.legend(frameon=False, fontsize=6.5, loc="upper right")
    ax2.set_ylim(-45, 190)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_validacion.{ext}"), dpi=300)


if __name__ == "__main__":
    main()
