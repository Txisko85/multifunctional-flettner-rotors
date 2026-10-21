#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figura de la base de coeficientes del rotor (A1, apartado 3.2): CL(k) y CD(k) de los conjuntos A (De Marco et al.
2016, AR 6, de/d 2) y B (Bordogna et al. 2019, Re >= 3,6e5, mas resistencia inducida de AR 6 sobre cubierta),
la mezcla central (w = 0,5), los ensayos de Ackeret recopilados por Seifert (2012) y los valores de nivel 0.
Datos digitalizados en datos_coef (vectores del PDF o segmentacion por color de la imagen).
19-10-2026: las leyendas citan autor y año en lugar del número de referencia, que cambia con la numeración
del manuscrito en LaTeX.
"""
import json
import math
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import coeficientes as cf
import ruta_modelo as rm
from matplotlib import font_manager

# Tipo de letra de las figuras del artículo (Inter); si no está instalado, el de Matplotlib por defecto.
if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))
DATOS = os.path.join(HERE, "datos_coef")
BLUE, ORANGE, INK, GREY = "#2a78d6", "#eb6834", "#0b0b0b", "#8a8984"


def main():
    p = rm.P()
    k = np.linspace(0.0, 3.0, 301)
    pa, pb = rm.P(w=1.0), rm.P(w=0.0)
    cl_a, cd_a = rm.rotor_coeffs(k, pa)
    cl_b, cd_b = rm.rotor_coeffs(k, pb)
    cl_c, cd_c = rm.rotor_coeffs(k, p)
    bo = json.load(open(os.path.join(DATOS, "bordogna2019_fig4_fig7.json")))
    ak = json.load(open(os.path.join(DATOS, "ackeret_seifert_fig27.json")))["series"]["placas de 140 mm"]
    sets = ("Re3.6e5", "Re5.5e5", "Re1.0e6")
    bcl = np.array([q for s in sets for q in bo["CL"][s]]); bcd = np.array([q for s in sets for q in bo["CD"][s]])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.48, 3.1), constrained_layout=True)
    for ax, lab in ((ax1, "(a)"), (ax2, "(b)")):
        ax.grid(True, color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
        ax.set_xlim(0, 3.2); ax.set_xlabel("Spin ratio k"); ax.set_title(lab, loc="left", fontsize=9)
    ax1.plot(bcl[:, 0], bcl[:, 1], "o", ms=3.2, mfc="none", mec=ORANGE, mew=0.8,
             label="Bordogna et al. (2019), Re 3.6×10⁵–10⁶")
    ax1.plot(cf.K_A[1:], cf.CL_A[1:], "o", ms=5, color=BLUE, label="De Marco et al. (2016), AR 6, $d_e/d$ 2")
    ax1.plot(ak["k_CL"][1:], ak["CL"][1:], "D", ms=4, mfc="none", mec=GREY, mew=0.9,
             label="Ackeret, via Seifert (2012)")
    ax1.plot(k, cl_a, color=BLUE, lw=1.2)
    ax1.plot(k, cl_b, color=ORANGE, lw=1.2)
    ax1.plot(k, cl_c, color=INK, lw=2.0, label="Central, $w$ = 0.5")
    ax1.plot(k, 8.0 * k / 3.0, "--", color=GREY, lw=1.0, label="Level 0")
    ax1.set_ylabel("Lift coefficient $C_L$"); ax1.set_ylim(0, 10.5)
    ax2.plot(bcd[:, 0], bcd[:, 1], "o", ms=3.2, mfc="none", mec=ORANGE, mew=0.8)
    ax2.plot(cf.K_A[1:], cf.CD_A[1:], "o", ms=5, color=BLUE)
    ax2.plot(ak["k_CD"][1:], ak["CD"][1:], "D", ms=4, mfc="none", mec=GREY, mew=0.9)
    ax2.plot(k, cd_a, color=BLUE, lw=1.2, label="Set A")
    ax2.plot(k, cd_b, color=ORANGE, lw=1.2, label="Set B (with induced drag)")
    ax2.plot(k, cd_c, color=INK, lw=2.0)
    ax2.plot(k, np.full_like(k, 2.0), "--", color=GREY, lw=1.0)
    ax2.set_ylabel("Drag coefficient $C_D$"); ax2.set_ylim(0, 4.5)
    h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1, l1, frameon=False, fontsize=6.5, loc="upper left")
    ax2.legend(h2, l2, frameon=False, fontsize=6.5, loc="upper left")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_coeficientes.{ext}"), dpi=300)
    print("k  CL_A CD_A  CL_B CD_B  CL_c CD_c")
    for kk in (1.0, 1.5, 2.0, 2.5, 3.0):
        i = int(round(kk * 100))
        print(f"{kk:3.1f} {cl_a[i]:5.2f} {cd_a[i]:5.2f} {cl_b[i]:5.2f} {cd_b[i]:5.2f} {cl_c[i]:5.2f} {cd_c[i]:5.2f}")


if __name__ == "__main__":
    main()
