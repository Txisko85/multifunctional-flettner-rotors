#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure of the rotor coefficient base (A1, section 3.2). Version 4 (22-10-2026): epistemic scenarios.
Version 5 (27-10-2026): band of scenario A-VG (C_L of A times 1.14-1.36, fitted to the Viking Grace; C_D unchanged).
CL(k) and CD(k) of sets A (De Marco et al. 2016, AR 6, de/d 2; reference scenario) and B (Bordogna et al. 2019,
Re >= 3.6e5, plus the induced drag of AR 6 on deck), B' (B with the lifting-line correction) and the blend AB of
version 3; extended range k = 3-4.6 (scenario A+ and sensitivity), with the points of Bordogna at Re 3.6e5 and the
shape of the tests of Ackeret with 120 mm plates. Points: De Marco (vectors of the PDF), Bordogna and Ackeret (120 mm
plates) compiled by Seifert (2012), digitised (coef_data). Writes rotor_coefficients.pdf/png.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import coefficients as cf
import route_model as rm
from matplotlib import font_manager

if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "coef_data")
BLUE, ORANGE, INK, GREY = "#2a78d6", "#eb6834", "#0b0b0b", "#8a8984"


def main():
    k = np.linspace(0.0, 3.0, 301)
    kx = np.linspace(3.0, 4.6, 161)
    pa, pb = rm.P(**cf.COEF_SCENARIOS["A"]), rm.P(**cf.COEF_SCENARIOS["B"])
    pbp, pab = rm.P(**cf.COEF_SCENARIOS["B'"]), rm.P(**cf.COEF_SCENARIOS["AB"])
    pax, pbx = rm.P(w=1.0, k_max=4.6), rm.P(w=0.0, k_max=4.6)
    cl_a, cd_a = rm.rotor_coeffs(k, pa)
    cl_b, cd_b = rm.rotor_coeffs(k, pb)
    cl_bp, cd_bp = rm.rotor_coeffs(k, pbp)
    cl_ab, cd_ab = rm.rotor_coeffs(k, pab)
    cl_ax, cd_ax = rm.rotor_coeffs(kx, pax)
    cl_bx, cd_bx = rm.rotor_coeffs(kx, pbx)
    bo = json.load(open(os.path.join(DATA, "bordogna2019_fig4_fig7.json")))
    ak = json.load(open(os.path.join(DATA, "ackeret_seifert_fig27.json")))["series"]["120 mm plates"]
    sets = ("Re3.6e5", "Re5.5e5", "Re1.0e6")
    bcl = np.array([q for s in sets for q in bo["CL"][s]])
    bcd = np.array([q for s in sets for q in bo["CD"][s]])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.48, 3.2), constrained_layout=True)
    for ax, lab in ((ax1, "(a)"), (ax2, "(b)")):
        ax.grid(True, color="#e6e5e0", lw=0.6)
        ax.set_axisbelow(True)
        ax.set_xlim(0, 4.8)
        ax.set_xlabel("Spin ratio k")
        ax.set_title(lab, loc="left", fontsize=9)
        ax.axvspan(3.0, 4.8, color="#f3f2ee", lw=0)
    ax1.text(3.9, 0.4, "Extended\n(scenario A+)", ha="center", fontsize=6.5, color=GREY)
    ax1.plot(bcl[:, 0], bcl[:, 1], "o", ms=3.2, mfc="none", mec=ORANGE, mew=0.8,
             label="Bordogna et al. (2019), Re 3.6×10⁵–10⁶")
    ax1.plot(cf.K_A[1:], cf.CL_A[1:], "o", ms=5, color=BLUE, label="De Marco et al. (2016), AR 6, $D_e/D_F$ = 2")
    ax1.plot(ak["k_CL"][1:], ak["CL"][1:], "D", ms=4, mfc="none", mec=GREY, mew=0.9,
             label="Ackeret, 120 mm endplates,\nvia Seifert (2012)")
    lo_vg, hi_vg = cf.CL_MULT_VG["207"], cf.CL_MULT_VG["282"]
    ax1.fill_between(k, lo_vg * cl_a, hi_vg * cl_a, color=BLUE, alpha=0.13, lw=0,
                     label=f"A-VG: A with $C_L$ ×{lo_vg:.2f}–{hi_vg:.2f}\n(fitted to the Viking Grace)")
    ax1.plot(k, cf.CL_MULT_VG["244.5"] * cl_a, color=BLUE, lw=0.9, ls=(0, (1, 1.2)))
    ax1.plot(k, cl_a, color=BLUE, lw=2.0, label="A (reference)")
    ax1.plot(kx, cl_ax, color=BLUE, lw=1.4, ls=(0, (4, 2)))
    ax1.plot(k, cl_ab, color=GREY, lw=1.0, label="AB (equal blend of A and B)")
    ax1.plot(k, cl_b, color=ORANGE, lw=1.4, label="B")
    ax1.plot(kx, cl_bx, color=ORANGE, lw=1.1, ls=(0, (4, 2)))
    ax1.plot(k, cl_bp, color=ORANGE, lw=1.0, ls=(0, (1.5, 1.5)), label="B′ (lifting-line correction)")
    ax1.set_ylabel("Lift coefficient $C_L$")
    ax1.set_ylim(0, 14)
    ax2.plot(bcd[:, 0], bcd[:, 1], "o", ms=3.2, mfc="none", mec=ORANGE, mew=0.8)
    ax2.plot(cf.K_A[1:], cf.CD_A[1:], "o", ms=5, color=BLUE)
    ax2.plot(ak["k_CD"][1:], ak["CD"][1:], "D", ms=4, mfc="none", mec=GREY, mew=0.9)
    ax2.plot(k, cd_a, color=BLUE, lw=2.0)
    ax2.plot(kx, cd_ax, color=BLUE, lw=1.4, ls=(0, (4, 2)))
    ax2.plot(k, cd_ab, color=GREY, lw=1.0)
    ax2.plot(k, cd_b, color=ORANGE, lw=1.4)
    ax2.plot(kx, cd_bx, color=ORANGE, lw=1.1, ls=(0, (4, 2)))
    ax2.plot(k, cd_bp, color=ORANGE, lw=1.0, ls=(0, (1.5, 1.5)))
    ax2.set_ylabel("Drag coefficient $C_D$ (B with induced drag)")
    ax2.set_ylim(0, 6.5)
    h1, l1 = ax1.get_legend_handles_labels()
    ax2.legend(h1, l1, frameon=False, fontsize=6.3, loc="upper left", labelspacing=0.3)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"rotor_coefficients.{ext}"), dpi=300)
    print("k  CL_A CD_A  CL_B CD_B  CL_B' CD_B'")
    for kk in (1.0, 1.5, 2.0, 2.5, 3.0):
        i = int(round(kk * 100))
        print(f"{kk:3.1f} {cl_a[i]:5.2f} {cd_a[i]:5.2f} {cl_b[i]:5.2f} {cd_b[i]:5.2f} {cl_bp[i]:5.2f} {cd_bp[i]:5.2f}")


if __name__ == "__main__":
    main()
