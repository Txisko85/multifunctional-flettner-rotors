#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map of where each added function pays (version 4; simulated review R3-6; version 5, 27-10-2026: shares weighted with
the wind climate of the crossings of the Galicia, simulated review v2, R3 minor 3): as a function of the ratio between
the true wind at 40 m and the ship speed, U_w/V_s, and of the true wind angle from the bow, gamma, route_model.bin_terms
(the same physics as the route, scenario A) gives the benefit per unit of the conventional rotor and the increment of
each function: extension (T3x), retractable turbine (T4) and Savonius conversion (T5, at sea). So that the map does
not depend on the hull of the Galicia, it is computed without induced resistance and without wind on the ship; the
case with the induced resistance of the Galicia is also saved. Two speeds: 18 kn (fast ferry, like the Galicia) and
12 kn (slow ship).
The shares "of the map" depend on the chosen axes; they are also given weighted with the sailing hours from April to
October (histograms conditioned on the timetable, each crossing at its speed).
Output: function_map.json and function_map.{pdf,png}.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import font_manager
import route_model as rm

if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))
RATIOS = np.round(np.arange(0.1, 2.501, 0.05), 3)
GAMMAS = np.arange(2.5, 180.0, 5.0)
RATIO_55_40 = 1.025
SPEEDS = (18.0, 12.0)
GREY, BLUE, ORANGE, PURPLE, INK = "#e6e5e0", "#2a78d6", "#eb6834", "#7b5fc4", "#0b0b0b"


def grid(v_kn, p):
    v = v_kn * rm.KN
    n_r, n_g = len(RATIOS), len(GAMMAS)
    r1, ext, turb, sav, run = (np.zeros((n_r, n_g)) for _ in range(5))
    for i, rr in enumerate(RATIOS):
        for j, g in enumerate(GAMMAS):
            t = rm.bin_terms(float(rr * v), RATIO_55_40, float(g), v, p)
            tt = rm.topologies(t)
            r1[i, j] = tt["R1"]
            ext[i, j] = tt["T3x"]
            turb[i, j] = tt["T4"]
            sav[i, j] = tt["T5"]
            run[i, j] = t["running"]
    return dict(R1=r1, T3x=ext, T4=turb, T5=sav, running=run)


def winner(g):
    stack = np.stack([np.zeros_like(g["T3x"]), g["T3x"], g["T4"], g["T5"]])
    w = np.argmax(stack, axis=0)
    w[stack.max(axis=0) < 0.05] = 0
    return w


NAMES = ("none", "extension", "turbine", "Savonius")


def climate_shares(p):
    """Share of the sailing hours (April-October, the four crossings) in which each added function contributes most,
    and share with the turbine > 1 kW per unit."""
    d, groups, nu, ng = rm.load_wind()
    acc = {k: 0.0 for k in NAMES}
    turb, tot = 0.0, 0.0
    for name, direction, hours in rm.crossings(p):
        v_kn, t_open = rm.crossing_speed(hours, p)
        v = v_kn * rm.KN
        for wg, F, ratio in rm.hist_entries(groups, 0, name):
            for iu in range(nu):
                for ig in range(ng):
                    f = F[iu, ig] * wg * t_open
                    if f == 0:
                        continue
                    gam = (ig + 0.5) * 15.0
                    if direction == "N":
                        gam = 180.0 - gam
                    tt = rm.topologies(rm.bin_terms(iu + 0.5, ratio[iu], gam, v, p))
                    vals = [0.0, tt["T3x"], tt["T4"], tt["T5"]]
                    i = int(np.argmax(vals)) if max(vals) >= 0.05 else 0
                    acc[NAMES[i]] += f
                    turb += f * (tt["T4"] > 1.0)
                    tot += f
    out = {k: v / tot for k, v in acc.items()}
    out["turbine > 1 kW"] = turb / tot
    return out


def main():
    out = {"ratios": RATIOS.tolist(), "gammas": GAMMAS.tolist(), "cases": {}}
    cases = {f"{v:g} kn, without induced resistance": (v, rm.P(sir=False, windage=False, r0_wind=False)) for v in SPEEDS}
    cases.update({f"{v:g} kn, with the induced resistance of the Galicia": (v, rm.P()) for v in SPEEDS})
    res = {}
    for lab, (v, p) in cases.items():
        g = grid(v, p)
        res[lab] = g
        w = winner(g)
        out["cases"][lab] = {k: np.round(a, 3).tolist() for k, a in g.items()}
        out["cases"][lab]["winner"] = w.tolist()
        share = {name: float((w == i).mean()) for i, name in enumerate(NAMES)}
        out["cases"][lab]["map_share"] = share
        print(f"{lab}: share of the map " + ", ".join(f"{k} {100 * s:.0f} %" for k, s in share.items()) +
              f"; turbine > 0 in {100 * float((g['T4'] > 0.05).mean()):.0f} %")
    out["galicia_climate"] = {"without induced resistance": climate_shares(rm.P(sir=False, windage=False, r0_wind=False)),
                              "with the induced resistance of the Galicia": climate_shares(rm.P())}
    for lab, v in out["galicia_climate"].items():
        print(f"Sailing hours of the Galicia, April-October, {lab}: " + ", ".join(f"{k} {100 * s:.0f} %"
                                                                                   for k, s in v.items()))
    json.dump(out, open(os.path.join(HERE, "function_map.json"), "w", encoding="utf-8"), ensure_ascii=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.48, 3.3), constrained_layout=True, sharey=True)
    cmap = ListedColormap([GREY, BLUE, ORANGE, PURPLE])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    for ax, v, lab in zip(axes, SPEEDS, ("(a)", "(b)")):
        g = res[f"{v:g} kn, without induced resistance"]
        w = winner(g)
        ax.pcolormesh(GAMMAS, RATIOS, w, cmap=cmap, norm=norm, shading="nearest")
        plt.rcParams["hatch.color"] = ORANGE
        plt.rcParams["hatch.linewidth"] = 0.5
        cs = ax.contour(GAMMAS, RATIOS, g["T3x"], levels=[5, 10, 20, 40], colors=INK, linewidths=0.6)
        ax.clabel(cs, fmt="%g kW", fontsize=6)
        ax.contour(GAMMAS, RATIOS, g["running"], levels=[0.5], colors=INK, linewidths=1.0, linestyles="--")
        ax.contourf(GAMMAS, RATIOS, g["T4"], levels=[1.0, 1e9], colors="none", hatches=["////"])
        ax.contour(GAMMAS, RATIOS, g["T4"], levels=[1.0], colors=ORANGE, linewidths=0.9)
        ax.set_title(f"{lab} $V_s$ = {v:g} kn", loc="left", fontsize=9)
        ax.set_xlabel("True wind angle from the bow, $\\gamma$ (°)")
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 45, 90, 135, 180])
    axes[0].set_ylabel("True wind at 40 m / ship speed, $U_w/V_s$")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (GREY, BLUE)]
    labels = ["No added function pays", "Extension adds most"]
    handles.append(plt.Rectangle((0, 0), 1, 1, facecolor="none", edgecolor=ORANGE, hatch="////", lw=0.9))
    labels.append("Turbine adds > 1 kW (less than the extension)")
    handles.append(plt.Line2D([0], [0], color=INK, lw=0.6))
    labels.append("Extension gain per unit (kW)")
    handles.append(plt.Line2D([0], [0], color=INK, lw=1.0, ls="--"))
    labels.append("Rotor stopped below or left of this line")
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False, fontsize=6.5)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"function_map.{ext}"), dpi=300)
    return out


if __name__ == "__main__":
    main()
