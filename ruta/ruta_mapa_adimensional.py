#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapa de dónde paga cada función añadida (versión 4; revisión simulada R3-6; versión 5, 27-10-2026: fracciones
ponderadas con el clima de viento de las travesías del Galicia, revisión simulada v2, R3 menor 3): en función de la razón entre el viento
verdadero a 40 m y la velocidad del buque, U_w/V_s, y del ángulo del viento verdadero desde la proa, gamma, se
calcula con ruta_modelo.bin_terms (la misma física que la ruta, escenario A) el beneficio por unidad del rotor
convencional y el incremento de cada función: prolongación (T3x), turbina retráctil (T4) y conversión en Savonius
(T5, en el mar). Para que el mapa no dependa del casco del Galicia, se calcula sin resistencia inducida ni viento
sobre el buque; el caso con la resistencia inducida del Galicia se guarda también. Dos velocidades: 18 kn (ferri
rápido, como el Galicia) y 12 kn (buque lento).
Las fracciones «del mapa» dependen de los ejes elegidos; se dan también ponderadas con las horas de navegación de abril
a octubre (histogramas condicionados al horario, cada travesía a su velocidad).
Salida: mapa_adimensional.json y fig_mapa_adimensional.{pdf,png}.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import font_manager
import ruta_modelo as rm

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
    return dict(R1=r1, T3x=ext, T4=turb, T5=sav, marcha=run)


def winner(g):
    stack = np.stack([np.zeros_like(g["T3x"]), g["T3x"], g["T4"], g["T5"]])
    w = np.argmax(stack, axis=0)
    w[stack.max(axis=0) < 0.05] = 0
    return w


NAMES = ("ninguna", "prolongación", "turbina", "Savonius")


def climate_shares(p):
    """Fracción de las horas de navegación (abril-octubre, las cuatro travesías) en que cada función añadida es la
    que más aporta, y fracción con la turbina > 1 kW por unidad."""
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
    out["turbina > 1 kW"] = turb / tot
    return out


def main():
    out = {"ratios": RATIOS.tolist(), "gammas": GAMMAS.tolist(), "casos": {}}
    cases = {f"{v:g} kn, sin resistencia inducida": (v, rm.P(sir=False, windage=False, r0_wind=False)) for v in SPEEDS}
    cases.update({f"{v:g} kn, con la resistencia inducida del Galicia": (v, rm.P()) for v in SPEEDS})
    res = {}
    for lab, (v, p) in cases.items():
        g = grid(v, p)
        res[lab] = g
        w = winner(g)
        out["casos"][lab] = {k: np.round(a, 3).tolist() for k, a in g.items()}
        out["casos"][lab]["ganadora"] = w.tolist()
        share = {name: float((w == i).mean()) for i, name in enumerate(("ninguna", "prolongación", "turbina", "Savonius"))}
        out["casos"][lab]["fraccion_del_mapa"] = share
        print(f"{lab}: fracción del mapa " + ", ".join(f"{k} {100 * s:.0f} %" for k, s in share.items()) +
              f"; turbina > 0 en {100 * float((g['T4'] > 0.05).mean()):.0f} %")
    out["clima_galicia"] = {"sin resistencia inducida": climate_shares(rm.P(sir=False, windage=False, r0_wind=False)),
                            "con la resistencia inducida del Galicia": climate_shares(rm.P())}
    for lab, v in out["clima_galicia"].items():
        print(f"Horas de navegación del Galicia, abril-octubre, {lab}: " + ", ".join(f"{k} {100 * s:.0f} %"
                                                                                       for k, s in v.items()))
    json.dump(out, open(os.path.join(HERE, "mapa_adimensional.json"), "w", encoding="utf-8"), ensure_ascii=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.48, 3.3), constrained_layout=True, sharey=True)
    cmap = ListedColormap([GREY, BLUE, ORANGE, PURPLE])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    for ax, v, lab in zip(axes, SPEEDS, ("(a)", "(b)")):
        g = res[f"{v:g} kn, sin resistencia inducida"]
        w = winner(g)
        ax.pcolormesh(GAMMAS, RATIOS, w, cmap=cmap, norm=norm, shading="nearest")
        plt.rcParams["hatch.color"] = ORANGE
        plt.rcParams["hatch.linewidth"] = 0.5
        cs = ax.contour(GAMMAS, RATIOS, g["T3x"], levels=[5, 10, 20, 40], colors=INK, linewidths=0.6)
        ax.clabel(cs, fmt="%g kW", fontsize=6)
        ax.contour(GAMMAS, RATIOS, g["marcha"], levels=[0.5], colors=INK, linewidths=1.0, linestyles="--")
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
        fig.savefig(os.path.join(HERE, f"fig_mapa_adimensional.{ext}"), dpi=300)
    return out


if __name__ == "__main__":
    main()
