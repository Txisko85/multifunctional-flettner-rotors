#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validación del modelo DMST de la turbina (A1, 21-10-2026).

1) Battisti et al. (2018), rotor H de tres palas NACA 0021 (c 85 mm, R 0,514 m, L 1,46 m, solidez 0,25) en el
   túnel del Politecnico di Milano a 400 rpm, U 6-16 m/s (Re de cuerda ~1-2e5): Cp, empuje Cx y fuerza lateral Cy
   (Tabla 8, datos/battisti2018_tabla8.json). Solo se modelan las palas: los brazos y el eje del ensayo no están
   publicados, y su resistencia está incluida en el empuje medido.
2) Kjellin et al. (2011), rotor H de 12 kW de la Universidad de Uppsala (tres palas NACA 0021, c 0,25 m, R 3 m,
   H 5 m; geometría en Apelfröjd et al. 2016) medido en campo: Cp máximo 0,29 a lambda 3,3. Se comparan los brazos
   con un coeficiente de resistencia de 0,02 y de 0,05 (cuerda 0,15 m; supuestos).
Variantes: Gormont-Berg con K1 = 0 al decrecer |alfa| (central, Shires 2013), Gormont-Berg original (K1 = 0,5) y
polares estáticas.
Salida: validacion_dmst.json y fig_validacion_turbina.pdf/png.
Uso: python validacion_dmst.py   (unos 7 minutos); python validacion_dmst.py --figura   (solo la figura, desde el JSON)
"""
import json
import math
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import dmst

if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))
BLUE, ORANGE, INK, GREY = "#2a78d6", "#eb6834", "#0b0b0b", "#8a8984"
VARIANTS = (("central", {"ds": True, "k1_dec": 0.0}, "Gormont–Berg, $K_1$ = 0 when |α| falls (central)"),
            ("gb05", {"ds": True, "k1_dec": 0.5}, "Gormont–Berg, $K_1$ = 0.5"),
            ("static", {"ds": False}, "Static polars"))


def air_nu(t_c, rho):
    mu = 1.716e-5 * ((t_c + 273.15) / 273.15) ** 1.5 * (273.15 + 110.4) / (t_c + 273.15 + 110.4)   # Sutherland
    return mu / rho


def battisti():
    d = json.load(open(os.path.join(HERE, "datos", "battisti2018_tabla8.json"), encoding="utf-8"))
    A, R = 1.5, 0.514
    tur = dmst.Turbine(n=3, R=R, H=1.46, c=0.085, tc=0.21, n_struts=2, c_strut=0.0, cd_strut=0.0, d_shaft=0.0)
    exp = []
    for r in d["filas"]:
        v, t_c, rho, rpm, q_t, uq, tx, ty = r[0], r[2], r[4], r[7], r[9], r[10], r[11], r[13]
        w = rpm * 2 * math.pi / 60
        q = 0.5 * rho * A * v * v
        exp.append(dict(U=v, lam=w * R / v, cp=q_t * w / (q * v), u_cp=uq * w / (q * v), cx=tx / q, cy=ty / q,
                        rho=rho, nu=air_nu(t_c, rho)))
    w = 400 * 2 * math.pi / 60
    model = {k: [] for k, _, _ in VARIANTS}
    us = np.linspace(5.8, 16.4, 23)
    for u in us:
        lam = w * R / u
        for k, o, _ in VARIANTS:
            s = dmst.solve(tur, float(u), lam, rho=1.156, nu=air_nu(25.0, 1.156), opts=o)
            model[k].append(dict(U=float(u), lam=lam, cp=s["cp"], cx=s["ct"], cy=abs(s["cs"])))
    at_exp = {}
    for k, o, _ in VARIANTS:
        at_exp[k] = [dmst.solve(tur, e["U"], e["lam"], rho=e["rho"], nu=e["nu"], opts=o) for e in exp]
    return exp, model, at_exp


def marsta():
    curves = {}
    lams = np.round(np.arange(2.0, 5.01, 0.25), 2)
    for cd_s in (0.02, 0.05):
        tur = dmst.Turbine(n=3, R=3.0, H=5.0, c=0.25, tc=0.21, n_struts=2, c_strut=0.15, cd_strut=cd_s, d_shaft=0.0)
        for k, o, _ in VARIANTS:
            curves[f"{k}_cds{cd_s}"] = [dict(lam=float(l), cp=dmst.solve(tur, 8.0, float(l), opts=o)["cp"]) for l in lams]
    return curves


def figure(exp, model, mar):
    fig, axs = plt.subplots(2, 2, figsize=(7.48, 5.4), constrained_layout=True)
    (a1, a2), (a3, a4) = axs
    styles = {"central": dict(color=BLUE, lw=2.0), "gb05": dict(color=ORANGE, lw=1.2),
              "static": dict(color=GREY, lw=1.0, ls="--")}
    for ax, key, ylab, lab in ((a1, "cp", "Power coefficient $C_p$", "(a) H-rotor in wind tunnel"),
                               (a2, "cx", "Streamwise-force coefficient $C_t$", "(b)"),
                               (a3, "cy", "Side-force coefficient $|C_s|$", "(c)")):
        ax.grid(True, color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
        for k, _, name in VARIANTS:
            xs = [r["lam"] for r in model[k]]
            ax.plot(xs, [r[key] for r in model[k]], label=name, **styles[k])
        if key == "cp":
            ax.errorbar([e["lam"] for e in exp], [e["cp"] for e in exp], yerr=[e["u_cp"] for e in exp], fmt="o",
                        ms=4.5, color=INK, elinewidth=0.8, capsize=2, label="Battisti et al. (2018), wind tunnel")
        else:
            ax.plot([e["lam"] for e in exp], [e[key] for e in exp], "o", ms=4.5, color=INK)
        ax.set_xlabel("Tip-speed ratio λ"); ax.set_ylabel(ylab); ax.set_title(lab, loc="left", fontsize=9)
        ax.set_xlim(1.2, 3.7)
    a1.set_ylim(-0.1, 0.6); a2.set_ylim(0, 1.1); a3.set_ylim(0, 0.3)
    a1.legend(frameon=False, fontsize=6.3, loc="upper left")
    a4.grid(True, color="#e6e5e0", lw=0.6); a4.set_axisbelow(True)
    for k, _, name in VARIANTS:
        for cds, alpha in (("0.02", 1.0), ("0.05", 0.55)):
            c = mar[f"{k}_cds{cds}"]
            a4.plot([q["lam"] for q in c], [q["cp"] for q in c], alpha=alpha,
                    label=(f"strut $C_D$ {cds}" if k == "central" else None), **styles[k])
    a4.plot([3.3], [0.29], "s", ms=6, color=INK, label="Kjellin et al. (2011), field peak")
    a4.set_xlabel("Tip-speed ratio λ"); a4.set_ylabel("Power coefficient $C_p$")
    a4.set_title("(d) 12 kW H-rotor, Uppsala", loc="left", fontsize=9)
    a4.set_ylim(0, 0.5); a4.set_xlim(2.0, 5.0)
    a4.legend(frameon=False, fontsize=6.3, loc="lower right")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"fig_validacion_turbina.{ext}"), dpi=300)


def main():
    if "--figura" in sys.argv:
        d = json.load(open(os.path.join(HERE, "validacion_dmst.json"), encoding="utf-8"))
        figure(d["battisti_exp"], d["battisti_model"], d["marsta"])
        return
    exp, model, at_exp = battisti()
    mar = marsta()
    # métricas de la variante central
    m = {}
    for k, _, _ in VARIANTS:
        rows = at_exp[k]
        i_pk = int(np.argmax([e["cp"] for e in exp]))
        cp_model_max = max(r["cp"] for r in model[k])
        near = [i for i, e in enumerate(exp) if 2.0 <= e["lam"] <= 2.8]
        m[k] = dict(
            cp_exp_max=exp[i_pk]["cp"], lam_exp_max=exp[i_pk]["lam"],
            cp_model_max=cp_model_max,
            lam_model_max=model[k][int(np.argmax([r["cp"] for r in model[k]]))]["lam"],
            cx_ratio_near_peak=float(np.mean([rows[i]["ct"] / exp[i]["cx"] for i in near])),
            cy_ratio_near_peak=float(np.mean([abs(rows[i]["cs"]) / exp[i]["cy"] for i in near])),
            cx_ratio_all=[rows[i]["ct"] / exp[i]["cx"] for i in range(len(exp))],
            cy_ratio_all=[abs(rows[i]["cs"]) / exp[i]["cy"] for i in range(len(exp))],
            cp_model_at_exp=[r["cp"] for r in rows],
            marsta_cp_max={c: max(p["cp"] for p in mar[f"{k}_cds{c}"]) for c in (0.02, 0.05)},
        )
        m[k]["marsta_ratio_meas_to_model"] = {c: 0.29 / v for c, v in m[k]["marsta_cp_max"].items()}
    json.dump(dict(battisti_exp=exp, battisti_model=model, marsta=mar, metricas=m), open(
        os.path.join(HERE, "validacion_dmst.json"), "w"), indent=1, ensure_ascii=False)

    figure(exp, model, mar)
    for k, _, _ in VARIANTS:
        x = m[k]
        print(f"{k:8s} Cp máx modelo {x['cp_model_max']:.3f} a λ {x['lam_model_max']:.2f} (ensayo {x['cp_exp_max']:.3f} "
              f"a {x['lam_exp_max']:.2f}); Cx/Cx_ens cerca del máximo {x['cx_ratio_near_peak']:.2f}; "
              f"|Cy|/Cy_ens {x['cy_ratio_near_peak']:.2f}; Uppsala: Cp máx {x['marsta_cp_max'][0.02]:.3f}/"
              f"{x['marsta_cp_max'][0.05]:.3f} (brazos 0,02/0,05), medido/modelo "
              f"{x['marsta_ratio_meas_to_model'][0.02]:.2f}/{x['marsta_ratio_meas_to_model'][0.05]:.2f}")


if __name__ == "__main__":
    main()
