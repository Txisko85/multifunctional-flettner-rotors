#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validación tridimensional del modelo de línea sustentadora (versión 5, 27-10-2026; revisión simulada v2, R2-3).

1. Dependencia de la sustentación con el alargamiento: De Marco et al. (2016, Fig. 10; URANS con el rotor sobre un
   plano de simetría, como en los ensayos de Badalamenti y Prince) dan C_L para AR = 2, 4, 6 y 8 y placas de extremo
   de/d = 1 (sin placa), 2 y 3, a k = 1,5, 2 y 3 (datos_coef/demarco2016_fig10.json). El modelo (cubierta como espejo
   perfecto, m_d = 1; placa con m_p = 0 sin placa y m_p = 0,5 o 1 con placa) predice C_L(AR)/C_L(AR 6) con los dos
   modelos de sección. Se compara el cociente, no el valor absoluto: el modelo da factores relativos (f_X, f30,
   eta_B), y el nivel absoluto lo fija el escenario de coeficientes.
2. Factor empírico del rotor de 30 m: con placa de/d = 2, f30 = [(AR'/6) r - 1]/((AR' - 6)/6), con r = C_L(AR')/
   C_L(AR 6) interpolado en AR' = 7,5 (30 m / 4 m).
3. Insensibilidad de los factores al nivel de C_l,2D: f_X y f30 con C_l,2D del conjunto B multiplicado por 1,7 (el
   cociente A/B a k = 3) y con la forma del conjunto A. La analogía de ala es lineal: sus factores no dependen ni del
   nivel ni de la forma; el modelo de giro no depende del nivel, pero sí de la forma a k >= 2,5.
4. Estimación de f_X: f30 empírico por el cociente f_X/f30 del modelo (casi constante en todas las variantes).
Salida: linea_sustentadora_validacion.json.
"""
import json
import os
import numpy as np
import coeficientes as cf
import linea_sustentadora as ls

HERE = os.path.dirname(os.path.abspath(__file__))
DM = json.load(open(os.path.join(HERE, "datos_coef", "demarco2016_fig10.json"), encoding="utf-8"))["datos"]
D = 4.0
KS = (1.5, 2.0, 3.0)
ARS = (2, 4, 6, 8)


def cl3d(ar, k, modelo, m_p):
    h = ar * D
    r = ls.Rotor([(0.0, h, D)], z_plate=h if m_p > 0 else None, m_d=1.0, m_p=m_p)
    s = r.solve(k, modelo)
    return r.tramo(s, 0, h)[0] / (0.5 * D * h)


def ratios_obs(de, k):
    return {ar: DM[de][f"AR{ar}"]["CL"][str(k)] / DM[de]["AR6"]["CL"][str(k)] for ar in ARS}


def validacion():
    out = {}
    for modelo in ("ala", "giro"):
        for de, m_ps in (("de1", (0.0,)), ("de2", (0.5, 1.0)), ("de3", (0.5, 1.0))):
            for m_p in m_ps:
                rows, err = {}, []
                for k in KS:
                    cl = {ar: cl3d(ar, k, modelo, m_p) for ar in ARS}
                    mod = {ar: cl[ar] / cl[6] for ar in ARS}
                    obs = ratios_obs(de, k)
                    rows[str(k)] = dict(modelo={str(a): mod[a] for a in ARS}, URANS={str(a): obs[a] for a in ARS},
                                        eta_AR6=cl[6] / float(ls.cl2d(k)))
                    err += [mod[a] - obs[a] for a in ARS if a != 6]
                out[f"{modelo}, {de}, m_p {m_p:g}"] = dict(k=rows, rms=float(np.sqrt(np.mean(np.square(err)))),
                                                         sesgo=float(np.mean(err)))
    return out


def f30_empirico():
    out = {}
    for de in ("de2", "de3"):
        for k in KS:
            r6, r8 = 1.0, ratios_obs(de, k)[8]
            r75 = r6 + (r8 - r6) * 0.75
            out[f"{de}, k {k:g}"] = ((7.5 / 6.0) * r75 - 1.0) / (1.5 / 6.0)
    return out


def insensibilidad():
    base_cl = ls.CL_2D.copy()
    k_a = np.array(cf.K_A_EXT)
    cl_a = np.array(cf.CL_A_EXT)
    variantes = {"B": base_cl, "B x 1,7": 1.7 * base_cl, "forma de A": np.interp(ls.K_2D, k_a, cl_a)}
    out = {}
    for lab, arr in variantes.items():
        ls.CL_2D = arr
        for modelo in ("ala", "giro"):
            for m_p in (0.5, 1.0):
                f = ls.factors(1.0, m_p, modelo)
                out[f"{lab}; {modelo}; m_d 1, m_p {m_p:g}"] = dict(f_X=f["f_X"], f30=f["f30"], eta0=f["eta0"], k=f["k"],
                                                                 f_X_f30=[a / b for a, b in zip(f["f_X"], f["f30"])])
    ls.CL_2D = base_cl
    return out


def main():
    res = dict(validacion=validacion(), f30_empirico=f30_empirico(), insensibilidad=insensibilidad())
    # f_X estimado: f30 empírico (sin el punto atípico de/d = 2, k = 2, AR 8) por el cociente f_X/f30 del modelo, que
    # apenas depende del modelo de sección, del nivel y de la forma de C_l,2D
    f30e = [v for k, v in res["f30_empirico"].items() if k != "de2, k 2"]
    ratio = [x for v in res["insensibilidad"].values() for x in v["f_X_f30"]]
    res["f_X_estimado"] = dict(f30_empirico=(min(f30e), max(f30e)), f_X_f30=(min(ratio), max(ratio)),
                               f_X=(min(f30e) * min(ratio), max(f30e) * max(ratio)))
    json.dump(res, open(os.path.join(HERE, "linea_sustentadora_validacion.json"), "w"), indent=1, ensure_ascii=False)
    print("C_L(AR)/C_L(AR 6): modelo frente a URANS (De Marco et al. 2016)")
    for key, v in res["validacion"].items():
        print(f"  {key:22s} RMS {v['rms']:.3f}  sesgo {v['sesgo']:+.3f}  " + "; ".join(
            f"k {k}: " + "/".join(f"{v['k'][k]['modelo'][a]:.2f}~{v['k'][k]['URANS'][a]:.2f}" for a in ("2", "4", "8"))
            for k in v["k"]))
    print("f30 empírico (De Marco, AR 7,5):", {k: round(x, 2) for k, x in res["f30_empirico"].items()})
    print("f_X estimado:", {k: tuple(round(x, 2) for x in v) for k, v in res["f_X_estimado"].items()})
    print("Factores con distintos C_l,2D:")
    for k, v in res["insensibilidad"].items():
        print(f"  {k:36s} f_X {np.round(v['f_X'], 3)}  f30 {np.round(v['f30'], 3)}  eta0 {np.round(v['eta0'], 3)}  "
              f"f_X/f30 {np.round(v['f_X_f30'], 3)}")
    return res


if __name__ == "__main__":
    main()
