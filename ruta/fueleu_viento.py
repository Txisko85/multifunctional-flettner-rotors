#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Factor de recompensa del viento de FuelEU Maritime (f_wind) con el método oficial (versión 5, 27-10-2026; revisión
simulada v2, R3-1).

Reglamento (UE) 2023/1805, anexo I: f_wind = 0,99, 0,97 y 0,95 para P_Wind/P_Prop = 0,05, 0,1 y >= 0,15, donde
P_Wind es la potencia efectiva disponible del sistema de propulsión eólica, f_eff P_eff, calculada según la
MEPC.1/Circ.896, y P_Prop es P_ME de las directrices del EEDI (75 % del MCR; más la de los motores de eje en modo PTI,
que aquí no se considera). El texto del anexo no prevé interpolación entre escalones: aquí f_wind se toma escalonado
(el valor del escalón alcanzado) y, como cota superior, interpolado linealmente.

MEPC.1/Circ.896, anexo 1, apartado 2 (categoría B-2):
  (f_eff P_eff) = (1 / sum_{k<=q} W_k) [ (0,5144 V_ref / eta_D) sum_{k<=q} F(V_ref)_k W_k - sum_{k<=q} P(V_ref)_k W_k ],
  con los elementos ordenados por F de mayor a menor y q el número de elementos con el que la probabilidad acumulada
  alcanza 1/2 (la media de la mitad más favorable de las condiciones de viento);
  W_k: matriz global de probabilidad de viento (apéndice 2, tabla 1; datos_coef/circ896_wk.json): viento a 10 m en
      intervalos de 1 m/s (centros 0,5 ... 25,5 m/s) y ángulo respecto a la proa en 72 sectores de 5 grados;
  viento a la altura de referencia (centro de cada rotor) con la ley potencial de exponente 1/9 (apartado 2.4.3.1);
  F(V_ref)_k: fuerza propulsora (kN) del sistema a V_ref con el mejor reglaje en cada condición, incluida la
      resistencia del rotor parado donde no compensa (2.4.3.1.7); P(V_ref)_k: potencia eléctrica que consume (kW);
  efectos secundarios (deriva, timón, escora) ignorados (2.1.2): sin resistencia inducida ni viento sobre el buque;
  eta_D = 0,7 (valor por defecto); V_ref: velocidad del EEDI al 75 % del MCR en aguas tranquilas (aquí, de la ley de
      potencia del modelo de ruta sin margen de mar: 22 (0,75/0,85)^(1/3) = 21,1 kn).
Reglaje: el rotor gira a la relación de giro que maximiza la potencia neta, F V_ref/eta_D - P (el mismo criterio que
el plan de control del modelo de ruta), o se para si eso rinde más.
Opciones: R1 (cuatro rotores), R1x5 (cinco), T3x (rotor con prolongación desplegable), R1-30 (rotor de 30 m) y R1t
(rotor abatible: abatido no tiene resistencia). Escenarios de coeficientes: A, A+, AB, B y A-VG (A ajustado al
Viking Grace, con los extremos 207 y 282 kW).
Para comparar se da también la potencia media de la ruta (empuje menos giro, la aproximación de la versión 4).
Valor del escalón: f_wind multiplica la intensidad de GEI de todo el periodo de notificación; un escalón de 0,01
mejora el balance de cumplimiento en 0,01 I E (I, intensidad real; E, energía en el ámbito de FuelEU) y la
penalización evitada es 0,01 E / 41 000 MJ x 2 400 EUR (anexo IV), sin depender de I.
Salida: fueleu_viento.json.
"""
import json
import math
import os
import numpy as np
import coeficientes as cf
import ruta_modelo as rm

HERE = os.path.dirname(os.path.abspath(__file__))
WK = json.load(open(os.path.join(HERE, "datos_coef", "circ896_wk.json"), encoding="utf-8"))
W = np.array(WK["W"], float)                       # 26 velocidades x 72 ángulos
U10 = np.array(WK["velocidad_10m_centro"], float)
ANG = np.array(WK["angulo_grados"], float)
ALPHA = 1.0 / 9.0
ETA_D_CIRC = 0.7
P_ME = 0.75 * rm.P_MCR
STEPS = ((0.15, 0.95), (0.10, 0.97), (0.05, 0.99))


def v_ref(p):
    return rm.V_DESIGN * (0.75 / rm.LOAD_DESIGN) ** (1.0 / p.n_prop)


def f_wind_step(r):
    for thr, f in STEPS:
        if r >= thr - 1e-12:
            return f
    return 1.0


def f_wind_interp(r):
    xs, fs = [0.0, 0.05, 0.10, 0.15], [1.0, 0.99, 0.97, 0.95]
    if r < 0.05:
        return 1.0
    return float(np.interp(r, xs, fs))


def unit_terms(u_rot, gamma, v_s, p, config):
    """Fuerza propulsora (N) y potencia eléctrica (W) de una unidad con el mejor reglaje (Circ.896)."""
    ratio = (55.0 / 40.0) ** ALPHA
    u_a, b = rm.apparent(u_rot, gamma, v_s)
    if u_a < 0.5:
        return 0.0, 0.0
    q = 0.5 * rm.RHO * u_a ** 2
    sb, cb = math.sin(b), math.cos(b)
    k_max = min(p.k_max, p.us_max / u_a)
    grid = rm.K_GRID if p.k_max <= rm.SR_NOM else np.round(np.arange(0.25, p.k_max + 1e-4, 0.25), 2)
    ks = grid[grid <= k_max + 1e-9]
    if ks.size == 0 or ks[-1] < k_max - 1e-9:
        ks = np.append(ks, k_max)
    cl, cd = rm.rotor_coeffs(ks, p)
    fx = q * rm.A_F * (cl * sb - cd * cb)
    pel = rm.spin_power(u_a, ks, p.cf) / (rm.ETA_MOT * rm.ETA_CONV)
    stop_f = -p.cd0 * q * rm.A_F * cb
    u_at, bt = rm.apparent(u_rot * ratio, gamma, v_s)
    qt = 0.5 * rm.RHO * u_at ** 2
    sbt, cbt = math.sin(bt), math.cos(bt)
    opts = [(stop_f, 0.0)] if config != "R1t" else [(0.0, 0.0)]
    if config in ("R1", "R1t", "T3x"):
        opts += list(zip(fx, pel))
    if config == "T3x":
        u_sx = ks * u_a * rm.D_X / rm.D_F
        clx, cdx = rm.rotor_coeffs(u_sx / u_at, p)
        fx_x = qt * rm.A_X * (p.cl_ext * clx * sbt - cdx * cbt)
        p_x = 0.5 * rm.RHO * p.cf * u_sx ** 3 * math.pi * rm.D_X * rm.H_X / (rm.ETA_MOT * rm.ETA_CONV)
        opts += list(zip(fx + fx_x, pel + p_x))
    if config == "R1-30":
        u_s30 = ks * u_a
        clu, cdu = rm.rotor_coeffs(u_s30 / u_at, p)
        fx_u = qt * rm.A_UP * (p.f30 * clu * sbt - cdu * cbt)
        p_u = 0.5 * rm.RHO * p.cf * u_s30 ** 3 * math.pi * rm.D_F * rm.H_X / (rm.ETA_MOT * rm.ETA_CONV)
        opts = [(stop_f - p.cd0 * qt * rm.A_UP * cbt, 0.0)] + list(zip(fx + fx_u, pel + p_u))
    best = max(opts, key=lambda o: o[0] * v_s / ETA_D_CIRC - o[1])
    return float(best[0]), float(best[1])


def p_wind_circ(p, config="R1", n_units=4):
    """P_Wind (kW) de n_units unidades con la Circ.896, y la fracción de probabilidad usada."""
    v_kn = v_ref(p)
    v = v_kn * rm.KN
    F, P = np.zeros_like(W), np.zeros_like(W)
    for i, u10 in enumerate(U10):
        u_rot = u10 * (40.0 / 10.0) ** ALPHA
        for j, ang in enumerate(ANG):
            g = ang if ang <= 180.0 else 360.0 - ang
            if j > 36:                                   # simetría babor-estribor
                F[i, j], P[i, j] = F[i, 72 - j], P[i, 72 - j]
                continue
            f, pw = unit_terms(u_rot, g, v, p, config)
            F[i, j], P[i, j] = n_units * f / 1e3, n_units * pw / 1e3      # kN, kW
    order = np.argsort(-F, axis=None, kind="stable")
    wf, ff, pf = W.ravel()[order], F.ravel()[order], P.ravel()[order]
    cum = np.cumsum(wf)
    q = int(np.searchsorted(cum, 0.5)) + 1
    ws = wf[:q].sum()
    p_wind = (0.5144 * v_kn / ETA_D_CIRC * (ff[:q] * wf[:q]).sum() - (pf[:q] * wf[:q]).sum()) / ws
    p_all = (0.5144 * v_kn / ETA_D_CIRC * (F * W).sum() - (P * W).sum()) / W.sum()
    return dict(p_wind_kw=float(p_wind), q=q, prob_usada=float(ws), v_ref_kn=v_kn, p_wind_toda_la_matriz_kw=float(p_all))


def p_wind_route(groups, nu, ng, p, n_units=4):
    """Aproximación de la versión 4: potencia media del empuje menos la de giro en la ruta (abril-octubre)."""
    acc = hrs = 0.0
    for s in rm.BASE_SEASONS:
        for name, direction, h in rm.crossings(p):
            v_kn, t_open = rm.crossing_speed(h, p)
            mb = rm.mean_benefits(groups, nu, ng, s, direction, v_kn, p, name=name)
            acc += (mb["_R1 empuje"] - mb["_R1 giro"]) * t_open
            hrs += t_open
    return n_units * acc / hrs


ESCENARIOS = ("A", "A+", "AB", "B", "A-VG 207", "A-VG", "A-VG 282")
OPCIONES = (("R1", "R1", 4), ("R1x5", "R1", 5), ("T3x", "T3x", 4), ("R1-30", "R1-30", 4), ("R1t", "R1t", 4))


def main():
    d, groups, nu, ng = rm.load_wind()
    out = dict(metodo="MEPC.1/Circ.896, anexo 1, apartado 2; Reglamento (UE) 2023/1805, anexo I", p_me_kw=P_ME,
               escalones=STEPS, v_ref_kn=v_ref(rm.P()), resultados={})
    for esc in ESCENARIOS:
        p = rm.P(**cf.ESCENARIOS_COEF[esc])
        res = {}
        for lab, cfg, n in OPCIONES:
            r = p_wind_circ(p, cfg, n)
            r["ratio"] = r["p_wind_kw"] / P_ME
            r["f_wind"] = f_wind_step(r["ratio"])
            r["f_wind_interpolado"] = f_wind_interp(r["ratio"])
            res[lab] = r
        res["R1"]["ruta_aprox_kw"] = p_wind_route(groups, nu, ng, p, 4)
        res["R1"]["ruta_aprox_ratio"] = res["R1"]["ruta_aprox_kw"] / P_ME
        out["resultados"][esc] = res
        print(f"{esc:9s} " + "  ".join(f"{lab} {v['p_wind_kw']:6.0f} kW ({v['ratio']:.3f}, f {v['f_wind']:.2f}/"
                                         f"{v['f_wind_interpolado']:.3f})" for lab, v in res.items())
              + f"  | ruta (v4) {res['R1']['ruta_aprox_kw']:.0f} kW ({res['R1']['ruta_aprox_ratio']:.3f})")
    json.dump(out, open(os.path.join(HERE, "fueleu_viento.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    return out


if __name__ == "__main__":
    main()
