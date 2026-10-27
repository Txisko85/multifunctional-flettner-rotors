#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contraste de la resistencia inducida por la fuerza lateral (versión 4; revisión simulada R1-3): equilibrio
estacionario de deriva y guiñada con timón, tipo MMG, frente al modelo de superficie sustentadora del modelo de
ruta, dR = c F_y^2 con c = 1/(pi e q_w (2T)^2).

Modelo (rumbo recto, r = 0, ángulos pequeños):
  casco (Kijima y Nakiri 2003, ec. 6 y 7, con kappa = 2d/L):
      Y'_beta = pi kappa/2 + 1,9257 (C_b B/L) sigma_a,
      N'_beta = kappa [15,0668 X^2 - 23,819 X + 1,802],  X = d (1 - C_b)/B e'_a K,
      e_a = (L/B)(1 - C_pa), e'_a = e_a / sqrt(1/4 + 1/(B/d)^2), sigma_a = (1 - C_wa)/(1 - C_pa),
      K = (1/e'_a + 1,5/(L/B) - 0,33)(0,95 sigma_a + 0,40);
      Y_H = -Q Y'_beta beta (hacia barlovento), N_H = -Q L N'_beta beta (momento de Munk, proa a barlovento),
      Q = rho_w L d U^2 / 2, beta = deriva hacia sotavento;
  timón (Yasukawa y Yoshimura 2015, ec. 18-24 y 38): F_N = rho_w A_R U_R^2 f_alpha alpha_R / 2,
      alpha_R = delta - gamma_R beta, f_alpha = 6,13 Lambda/(Lambda + 2,25), Y_R = -(1 + a_H) F_N,
      N_R = (|x_R| + a_H |x_H|) F_N, X_R = -(1 - t_R) F_N sin(delta);
  viento (rotores y buque): Y_A = F, N_A = F x_A (x_A > 0 a proa de la maestra).
  Equilibrio: Y_H + Y_R + Y_A = 0 y N_H + N_R + N_A = 0, lineal en beta y F_N.
  Resistencia añadida: dR = (1 - t_R) F_N sin(delta) + c_H Y_H beta / 2. El término del casco va de 0 (forma de
  X'_H de Kijima y Nakiri, sin término en beta^2) a la resistencia inducida de un cuerpo esbelto, Y beta/2 (c_H = 1).
  Eficacia equivalente: e_eq = F^2 / (pi q_w (2d)^2 dR), comparable con e = 0,4-0,8 del modelo de ruta (inferida de
  Bickert et al. 2025 con la misma expresión; equivale al doble de la eficacia de envergadura clásica del ala
  reflejada, porque la superficie libre actúa como plano de simetría).
Datos del buque no públicos (supuestos, con rangos): eslora entre perpendiculares 200 m (0,93 de la eslora total),
C_b 0,55-0,65, C_pa 0,65-0,75, C_wa 0,80-0,90; dos timones de área total L d/70 a L d/50 y alargamiento 1,4-2,0;
a_H 0,2-0,35, x'_H = -0,45, x'_R = -0,5, t_R 0,2-0,39, gamma_R 0,4-0,64 y u_R/U 1,0-1,3 (rangos alrededor de los
valores del KVLCC2 de Yasukawa y Yoshimura 2015: a_H 0,312, x'_H -0,464, t_R 0,387, gamma_R 0,395-0,640).
Disposición de las cuatro unidades declarada: x/L = +0,30, +0,10, -0,10, -0,30 (centroide en la maestra); como
sensibilidad, centroide a +-0,2 L. Centro de la fuerza lateral del viento sobre el buque: x_s/L 0 a +0,1 (supuesto).
Salida: maniobra.json.
"""
import json
import math
import numpy as np

RHO_W = 1025.0
L, B, D = 200.0, 27.8, 6.4
KN = 0.514444
SEED = 20261022

CENTRAL = dict(cb=0.60, cpa=0.70, cwa=0.85, ar_ld=1 / 60.0, lam=1.6, a_h=0.30, x_h=-0.45, x_r=-0.5, t_r=0.30,
               gamma_r=0.5, ur=1.15, c_h=0.5)
RANGES = dict(cb=(0.55, 0.65), cpa=(0.65, 0.75), cwa=(0.80, 0.90), ar_ld=(1 / 70.0, 1 / 50.0), lam=(1.4, 2.0),
              a_h=(0.2, 0.35), t_r=(0.2, 0.39), gamma_r=(0.4, 0.64), ur=(1.0, 1.3), c_h=(0.0, 1.0))


def hull_derivatives(cb, cpa, cwa):
    kappa = 2.0 * D / L
    e_a = (L / B) * (1.0 - cpa)
    e_ap = e_a / math.sqrt(0.25 + 1.0 / (B / D) ** 2)
    sigma_a = (1.0 - cwa) / (1.0 - cpa)
    k_form = (1.0 / e_ap + 1.5 / (L / B) - 0.33) * (0.95 * sigma_a + 0.40)
    y_b = 0.5 * math.pi * kappa + 1.9257 * (cb * B / L) * sigma_a
    x = D * (1.0 - cb) / B * e_ap * k_form
    n_b = kappa * (15.0668 * x * x - 23.819 * x + 1.802)
    return y_b, n_b


def equilibrium(force, x_a, v_kn, prm):
    """Deriva (rad), ángulo de timón (rad), fuerza normal del timón (N) y resistencia añadida (N) para una fuerza
    lateral 'force' (N) aplicada en x_a (m a proa de la maestra)."""
    u = v_kn * KN
    q = 0.5 * RHO_W * L * D * u * u
    y_b, n_b = hull_derivatives(prm["cb"], prm["cpa"], prm["cwa"])
    a_r = prm["ar_ld"] * L * D
    f_alpha = 6.13 * prm["lam"] / (prm["lam"] + 2.25)
    k_r = 0.5 * RHO_W * a_r * (prm["ur"] * u) ** 2 * f_alpha
    l_r = (abs(prm["x_r"]) + prm["a_h"] * abs(prm["x_h"])) * L
    beta = force * (1.0 + (1.0 + prm["a_h"]) * x_a / l_r) / (q * (y_b + (1.0 + prm["a_h"]) * L * n_b / l_r))
    f_n = (q * L * n_b * beta - force * x_a) / l_r
    alpha_r = f_n / k_r
    delta = alpha_r + prm["gamma_r"] * beta
    y_h = q * y_b * beta
    d_r = (1.0 - prm["t_r"]) * f_n * math.sin(delta) + prm["c_h"] * 0.5 * y_h * beta
    return dict(beta=beta, delta=delta, f_n=f_n, y_hull=y_h, y_rudder=(1.0 + prm["a_h"]) * f_n, dR=d_r)


def quad_coeffs(v_kn, prm):
    """Forma cuadrática de la resistencia añadida (versión 5): dR = a Y^2 + b Y N + c N^2 (N; Y en N, fuerza lateral
    total hacia sotavento; N en N m, su momento de guiñada respecto a la maestra, positivo con la fuerza a proa), con
    ángulos pequeños (sin delta = delta). Es la misma física que equilibrium(), lineal en (Y, N); los coeficientes
    escalan con 1/U^2."""
    u = v_kn * KN
    q = 0.5 * RHO_W * L * D * u * u
    y_b, n_b = hull_derivatives(prm["cb"], prm["cpa"], prm["cwa"])
    a_r = prm["ar_ld"] * L * D
    f_alpha = 6.13 * prm["lam"] / (prm["lam"] + 2.25)
    k_r = 0.5 * RHO_W * a_r * (prm["ur"] * u) ** 2 * f_alpha
    l_r = (abs(prm["x_r"]) + prm["a_h"] * abs(prm["x_h"])) * L
    den = q * (y_b + (1.0 + prm["a_h"]) * L * n_b / l_r)
    b_y, b_n = 1.0 / den, (1.0 + prm["a_h"]) / (l_r * den)
    f_y, f_n = q * L * n_b * b_y / l_r, (q * L * n_b * b_n - 1.0) / l_r
    d_y, d_n = f_y / k_r + prm["gamma_r"] * b_y, f_n / k_r + prm["gamma_r"] * b_n
    tr, ch = 1.0 - prm["t_r"], prm["c_h"] * q * y_b
    return (tr * f_y * d_y + 0.5 * ch * b_y ** 2, tr * (f_y * d_n + f_n * d_y) + ch * b_y * b_n,
            tr * f_n * d_n + 0.5 * ch * b_n ** 2)


def e_equivalent(x_a, v_kn, prm, force=1e5):
    u = v_kn * KN
    q_w = 0.5 * RHO_W * u * u
    eq = equilibrium(force, x_a, v_kn, prm)
    return force ** 2 / (math.pi * q_w * (2.0 * D) ** 2 * eq["dR"]), eq


def main():
    out = dict(central=CENTRAL, rangos=RANGES, derivadas=dict(zip(("Y'_beta", "N'_beta"), hull_derivatives(
        CENTRAL["cb"], CENTRAL["cpa"], CENTRAL["cwa"]))), casos={})
    yb, nb = out["derivadas"]["Y'_beta"], out["derivadas"]["N'_beta"]
    print(f"Derivadas del casco (Kijima y Nakiri 2003), valores centrales: Y'_beta = {yb:.4f}, N'_beta = {nb:.4f}; "
          f"centro de presión a {nb / yb:.2f} L a proa de la maestra")
    v = 18.0
    rng = np.random.default_rng(SEED)
    samples = [{k: (float(rng.uniform(*RANGES[k])) if k in RANGES else CENTRAL[k]) for k in CENTRAL}
               for _ in range(2000)]
    print(f"\nEficacia equivalente e_eq a {v:.0f} kn (modelo de ruta: e = 0,4-0,8, central 0,6)")
    for lab, xa in (("rotores con el centroide a -0,2 L", -0.2 * L), ("rotores con el centroide en la maestra", 0.0),
                    ("rotores con el centroide a +0,2 L", 0.2 * L), ("viento sobre el buque, x_s = +0,05 L", 0.05 * L),
                    ("viento sobre el buque, x_s = +0,1 L", 0.1 * L)):
        e_c, eq = e_equivalent(xa, v, CENTRAL)
        es = np.array([e_equivalent(xa, v, s)[0] for s in samples])
        es = es[np.isfinite(es) & (es > 0)]
        frac = eq["y_rudder"] / 1e5
        out["casos"][lab] = dict(e_central=e_c, e_p5=float(np.percentile(es, 5)), e_p50=float(np.percentile(es, 50)),
                                 e_p95=float(np.percentile(es, 95)), fraccion_timon=frac,
                                 beta_grados_100kN=math.degrees(eq["beta"]), delta_grados_100kN=math.degrees(eq["delta"]))
        print(f"  {lab:40s} e_eq central {e_c:5.2f}  (P5-P95 {np.percentile(es, 5):.2f}-{np.percentile(es, 95):.2f}); "
              f"timón {100 * frac:5.1f} % de la fuerza; con 100 kN: deriva {math.degrees(eq['beta']):.3f} grados, "
              f"timón {math.degrees(eq['delta']):.2f} grados")
    # caso fuerte (versión 5, revisión simulada v2, menor R1-5): U_a = 20 m/s a 18 kn, viento sobre el buque de
    # través centrado en 0,05 L. (a) Cota: la resultante de las cuatro unidades (|C| q A_F) como si fuera toda
    # lateral, lo que ninguna dirección del viento produce a la vez; (b) viento aparente de través: la fuerza lateral
    # de los rotores es su resistencia, C_D q A_F (la sustentación empuja hacia proa).
    import ruta_modelo as rm
    p = rm.P()
    u_a = 20.0
    k = min(p.k_max, p.us_max / u_a)
    cl, cd = (float(x[0]) for x in rm.rotor_coeffs(np.array([k]), p))
    q = 0.5 * rm.RHO * u_a ** 2
    f_ship = q * p.ayc
    out["caso_fuerte"] = {}
    for lab, f_rot in (("cota_resultante", 4 * q * rm.A_F * math.hypot(cl, cd)), ("traves", 4 * q * rm.A_F * cd)):
        eq = equilibrium(f_rot + f_ship, (f_ship * 0.05 * L) / (f_rot + f_ship), v, CENTRAL)
        out["caso_fuerte"][lab] = dict(F_rotores_kN=f_rot / 1e3, F_buque_kN=f_ship / 1e3, F_kN=(f_rot + f_ship) / 1e3,
                                       beta_grados=math.degrees(eq["beta"]), delta_grados=math.degrees(eq["delta"]),
                                       dR_kN=eq["dR"] / 1e3, k=k, CL=cl, CD=cd)
        print(f"\nCaso fuerte, {lab} (U_a = 20 m/s, 18 kn): rotores {f_rot / 1e3:.0f} kN y buque {f_ship / 1e3:.0f} kN: "
              f"deriva {math.degrees(eq['beta']):.2f} grados, timón {math.degrees(eq['delta']):.2f} grados, "
              f"resistencia añadida {eq['dR'] / 1e3:.1f} kN")
    json.dump(out, open("maniobra.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    return out


if __name__ == "__main__":
    main()
