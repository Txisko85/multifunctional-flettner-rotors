#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check of the resistance induced by the side force (version 4; simulated review R1-3): steady sway-yaw balance with
rudder, MMG type, against the lifting-surface model of the route model, dR = c F_y^2 with c = 1/(pi e q_w (2T)^2).

Model (straight course, r = 0, small angles):
  hull (Kijima and Nakiri 2003, eqs. 6 and 7, with kappa = 2d/L):
      Y'_beta = pi kappa/2 + 1.9257 (C_b B/L) sigma_a,
      N'_beta = kappa [15.0668 X^2 - 23.819 X + 1.802],  X = d (1 - C_b)/B e'_a K,
      e_a = (L/B)(1 - C_pa), e'_a = e_a / sqrt(1/4 + 1/(B/d)^2), sigma_a = (1 - C_wa)/(1 - C_pa),
      K = (1/e'_a + 1.5/(L/B) - 0.33)(0.95 sigma_a + 0.40);
      Y_H = -Q Y'_beta beta (to windward), N_H = -Q L N'_beta beta (Munk moment, bow to windward),
      Q = rho_w L d U^2 / 2, beta = drift to leeward;
  rudder (Yasukawa and Yoshimura 2015, eqs. 18-24 and 38): F_N = rho_w A_R U_R^2 f_alpha alpha_R / 2,
      alpha_R = delta - gamma_R beta, f_alpha = 6.13 Lambda/(Lambda + 2.25), Y_R = -(1 + a_H) F_N,
      N_R = (|x_R| + a_H |x_H|) F_N, X_R = -(1 - t_R) F_N sin(delta);
  wind (rotors and ship): Y_A = F, N_A = F x_A (x_A > 0 forward of midships).
  Balance: Y_H + Y_R + Y_A = 0 and N_H + N_R + N_A = 0, linear in beta and F_N.
  Added resistance: dR = (1 - t_R) F_N sin(delta) + c_H Y_H beta / 2. The hull term ranges from 0 (form of X'_H of
  Kijima and Nakiri, without a beta^2 term) to the induced resistance of a slender body, Y beta/2 (c_H = 1).
  Equivalent efficiency: e_eq = F^2 / (pi q_w (2d)^2 dR), comparable with e = 0.4-0.8 of the route model (inferred
  from Bickert et al. 2025 with the same expression; it equals twice the classical span efficiency of the mirrored
  wing, because the free surface acts as a symmetry plane).
Non-public ship data (assumed, with ranges): length between perpendiculars 200 m (0.93 of the length overall),
C_b 0.55-0.65, C_pa 0.65-0.75, C_wa 0.80-0.90; two rudders of total area L d/70 to L d/50 and aspect ratio 1.4-2.0;
a_H 0.2-0.35, x'_H = -0.45, x'_R = -0.5, t_R 0.2-0.39, gamma_R 0.4-0.64 and u_R/U 1.0-1.3 (ranges around the values
of the KVLCC2 in Yasukawa and Yoshimura 2015: a_H 0.312, x'_H -0.464, t_R 0.387, gamma_R 0.395-0.640).
Declared layout of the four units: x/L = +0.30, +0.10, -0.10, -0.30 (centroid at midships); as a sensitivity,
centroid at +-0.2 L. Centre of the side force of the wind on the ship: x_s/L 0 to +0.1 (assumed).
Output: sway_yaw.json.
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
    """Drift (rad), rudder angle (rad), rudder normal force (N) and added resistance (N) for a side force 'force' (N)
    applied at x_a (m forward of midships)."""
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
    """Quadratic form of the added resistance (version 5): dR = a Y^2 + b Y N + c N^2 (N; Y in N, total side force to
    leeward; N in N m, its yaw moment about midships, positive with the force forward), with small angles
    (sin delta = delta). It is the same physics as equilibrium(), linear in (Y, N); the coefficients scale with
    1/U^2."""
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
    out = dict(central=CENTRAL, ranges=RANGES, derivatives=dict(zip(("Y'_beta", "N'_beta"), hull_derivatives(
        CENTRAL["cb"], CENTRAL["cpa"], CENTRAL["cwa"]))), cases={})
    yb, nb = out["derivatives"]["Y'_beta"], out["derivatives"]["N'_beta"]
    print(f"Hull derivatives (Kijima and Nakiri 2003), central values: Y'_beta = {yb:.4f}, N'_beta = {nb:.4f}; "
          f"centre of pressure {nb / yb:.2f} L forward of midships")
    v = 18.0
    rng = np.random.default_rng(SEED)
    samples = [{k: (float(rng.uniform(*RANGES[k])) if k in RANGES else CENTRAL[k]) for k in CENTRAL}
               for _ in range(2000)]
    print(f"\nEquivalent efficiency e_eq at {v:.0f} kn (route model: e = 0.4-0.8, central 0.6)")
    for lab, xa in (("rotors with the centroid at -0.2 L", -0.2 * L), ("rotors with the centroid at midships", 0.0),
                    ("rotors with the centroid at +0.2 L", 0.2 * L), ("wind on the ship, x_s = +0.05 L", 0.05 * L),
                    ("wind on the ship, x_s = +0.1 L", 0.1 * L)):
        e_c, eq = e_equivalent(xa, v, CENTRAL)
        es = np.array([e_equivalent(xa, v, s)[0] for s in samples])
        es = es[np.isfinite(es) & (es > 0)]
        frac = eq["y_rudder"] / 1e5
        out["cases"][lab] = dict(e_central=e_c, e_p5=float(np.percentile(es, 5)), e_p50=float(np.percentile(es, 50)),
                                 e_p95=float(np.percentile(es, 95)), rudder_fraction=frac,
                                 beta_deg_100kN=math.degrees(eq["beta"]), delta_deg_100kN=math.degrees(eq["delta"]))
        print(f"  {lab:40s} e_eq central {e_c:5.2f}  (P5-P95 {np.percentile(es, 5):.2f}-{np.percentile(es, 95):.2f}); "
              f"rudder {100 * frac:5.1f} % of the force; with 100 kN: drift {math.degrees(eq['beta']):.3f} degrees, "
              f"rudder {math.degrees(eq['delta']):.2f} degrees")
    # strong case (version 5, simulated review v2, minor R1-5): U_a = 20 m/s at 18 kn, wind on the ship from abeam
    # centred at 0.05 L. (a) Bound: the resultant of the four units (|C| q A_F) as if it were all lateral, which no
    # wind direction produces at once; (b) apparent wind from abeam: the side force of the rotors is their drag,
    # C_D q A_F (the lift pushes forward).
    import route_model as rm
    p = rm.P()
    u_a = 20.0
    k = min(p.k_max, p.us_max / u_a)
    cl, cd = (float(x[0]) for x in rm.rotor_coeffs(np.array([k]), p))
    q = 0.5 * rm.RHO * u_a ** 2
    f_ship = q * p.ayc
    out["strong_case"] = {}
    for lab, f_rot in (("resultant_bound", 4 * q * rm.A_F * math.hypot(cl, cd)), ("abeam", 4 * q * rm.A_F * cd)):
        eq = equilibrium(f_rot + f_ship, (f_ship * 0.05 * L) / (f_rot + f_ship), v, CENTRAL)
        out["strong_case"][lab] = dict(F_rotors_kN=f_rot / 1e3, F_ship_kN=f_ship / 1e3, F_kN=(f_rot + f_ship) / 1e3,
                                       beta_deg=math.degrees(eq["beta"]), delta_deg=math.degrees(eq["delta"]),
                                       dR_kN=eq["dR"] / 1e3, k=k, CL=cl, CD=cd)
        print(f"\nStrong case, {lab} (U_a = 20 m/s, 18 kn): rotors {f_rot / 1e3:.0f} kN and ship {f_ship / 1e3:.0f} kN: "
              f"drift {math.degrees(eq['beta']):.2f} degrees, rudder {math.degrees(eq['delta']):.2f} degrees, "
              f"added resistance {eq['dR'] / 1e3:.1f} kN")
    json.dump(out, open("sway_yaw.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    return out


if __name__ == "__main__":
    main()
