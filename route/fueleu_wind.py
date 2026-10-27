#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wind reward factor of FuelEU Maritime (f_wind) with the official method (version 5, 27-10-2026; simulated review
v2, R3-1).

Regulation (EU) 2023/1805, annex I: f_wind = 0.99, 0.97 and 0.95 for P_Wind/P_Prop = 0.05, 0.1 and >= 0.15, where
P_Wind is the available effective power of the wind-assisted propulsion system, f_eff P_eff, computed according to
MEPC.1/Circ.896, and P_Prop is P_ME of the EEDI guidelines (75 % of MCR; plus that of the shaft motors in PTI mode,
not considered here). The text of the annex provides no interpolation between steps: here f_wind is taken stepwise
(the value of the step reached) and, as an upper bound, linearly interpolated.

MEPC.1/Circ.896, annex 1, section 2 (category B-2):
  (f_eff P_eff) = (1 / sum_{k<=q} W_k) [ (0.5144 V_ref / eta_D) sum_{k<=q} F(V_ref)_k W_k - sum_{k<=q} P(V_ref)_k W_k ],
  with the elements sorted by F from largest to smallest and q the number of elements with which the cumulative
  probability reaches 1/2 (the mean of the most favourable half of the wind conditions);
  W_k: global wind probability matrix (appendix 2, table 1; coef_data/circ896_wk.json): wind at 10 m in bins of 1 m/s
      (centres 0.5 ... 25.5 m/s) and angle relative to the bow in 72 sectors of 5 degrees;
  wind at the reference height (centre of each rotor) with the power law of exponent 1/9 (section 2.4.3.1);
  F(V_ref)_k: propulsive force (kN) of the system at V_ref with the best setting in each condition, including the drag
      of the stopped rotor where spinning does not pay (2.4.3.1.7); P(V_ref)_k: electrical power it consumes (kW);
  secondary effects (drift, rudder, heel) ignored (2.1.2): no induced resistance and no wind on the ship;
  eta_D = 0.7 (default value); V_ref: EEDI speed at 75 % of MCR in calm water (here, from the power law of the route
      model without sea margin: 22 (0.75/0.85)^(1/3) = 21.1 kn).
Setting: the rotor spins at the spin ratio that maximises the net power, F V_ref/eta_D - P (the same criterion as the
control plan of the route model), or stops if that yields more.
Options: R1 (four rotors), R1x5 (five), T3x (rotor with deployable extension), R1-30 (30 m rotor) and R1t (tiltable
rotor: tilted, it has no drag). Coefficient scenarios: A, A+, AB, B and A-VG (A fitted to the Viking Grace, with the
ends 207 and 282 kW).
For comparison, the mean power on the route is also given (thrust minus spin, the approximation of version 4).
Value of a step: f_wind multiplies the GHG intensity of the whole reporting period; a step of 0.01 improves the
compliance balance by 0.01 I E (I, actual intensity; E, energy within the scope of FuelEU) and the penalty avoided is
0.01 E / 41,000 MJ x 2,400 EUR (annex IV), independent of I.
Output: fueleu_wind.json.
"""
import json
import math
import os
import numpy as np
import coefficients as cf
import route_model as rm

HERE = os.path.dirname(os.path.abspath(__file__))
WK = json.load(open(os.path.join(HERE, "coef_data", "circ896_wk.json"), encoding="utf-8"))
W = np.array(WK["W"], float)                       # 26 speeds x 72 angles
U10 = np.array(WK["wind_speed_10m_centre"], float)
ANG = np.array(WK["angle_deg"], float)
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
    """Propulsive force (N) and electrical power (W) of one unit with the best setting (Circ.896)."""
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
    """P_Wind (kW) of n_units units with Circ.896, and the fraction of probability used."""
    v_kn = v_ref(p)
    v = v_kn * rm.KN
    F, P = np.zeros_like(W), np.zeros_like(W)
    for i, u10 in enumerate(U10):
        u_rot = u10 * (40.0 / 10.0) ** ALPHA
        for j, ang in enumerate(ANG):
            g = ang if ang <= 180.0 else 360.0 - ang
            if j > 36:                                   # port-starboard symmetry
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
    return dict(p_wind_kw=float(p_wind), q=q, prob_used=float(ws), v_ref_kn=v_kn, p_wind_whole_matrix_kw=float(p_all))


def p_wind_route(groups, nu, ng, p, n_units=4):
    """Approximation of version 4: mean power of the thrust minus that of the spin on the route (April-October)."""
    acc = hrs = 0.0
    for s in rm.BASE_SEASONS:
        for name, direction, h in rm.crossings(p):
            v_kn, t_open = rm.crossing_speed(h, p)
            mb = rm.mean_benefits(groups, nu, ng, s, direction, v_kn, p, name=name)
            acc += (mb["_R1 thrust"] - mb["_R1 spin"]) * t_open
            hrs += t_open
    return n_units * acc / hrs


SCENARIOS = ("A", "A+", "AB", "B", "A-VG 207", "A-VG", "A-VG 282")
OPTIONS = (("R1", "R1", 4), ("R1x5", "R1", 5), ("T3x", "T3x", 4), ("R1-30", "R1-30", 4), ("R1t", "R1t", 4))


def main():
    d, groups, nu, ng = rm.load_wind()
    out = dict(method="MEPC.1/Circ.896, annex 1, section 2; Regulation (EU) 2023/1805, annex I", p_me_kw=P_ME,
               steps=STEPS, v_ref_kn=v_ref(rm.P()), results={})
    for esc in SCENARIOS:
        p = rm.P(**cf.COEF_SCENARIOS[esc])
        res = {}
        for lab, cfg, n in OPTIONS:
            r = p_wind_circ(p, cfg, n)
            r["ratio"] = r["p_wind_kw"] / P_ME
            r["f_wind"] = f_wind_step(r["ratio"])
            r["f_wind_interpolated"] = f_wind_interp(r["ratio"])
            res[lab] = r
        res["R1"]["route_approx_kw"] = p_wind_route(groups, nu, ng, p, 4)
        res["R1"]["route_approx_ratio"] = res["R1"]["route_approx_kw"] / P_ME
        out["results"][esc] = res
        print(f"{esc:9s} " + "  ".join(f"{lab} {v['p_wind_kw']:6.0f} kW ({v['ratio']:.3f}, f {v['f_wind']:.2f}/"
                                         f"{v['f_wind_interpolated']:.3f})" for lab, v in res.items())
              + f"  | route (v4) {res['R1']['route_approx_kw']:.0f} kW ({res['R1']['route_approx_ratio']:.3f})")
    json.dump(out, open(os.path.join(HERE, "fueleu_wind.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    return out


if __name__ == "__main__":
    main()
