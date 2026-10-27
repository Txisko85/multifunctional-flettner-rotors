#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Route analysis, version 5 (27-10-2026), after the second simulated review: reference case (April-October season of
the 2026 timetable, histograms conditioned on the timetable, induced resistance from the sway-yaw balance, coefficients
of scenario A), scenario A-VG fitted to the Viking Grace, non-patented comparators (tiltable rotor, 30 m rotor, rotor of
5 m diameter, fifth unit and tiltable rotor with extension), common structural envelope, breakdown by crossing,
step by step from version 4, scenarios, one-factor sensitivities, interaction threshold, rules of the Savonius
conversion and its break-even C_t,S (route_results.json).
Version 4 (22-10-2026): April-November service with the histograms of all hours (P(v4=True) and era5_route.json
reproduce it)."""
import json
import coefficients as cf
import route_model as rm
import moments as mom
from route_model import P, annual, table, break_even, load_wind, mean_benefits, crossing_speed

d, groups, nu, ng = load_wind()
base = annual(d, groups, nu, ng)
KEYS_SEA = ("R1", "T3x", "T3", "T4", "T1", "T2", "T5", "R1t", "R1-30", "R1-D5", "R1t+X", "_pen R1", "_pen own",
            "_pen cross",
            "_R1 thrust", "_R1 spin", "_R1 stopped", "_R0 wind", "T5 claim1", "T5 intended", "T5 astern",
            "T2 claim1", "T2 intended", "T3 claim1", "T3 intended")


def summ(r, extra=None):
    t = {row[0]: row[3] for row in table(r)}
    out = {"R0": round(r["fuel_total"], 1), "R0_sea": round(r["fuel_sea"], 1), "R0_berth": round(r["fuel_berth"], 1),
           "R1": round(r["sea"]["R1"], 2), "R1_pct": round(100 * r["sea"]["R1"] / r["fuel_total"], 3),
           "pen_R1": round(r["sea"]["_pen R1"], 2)}
    for k in ("T3x", "T3", "T4", "T1", "T2", "T5"):
        out[k] = round(t[k], 2)
    for k in ("R1t", "R1-30", "R1-D5", "R1t+X", "T5 claim1", "T5 intended", "T5 astern"):
        out[k] = round(r["sea"][k], 2)
    out["T5_sea"] = round(r["sea"]["T5"], 2)
    out["thrust"], out["spin"], out["stopped"] = (round(r["sea"][k], 2) for k in ("_R1 thrust", "_R1 spin",
                                                                                  "_R1 stopped"))
    out["pen_own"], out["pen_cross"] = round(r["sea"]["_pen own"], 2), round(r["sea"]["_pen cross"], 2)
    if extra:
        out.update(extra)
    return out


# ---------------- reference case ----------------
r5 = annual(d, groups, nu, ng, P(n_units=5))
res = {"version": 5,
       "base": {"fuel_total": base["fuel_total"], "fuel_sea": base["fuel_sea"], "fuel_berth": base["fuel_berth"],
                "fuel_me": base["fuel_me"], "load_mean": base["load_mean"], "seasons": base["seasons"],
                "sea": {k: base["sea"][k] for k in KEYS_SEA}, "port": base["port"], "rows": table(base),
                "port_kw": base["port_kw"], "hours_sea": base["hours_sea"], "hours_port": base["hours_port"],
                "R1_fifth_unit": r5["sea"]["R1"] - base["sea"]["R1"], "R1_five_units": r5["sea"]["R1"],
                "parameters": {k: v for k, v in vars(P()).items() if isinstance(v, (int, float, bool)) or v is None}}}

# ---------------- breakdown by crossing (mean kW per unit, open sea) ----------------
per = []
for name, direction, hours in rm.crossings(P()):
    v, t_open = crossing_speed(hours, P())
    row = {"crossing": name, "hours": hours, "v_kn": round(v, 2), "h_sea": round(t_open, 2)}
    for s_, lab in ((0, "Apr-Oct"), (2, "Nov"), (1, "Dec-Mar")):
        mb = mean_benefits(groups, nu, ng, s_, direction, v, P(), name=name)
        row[lab] = {k: round(mb[k], 2) for k in ("R1", "_R1 thrust", "_R1 spin", "_R1 stopped", "_pen R1",
                                                    "_R1 running", "T3x", "T4", "T1", "T5", "R1t", "R1-30")}
        row[lab]["_R0 wind kW ship"] = round(mb["_R0 wind"], 1)
    per.append(row)
res["per_crossing_kW_unit"] = per

# ---------------- interaction threshold (T4 = T3x) ----------------
res["eps_star"] = break_even(d, groups, nu, ng, base)

# ---------------- step by step from version 4 ----------------
d4, g4, _, _ = load_wind("era5_route.json")
steps = [
    ("version 4 (26-10-2026): April-November, histograms of all hours", (d4, g4), P(v4=True), (0,)),
    ("+ histograms conditioned on the timetable (April-November)", (d, groups), P(v4=True), rm.NOV_SEASONS),
    ("+ April-October season", (d, groups), P(v4=True), None),
    ("+ induced resistance from the sway-yaw balance", (d, groups), P(v4=True, sir_model="mmg"), None),
    ("+ recalibrated sea margin (0.077)", (d, groups),
     P(v4=True, sir_model="mmg", sea_margin=rm.SEA_MARGIN[3.0]), None),
    ("+ strict Savonius bound (central C_t,S 0.62)", (d, groups),
     P(v4=True, sir_model="mmg", sea_margin=rm.SEA_MARGIN[3.0], ct_s=cf.CENTRAL_V4["ct_s"]), None),
    ("+ f30 checked against De Marco et al. (central 1.025): version 5", (d, groups), P(), None),
]
res["step_by_step"] = {lab: summ(annual(dd[0], dd[1], nu, ng, p, seasons=ss)) for lab, dd, p, ss in steps}

# ---------------- epistemic scenarios ----------------
UA_LIM = mom.ext_limit(P())[0]
res["ua_ext_max_structural"] = UA_LIM
scen = {f"coefficients {k}": (P(**v), None) for k, v in cf.COEF_SCENARIOS.items()}
scen.update({
    "April-November (with November)": (P(), rm.NOV_SEASONS),
    "whole year (two rotations all year)": (P(), rm.ALL_SEASONS),
    "November only": (P(), (2,)),
    "December-March only": (P(), (1,)),
    "exponent 4 of the power-speed curve": (P(n_prop=4.0), None),
    "five units": (P(n_units=5), None),
    "speed -2 kn": (P(dv=-2.0), None),
    "speed -4 kn": (P(dv=-4.0), None),
    "speed -6 kn": (P(dv=-6.0), None),
    "tandem interference: -10 % lift downstream (|beta| < 30°)": (P(int_frac=0.10), None),
    "tandem interference: -20 % lift downstream (|beta| < 30°)": (P(int_frac=0.20), None),
    "tandem interference: -10 % downstream (|beta| < 60°)": (P(int_frac=0.10, int_beta=60.0), None),
    "tandem interference: -20 % downstream (|beta| < 60°)": (P(int_frac=0.20, int_beta=60.0), None),
    "induced: low-aspect-ratio wing, e = 0.6 (version 4)": (P(sir_model="wing", e_span=0.6), None),
    "induced: low-aspect-ratio wing, e = 0.4": (P(sir_model="wing", e_span=0.4), None),
    "induced: units with the centroid 0.2 L aft": (P(x_u=-0.2), None),
    "induced: units with the centroid 0.2 L forward": (P(x_u=0.2), None),
    "induced: c_H = 0": (P(c_h=0.0), None),
    "induced: c_H = 0.25 (KVLCC2, Yasukawa and Yoshimura 2015)": (P(c_h=0.25), None),
    "induced: c_H = 1": (P(c_h=1.0), None),
    "induced: wind on the ship centred at x_s = 0.05 L": (P(x_s=0.05), None),
    "induced: wind on the ship centred at x_s = 0.10 L": (P(x_s=0.10), None),
    "f_X 0.5 (endplate not moved)": (P(cl_ext=0.5), None),
    "f_X 0.90": (P(cl_ext=0.90), None),
    "f_X 1.30": (P(cl_ext=1.30), None),
    "f30 0.90": (P(f30=0.90), None),
    "f30 1.15": (P(f30=1.15), None),
    f"structural limit of the extension Ua <= {UA_LIM:.2f} m/s (version 4)": (P(ua_ext_max=UA_LIM), None),
    "common envelope: moment of the reference rotor at 20 m/s": (P(m_cap="ref"), None),
    "common envelope: moment of the 30 m rotor at 20 m/s (reinforced foundation)": (P(m_cap="30m"), None),
    "A-VG, common envelope of the reference rotor": (P(m_cap="ref", **cf.COEF_SCENARIOS["A-VG"]), None),
    "A-VG, common envelope of the 30 m rotor": (P(m_cap="30m", **cf.COEF_SCENARIOS["A-VG"]), None),
    "interaction eps = -0.01": (P(eps=-0.01), None),
    "interaction eps = +0.01": (P(eps=0.01), None),
    "side force of the Savonius C_y,S = +0.5": (P(cy_s=0.5), None),
    "side force of the Savonius C_y,S = -0.5": (P(cy_s=-0.5), None),
    "wind at the turbine x0.9 (slowed by the rotor)": (P(u_turb=0.9), None),
    "wind at the turbine x1.17 (lifting line)": (P(u_turb=1.17), None),
    "flow at the turbine turned -4°": (P(turb_dang=-4.0), None),
    "flow at the turbine turned +19°": (P(turb_dang=19.0), None),
    "wind at the turbine x1.17 and turned +19°": (P(u_turb=1.17, turb_dang=19.0), None),
    "kappa 1.2 (auxiliary generators)": (P(kappa=1.2), None),
    "no drag of the stopped rotor": (P(stop_drag=False), None),
    "no induced resistance": (P(sir=False), None),
    "fixed spin ratio k = min(k_max, U_s,max/Ua)": (P(k_opt=False), None),
    "Savonius not widened": (P(ds=1.0), None),
    "assumed turbine (version 2: Cp 0.25, Cp/Ct 0.4)": (P(turb="assumed"), None),
    "turbine turning opposite to the rotor": (P(corot=False), None),
})
res["scenarios"] = {k: summ(annual(d, groups, nu, ng, p, seasons=ss)) for k, (p, ss) in scen.items()}

# ---------------- Savonius conversion: break-even C_t,S (version 5, R2-1) ----------------
# C_t,S at which the rule of claim 1 (T5 claim1) stops losing at sea, and C_t,S below which opening the rotor at sea
# optimally saves at least 1 t/yr, for C_p,S = 0.14, 0.19 and 0.24; and the effective radius of the buckets that gives
# that C_t,S (C_t,S = C_p,S R_S / (lambda r_ef)).


def _bisect(fun, lo, hi, n=30):
    flo = fun(lo)
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if (fun(mid) > 0) == (flo > 0):
            lo, flo = mid, fun(mid)
        else:
            hi = mid
    return 0.5 * (lo + hi)


sav = {}
for cp in (0.14, 0.19, 0.24):
    lit = lambda ct: annual(d, groups, nu, ng, P(cp_s=cp, ct_s=ct))["sea"]["T5 claim1"]
    opt = lambda ct: annual(d, groups, nu, ng, P(cp_s=cp, ct_s=ct))["sea"]["T5"] - 1.0
    ct_lit, ct_opt = _bisect(lit, 0.05, 1.2), _bisect(opt, 0.05, 1.2)
    curve = {f"{ct:.2f}": {k: round(annual(d, groups, nu, ng, P(cp_s=cp, ct_s=ct))["sea"][k], 2)
                          for k in ("T5", "T5 claim1", "T2 claim1", "T3 claim1")}
             for ct in (cp / cf.LAMBDA_S, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0)}
    sav[f"{cp:.2f}"] = dict(ct_break_even_claim1=ct_lit, ct_break_even_optimal=ct_opt,
                            strict_bound=cp / cf.LAMBDA_S, centroid_estimate=cf.KAPPA_S_CENTROID * cp,
                            effective_radius_claim1=cp * cf.R_S / (cf.LAMBDA_S * ct_lit),
                            effective_radius_optimal=cp * cf.R_S / (cf.LAMBDA_S * ct_opt), sea_curve=curve)
res["savonius_break_even"] = sav

# ---------------- one-factor sensitivities (ends of the sampled ranges) ----------------
sens = {}
for k, (lo, hi) in cf.RANGES_V5_MC.items():
    for lab, val in (("min", lo), ("max", hi)):
        sens[f"{k} = {val:g}"] = summ(annual(d, groups, nu, ng, P(**{k: val})))
res["sensitivity"] = sens
json.dump(res, open("route_results.json", "w"), indent=1, ensure_ascii=False)

# ---------------- output ----------------
b = res["base"]
print(f"Reference case (scenario A, April-October): R0 {b['fuel_total']:.0f} t (sea {b['fuel_sea']:.0f}, "
      f"berth {b['fuel_berth']:.0f}); mean load {b['load_mean']:.2f}")
print(f"R1 {b['sea']['R1']:.1f} t ({100 * b['sea']['R1'] / b['fuel_total']:.2f} %): thrust {b['sea']['_R1 thrust']:.1f}, "
      f"spin {b['sea']['_R1 spin']:.1f}, stopped {b['sea']['_R1 stopped']:.1f}, induced {b['sea']['_pen R1']:.1f} "
      f"(own {b['sea']['_pen own']:.1f}, cross {b['sea']['_pen cross']:.1f}); wind on the ship "
      f"{b['sea']['_R0 wind']:.0f} t")
for r in table(base):
    print(f"  {r[0]:4s} sea {r[1]:7.1f}  port {r[2]:5.1f}  total {r[3]:7.1f}  ({r[4]:.2f} %)")
print(f"  comparators: R1t {b['sea']['R1t']:.1f}, R1-30 {b['sea']['R1-30']:.1f}, R1-D5 {b['sea']['R1-D5']:.1f}, "
      f"R1t+X {b['sea']['R1t+X']:.1f}, fifth unit {b['R1_fifth_unit']:.1f} t")
print("  rules:", {k: round(b['sea'][k], 1) for k in ("T5 claim1", "T5 intended", "T5 astern", "T2 claim1",
                                                      "T2 intended", "T3 claim1", "T3 intended")})
print("eps* =", round(res["eps_star"], 3))
print("\nSTEP BY STEP")
for k, v in res["step_by_step"].items():
    print(f"  {k[:70]:70s} R0 {v['R0']:7.0f}  R1 {v['R1']:6.1f} ({v['R1_pct']:.2f} %)  T3x {v['T3x']:5.1f}  T3 {v['T3']:5.1f}  "
          f"T5 {v['T5']:5.1f}  T4 {v['T4']:4.1f}  T1 {v['T1']:6.1f}  R1t {v['R1t']:6.1f}  R1-30 {v['R1-30']:6.1f}")
print("\nBY CROSSING (kW per unit, April-October)")
for row in per:
    print(f"  {row['crossing']:28s} {row['v_kn']:5.2f} kn  {row['Apr-Oct']}")
print("\nSAVONIUS: break-even C_t,S")
for cp, v in res["savonius_break_even"].items():
    print(f"  C_p,S {cp}: claim 1 {v['ct_break_even_claim1']:.3f} (r_ef {v['effective_radius_claim1']:.2f} m), "
          f"optimal {v['ct_break_even_optimal']:.3f} (r_ef {v['effective_radius_optimal']:.2f} m); bound "
          f"{v['strict_bound']:.3f}, centroid {v['centroid_estimate']:.3f}")
print("\nSCENARIOS")
for k, v in res["scenarios"].items():
    print(f"  {k[:58]:58s} R0 {v['R0']:7.0f} R1 {v['R1']:6.1f} ({v['R1_pct']:.2f} %) T3x {v['T3x']:5.1f} T3 {v['T3']:5.1f} "
          f"T5 {v['T5']:5.1f} T4 {v['T4']:5.1f} T1 {v['T1']:6.1f} T2 {v['T2']:6.1f} R1t {v['R1t']:6.1f} R1-30 {v['R1-30']:6.1f}")
print("\nONE-FACTOR SENSITIVITY")
for k, v in sens.items():
    print(f"  {k:22s} R1 {v['R1']:6.1f} T3x {v['T3x']:5.1f} T3 {v['T3']:5.1f} T5 {v['T5']:5.1f} T4 {v['T4']:5.1f} "
          f"T1 {v['T1']:6.1f} R1t {v['R1t']:6.1f} R1-30 {v['R1-30']:6.1f}")
