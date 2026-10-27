#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structural indicators of A1 (section 4.6). Version 4 (22-10-2026): spin limits of the parameter p (k_max, U_s,max),
30 m rotor (R1-30), absolute air draught and heeling moment of the ship.

Bending moment at the base of one unit, in operation, with the same apparent wind over the whole height:
- rotor spinning at k = min(k_max, U_s,max/U_a): resultant q A_F (CL^2 + CD^2)^0.5 applied at H_F/2;
- extension (T3, T3x) at k_X = k D_X/D_F, with the factor f_X on the lift, applied at H_F + H_X/2;
- upper segment of the 30 m rotor (R1-30), with the diameter of the rotor and the factor f30, at H_F + H_X/2;
- VAWT (T1, T2, T4): resultant of its thrust and its side force (DMST model of ../turbine) in the state of largest
  force allowed at that wind, applied at H_F + H_X/2;
- rotor opened as a Savonius (T2, T3, T5): T_S = C_t,S q D_S H_F at H_F/2; its lift is unknown.
The forces of rotor and device are added in magnitude (upper bound). Reference: moment of the rotor alone at the
assumed operating limit, U_a = 20 m/s. The largest apparent wind at which the extension fits within that envelope is
also computed. Survival with the rotor stopped in a storm is not assessed.
Air draught: base of the rotors at Z_DECK = Z_ROT - H_F/2 = 28 m above the waterline (assumed, the same that puts the
rotor centre at 40 m): 52 m at the top of the rotor and 58 m with the extension, the VAWT or the 30 m rotor. The port
profile of Portsmouth gives no air-draught limit (Cruise Britain) and the berth at Bilbao is in the outer harbour, at
Santurtzi, seaward of the bridges over the estuary (Bilbao Port Authority).
Heeling moment: side force of the four units with the apparent wind abeam at U_a = 20 m/s, applied at the height of
their centre relative to the centre of lateral resistance (T/2 below the waterline); the static heel is also given
for an assumed displacement and metacentric height (no public stability data).
"""
import json
import math
import numpy as np
import route_model as rm

UA_DESIGN = 20.0
Z_DECK = rm.Z_ROT - rm.H_F / 2.0            # 28 m: base of the rotors above the waterline (assumed)
DISPL_T = (20000.0, 24000.0)               # t, assumed: L_pp ~ 200 m x 27.8 m x 6.4 m x C_B 0.55-0.65
GM_M = (1.5, 2.5)                          # m, assumed metacentric height


def coeffs(k, p):
    cl, cd = rm.rotor_coeffs(np.array([k]), p)
    return float(cl[0]), float(cd[0])


def k_op(u_a, p):
    return min(p.k_max, p.us_max / u_a)


def moments(u_a, p):
    """Moments at the base (N m) per unit and side force of the rotor (N): rotor alone, extension, VAWT, opened
    Savonius and upper segment of the 30 m rotor."""
    q = 0.5 * rm.RHO * u_a ** 2
    k = k_op(u_a, p)
    cl, cd = coeffs(k, p)
    m_rot = q * rm.A_F * math.hypot(cl, cd) * rm.H_F / 2
    clx, cdx = coeffs(k * rm.D_X / rm.D_F, p)
    m_ext = q * rm.A_X * math.hypot(p.cl_ext * clx, cdx) * (rm.H_F + rm.H_X / 2)
    m_30 = q * rm.A_UP * math.hypot(p.f30 * cl, cd) * (rm.H_F + rm.H_X / 2)
    if p.turb == "dmst":
        _, th, sf = rm.vawt_states(u_a, p)
        f_vawt = float(np.max(np.hypot(th, sf)))
    else:
        f_vawt = q * rm.A_X * p.cp_t / p.cpct_t
    m_vawt = f_vawt * (rm.H_F + rm.H_X / 2)
    m_sav = p.ct_s * q * p.ds * rm.A_F * rm.H_F / 2
    return dict(k=k, rot=m_rot, ext=m_ext, r30=m_30, vawt=m_vawt, sav=m_sav,
                f_rot=q * rm.A_F * math.hypot(cl, cd), f_ext=q * rm.A_X * math.hypot(p.cl_ext * clx, cdx),
                f_30=q * rm.A_UP * math.hypot(p.f30 * cl, cd))


def ext_limit(p):
    m_design = moments(UA_DESIGN, p)["rot"]
    for u_a in np.arange(1.0, UA_DESIGN, 0.01):
        m = moments(u_a, p)
        if m["rot"] + m["ext"] > m_design:
            return round(float(u_a) - 0.01, 2), m_design
    return UA_DESIGN, m_design


def heeling(p, u_a=UA_DESIGN):
    """Heeling moment (MN m) of the units at U_a abeam and static heel (degrees) for the assumed range of displacement
    and GM: rotor alone (R1) and with the extension or the 30 m rotor."""
    m = moments(u_a, p)
    arm_rot = rm.Z_ROT + p.t_draught / 2.0
    arm_top = Z_DECK + rm.H_F + rm.H_X / 2.0 + p.t_draught / 2.0
    # version 5: (a) bound with the resultant of each unit; (b) apparent wind abeam, with the side force of the
    # rotors equal to their drag (C_D q A)
    q = 0.5 * rm.RHO * u_a ** 2
    k = k_op(u_a, p)
    _, cd = coeffs(k, p)
    _, cdx = coeffs(k * rm.D_X / rm.D_F, p)
    lat = dict(rot=q * rm.A_F * cd, ext=q * rm.A_X * cdx, r30=q * rm.A_UP * cd)
    out = {}
    for lab, mom, mom_t in (("R1", m["f_rot"] * arm_rot, lat["rot"] * arm_rot),
                            ("T3x", m["f_rot"] * arm_rot + m["f_ext"] * arm_top, lat["rot"] * arm_rot + lat["ext"] * arm_top),
                            ("R1-30", m["f_rot"] * arm_rot + m["f_30"] * arm_top, lat["rot"] * arm_rot + lat["r30"] * arm_top)):
        mh, mt = p.n_units * mom, p.n_units * mom_t
        heel = [math.degrees(math.asin(min(1.0, mh / (dsp * 1e3 * 9.81 * gm)))) for dsp in DISPL_T for gm in GM_M]
        heel_t = [math.degrees(math.asin(min(1.0, mt / (dsp * 1e3 * 9.81 * gm)))) for dsp in DISPL_T for gm in GM_M]
        out[lab] = dict(M_MNm=mh / 1e6, heel_deg=(min(heel), max(heel)), M_abeam_MNm=mt / 1e6,
                        heel_abeam_deg=(min(heel_t), max(heel_t)))
    return out


def report(p, label):
    k_lim = p.us_max / p.k_max                # apparent wind above which U_s,max governs
    m_k = moments(k_lim, p)
    m20 = moments(UA_DESIGN, p)
    lim, m_design = ext_limit(p)
    m_sav_cut = moments(rm.UA_CUT_S, p)["sav"]
    return dict(case=label, ua_kmax=k_lim,
                M_rotor_kNm_kmax=m_k["rot"] / 1e3, M_rotor_kNm_20=m20["rot"] / 1e3,
                ext_pct_kmax=100 * m_k["ext"] / m_k["rot"], vawt_pct_kmax=100 * m_k["vawt"] / m_k["rot"],
                r30_pct_kmax=100 * m_k["r30"] / m_k["rot"],
                ext_pct_20=100 * m20["ext"] / m20["rot"], vawt_pct_20=100 * m20["vawt"] / m20["rot"],
                r30_pct_20=100 * m20["r30"] / m20["rot"],
                sav_pct_of_design=100 * m_sav_cut / m_design, ua_ext_max=lim,
                heeling=heeling(p))


def air_draft():
    return {"R1": Z_DECK + rm.H_F, "T3x, T3, T1, T2, T4, R1-30": Z_DECK + rm.H_F + rm.H_X}


if __name__ == "__main__":
    import coefficients as cf
    cases = [("A (reference)", rm.P()), ("A+ (k <= 4, U_s <= 46 m/s)", rm.P(**cf.COEF_SCENARIOS["A+"])),
             ("A-VG (fitted to the Viking Grace)", rm.P(**cf.COEF_SCENARIOS["A-VG"])),
             ("AB (version 3)", rm.P(**cf.COEF_SCENARIOS["AB"])), ("B", rm.P(**cf.COEF_SCENARIOS["B"])),
             ("f_X = 1.30", rm.P(cl_ext=1.30)), ("f_X = 0.90", rm.P(cl_ext=0.90)), ("f_X = 0.5 (plate below)", rm.P(cl_ext=0.5)),
             ("Savonius Ct 1.0", rm.P(ct_s=1.0)), ("Savonius Ct at the bound", rm.P(ct_s=cf.KAPPA_S * cf.CENTRAL["cp_s"]))]
    res = [report(p, lab) for lab, p in cases]
    for r in res:
        print(f"{r['case']:26s} M rotor {r['M_rotor_kNm_kmax']:6.0f} kN m (Ua {r['ua_kmax']:.0f} m/s, k_max), "
              f"{r['M_rotor_kNm_20']:6.0f} (20 m/s); extension +{r['ext_pct_kmax']:4.1f} % / +{r['ext_pct_20']:4.1f} %; "
              f"30 m +{r['r30_pct_kmax']:4.1f} % / +{r['r30_pct_20']:4.1f} %; VAWT +{r['vawt_pct_kmax']:4.1f} % / "
              f"+{r['vawt_pct_20']:4.1f} %; Savonius {r['sav_pct_of_design']:5.1f} % of the design moment; extension "
              f"up to Ua {r['ua_ext_max']:.2f} m/s")
        e = r["heeling"]
        print("    heeling moment at 20 m/s (bound with the resultant | abeam): " + "; ".join(
            f"{k} {v['M_MNm']:.1f} MN m ({v['heel_deg'][0]:.1f}-{v['heel_deg'][1]:.1f} degrees) | "
            f"{v['M_abeam_MNm']:.1f} MN m ({v['heel_abeam_deg'][0]:.1f}-{v['heel_abeam_deg'][1]:.1f})"
            for k, v in e.items()))
    print("air draught (m above the waterline):", air_draft())
    json.dump(dict(cases=res, air_draught=air_draft(), z_deck=Z_DECK, displacement_t=DISPL_T, gm_m=GM_M),
              open("moments.json", "w"), indent=1, ensure_ascii=False)
