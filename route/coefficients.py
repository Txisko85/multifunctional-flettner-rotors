#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coefficient base of the route model.

Spinning rotor, CL(k) and CD(k), with k = U_surface / U_apparent:
  Set A  De Marco et al. (2016), Fig. 10, URANS: AR = 6 and d_e/d = 2, the geometry of the case study (digitised from
         the vectors of the PDF). It agrees with the tests of Badalamenti (AR 5.1, d_e/d = 2; Seifert 2012, Figs.
         21-24) and of Ackeret (Seifert 2012, Fig. 27). Model-scale Reynolds number.
  Set B  Bordogna et al. (2019), Figs. 4 and 7: wall-to-wall tests at Re 3.6e5-1e6 (almost 2D), mean of the
         supercritical points, plus the induced drag of a rotor of AR = 6 on a deck (effective AR 12 by reflection in
         the deck): CD_i = CL^2 / (pi * 12).
  The model blends both with a weight w (w = 1: A; w = 0: B).
Stopped rotor: CD0 0.4-0.7 (Bordogna, k = 0, Re >= 3.6e5: 0.44-0.56; Ackeret, low Re: 0.84).
Spin power: friction of the cylinder P = Cf * rho/2 * U_s^3 * pi * D * H, with Cf between 0.002 (Theodorsen-Regier,
  smooth cylinder at Re ~1e7) and 0.007 (Bordogna, test with tip and mechanical losses), plus the disc of the endplate
  (free disc, CM = 0.146 Re^-0.2).
Savonius: Cp 0.14-0.24 in a wind tunnel with semicircular buckets (Akwa et al. 2012, Table 1). Thrust C_t,S between
  its kinematic bound and 1.0 (see below); Chauvin and Benghrib (1989) give 0.14-0.20 at lambda 0.6-1.0 from six
  pressure taps at mid-height; upper bound, the drag of a bluff body of the same area (assumption).
VAWT: double-multiple streamtube model of ../turbine (turbine_map.json) with three factors from its validation
  (dmst_validation.py). f_cp_t 0.7-0.9: the peak Cp measured in the field on the 12 kW H-rotor of Uppsala is 0.68-0.73
  times the model's, and in the wind tunnel of Battisti et al. (2018) 0.92 times. f_ct_t 1.0-1.1: the model gives 89 %
  of the thrust measured in the tunnel near the peak, which includes the shaft and arms of the rig, not modelled.
  f_cs_t 0.6-1.0: the model overestimates the side force of the tunnel by 35 % near the peak (0-80 % over the test).
  Earlier versions (1 and 2): Cp 0.20-0.30 and Cp/Ct 0.3-0.5, assumed (RANGES_V2).
Ship: propulsive efficiency eta_D 0.65-0.75 and SFOC of the main engines 170-190 g/kWh.
Induced resistance (versions 2-4; a sensitivity from version 5): dR = Fy^2 / (pi e q_w (2T)^2), a wing of low aspect
  ratio with effective span 2T (reflection in the free surface). Bickert et al. (MARINE 2025) give dR/R0 with Fy = 2 R0
  at full scale: KVLCC2 at 15.5 kn 6.7 % (3 DOF), KCS at 20 kn 10.0 % (3+1 DOF) and 21.1 % (4 DOF, heel of 15
  degrees). With R0 estimated by ITTC-78 (+-10 %), they are e = 0.51-0.62, 0.43-0.52 and 0.20-0.25 (infer_eh.py).
Side force of the wind on hull and superstructure: F_s = q_s (A_L C_Y,max) sin(beta_s), with A_L C_Y,max 3,000-6,000
  m2 (assumption: A_L ~ 5,000 m2 for 214.5 m length overall and ~25 m mean height above the waterline, C_Y,max
  0.6-1.2).
Version 4 (22-10-2026), after the simulated review:
  - Wind resistance of the ship without rotors (R0), explicit: longitudinal X_s = q_s (A_T C_X,max) cos(beta_s), with
    A_T C_X,max 300-600 m2 (assumption: frontal area ~ 27.8 m beam x ~ 27 m height, C_X,max 0.4-0.8, as the ferries of
    Blendermann 1994), and induced resistance from its own side force.
  - Rotor coefficient sets as epistemic scenarios (COEF_SCENARIOS), not a random variable: the weight w is no longer
    sampled. The reference scenario is A: with the same code, on the route of the Viking Grace
    (viking_grace_check.py), A gives 141 kW of propulsion power equivalent to the thrust of the 24 m x 4 m rotor,
    against the 207-282 kW reported (Norsepower 2019); the A/B blend of version 3 gives 86 kW; B, 37 kW; B', 22 kW.
Version 5 (27-10-2026), after the second simulated review: strict Savonius bound, f_X and f30 checked against URANS,
  scenario A-VG, sampled hull and windage parameters of the sway-yaw balance (see below).
"""
import math
import numpy as np

K_A = [0.0, 1.0, 1.5, 2.0, 3.0]
CL_A = [0.0, 1.807, 3.298, 5.451, 9.729]
CD_A = [None, 0.891, 0.869, 1.465, 3.977]          # k = 0: CD0 of the stopped rotor

K_B = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
CL_B = [0.0, 1.12, 2.20, 3.05, 3.85, 5.04, 5.58]
CD_B2D = [0.51, 0.77, 1.03, 1.33, 1.61, 1.87, 2.16]
AR_EFF_B = 12.0

# Version 4: extension beyond k = 3, only for the sensitivity with k_max = 4 (and the Viking Grace check). Set B:
# points of Bordogna et al. (2019) at Re = 3.6e5 (CD, mean of the two series of Fig. 7). Set A: De Marco et al. (2016)
# stop at k = 3; extended with the relative shape of Ackeret's tests with 120 mm plates (Seifert 2012, Fig. 27), which
# saturates: CL x1.08 / 1.10 / 1.12 and CD x1.23 / 1.34 / 1.45 at k = 3.5 / 4.0 / 4.5 relative to k = 3.
K_B_EXT = K_B + [3.65, 4.13, 4.60]
CL_B_EXT = CL_B + [6.50, 7.31, 8.46]
CD_B2D_EXT = CD_B2D + [2.37, 2.50, 2.68]
K_A_EXT = K_A + [3.5, 4.0, 4.5]
CL_A_EXT = CL_A + [round(9.729 * r, 3) for r in (1.0822, 1.1032, 1.1244)]
CD_A_EXT = CD_A + [round(3.977 * r, 3) for r in (1.2337, 1.3405, 1.4472)]

# Three-dimensional correction of set B (version 4, lifting_line.py): lift factor eta_b and effective aspect ratio of
# the induced drag. By default, no correction (1.0 and 12, as up to version 3); scenario B' with the wing analogy and
# the deck as a symmetry plane (eta 0.83, AR_eff 12) and B'' without reflection in the deck (eta 0.725, AR_eff 6).
ETA_B_SCEN = {"B": (1.0, 12.0), "B'": (0.83, 12.0), "B''": (0.725, 6.0)}

# Central value = midpoint of each range (the uncertainty analysis samples the ranges uniformly).
# Savonius geometry: bucket centroid r_c = (d - e)/2 = 1.6 m, R = D_S/2 = 3.6 m, tip-speed ratio lambda = 0.8.
R_S, RC_S, LAMBDA_S = 3.6, 1.6, 0.8
# Version 5 (second simulated review, R2-1): the strict bound uses the largest speed of the buckets, that of the outer
# edge (omega R_S): C_t,S >= C_p,S / lambda = 1.25 C_p,S. Evaluated at the centroid of the bucket projection (version
# 4) it gave 2.81 C_p,S, an estimate with an assumed effective radius, not a bound.
KAPPA_S = 1.0 / LAMBDA_S
KAPPA_S_CENTROID = R_S / (LAMBDA_S * RC_S)          # 2.8125: estimate of version 4
CT_S_MAX = 1.0

RANGES = dict(w=(0.0, 1.0), cd_mult=(0.8, 1.2), cd0=(0.4, 0.7), cf=(0.002, 0.007), cp_s=(0.14, 0.24),
              ct_s=(0.4, 1.0), f_cp_t=(0.7, 0.9), f_ct_t=(1.0, 1.1), f_cs_t=(0.6, 1.0), cl_ext=(0.6, 1.0),
              eta_d=(0.65, 0.75), sfoc_me=(170.0, 190.0), e_span=(0.4, 0.8), ayc=(3000.0, 6000.0))
CENTRAL = {k: 0.5 * (a + b) for k, (a, b) in RANGES.items()}
# Version 4: frontal wind resistance of the ship, C_t,S bounded by its kinematic bound, and lift factors of the
# extension and of the 30 m rotor from the lifting-line model (lifting_line.py).
# Version 5 (second simulated review, R2-3): f_X and f30 checked against the aspect-ratio dependence of C_L in De Marco
# et al. (2016, URANS with endplates d_e/d = 2 and 3): empirical f30 0.89-1.13 and f_X = f30 x (1.01-1.18), the ratio
# from the lifting-line model with both section models (lifting_line_validation.py): f_X 0.90-1.30.
RANGES_V4 = dict(RANGES, ct_s=(KAPPA_S * RANGES["cp_s"][0], CT_S_MAX), axc=(300.0, 600.0), cl_ext=(0.90, 1.30),
                 f30=(0.90, 1.15))
CENTRAL_V4 = dict(CENTRAL, ct_s=0.5 * (KAPPA_S * CENTRAL["cp_s"] + CT_S_MAX), axc=450.0, cl_ext=1.10, f30=1.025,
                  w=1.0)
RANGES_V4_ORIG = dict(RANGES, ct_s=(KAPPA_S_CENTROID * RANGES["cp_s"][0], CT_S_MAX), axc=(300.0, 600.0),
                      cl_ext=(0.95, 1.25), f30=(0.95, 1.15))           # version 4 (26-10-2026)
CENTRAL_V4_ORIG = dict(CENTRAL, ct_s=0.5 * (KAPPA_S_CENTROID * CENTRAL["cp_s"] + CT_S_MAX), axc=450.0, cl_ext=1.10,
                       f30=1.05, w=1.0)
# Epistemic scenarios of the rotor coefficients (version 4). A: reference (the closest to the full-scale figure); A+:
# A extended up to k = 4 with U_s <= 46 m/s; AB: the blend of version 3; B: Bordogna with AR_eff 12; B': B corrected
# with the lifting line (wing analogy, deck as a symmetry plane); B'': without reflection.
# Version 5 (second simulated review, R1-1 and R2-2): scenario A-VG, A with C_L multiplied so that the thrust of the
# model on the route of the Viking Grace, with the ERA5 cells over water, equals 207-282 kW (centre 244.5 kW);
# viking_grace_check.py checks these values.
CL_MULT_VG = {"207": 1.141, "244.5": 1.25, "282": 1.359}
COEF_SCENARIOS = {"A": dict(w=1.0), "A+": dict(w=1.0, k_max=4.0, us_max=46.0), "AB": dict(w=0.5),
                  "B": dict(w=0.0), "B'": dict(w=0.0, eta_b=ETA_B_SCEN["B'"][0], ar_b=ETA_B_SCEN["B'"][1]),
                  "B''": dict(w=0.0, eta_b=ETA_B_SCEN["B''"][0], ar_b=ETA_B_SCEN["B''"][1]),
                  "A-VG": dict(w=1.0, cl_mult=CL_MULT_VG["244.5"]),
                  "A-VG 207": dict(w=1.0, cl_mult=CL_MULT_VG["207"]),
                  "A-VG 282": dict(w=1.0, cl_mult=CL_MULT_VG["282"])}
# ranges sampled in version 4 (scenario fixed): all but w
RANGES_V4_MC = {k: v for k, v in RANGES_V4.items() if k != "w"}
# Version 5: the induced resistance comes from the sway-yaw balance (sway_yaw.py); instead of e_span, the hull term c_H
# (0-1; the KVLCC2 of Yasukawa and Yoshimura 2015, X'_vv = -0.040 and Y'_v = -0.315, gives c_H = 2 x 0.040 / 0.315 =
# 0.25), the aft-body coefficients C_pa and C_wa and the centre of the side force of the wind on the ship x_s (fraction
# of Lpp forward of midships) are sampled.
RANGES_V5_MC = {k: v for k, v in RANGES_V4.items() if k not in ("w", "e_span")}
RANGES_V5_MC.update(c_h=(0.0, 1.0), cpa=(0.65, 0.75), cwa=(0.80, 0.90), x_s=(0.05, 0.10))


def ct_s_range(cp_s):
    """Range of C_t,S compatible with the kinematic bound for a given C_p,S."""
    return KAPPA_S * cp_s, CT_S_MAX


# ranges of earlier versions, with the assumed turbine (to reproduce them)
RANGES_V2 = dict(w=(0.0, 1.0), cd_mult=(0.8, 1.2), cd0=(0.4, 0.7), cf=(0.002, 0.007), cp_s=(0.14, 0.24),
                 ct_s=(0.4, 1.0), cp_t=(0.20, 0.30), cpct_t=(0.3, 0.5), cl_ext=(0.6, 1.0),
                 eta_d=(0.65, 0.75), sfoc_me=(170.0, 190.0), e_span=(0.4, 0.8), ayc=(3000.0, 6000.0))   # 19-10-2026
RANGES_V1 = {k: v for k, v in RANGES_V2.items() if k not in ("e_span", "ayc")}                       # 18-10-2026


def _interp(k, ks, vs):
    return float(np.interp(k, ks, vs))


def rotor_cl_cd(k, w=0.5, cd0=0.55, cd_mult=1.0):
    """CL and CD of the spinning rotor at spin ratio k (0 <= k <= 3)."""
    k = min(max(k, 0.0), 3.0)
    cd_a = [cd0] + CD_A[1:]
    cl_a, cd_a_k = _interp(k, K_A, CL_A), _interp(k, K_A, cd_a)
    cl_b = _interp(k, K_B, CL_B)
    cd_b = _interp(k, K_B, [cd0] + CD_B2D[1:]) + cl_b ** 2 / (math.pi * AR_EFF_B)
    cl = w * cl_a + (1.0 - w) * cl_b
    cd = w * cd_a_k + (1.0 - w) * cd_b
    return cl, cd0 + (cd - cd0) * cd_mult if k > 0 else cd0


if __name__ == "__main__":
    print(" k    CL_A  CD_A   CL_B  CD_B   blend CL  CD   L/D")
    for k in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        a = rotor_cl_cd(k, w=1.0); b = rotor_cl_cd(k, w=0.0); c = rotor_cl_cd(k)
        print(f"{k:3.1f}  {a[0]:5.2f} {a[1]:5.2f}  {b[0]:5.2f} {b[1]:5.2f}   {c[0]:5.2f} {c[1]:5.2f} {c[0]/c[1]:5.2f}")
