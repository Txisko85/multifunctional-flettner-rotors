#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Route model (level 1): ro-pax Galicia (Brittany Ferries), Portsmouth-Bilbao.

Version 5 (27-10-2026, after the second simulated review of Applied Energy):
1. Service: the base case is April-October (BASE_SEASONS), with two rotations a week; November and the whole year
   are scenarios. Wind histograms conditioned on the timetable (era5_conditioned.py): by crossing and month, only the
   hours and points the ship sails through, and at berth only the berth hours.
2. Induced resistance from a sway-yaw balance with rudder (sway_yaw.py), a quadratic form of the total side force
   and of its yaw moment about midships; the wing form of version 4 (dR = c Fy^2 with e_span) is kept as a
   sensitivity (sir_model = "wing").
3. Strict Savonius bound (coefficients.py), scenario A-VG (cl_mult), comparator R1-D5 (wider rotor) and R1t+X
   (tiltable rotor with extension), common structural envelope (m_cap), turbine inflow sensitivities (turb_dang).

Version 4 (22-10-2026, after the first simulated review; reproduced with P(v4=True) and era5_route.json):
1. Service: April-November with two rotations a week; the whole year as a scenario.
2. The ship without rotors (R0) carries the wind resistance of hull and superstructure explicitly: longitudinal,
   X_s = q_s (A_T C_X,max) cos(beta_s), and induced by its own side force. R0, R1 and the topologies use the same model.
3. Power-speed curve: power law (exponent 3; sensitivity 4) through 85 % MCR at 22 kn, multiplied by (1 + sea margin)
   calibrated against EU MRV at the same mean speed (mrv_calibration.py; SEA_MARGIN).
4. Specific fuel consumption with engine load (Fourth IMO GHG Study, 2020: SFOC = SFOC_base (0.455 L^2 - 0.710 L +
   1.280)). The load includes the power of the shaft generators. Savings, small compared with the ship's power, are
   converted with the marginal consumption of each crossing, d(SFOC P)/dP.
5. Non-patented comparators: tiltable rotor (R1t: tilted, it has no drag and no side force), 30 m rotor (R1-30: the
   same rotor with 6 m more of the same diameter, always running or stopped) and five units (n_units = 5).
6. Savonius: side force of the opened rotor as a sensitivity (cy_s).
7. Open-sea distance 539.8 nm (208.4 + 331.4).
8. Ship particulars from primary sources: draught 6.4 m. Hotel load at berth 2,000 kW: EU MRV 2021-2022 gives at least
   1,700-1,850 kW.
9. Rotor coefficients as epistemic scenarios (coefficients.COEF_SCENARIOS); the reference is A, the closest to the
   full-scale figure of the Viking Grace (viking_grace_check.py).
10. Sensitivities: interference between units in tandem (int_frac, int_beta), wind at the rotor (f_u), wind at the
   turbine (u_turb), k > 3 (k_max, us_max) and the term 'rot_run' for the start and stop cost of R1 in the hourly
   model.
Version 2 (19-10-2026; version 1 of 18-10-2026 is reproduced with P(v1=True)): timetable in real durations (Spain is
one hour ahead of the UK); induced resistance of the total side force; the interaction factor eps scales the
aerodynamic force of the rotor (both components) with the VAWT deployed. Version 3 (21-10-2026): turbine with the DMST
model (../turbine).

Public data:
- Ship (primary sources): 214.5 m length overall, 27.8 m beam, 6.4 m draught, 2 x 12,600 kW and 22 kn design speed
  (Brittany Ferries, press release of 24-08-2020; Stena RoRo, 03-09-2020); two shaft lines with one engine and one
  scrubber each. Shaft generators assumed to supply the electricity at sea.
- 2026 Portsmouth-Bilbao timetable: two rotations a week from April to November (Brittany Ferries newsroom).
- Wind: hourly ERA5 2015-2024 at 10 and 100 m (Open-Meteo, models=era5), at 10 route points and the two ports;
  extrapolated with a power law to 40 m (rotor centre) and 55 m (VAWT or extension).
Coefficients: coefficients.py. Assumptions: efficiencies, power-speed curve, SFOC, hotel loads, deck height, wind
areas of the ship and turbine factors.
"""
import json
import math
import os
import numpy as np
import coefficients as cf

RHO, NU_AIR, RHO_W = 1.225, 1.48e-5, 1025.0
KN = 0.514444

# ---------------- Rotor and devices ----------------
D_F, H_F, R_E = 4.0, 24.0, 4.0           # R_E: endplate radius (De = 2D)
D_X, H_X = 3.4, 6.0                      # VAWT or extension
A_F, A_X = D_F * H_F, D_X * H_X
A_UP = D_F * H_X                         # upper segment of the 30 m rotor (R1-30)
CL3, CD0_ROT, SR_NOM, US_MAX, CF_SKIN = 8.0, 2.0, 3.0, 36.0, 0.0025   # level 0 (only P(level0=True))
CP_T, CPCT = 0.30, 0.40                  # level-0 turbine (P(level0=True))
CT_T = CP_T / CPCT
CP_T_V2, CPCT_V2 = 0.25, 0.40            # assumed turbine of versions 1 and 2 (centre of their ranges)
UA_CUT_T = 20.0
# Turbine with the DMST model (21-10-2026): maps of ../turbine/turbine_map.json, validation factors f_cp_t, f_ct_t and
# f_cs_t (coefficients.py), maximum tip speed and rated and cut-in winds of the 12 kW H-rotor of Uppsala (Apelfröjd et
# al. 2016: 40 m/s and 12 m/s; cut-in 3 m/s, assumed).
TURBINE_MAP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "turbine", "turbine_map.json")
V_TIP_T, U_RATED_T, U_IN_T = 40.0, 12.0, 3.0
DS_FACTOR, CP_S, CT_S, CD_STOP, UA_CUT_S = 1.8, 0.15, 0.9, 0.6, 20.0
ETA_M, ETA_MOT, ETA_CONV, ETA_GEN = 0.97, 0.95, 0.97, 0.94
N_UNITS = 4

# ---------------- Ship (public) and assumptions ----------------
P_MCR, V_DESIGN = 25200.0, 22.0
LOAD_DESIGN, N_PROP, ETA_D = 0.85, 3.0, 0.70       # assumptions
SFOC_ME, SFOC_AE, ETA_PTO = 180.0, 200.0, 0.95       # assumptions
P_HOTEL_SEA, P_HOTEL_BERTH = 2000.0, 2000.0          # kW, assumptions. Berth: EU MRV 2021-2022 requires >= 1,700-1,850 kW
P_HOTEL_BERTH_V3 = 1500.0                            # kW, up to version 3 (below the EU MRV bound)
T_DRAUGHT = 6.4                                      # m, draught (Brittany Ferries, Stena RoRo)
T_DRAUGHT_V3 = 6.7                                   # m, draught from a secondary source (up to version 3)
Z_SHIP, Z_ROT, ALPHA_SHIP = 20.0, 40.0, 0.08         # centre of the lateral area, rotor centre, mean exponent
CF_FUEL = 3.15    # t CO2/t fuel: CO2/fuel ratio of the Galicia in EU MRV 2021-2025
                  # (3.1518; annual 3.13-3.20). Annex I of Regulation (EU) 2015/757: HFO 3.114, LFO 3.151, MGO 3.206.
# Sea margin (waves, fouling, off-design operation) on the calm-water power, calibrated against EU MRV 2021-2023 at the
# same mean speed (mrv_calibration.py), by exponent of the curve.
# Version 5: recalibrated with the timetable-conditioned histograms and the induced resistance of the sway-yaw balance
# (version 4: 0.073 and 0.348).
SEA_MARGIN = {3.0: 0.077, 4.0: 0.353}
SEA_MARGIN_V4 = {3.0: 0.073, 4.0: 0.348}

# ---------------- Route and timetable (public) ----------------
NM_OPEN = {"C": 208.4, "B": 331.4}                   # open-sea nm: Channel and Bay of Biscay
NM_OPEN_TOTAL = sum(NM_OPEN.values())                # 539.8 nm
NM_APPROACH, V_APPROACH = 32.0, 10.0                 # port approaches (assumed 10 kn)
CROSSINGS = [("Portsmouth-Bilbao, Thursday", "S", 36.0), ("Portsmouth-Bilbao, Sunday", "S", 33.0),
             ("Bilbao-Portsmouth, Tuesday", "N", 32.75), ("Bilbao-Portsmouth, Saturday", "N", 29.5)]
CROSSINGS_V1 = [("Portsmouth-Bilbao, Thursday", "S", 37.0), ("Portsmouth-Bilbao, Sunday", "S", 34.0),
                ("Bilbao-Portsmouth, Tuesday", "N", 31.75), ("Bilbao-Portsmouth, Saturday", "N", 28.5)]
PORT_HOURS_WEEK = {"BILBAO": 10.0, "PORTSMOUTH": 26.75}
# Version 5: season of two rotations a week from April to October. The published 2026 timetable ends on Sunday 1
# November 2026, and in 2025 the winter frequency (one rotation a week) started in November (Brittany Ferries, press
# release); November becomes a sensitivity. Seasons: 0 April-October, 1 December-March, 2 November. Up to version 4:
# 0 April-November (244 days) and 1 December-March (DAYS_V4).
MONTH_DAYS = {1: 31.0, 2: 28.25, 3: 31.0, 4: 30.0, 5: 31.0, 6: 30.0, 7: 31.0, 8: 31.0, 9: 30.0, 10: 31.0, 11: 30.0,
              12: 31.0}
SEASON_MONTHS = {0: (4, 5, 6, 7, 8, 9, 10), 1: (12, 1, 2, 3), 2: (11,)}
DAYS = {s: sum(MONTH_DAYS[m] for m in ms) for s, ms in SEASON_MONTHS.items()}   # 214; 121.25; 30
DAYS_V4 = {0: 244.0, 1: 121.25}
BASE_SEASONS = (0,)                                  # version 5: April-October (version 4: April-November)
ALL_SEASONS = (0, 1, 2)                              # whole year
NOV_SEASONS = (0, 2)                                 # April-November (sensitivity)
WIND_FILE = "era5_route_conditioned.json"            # version 5: histograms conditioned on the timetable


def season_of_month(m):
    return next(s for s, ms in SEASON_MONTHS.items() if m in ms)


def apparent(u_w, gamma_deg, v_s):
    g = math.radians(gamma_deg)
    ux, uy = u_w * math.cos(g) + v_s, u_w * math.sin(g)
    return math.hypot(ux, uy), math.atan2(uy, ux)


def spin_power(u_a, k, cf_skin=CF_SKIN, height=H_F, d=D_F):
    """Spin power (W) of the rotor: friction of the cylinder and disc of the endplate (diameter 2d); k scalar or
    vector."""
    u_s = np.asarray(k, float) * u_a
    r_e = d                                              # endplate De = 2D
    p_cyl = 0.5 * RHO * cf_skin * u_s ** 3 * math.pi * d * height
    w = u_s / (d / 2)
    re_d = np.maximum(w * r_e ** 2 / NU_AIR, 1.0)
    p_disc = 0.5 * 0.146 * re_d ** -0.2 * RHO * w ** 3 * r_e ** 5
    return p_cyl + p_disc


def sfoc_factor(load):
    """SFOC/SFOC_base of the main engine at load 'load' (fraction of MCR): Fourth IMO GHG Study (2020), after Jalkanen
    et al. (2012)."""
    return 0.455 * load ** 2 - 0.710 * load + 1.280


def sfoc_marginal_factor(load):
    """d(SFOC P)/dP / SFOC_base: consumption of one kW more or less at load 'load'."""
    return 1.365 * load ** 2 - 1.420 * load + 1.280


K_GRID = np.round(np.arange(0.25, 3.0001, 0.25), 2)


class P:
    """Parameters varied in the sensitivities and the uncertainty analysis."""
    def __init__(self, **kw):
        self.v1 = False        # True: version of 18-10-2026 (uncorrected timetable, no induced resistance)
        self.v3 = False        # True: version of 21-10-2026 (without the changes of version 4)
        self.v4 = False        # True: version 4 (induced resistance of a low-aspect-ratio wing with e_span)
        self.level0 = False    # True: level-0 coefficients and fixed spin ratio
        self.cd = CD0_ROT      # level 0 only
        self.cl_ext = 1.0
        self.kappa = 1.0 / ETA_PTO
        self.ds = DS_FACTOR
        self.eps = 0.0
        self.ua_ext_max = None
        self.dv = 0.0          # change of ship speed (kn)
        self.stop_drag = True  # drag of the stopped rotor in R1
        self.k_opt = True      # optimal spin ratio in each state
        self.tz_fix = True     # timetable in real durations
        self.sir = True        # induced resistance of the side force
        self.windage = True    # side force of the wind on the ship (cross term)
        self.eta_d = ETA_D
        self.sfoc_me = SFOC_ME
        self.turb = "dmst"     # "dmst": turbine model (21-10-2026); "assumed": assumed Cp and Cp/Ct
        self.cp_t, self.cpct_t = CP_T_V2, CPCT_V2      # only with turb = "assumed"
        self.v_tip_t = V_TIP_T
        self.corot = True      # the turbine turns with the rotor when it assists it (its side force adds to the Magnus)
        # version 4
        self.n_units = N_UNITS
        self.r0_wind = True    # wind resistance of the ship explicit in R0 (and in every case)
        self.sfoc_curve = True # specific consumption with load
        self.n_prop = N_PROP   # exponent of the power-speed curve
        self.sea_margin = None # None: the value calibrated against EU MRV for n_prop (SEA_MARGIN)
        self.cy_s = 0.0        # side-force coefficient of the opened Savonius (sensitivity)
        self.f30 = 1.0         # lift factor of the upper segment of the 30 m rotor (lifting line)
        self.nm_open = NM_OPEN_TOTAL
        self.eta_b, self.ar_b = 1.0, cf.AR_EFF_B   # 3D correction of set B (scenarios B', B'')
        self.k_max, self.us_max = SR_NOM, US_MAX    # limits of the spin ratio and of the surface speed
        self.u_turb = 1.0      # factor on the wind seen by the turbine (flow induced by the rotor)
        self.t_draught = T_DRAUGHT
        self.hotel_berth = P_HOTEL_BERTH
        self.f_u = 1.0         # factor on the wind at the rotor (acceleration over the superstructure; validation)
        self.r1_only = False   # True: bin_terms computes only the conventional rotor (Viking Grace check)
        self.int_frac, self.int_beta = 0.0, 30.0   # interference in tandem (sensitivity): fraction and angle (degrees)
        self.cl_mult = 1.0     # version 5: multiplier of C_L (scenario fitted to the Viking Grace)
        # version 5: induced resistance from the sway-yaw balance with rudder (sway_yaw.py) as a quadratic form of the
        # total side force and of its yaw moment; "wing": dR = c Fy^2 with e_span (version 4)
        self.sir_model = "mmg"
        self.x_s = 0.075       # centre of the side force of the wind on the ship (fraction of Lpp forward)
        self.x_u = 0.0         # centroid of the units (fraction of Lpp forward; declared: at midships)
        self.c_h, self.cpa, self.cwa = 0.5, 0.70, 0.85   # hull term and aft-body form (sway_yaw.py, central values)
        self.m_cap = None      # common structural envelope (m_cap_value)
        self.d_alt = 5.0       # diameter of the wider comparator (R1-D5), same height and De = 2D
        self.turb_dang = 0.0   # rotation (degrees) of the flow reaching the turbine towards the lift side
        self.__dict__.update(cf.CENTRAL_V4)
        self.__dict__.update(kw)
        if self.v1 or self.v3 or self.level0 or self.v4:
            self.sir_model = kw.get("sir_model", "wing")
        if self.v4:            # values of version 4 (26-10-2026) that change in version 5
            for k_, v_ in cf.CENTRAL_V4_ORIG.items():
                if k_ not in kw:
                    setattr(self, k_, v_)
            if "sea_margin" not in kw:
                self.sea_margin = SEA_MARGIN_V4.get(float(self.n_prop), 0.0)
        if self.v1 or self.v3 or self.level0:
            self.r0_wind, self.sfoc_curve, self.sea_margin, self.nm_open = False, False, 0.0, 540.0
            self.ct_s, self.cl_ext, self.t_draught = cf.CENTRAL["ct_s"], cf.CENTRAL["cl_ext"], T_DRAUGHT_V3
            self.hotel_berth, self.w = P_HOTEL_BERTH_V3, cf.CENTRAL["w"]
            self.__dict__.update({k: v for k, v in kw.items()
                                  if k in ("ct_s", "cl_ext", "sea_margin", "r0_wind", "sfoc_curve", "t_draught",
                                           "hotel_berth", "w")})
        if self.v1:
            self.tz_fix, self.sir, self.turb = False, False, "assumed"
            self.__dict__.update({k: v for k, v in kw.items() if k not in ("v1", "tz_fix", "sir")})
        if self.level0:
            self.k_opt, self.turb = False, "assumed"
            self.cd0, self.cf, self.cp_s, self.ct_s, self.cp_t, self.cpct_t = CD_STOP, CF_SKIN, CP_S, CT_S, CP_T, CPCT
            self.cl_ext, self.eta_d, self.sfoc_me = 1.0, ETA_D, SFOC_ME
            self.__dict__.update({k: v for k, v in kw.items() if k != "level0"})
        if self.sea_margin is None:
            self.sea_margin = SEA_MARGIN.get(float(self.n_prop), 0.0)


def crossings(p):
    return CROSSINGS if p.tz_fix else CROSSINGS_V1


def rotor_coeffs(k, p):
    """CL and CD (vectors) of the rotor at spin ratios k."""
    k = np.asarray(k, float)
    if p.level0:
        return CL3 * np.minimum(k, SR_NOM) / SR_NOM, np.full_like(k, p.cd)
    if p.k_max > 3.0:      # version 4, sensitivity: data extended beyond k = 3
        kc = np.clip(k, 0.0, p.k_max)
        ka, cla, cda, kb, clb, cdb = cf.K_A_EXT, cf.CL_A_EXT, cf.CD_A_EXT, cf.K_B_EXT, cf.CL_B_EXT, cf.CD_B2D_EXT
    else:
        kc = np.clip(k, 0.0, 3.0)
        ka, cla, cda, kb, clb, cdb = cf.K_A, cf.CL_A, cf.CD_A, cf.K_B, cf.CL_B, cf.CD_B2D
    cl_a = np.interp(kc, ka, cla)
    cd_a = np.interp(kc, ka, [p.cd0] + cda[1:])
    cl_b = p.eta_b * np.interp(kc, kb, clb)
    cd_b = np.interp(kc, kb, [p.cd0] + cdb[1:]) + cl_b ** 2 / (math.pi * p.ar_b)
    cl = (p.w * cl_a + (1.0 - p.w) * cl_b) * p.cl_mult
    cd = p.w * cd_a + (1.0 - p.w) * cd_b
    return cl, p.cd0 + (cd - p.cd0) * p.cd_mult


_QC = {}
LPP = 200.0


def _quad(v_s, p):
    """Coefficients (a, b, c) of dR = a Y^2 + b Y N + c N^2 at speed v_s (m/s), with N in N m."""
    key = (p.c_h, p.cpa, p.cwa)
    if key not in _QC:
        import sway_yaw
        prm = dict(sway_yaw.CENTRAL, c_h=p.c_h, cpa=p.cpa, cwa=p.cwa)
        _QC[key] = sway_yaw.quad_coeffs(1.0 / KN, prm)       # at 1 m/s; they scale with 1/U^2
    a, b, c = _QC[key]
    return a / v_s ** 2, b / v_s ** 2, c / v_s ** 2


def induced(y, n_m, v_s, p):
    """Induced resistance (N) of a total side force y (N) with yaw moment n_m (N m)."""
    if not p.sir or v_s <= 0.0:
        return 0.0 * np.asarray(y, float)
    if p.sir_model == "mmg":
        a, b, c = _quad(v_s, p)
        return a * y * y + b * y * n_m + c * n_m * n_m
    return sir_coef(v_s, p) * y * y


def m_cap_value(p):
    """Largest base moment (N m) of the common envelope: None (no limit), a number, 'ref' (that of the reference rotor
    at U_a = 20 m/s with k = min(k_max, U_s,max/U_a)) or '30m' (that of the 30 m rotor in the same state, with a
    foundation reinforced for it)."""
    cap = getattr(p, "m_cap", None)
    if cap is None or isinstance(cap, (int, float)):
        return cap
    key = ("_mcap", cap)
    if key not in p.__dict__:
        u_a = 20.0
        k = np.array([min(p.k_max, p.us_max / u_a)])
        cl, cd = rotor_coeffs(k, p)
        q = 0.5 * RHO * u_a ** 2
        m = float(q * A_F * np.hypot(cl, cd)[0] * H_F / 2.0)
        if cap == "30m":
            m += float(q * A_UP * np.hypot(p.f30 * cl, cd)[0] * (H_F + H_X / 2.0))
        p.__dict__[key] = m
    return p.__dict__[key]


def sir_coef(v_s, p):
    """Coefficient c (1/N) of the induced resistance, dR = c Fy^2, at speed v_s (m/s) (wing form)."""
    if not p.sir or v_s <= 0.0:
        return 0.0
    return 1.0 / (math.pi * p.e_span * 0.5 * RHO_W * v_s ** 2 * (2.0 * p.t_draught) ** 2)


def ship_wind(u_rot, gamma, v_s, p):
    """Wind on hull and superstructure (at Z_SHIP): longitudinal resistance (N, positive aft) and side force (N, to
    leeward)."""
    u_sh, b_sh = apparent(u_rot * (Z_SHIP / Z_ROT) ** ALPHA_SHIP, gamma, v_s)
    q_sh = 0.5 * RHO * u_sh ** 2
    return q_sh * getattr(p, "axc", 0.0) * math.cos(b_sh), q_sh * p.ayc * math.sin(b_sh)


def ship_side_force(u_rot, gamma, v_s, p):
    """Side force (N, to leeward) of the wind on hull and superstructure (cross term)."""
    if not (p.sir and p.windage):
        return 0.0
    return ship_wind(u_rot, gamma, v_s, p)[1]


class _TurbineMap:
    """Cp, Ct and |Cs| of the turbine (DMST model, without factors) as functions of lambda and of the wind U."""

    def __init__(self, path=TURBINE_MAP):
        d = json.load(open(path, encoding="utf-8"))
        self.lam, self.u = np.array(d["lam"], float), np.array(d["U"], float)
        self.cp, self.ct, self.cs = (np.array(d[k], float) for k in ("cp", "ct", "cs"))
        self.ct_park = np.array(d["ct_parked"], float)
        self._cache = {}

    LAM_FINE = np.round(np.arange(0.5, 5.0001, 0.05), 2)

    def row(self, u):
        """Cp, Ct and |Cs| on the fine lambda grid (linear interpolation of the map) and parked Ct, at wind u."""
        key = round(float(u), 3)
        if key not in self._cache:
            uc = min(max(float(u), self.u[0]), self.u[-1])
            j = int(np.clip(np.searchsorted(self.u, uc) - 1, 0, len(self.u) - 2))
            t = (uc - self.u[j]) / (self.u[j + 1] - self.u[j])
            rows = [a[j] * (1 - t) + a[j + 1] * t for a in (self.cp, self.ct, self.cs)]
            self._cache[key] = tuple(np.interp(self.LAM_FINE, self.lam, r) for r in rows) + \
                (float(self.ct_park[j] * (1 - t) + self.ct_park[j + 1] * t),)
        return self._cache[key]


_TM = None


def turbine_map():
    global _TM
    if _TM is None:
        _TM = _TurbineMap()
    return _TM


def rated_power_t(p):
    """Rated power (W, shaft) of the turbine: the largest at U_RATED_T within the allowed tip speed."""
    mt = turbine_map()
    cp, _, _, _ = mt.row(U_RATED_T)
    ok = mt.LAM_FINE <= p.v_tip_t / U_RATED_T + 1e-9
    return p.f_cp_t * float(cp[ok].max()) * 0.5 * RHO * U_RATED_T ** 3 * A_X


def vawt_states(u, p):
    """Possible states of the turbine with wind u (m/s) at its height: shaft power (W), thrust along the wind (N) and
    magnitude of the side force (N). The first is the parked turbine; the others, running at each lambda of the fine
    grid (steps of 0.05) allowed by the maximum tip speed, the rated power and the cut-in and cut-out winds."""
    mt = turbine_map()
    q = 0.5 * RHO * u * u
    cp, ct, cs, ct_park = mt.row(u)
    pw, th, sf = [0.0], [ct_park * q * A_X], [0.0]
    if U_IN_T <= u <= UA_CUT_T:
        lam = mt.LAM_FINE
        l_max = min(p.v_tip_t / u, lam[-1])
        p_r = rated_power_t(p)
        pm = p.f_cp_t * cp * q * u * A_X
        # candidates: the grid up to l_max, l_max and the crossings of the rated power (continuous lambda control)
        cand = list(lam[lam <= l_max + 1e-9]) + [l_max]
        for i in range(len(lam) - 1):
            if lam[i] < l_max and (pm[i] - p_r) * (pm[i + 1] - p_r) < 0.0:
                lc = lam[i] + (p_r - pm[i]) / (pm[i + 1] - pm[i]) * (lam[i + 1] - lam[i])
                if lc <= l_max:
                    cand.append(lc)
        cand = np.array(cand)
        pmc = p.f_cp_t * np.interp(cand, lam, cp) * q * u * A_X
        ok = (pmc > 0.0) & (pmc <= p_r * (1.0 + 1e-6))
        pw += list(pmc[ok])
        th += list(p.f_ct_t * np.interp(cand[ok], lam, ct) * q * A_X)
        sf += list(p.f_cs_t * np.interp(cand[ok], lam, cs) * q * A_X)
    return np.array(pw), np.array(th), np.array(sf)


def vawt_power_curve(u, p):
    """Electrical power (kW) of the turbine generating (at berth), for a vector of winds u (m/s)."""
    if p.turb != "dmst":
        return power_curve(A_X, p.cp_t, np.asarray(u, float))
    return np.array([float(vawt_states(float(x), p)[0].max()) for x in np.atleast_1d(u)]) * ETA_GEN * ETA_CONV / 1e3


def bin_terms(u_rot, ratio, gamma, v_s, p):
    """Terms per unit (equivalent kW of the main engines) in one wind state, with the induced resistance of the total
    side force already deducted. 'r0_wind' is the power (kW, whole ship) that the wind on hull and superstructure adds
    to R0: its longitudinal resistance and the resistance induced by its side force."""
    u_a, b = apparent(u_rot * p.f_u, gamma, v_s)
    out = dict(rot=0.0, rot_fx=0.0, ext=0.0, vawt_opt=0.0, vawt_raw=0.0, sav=None, sav3=None,
               running=False, cosb=1.0, k=0.0, pen_r1=0.0, fy_r1=0.0, f_s=0.0,
               b_stop=0.0, ext_raw=None, vawt_ok=False, r0_wind=0.0, pen_self=0.0, pen_cross=0.0,
               r30=0.0, b_stop30=0.0, spin_r1=0.0, rot_run=0.0, r_d=0.0)
    eta_d = p.eta_d
    xs_m, xu_m = p.x_s * LPP, p.x_u * LPP
    if p.r0_wind:
        x_s, fy_s0 = ship_wind(u_rot, gamma, v_s, p)
        out["r0_wind"] = (x_s + float(induced(fy_s0, fy_s0 * xs_m, v_s, p))) * v_s / eta_d / 1e3
    if u_a < 0.5:
        return out
    q = 0.5 * RHO * u_a ** 2
    sb, cb = math.sin(b), math.cos(b)
    out["cosb"] = cb
    n_u = p.n_units
    # induced resistance: penalty per unit (kW) for a side force per unit fy (N)
    f_s = ship_side_force(u_rot, gamma, v_s, p)
    out["f_s"] = f_s
    k_pen = v_s / eta_d / 1e3 / n_u
    r_s0 = float(induced(f_s, f_s * xs_m, v_s, p))

    def pen(fy_unit):
        fr = n_u * np.asarray(fy_unit, float)
        return k_pen * (induced(f_s + fr, f_s * xs_m + fr * xu_m, v_s, p) - r_s0)

    k_max = min(p.k_max, p.us_max / u_a)
    grid = K_GRID if p.k_max <= SR_NOM else np.round(np.arange(0.25, p.k_max + 1e-4, 0.25), 2)
    if p.k_opt:
        ks = grid[grid <= k_max + 1e-9]
        if ks.size == 0 or ks[-1] < k_max - 1e-9:
            ks = np.append(ks, k_max)
    else:
        ks = np.array([k_max])
    cl, cd = rotor_coeffs(ks, p)
    # interference between units in tandem (sensitivity, R1-6): with the apparent wind within int_beta of the
    # centreline, the n - 1 downstream units lose the fraction int_frac of their lift
    f_int = 1.0
    if p.int_frac > 0.0:
        bd = math.degrees(abs(b))
        if bd < p.int_beta or bd > 180.0 - p.int_beta:
            f_int = 1.0 - p.int_frac * (n_u - 1) / n_u
    cl = cl * f_int
    fx = q * A_F * (cl * sb - cd * cb)
    fy = q * A_F * (cl * cb + cd * sb)
    p_spin = spin_power(u_a, ks, p.cf)
    br_all = (fx * v_s / eta_d - p.kappa * p_spin / (ETA_MOT * ETA_CONV)) / 1e3
    run_all = br_all - pen(fy)
    # version 5 (second simulated review, R1-3): common structural envelope. With m_cap, each configuration may only
    # spin at the spin ratios whose base bending moment (forces added in magnitude, as in moments.py) does not exceed
    # m_cap; if none fits, the rotor stops (the extension is retracted).
    cap = m_cap_value(p)
    m_rot = q * A_F * np.hypot(cl, cd) * H_F / 2.0
    if cap is not None:
        run_all = np.where(m_rot <= cap * (1 + 1e-9), run_all, -1e30)
    i = int(np.argmax(run_all))
    br = float(run_all[i])
    # stopped rotor: drag of the cylinder (counts if p.stop_drag)
    b_stop = -p.cd0 * q * A_F * cb * v_s / eta_d / 1e3
    fy_stop = p.cd0 * q * A_F * sb
    b_ref, fy_ref = (b_stop, fy_stop) if p.stop_drag else (0.0, 0.0)
    ref = b_ref - float(pen(fy_ref))
    running = br > ref
    out["running"] = running
    out["rot_run"] = br                       # rotor running at its best k (even if stopping it pays)
    out["rot"] = max(br, ref)
    out["b_stop"] = ref                       # stopped rotor (with its induced resistance), for the switching cost
    fy_r1 = float(fy[i]) if running else fy_ref
    out["fy_r1"] = fy_r1
    out["pen_r1"] = float(pen(fy_r1))
    # split of the induced resistance of R1 (per unit): own term of the units and cross term with the ship
    out["pen_self"] = k_pen * float(induced(n_u * fy_r1, n_u * fy_r1 * xu_m, v_s, p))
    out["pen_cross"] = out["pen_r1"] - out["pen_self"]
    if running:
        out["k"] = float(ks[i])
        out["rot_fx"] = float(fx[i]) * v_s / eta_d / 1e3
        out["spin_r1"] = p.kappa * float(p_spin[i]) / (ETA_MOT * ETA_CONV) / 1e3
    if p.r1_only:
        return out
    # upper end: wind at 55 m
    u_at, bt_ang = apparent(u_rot * p.f_u * ratio, gamma, v_s)
    qt = 0.5 * RHO * u_at ** 2
    sbt, cbt = math.sin(bt_ang), math.cos(bt_ang)
    # extension fixed to the rotor: same angular velocity as the rotor; k chosen for the whole
    if u_at > 0.5 and (p.ua_ext_max is None or u_a <= p.ua_ext_max) and (running or not p.level0):
        u_sx = ks * u_a * D_X / D_F
        cl_x, cd_x = rotor_coeffs(u_sx / u_at, p)
        cl_x = cl_x * f_int
        fx_x = qt * A_X * (p.cl_ext * cl_x * sbt - cd_x * cbt)
        fy_x = qt * A_X * (p.cl_ext * cl_x * cbt + cd_x * sbt)
        p_x = 0.5 * RHO * p.cf * u_sx ** 3 * math.pi * D_X * H_X
        tot = br_all + (fx_x * v_s / eta_d - p.kappa * p_x / (ETA_MOT * ETA_CONV)) / 1e3 - pen(fy + fy_x)
        if cap is not None:
            m_x = m_rot + qt * A_X * np.hypot(p.cl_ext * cl_x, cd_x) * (H_F + H_X / 2.0)
            tot = np.where(m_x <= cap * (1 + 1e-9), tot, -1e30)
        out["ext_raw"] = float(tot.max()) - out["rot"]          # increment without clipping (hourly model)
        out["ext"] = max(out["ext_raw"], 0.0)
    # 30 m rotor (comparator R1-30): upper segment of 6 m with the rotor's diameter and the same angular velocity,
    # always mounted; running or stopped
    if u_at > 0.5 and not p.level0:
        u_s30 = ks * u_a
        cl_u, cd_u = rotor_coeffs(u_s30 / u_at, p)
        cl_u = cl_u * f_int
        fx_u = qt * A_UP * (p.f30 * cl_u * sbt - cd_u * cbt)
        fy_u = qt * A_UP * (p.f30 * cl_u * cbt + cd_u * sbt)
        p_u = 0.5 * RHO * p.cf * u_s30 ** 3 * math.pi * D_F * H_X
        run30_all = br_all + (fx_u * v_s / eta_d - p.kappa * p_u / (ETA_MOT * ETA_CONV)) / 1e3 - pen(fy + fy_u)
        if cap is not None:
            m_30 = m_rot + qt * A_UP * np.hypot(p.f30 * cl_u, cd_u) * (H_F + H_X / 2.0)
            run30_all = np.where(m_30 <= cap * (1 + 1e-9), run30_all, -1e30)
        run30 = float(np.max(run30_all))
        bs30 = -p.cd0 * qt * A_UP * cbt * v_s / eta_d / 1e3
        fys30 = p.cd0 * qt * A_UP * sbt
        if p.stop_drag:
            stop30 = b_stop + bs30 - float(pen(fy_stop + fys30))
        else:
            stop30 = 0.0
        out["r30"], out["b_stop30"] = max(run30, stop30), stop30
    # wider rotor of the same height (comparator R1-D5, version 5): same coefficients (with an endplate, C_L hardly
    # depends on the aspect ratio in De Marco et al. 2016), area and spin power with its diameter, the same largest
    # surface speed
    if not p.level0 and not p.r1_only:
        d_a = p.d_alt
        a_d = d_a * H_F
        ks_d = grid[grid <= k_max + 1e-9]
        if ks_d.size == 0 or ks_d[-1] < k_max - 1e-9:
            ks_d = np.append(ks_d, k_max)
        cl_d, cd_d = rotor_coeffs(ks_d, p)
        cl_d = cl_d * f_int
        fx_d = q * a_d * (cl_d * sb - cd_d * cb)
        fy_d = q * a_d * (cl_d * cb + cd_d * sb)
        p_d = spin_power(u_a, ks_d, p.cf, d=d_a)
        run_d = (fx_d * v_s / eta_d - p.kappa * p_d / (ETA_MOT * ETA_CONV)) / 1e3 - pen(fy_d)
        if cap is not None:
            run_d = np.where(q * a_d * np.hypot(cl_d, cd_d) * H_F / 2.0 <= cap * (1 + 1e-9), run_d, -1e30)
        b_stop_d = -p.cd0 * q * a_d * cb * v_s / eta_d / 1e3
        fy_stop_d = p.cd0 * q * a_d * sb
        stop_d = (b_stop_d - float(pen(fy_stop_d))) if p.stop_drag else 0.0
        out["r_d"] = max(float(np.max(run_d)), stop_d)
    # VAWT at the top
    if p.turb != "dmst" and u_at <= UA_CUT_T:          # versions 1 and 2: assumed Cp and Cp/Ct, no parked turbine
        p_t = qt * u_at * A_X * p.cp_t
        t_t = qt * A_X * p.cp_t / p.cpct_t
        fy_t = t_t * sbt                       # side force of the turbine itself (C_s) not modelled
        if running:   # mechanical assistance up to the spin power; the excess, through the generator
            ps = float(p_spin[i])
            mech = min(ETA_M * p_t, ps)
            saved_e = mech / (ETA_MOT * ETA_CONV) + max(ETA_M * p_t - ps, 0.0) * ETA_GEN * ETA_CONV
            inter = p.eps * float(fx[i]) * v_s / eta_d
            fy_inter = p.eps * float(fy[i])
        else:
            saved_e = p_t * ETA_GEN * ETA_CONV
            inter, fy_inter = 0.0, 0.0
        bt = (-t_t * cbt * v_s / eta_d + p.kappa * saved_e + inter) / 1e3 \
            - float(pen(fy_r1 + fy_t + fy_inter) - pen(fy_r1))
        out["vawt_raw"], out["vawt_opt"] = bt, max(bt, 0.0)
        out["vawt_ok"] = True                                   # above cut-out it is not deployed
    elif p.turb == "dmst" and u_at > 0.5:
        # DMST model: the turbine chooses the state (parked or lambda) that saves most; its thrust follows the apparent
        # wind and its side force is perpendicular, on the side of the rotor's lift if it turns with it
        pw, th, sf = vawt_states(u_at * p.u_turb, p)
        if running:   # mechanical assistance up to the spin power; the excess, through the generator
            ps = float(p_spin[i])
            saved_e = np.minimum(ETA_M * pw, ps) / (ETA_MOT * ETA_CONV) \
                + np.maximum(ETA_M * pw - ps, 0.0) * ETA_GEN * ETA_CONV
            inter = p.eps * float(fx[i]) * v_s / eta_d
            fy_inter = p.eps * float(fy[i])
            signs = (1.0,) if p.corot else (-1.0,)
        else:         # stopped rotor: the turbine starts in the direction that pays most
            saved_e = pw * ETA_GEN * ETA_CONV
            inter, fy_inter = 0.0, 0.0
            signs = (1.0, -1.0)
        best = -np.inf
        bt_t = bt_ang + math.radians(p.turb_dang)        # flow turned by the rotor (sensitivity)
        sbt_t, cbt_t = math.sin(bt_t), math.cos(bt_t)
        for sg in signs:
            fx_t = sg * sf * sbt_t - th * cbt_t
            fy_t = sg * sf * cbt_t + th * sbt_t
            bt = (fx_t * v_s / eta_d + p.kappa * saved_e + inter) / 1e3 - (pen(fy_r1 + fy_t + fy_inter) - pen(fy_r1))
            best = max(best, float(np.max(bt)))
        out["vawt_raw"], out["vawt_opt"] = best, max(best, 0.0)
        out["vawt_ok"] = u_at <= UA_CUT_T                       # above cut-out it is not deployed (T4)
    # rotor opened as a Savonius (at rotor height): absolute benefit, with the same criterion as R1 (if the drag of
    # the stopped rotor is not counted in R1, it is measured against the stopped rotor with drag). Side force of the
    # opened rotor itself (perpendicular to the wind) only in the sensitivity cy_s.
    if u_a <= UA_CUT_S:
        for key, area_c in (("sav", A_F), ("sav3", A_F + A_X)):
            a_s = p.ds * area_c
            p_e = ETA_GEN * ETA_CONV * p.cp_s * q * u_a * a_s
            fx_sav = -p.ct_s * a_s * q * cb + p.cy_s * a_s * q * sb
            s_abs = (fx_sav * v_s / eta_d + p.kappa * p_e) / 1e3
            fy_sav = p.ct_s * a_s * q * sb + p.cy_s * a_s * q * cb
            out[key] = (s_abs - float(pen(fy_sav))) - (b_stop - float(pen(fy_stop))) + ref
    return out


def topologies(t):
    """Incremental benefit per unit over R1 (kW) at sea, by topology; t = bin_terms.
    Rules of the Savonius conversion:
      claim1   : opens whenever the rotor does not pay (text of claim 1 of the patent, headwinds included);
      intended : opens only if the rotor does not pay and the apparent wind pushes (from abaft the beam);
      astern   : opens whenever the apparent wind pushes, even if the rotor pays.
    Comparators: R1t, tiltable rotor (tilted, it has no drag and no side force); R1-30, 30 m rotor; R1-D5, wider rotor;
    R1t+X, tiltable rotor with extension. Keys starting with '_' are diagnostics per unit."""
    r1, be, run = t["rot"], t["ext"], t["running"]
    s, s3 = t["sav"], t["sav3"]
    push = t["cosb"] < 0.0
    sav_opt = max(0.0, s - r1) if s is not None else 0.0
    sav_c1 = (s - r1) if (s is not None and not run) else 0.0
    sav_int = (s - r1) if (s is not None and not run and push) else 0.0
    sav_ast = (s - r1) if (s is not None and push) else 0.0
    s_best = max(v for v in (s, s3) if v is not None) if s is not None else None
    sav3_opt = max(0.0, s_best - (r1 + be)) if s_best is not None else 0.0
    sav3_c1 = (s3 - r1) if (s3 is not None and not run) else 0.0
    return {
        "R1": r1,
        "T1": t["vawt_raw"],
        "T2": t["vawt_raw"] + sav_opt,
        "T2 claim1": t["vawt_raw"] + sav_c1,
        "T2 intended": t["vawt_raw"] + sav_int,
        "T3": be + sav3_opt,
        "T3 claim1": be + sav3_c1,
        "T3 intended": be + sav_int,
        "T3x": be,
        "T4": t["vawt_opt"],
        "T5": sav_opt,
        "T5 claim1": sav_c1,
        "T5 intended": sav_int,
        "T5 astern": sav_ast,
        "R1t": max(0.0, -r1),
        "R1-30": t["r30"] - r1,
        "R1-D5": t["r_d"] - r1,
        "R1t+X": max(0.0, t["rot_run"], r1 + be) - r1,
        "_pen R1": t["pen_r1"],          # diagnostic: induced resistance paid by R1 (kW per unit)
        "_pen own": t["pen_self"],
        "_pen cross": t["pen_cross"],
        "_R1 running": 1.0 if run else 0.0,
        "_R1 thrust": t["rot_fx"],
        "_R1 spin": t["spin_r1"],
        "_R1 stopped": -t["b_stop"] if not run else 0.0,
        "_R0 wind": t["r0_wind"],        # kW of the whole ship (not per unit)
    }


def load_wind(path=None):
    """Wind histograms. By default (version 5), those conditioned on the timetable (era5_conditioned.py): by crossing
    and season, and at berth only in the berth hours. With path = 'era5_route.json', those of version 4 (all hours of
    each point; seasons April-November and December-March)."""
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, WIND_FILE) if path is None else path
    if not os.path.isabs(path) and not os.path.exists(path):
        path = os.path.join(here, path)
    d = json.load(open(path, encoding="utf-8"))
    if "crossings" in d:
        nu, ng = d["meta"]["NU"], d["meta"]["NG"]
        groups = {"_cond": True, "cross": {}}
        for name, months in d["crossings"].items():
            groups["cross"][name] = {}
            for s, ms in SEASON_MONTHS.items():
                F, rs, rc = np.zeros((nu, ng)), np.zeros(nu), np.zeros(nu)
                for m in ms:
                    e = months.get(str(m))
                    if e:
                        F += np.array(e["F"], float)
                        rs += np.array(e["rs"], float)
                        rc += np.array(e["rc"], float)
                ratio = np.where(rc > 0, rs / np.where(rc > 0, rc, 1.0), np.nan)
                last = 1.0
                for i in range(nu):                     # bins without hours: the ratio of the previous row
                    if np.isfinite(ratio[i]):
                        last = ratio[i]
                    ratio[i] = last
                groups["cross"][name][s] = dict(F=F / F.sum() if F.sum() > 0 else F, ratio=ratio)
        port_months = d["port_months"]
        d["ports"] = {}
        for prt, months in port_months.items():
            hr, ht, ns = {}, {}, {}
            for s, ms in SEASON_MONTHS.items():
                es = [months[str(m)] for m in ms if str(m) in months]
                hr[s] = np.sum([e["hr"] for e in es], axis=0)
                ht[s] = np.sum([e["ht"] for e in es], axis=0)
                ns[s] = float(sum(e["h"] for e in es))
            d["ports"][prt] = dict(hr=hr, ht=ht, ns=ns)
        return d, groups, nu, ng
    nu, ng = d["meta"]["NU"], d["meta"]["NG"]

    def ffill(vals):
        out, last = [], vals[0]
        for v in vals:
            last = v if v is not None else last
            out.append(last)
        return np.array(out, float)
    groups = {}
    for g, gd in d["groups"].items():
        groups[g] = dict(F=[np.array(gd["F"][s], float).reshape(nu, ng) / 1e6 for s in (0, 1)],
                         ratio=[ffill(gd["ratio"][s]) for s in (0, 1)])
    return d, groups, nu, ng


def is_cond(groups):
    return isinstance(groups, dict) and groups.get("_cond", False)


def season_days(groups):
    return DAYS if is_cond(groups) else DAYS_V4


def hist_entries(groups, season, name):
    """(weight, F, ratio) of the histograms that represent a crossing in a season."""
    if is_cond(groups):
        e = groups["cross"][name][season]
        return [(1.0, e["F"], e["ratio"])]
    w_tot = sum(NM_OPEN.values())
    return [(NM_OPEN[g] / w_tot, gd["F"][season], gd["ratio"][season]) for g, gd in groups.items()]


def crossing_speed(hours, p):
    t_open = hours - NM_APPROACH / V_APPROACH
    return p.nm_open / t_open + p.dv, t_open


def mean_benefits(groups, nu, ng, season, direction, v_kn, p, name=None):
    """Time mean per unit (kW) of each topology on one crossing (name: crossing, needed with the histograms
    conditioned on the timetable)."""
    v = v_kn * KN
    acc = {}
    for wg, F, ratio in hist_entries(groups, season, name):
        for iu in range(nu):
            u = iu + 0.5
            for ig in range(ng):
                f = F[iu, ig]
                if f == 0:
                    continue
                gam = (ig + 0.5) * 15.0
                if direction == "N":
                    gam = 180.0 - gam
                tt = topologies(bin_terms(u, ratio[iu], gam, v, p))
                for key, val in tt.items():
                    acc[key] = acc.get(key, 0.0) + wg * f * val
    return acc


def r0_wind_power(groups, nu, ng, season, direction, v_kn, p, name=None):
    """Mean power (kW, whole ship) that the wind on hull and superstructure adds to R0 on one crossing: the same as the
    key '_R0 wind' of mean_benefits, without computing the topologies."""
    if not p.r0_wind:
        return 0.0
    v = v_kn * KN
    acc = 0.0
    for wg, F, _ in hist_entries(groups, season, name):
        for iu in range(nu):
            for ig in range(ng):
                f = F[iu, ig]
                if f == 0:
                    continue
                gam = (ig + 0.5) * 15.0
                if direction == "N":
                    gam = 180.0 - gam
                x_s, fy_s = ship_wind(iu + 0.5, gam, v, p)
                acc += wg * f * (x_s + float(induced(fy_s, fy_s * p.x_s * LPP, v, p))) * v / p.eta_d / 1e3
    return acc


def power_curve(area, cp, u):
    """Electrical power (kW) of a wind generator: cut-in 3 m/s, rated 12 m/s, cut-out 20 m/s."""
    p_r = 0.5 * RHO * area * cp * 12.0 ** 3
    pw = np.where((u >= 3.0) & (u < 20.0), np.minimum(0.5 * RHO * area * cp * u ** 3, p_r), 0.0)
    return pw * ETA_GEN * ETA_CONV / 1e3


def port_power(d, nu, p=None):
    """Mean electrical power per unit at berth (kW): VAWT at 55 m; Savonius at 40 m."""
    p = p or P()
    res = {}
    u = np.arange(nu) + 0.5
    for port, pd in d["ports"].items():
        res[port] = []
        for s in pd["ns"] if isinstance(pd["ns"], dict) else (0, 1):
            hr, ht, n = np.array(pd["hr"][s], float), np.array(pd["ht"][s], float), pd["ns"][s]
            vawt = float((vawt_power_curve(u, p) * ht).sum() / n)
            sav = float((power_curve(p.ds * A_F, p.cp_s, u) * hr).sum() / n)
            sav3 = float((power_curve(p.ds * (A_F + A_X), p.cp_s, u) * hr).sum() / n)
            if isinstance(res[port], list):
                res[port] = {}
            res[port][s] = dict(vawt=vawt, sav=sav, sav3=sav3)
    return res


def calm_power(v_kn, p):
    """Brake power (kW) in calm water with the sea margin: power law through 85 % MCR at 22 kn."""
    return (1.0 + p.sea_margin) * LOAD_DESIGN * P_MCR * (v_kn / V_DESIGN) ** p.n_prop


def crossing_fuel(v_kn, t_open, p_wind_kw, p):
    """Fuel of R0 on one crossing (t) and factor of the marginal consumption in open sea. The engine load includes the
    shaft generators; port approaches have no wind and no rotors."""
    p_pto = P_HOTEL_SEA / ETA_PTO
    p_open = calm_power(v_kn, p) + p_wind_kw
    p_app = calm_power(V_APPROACH, p)
    t_app = NM_APPROACH / V_APPROACH
    if p.sfoc_curve:
        l_open, l_app = (p_open + p_pto) / P_MCR, (p_app + p_pto) / P_MCR
        f_open, f_app, f_marg = sfoc_factor(l_open), sfoc_factor(l_app), sfoc_marginal_factor(l_open)
    else:
        l_open, f_open, f_app, f_marg = (p_open + p_pto) / P_MCR, 1.0, 1.0, 1.0
    fuel = ((p_open + p_pto) * t_open * f_open + (p_app + p_pto) * t_app * f_app) * p.sfoc_me / 1e6
    return fuel, f_marg, l_open


def annual(d, groups, nu, ng, p=None, seasons=None):
    """Annual fuel (t) of R0 and savings of R1 and of each topology, in the given seasons (by default, the base
    service). Savings use the marginal consumption of each crossing."""
    p = p or P()
    seasons = BASE_SEASONS if seasons is None else seasons
    sea = {}
    fuel_sea = 0.0
    fuel_me = 0.0
    hours_sea = 0.0
    loads = []
    days = season_days(groups)
    for s in seasons:
        weeks = days[s] / 7.0
        for name, direction, hours in crossings(p):
            v_kn, t_open = crossing_speed(hours, p)
            mb = mean_benefits(groups, nu, ng, s, direction, v_kn, p, name=name)
            p_wind = mb.get("_R0 wind", 0.0)
            f_cross, f_marg, load = crossing_fuel(v_kn, t_open, p_wind, p)
            loads.append(load)
            for key, val in mb.items():
                n = 1 if key == "_R0 wind" else p.n_units
                sea[key] = sea.get(key, 0.0) + n * val * t_open * weeks * p.sfoc_me * f_marg / 1e6
            fuel_sea += f_cross * weeks
            fuel_me += (calm_power(v_kn, p) + p_wind) * t_open * weeks * p.sfoc_me / 1e6
            hours_sea += hours * weeks
    pp = port_power(d, nu, p)
    port = {"T1": 0.0, "T2": 0.0, "T3": 0.0, "T4": 0.0, "T5": 0.0}
    hours_port = 0.0
    for s in seasons:
        weeks = days[s] / 7.0
        for prt, h in PORT_HOURS_WEEK.items():
            e = p.n_units * h * weeks * SFOC_AE / 1e6
            port["T1"] += e * pp[prt][s]["vawt"]
            port["T4"] += e * pp[prt][s]["vawt"]
            port["T2"] += e * (pp[prt][s]["vawt"] + pp[prt][s]["sav"])
            port["T5"] += e * pp[prt][s]["sav"]
            port["T3"] += e * pp[prt][s]["sav3"]
            hours_port += h * weeks
    fuel_berth = p.hotel_berth * hours_port * SFOC_AE / 1e6
    return dict(sea=sea, port=port, fuel_total=fuel_sea + fuel_berth, fuel_sea=fuel_sea, fuel_berth=fuel_berth,
                fuel_me=fuel_me, hours_sea=hours_sea, hours_port=hours_port, port_kw=pp,
                load_mean=float(np.mean(loads)) if loads else 0.0, seasons=list(seasons))


def table(res):
    rows = []
    for key in ("T3", "T3x", "T2", "T5", "T4", "T1"):
        s = res["sea"][key]
        pt = res["port"].get(key.replace("x", ""), 0.0) if key != "T3x" else 0.0
        rows.append((key, s, pt, s + pt, 100 * (s + pt) / res["fuel_total"]))
    return rows


def comparators(res):
    """Non-patented comparators relative to R1 (t/yr, at sea): tiltable rotor and 30 m rotor."""
    return {k: res["sea"][k] for k in ("R1t", "R1-30")}


def break_even(d, groups, nu, ng, base, lo=0.0, hi=0.6, seasons=None, **kw):
    """Interaction factor eps* at which T4 equals T3x at sea (same parameters kw)."""
    target = base["sea"]["T3x"]
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if annual(d, groups, nu, ng, P(eps=mid, **kw), seasons=seasons)["sea"]["T4"] < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    d, groups, nu, ng = load_wind()
    for label, p, ss in (("version 5, April-October", P(), None), ("version 5, whole year", P(), ALL_SEASONS)):
        base = annual(d, groups, nu, ng, p, seasons=ss)
        print(f"== {label}")
        print(f"Fuel R0: {base['fuel_total']:.0f} t (sea {base['fuel_sea']:.0f}); mean load "
              f"{base['load_mean']:.2f}; wind on the ship {base['sea'].get('_R0 wind', 0):.0f} t")
        print(f"R1, {p.n_units} rotors: {base['sea']['R1']:.1f} t ({100 * base['sea']['R1'] / base['fuel_total']:.2f} %); "
              f"induced resistance paid by R1 {base['sea']['_pen R1']:.1f} t (own {base['sea']['_pen own']:.1f}, "
              f"cross {base['sea']['_pen cross']:.1f})")
        for r in table(base):
            print(f"  {r[0]:4s} sea {r[1]:8.1f}  berth {r[2]:6.1f}  total {r[3]:8.1f}  ({r[4]:.2f} %)")
        print("  comparators:", {k: round(v, 1) for k, v in comparators(base).items()})
        for key in ("T5 claim1", "T5 intended", "T5 astern", "T2 claim1", "T2 intended", "T3 claim1", "T3 intended"):
            print(f"  {key:12s} sea {base['sea'][key]:8.1f}")
