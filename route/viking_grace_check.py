#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full-scale validation case of the conventional rotor (R1): Viking Grace (Viking Line), with a Norsepower rotor of
24 m x 4 m (Lloyd's Register, press release of 30-04-2018), the same size as that of this article, on the route
Turku-Mariehamn-Stockholm. Version 4 (22-10-2026), in response to the simulated review (R1-1, R2-6).

Reference figure (Norsepower, release of 16-05-2019): with a regression analysis by NAPA, the breakdown of the
propulsion power with on-board data from ABB and strain gauges on the rotor, the verified saving is 231-315 t of LNG
per year, «equalling an average propulsion power between 207kW and 282kW», at a service speed of 21 kn. The release
does not say whether the figure is net of the electrical consumption of the rotor; from how it was obtained
(propulsion power and measured thrust), it is treated here as gross of that consumption. Norsepower (presentation of
2018) gives for its rotors a mean electrical consumption of 15-35 kW and a maximum spin speed of 250 rpm.

Method: the same rotor code as in Portsmouth-Bilbao (route_model.bin_terms with r1_only = True): coefficients, optimal
spin ratio limited by k_max and by the maximum surface speed, spin power and stopped rotor.
Wind: hourly ERA5 2015-2024 at 10 points of the route (download_era5_viking_grace.py; Open-Meteo, models=era5),
extrapolated to 40 m with the hourly exponent between 10 and 100 m (as on the route of the Galicia), in histograms of
1 m/s and 15 degrees per point, leg and direction. Each point represents half of the legs around it; the ports of
Turku and Stockholm are assigned to the nearest point. The coordinates approximate the route, not the actual track.
Ship (Wärtsilä, Viking Grace reference sheet): electric propulsion with four dual-fuel 8L50DF engines (30,400 kW),
21.8 kn service speed, 6.80 m design draught. The rotor is fed from the ship's grid: kappa = 1.
Quantities (kW, means over all sailing hours, in both directions):
  thrust      F_x V / eta_D with the rotor running (propulsion power equivalent to the thrust, gross);
  thrust_eff  F_x V (effective power of the thrust);
  spin        electrical consumption of the rotor; stopped: drag of the stopped rotor / eta_D;
  induced     resistance induced by the side force of the rotor (only its own term: that of the wind on the Viking
              Grace is unknown), with the draught of 6.80 m;
  net         thrust - spin - induced - stopped (what the route model counts as the saving of R1).
Output: viking_grace_check.json.
"""
import datetime as dt
import gzip
import json
import math
import os
import numpy as np
import route_model as rm
import era5_hourly as eh
import coefficients as cf

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "era5_hourly_viking_grace.csv.gz")
N_HOURS = 87672
T0 = dt.datetime(2015, 1, 1, tzinfo=dt.timezone.utc)
NU, NG = 30, 12
NM_KM = 1.852

# ERA5 points (id, latitude, longitude) and ports
POINTS = {"VG1": (60.380, 22.000), "VG2": (60.240, 21.620), "VG3": (60.150, 21.150), "VG4": (60.050, 20.650),
          "VG5": (60.050, 20.200), "MARIEHAMN": (60.095, 19.935), "VG6": (59.950, 19.500), "VG7": (59.800, 19.150),
          "VG8": (59.620, 18.850), "VG9": (59.420, 18.400)}
PORTS = {"TURKU": (60.434, 22.222), "STOCKHOLM": (59.316, 18.090)}
ROUTE = ["TURKU", "VG1", "VG2", "VG3", "VG4", "VG5", "MARIEHAMN", "VG6", "VG7", "VG8", "VG9", "STOCKHOLM"]

T_VG = 6.80                 # m, design draught (Wärtsilä)
V_SERVICE = 21.0            # kn, service speed of the Norsepower release (21.8 kn according to Wärtsilä)
REF_KW = (207.0, 282.0)     # Norsepower (2019): verified mean propulsion power
REF_T = (231.0, 315.0)      # t of LNG per year
TEST_PERIOD = (dt.datetime(2018, 5, 1, tzinfo=dt.timezone.utc), dt.datetime(2019, 5, 1, tzinfo=dt.timezone.utc))


def load_vg(path=DATA):
    data = {p: {k: np.full(N_HOURS, np.nan) for k in ("u10", "u100", "d100")} for p in POINTS}
    meta = []
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                meta.append(line.strip())
                continue
            if line.startswith("point,"):
                continue
            pid, h, a, b, c = line.rstrip("\n").split(",")
            h = int(h)
            if pid in data and 0 <= h < N_HOURS:
                dd = data[pid]
                if a:
                    dd["u10"][h] = float(a)
                if b:
                    dd["u100"][h] = float(b)
                if c:
                    dd["d100"][h] = float(c)
    return data, meta


def land_cells(meta):
    """Points whose ERA5 cell has elevation > 0 in the header (cell with land)."""
    out = set()
    for m in meta:
        if m.startswith("# VG") or m.startswith("# MARIEHAMN"):
            pid = m[2:].split(":")[0]
            elev = float(m.split("elevation")[1].split(";")[0])
            if elev > 0.0:
                out.add(pid)
    return out


def coord(n):
    return POINTS.get(n) or PORTS[n]


def dist_bearing(a, b):
    """Distance (nm) and initial heading (degrees) from a to b (great circle)."""
    la1, lo1, la2, lo2 = map(math.radians, (a[0], a[1], b[0], b[1]))
    dlo = lo2 - lo1
    h = math.sin((la2 - la1) / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin(dlo / 2) ** 2
    d_km = 2 * 6371.0 * math.asin(math.sqrt(h))
    y = math.sin(dlo) * math.cos(la2)
    x = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(dlo)
    return d_km / NM_KM, (math.degrees(math.atan2(y, x)) + 360.0) % 360.0


def half_legs():
    """(ERA5 point, heading towards Stockholm, miles): each leg is split into two halves, each with its point."""
    out = []
    for a, b in zip(ROUTE[:-1], ROUTE[1:]):
        d, psi = dist_bearing(coord(a), coord(b))
        first = a if a in POINTS else b
        second = b if b in POINTS else a
        out += [(first, psi, d / 2.0), (second, psi, d / 2.0)]
    return out


def histograms(data, mask_hours=None, exclude=()):
    """Normalised histograms (NU x NG) per half-leg and direction, with the weight of their miles."""
    hist = []
    sub = [s for s in half_legs() if s[0] not in exclude]
    w_tot = sum(s[2] for s in sub)
    for pid, psi, nm in sub:
        dd = data[pid]
        ok = ~(np.isnan(dd["u10"]) | np.isnan(dd["u100"]) | np.isnan(dd["d100"]))
        if mask_hours is not None:
            ok &= mask_hours
        ur, _, _ = eh.heights(dd)
        iu = np.minimum(NU - 1, np.floor(np.nan_to_num(ur)).astype(int))
        for direction, heading in (("W", psi), ("E", (psi + 180.0) % 360.0)):
            ig = np.minimum(NG - 1, np.floor(eh.fold(dd["d100"] - heading) / 15.0).astype(int))
            cnt = np.bincount((iu * NG + ig)[ok], minlength=NU * NG).reshape(NU, NG).astype(float)
            hist.append((pid, direction, 0.5 * nm / w_tot, cnt / cnt.sum()))
    return hist, w_tot


def evaluate(hist, v_kn, p):
    """Means (kW) of the quantities of the rotor over all sailing hours."""
    v = v_kn * rm.KN
    acc = dict(thrust=0.0, thrust_eff=0.0, spin=0.0, stopped=0.0, induced=0.0, net=0.0, running=0.0, k_mean=0.0)
    cache = {}
    for pid, direction, w, F in hist:
        for iu in range(NU):
            for ig in range(NG):
                f = F[iu, ig]
                if f == 0.0:
                    continue
                key = (iu, ig)
                if key not in cache:
                    t = rm.bin_terms(iu + 0.5, 1.0, (ig + 0.5) * 15.0, v, p)
                    run = t["running"]
                    cache[key] = dict(thrust=t["rot_fx"], thrust_eff=t["rot_fx"] * p.eta_d, spin=t["spin_r1"],
                                      stopped=(-t["b_stop"] if not run else 0.0),
                                      induced=(t["pen_r1"] if run else 0.0), net=t["rot"],
                                      running=1.0 if run else 0.0, k_mean=t["k"])
                for k2, val in cache[key].items():
                    acc[k2] += w * f * val
    acc["k_mean"] = acc["k_mean"] / acc["running"] if acc["running"] > 0 else 0.0
    return acc


def p_vg(**kw):
    base = dict(n_units=1, r1_only=True, t_draught=T_VG, kappa=1.0, windage=False, r0_wind=False, sir_model="wing")
    base.update(kw)
    return rm.P(**base)


SCENARIOS = [  # epistemic coefficient scenarios (coefficients.COEF_SCENARIOS)
    ("A: De Marco et al. (2016), k <= 3, U_s <= 36 m/s", "A"),
    ("A+: A extended, k <= 4, U_s <= 46 m/s", "A+"),
    ("AB: equal mix of A and B (version 3)", "AB"),
    ("B: Bordogna et al. (2019), AR_eff 12", "B"),
    ("B': B with the lifting line (eta 0.83)", "B'"),
    ("B'': B without reflection in the deck (eta 0.725, AR_eff 6)", "B''"),
    ("A-VG: A with C_L fitted to the Viking Grace (sea cells, 244.5 kW)", "A-VG"),
]


def main():
    data, meta = load_vg()
    land = land_cells(meta)
    sub = half_legs()
    tot = sum(s[2] for s in sub)
    print("Viking Grace: approximate Turku-Mariehamn-Stockholm route")
    print(f"  {tot:.1f} nm; ERA5 cells with land: {sorted(land)}")
    hours = np.arange(N_HOURS)
    t_ms = np.array([(T0 + dt.timedelta(hours=int(h))) for h in hours])
    test = np.array([TEST_PERIOD[0] <= t < TEST_PERIOD[1] for t in t_ms])
    h_all, _ = histograms(data)
    h_test, _ = histograms(data, mask_hours=test)
    h_sea, _ = histograms(data, exclude=land)
    u40 = {pid: float(np.nanmean(eh.heights(data[pid])[0])) for pid in POINTS}
    print("  mean wind at 40 m (m/s): " + ", ".join(f"{k} {v:.1f}" for k, v in u40.items()))
    out = dict(source="Norsepower (2019); LR (2018); Wärtsilä", ref_kw=REF_KW, ref_t=REF_T, miles=tot,
               land_cells=sorted(land), wind_40m=u40, scenarios={}, sensitivity={})
    print(f"\nReference: {REF_KW[0]:.0f}-{REF_KW[1]:.0f} kW of mean propulsion power (Norsepower, 2019)")
    hdr = f"{'scenario':58s} {'thrust':>7s} {'eff.':>6s} {'spin':>6s} {'stop':>6s} {'induc':>6s} {'net':>6s} {'run':>6s} {'k':>5s}"
    print(f"\n21 kn, ERA5 2015-2024 (kW per rotor)\n{hdr}")
    for label, esc in SCENARIOS:
        r = evaluate(h_all, V_SERVICE, p_vg(**cf.COEF_SCENARIOS[esc]))
        out["scenarios"][label] = r
        print(f"{label:58s} {r['thrust']:7.1f} {r['thrust_eff']:6.1f} {r['spin']:6.1f} {r['stopped']:6.1f} "
              f"{r['induced']:6.1f} {r['net']:6.1f} {100 * r['running']:5.0f}% {r['k_mean']:5.2f}")
    print(f"\nSensitivity (scenario A)\n{hdr}")
    sens = [("speed 21.8 kn", h_all, 21.8, {}), ("speed 18 kn", h_all, 18.0, {}),
            ("speed 16 kn", h_all, 16.0, {}), ("period of the test (05-2018 to 04-2019)", h_test, V_SERVICE, {}),
            ("without land cells", h_sea, V_SERVICE, {}), ("eta_D 0.65", h_all, V_SERVICE, dict(eta_d=0.65)),
            ("eta_D 0.75", h_all, V_SERVICE, dict(eta_d=0.75)), ("Cf 0.002 (spin)", h_all, V_SERVICE, dict(cf=0.002)),
            ("Cf 0.007 (spin)", h_all, V_SERVICE, dict(cf=0.007)),
            ("A+, without land cells", h_sea, V_SERVICE, dict(k_max=4.0, us_max=46.0)),
            ("A+, without land cells, wind x1.1 at the rotor", h_sea, V_SERVICE, dict(k_max=4.0, us_max=46.0, f_u=1.1))]
    for label, hh, v, kw in sens:
        r = evaluate(hh, v, p_vg(**kw))
        out["sensitivity"][label] = r
        print(f"{label:58s} {r['thrust']:7.1f} {r['thrust_eff']:6.1f} {r['spin']:6.1f} {r['stopped']:6.1f} "
              f"{r['induced']:6.1f} {r['net']:6.1f} {100 * r['running']:5.0f}% {r['k_mean']:5.2f}")
    # scenario fitted to the measurement (version 5; simulated review v2, R1-1 and R2-2): multiplier of C_L of
    # scenario A with which the thrust of the model, with the sea cells, equals 207, 244.5 (centre) and 282 kW
    out["cl_fit"] = {}
    for target in (REF_KW[0], 0.5 * (REF_KW[0] + REF_KW[1]), REF_KW[1]):
        lo, hi = 0.8, 3.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if evaluate(h_sea, V_SERVICE, p_vg(cl_mult=mid))["thrust"] < target:
                lo = mid
            else:
                hi = mid
        m = 0.5 * (lo + hi)
        r = evaluate(h_sea, V_SERVICE, p_vg(cl_mult=m))
        r_all = evaluate(h_all, V_SERVICE, p_vg(cl_mult=m))
        out["cl_fit"][f"{target:g}"] = dict(cl_mult=m, sea_cells=r, all_cells=r_all)
        print(f"Fit to {target:.1f} kW (sea cells): C_L x {m:.3f}; with all cells, thrust "
              f"{r_all['thrust']:.1f} kW, net {r_all['net']:.1f} kW")
        chk = cf.CL_MULT_VG.get(f"{target:g}")
        if chk is None or abs(chk - round(m, 3)) > 5e-4:
            print(f"  WARNING: coefficients.CL_MULT_VG['{target:g}'] = {chk}; the fitted value is {m:.3f}")
    # the same calculation for the Galicia (per rotor, April-October, crossings of the timetable), for comparison
    d, groups, nu, ng = rm.load_wind()
    gal = {}
    for label, esc in SCENARIOS:
        p = rm.P(**cf.COEF_SCENARIOS[esc])
        acc, hrs = {}, 0.0
        for name, direction, h in rm.crossings(p):
            v_kn, t_open = rm.crossing_speed(h, p)
            mb = rm.mean_benefits(groups, nu, ng, 0, direction, v_kn, p, name=name)
            for k2 in ("R1", "_R1 thrust", "_R1 spin", "_R1 stopped", "_pen R1", "_R1 running"):
                acc[k2] = acc.get(k2, 0.0) + mb[k2] * t_open
            hrs += t_open
        gal[label] = {k2: v2 / hrs for k2, v2 in acc.items()}
    out["galicia_apr_oct"] = gal
    print("\nGalicia, April-October, open sea (kW per rotor, same code)")
    for label, r in gal.items():
        print(f"{label:58s} thrust {r['_R1 thrust']:6.1f}  spin {r['_R1 spin']:5.1f}  stopped {r['_R1 stopped']:5.1f}  "
              f"induced {r['_pen R1']:5.1f}  net {r['R1']:6.1f}  running {100 * r['_R1 running']:3.0f}%")
    json.dump(out, open(os.path.join(HERE, "viking_grace_check.json"), "w", encoding="utf-8"), indent=1,
              ensure_ascii=False)
    return out


if __name__ == "__main__":
    main()
