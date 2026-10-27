#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A1 - Hourly voyage model. Version 5 (27-10-2026): April-October season (route_model.SEASON_MONTHS, by the month of
departure), induced resistance from the sway-yaw balance and new comparators (R1-D5).
Version 4 (22-10-2026; version 3 is in the history of the repository).
Simulates each crossing of the 2026 timetable over the hourly ERA5 wind of 2015-2024 (era5_hourly_route.csv.gz):
position, heading and speed of the ship in each hour, true and apparent wind, and benefit of each mode of each
topology with route_model.bin_terms (the same physics as the histogram model).

Changes of version 4 (simulated review R3-5 and R3-m4):
- Base case: documented service (route_model.BASE_SEASONS; April-November in version 4, April-October from version 5);
  the rest of the year as a scenario.
- R0 with the model of version 4 (wind resistance of the ship, SFOC curve, sea margin) and savings converted with the
  marginal consumption of each crossing, as in route_model.annual.
- Start and stop cost also for the conventional rotor: R1 has two modes (stopped, running) and each change lasts
  tau_R1 (0.05-0.1 h of acceleration or braking).
- Tiltable rotor (R1t, commercial comparator): modes stopped, running and tilted (no drag and no side force); each
  tilting or raising lasts tau_t (0.25-0.5 h, assumed).
- Control with imperfect forecast (model predictive control with horizon H): in each hour the current wind is known
  and, for the following hours, a forecast issued at departure with AR(1) errors of speed and direction that grow with
  the lead time; the optimal sequence over the horizon is solved, the first step is applied and the process is
  repeated. The result is evaluated with the actual wind. Two error levels (assumed, of the order of those of a
  short-range numerical model at sea): good, 1.0 -> 2.0 m/s and 10 -> 20 degrees in 24 h; poor, 2.0 -> 3.5 m/s and
  20 -> 40 degrees.

Strategies, per crossing:
  inst      the best mode in each hour, with no switching cost (equivalent to the histogram model);
  greedy    the best mode in each hour, paying each change;
  optimal   optimal sequence with switching cost by dynamic programming, knowing the wind of the whole crossing
            (upper bound of any real control);
  causal    in each hour the current wind is assumed to persist L hours (persistence);
  mpc       model predictive control with the imperfect forecast (horizon H);
  rules     rules of the patent for the Savonius (claim-1 and intended readings), paying each change.
Switching cost: each change of mode lasts tau hours, during which the unit performs as a stopped rotor; the
difference with the better of the two modes is charged (a prudent bound). The patented topologies start from the best
state of R1 in each hour (without a switching cost of R1 itself); the start and stop cost of R1 is given separately.
Generation in port is computed in the actual berth hours of the timetable.

Usage: python3 hourly_voyages.py [../data/era5_hourly_route.csv.gz]   ->  hourly_results.json
"""
import json
import math
import sys
import datetime as dt
from zoneinfo import ZoneInfo
import numpy as np
import route_model as rm
import era5_hourly as eh

LON, MAD = ZoneInfo("Europe/London"), ZoneInfo("Europe/Madrid")
T0, N_HOURS = eh.T0, eh.N_HOURS
ROUTE = [p for p in eh.POINTS if p[1] in ("C", "B")]         # southbound: C1 ... C4, B1 ... B6
NM_OPEN = rm.NM_OPEN_TOTAL
EDGES = np.concatenate([[0.0], np.cumsum([p[5] for p in ROUTE])]) * NM_OPEN / sum(p[5] for p in ROUTE)
TAUS = (0.1, 0.25, 0.5)
TAU_R1 = (0.05, 0.1)
TAU_TILT = (0.25, 0.5)
LS = (1.0, 2.0, 3.0, 6.0, 12.0)
HS = (3, 6, 12)
FORECAST = {"good": dict(su0=1.0, su1=2.0, sd0=10.0, sd1=20.0), "poor": dict(su0=2.0, su1=3.5, sd0=20.0, sd1=40.0)}
RHO_AR = 0.9
SEED = 20261022
COLS = ("w", "rot", "b_stop", "ext", "vawt", "vawt_ok", "sav", "sav3", "running", "cosb", "pen_r1", "rot_run",
        "r0_wind")


def gc_nm(a, b):
    p1, p2, dl = math.radians(a[0]), math.radians(b[0]), math.radians(b[1] - a[1])
    d = 2 * math.asin(math.sqrt(math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2))
    return d * 3440.065


# Approaches: Portsmouth-Nab-St Catherine and approach-Bilbao (route_geometry.py), scaled to 32 nm at 10 kn
_NM_P = gc_nm((50.790, -1.110), (50.670, -0.950)) + gc_nm((50.670, -0.950), (50.500, -1.300))
_NM_B = gc_nm((43.450, -3.080), (43.360, -3.030))
T_APP = {"P": _NM_P * rm.NM_APPROACH / (_NM_P + _NM_B) / rm.V_APPROACH,
         "B": _NM_B * rm.NM_APPROACH / (_NM_P + _NM_B) / rm.V_APPROACH}

# 2026 timetable in local hours: (name, direction, departure (ISO day, h, min, zone), arrival (days later, h, min, zone))
TIMETABLE = [("Portsmouth-Bilbao, Thursday", "S", (4, 19, 0, LON), (2, 8, 0, MAD)),
             ("Portsmouth-Bilbao, Sunday", "S", (7, 22, 0, LON), (2, 8, 0, MAD)),
             ("Bilbao-Portsmouth, Tuesday", "N", (2, 13, 0, MAD), (1, 20, 45, LON)),
             ("Bilbao-Portsmouth, Saturday", "N", (6, 13, 0, MAD), (1, 17, 30, LON))]
BERTHS = [("BILBAO", "Saturday", (6, 8, 0, MAD), (0, 13, 0, MAD)), ("BILBAO", "Tuesday", (2, 8, 0, MAD), (0, 13, 0, MAD)),
          ("PORTSMOUTH", "Sunday", (7, 17, 30, LON), (0, 22, 0, LON)),
          ("PORTSMOUTH", "Wednesday-Thursday", (3, 20, 45, LON), (1, 19, 0, LON))]


def _hours(d, hh, mm, tz):
    return (dt.datetime(d.year, d.month, d.day, hh, mm, tzinfo=tz) - T0).total_seconds() / 3600.0


def schedule(table, shift=0.0):
    """Intervals (hours since 2015-01-01 UTC) of each crossing or berth stay that falls entirely within the data.
    shift: shift (h) of the whole timetable, for the time-of-day diagnostic."""
    out, d = [], dt.date(2014, 12, 20)
    while d <= dt.date(2025, 1, 10):
        for name, tag, (wd, hh, mm, tz), (nd, ah, am, atz) in table:
            if d.isoweekday() == wd:
                t1 = _hours(d, hh, mm, tz) + shift
                t2 = _hours(d + dt.timedelta(days=nd), ah, am, atz) + shift
                if t1 >= 0 and t2 <= N_HOURS - 1:
                    out.append(dict(name=name, tag=tag, t1=t1, t2=t2, year=d.year, month=d.month,
                                    season=rm.season_of_month(d.month)))
        d += dt.timedelta(days=1)
    return out


def wind_fields(data):
    f = {}
    for pid, dd in data.items():
        u40, u55, _ = eh.heights(dd)
        bad = np.isnan(dd["u10"]) | np.isnan(dd["u100"]) | np.isnan(dd["d100"])   # as in the aggregation
        u40, u55 = np.where(bad, np.nan, u40), np.where(bad, np.nan, u55)
        f[pid] = dict(u40=u40, ratio=np.where(u40 > 0, u55 / np.where(u40 > 0, u40, 1.0), 1.0), u55=u55, d100=dd["d100"])
    return f


def era5_slots(t1, t2):
    """Split of the interval [t1, t2] (hours since T0) among hourly ERA5 values: each value is instantaneous, on the
    hour, and represents the hour centred on it. Returns (ERA5 hour, start, end) of each piece."""
    out = []
    for h in range(int(math.floor(t1 + 0.5)), int(math.floor(t2 + 0.5)) + 1):
        lo, hi = max(t1, h - 0.5), min(t2, h + 0.5)
        if hi > lo:
            out.append((min(N_HOURS - 1, max(0, h)), lo, hi))
    return out


def crossing_hours(c):
    """(weight h, point index, ERA5 hour) of each open-sea piece, and speed (kn)."""
    o, de = ("P", "B") if c["tag"] == "S" else ("B", "P")
    a, b = c["t1"] + T_APP[o], c["t2"] - T_APP[de]
    v_kn = NM_OPEN / (b - a)
    hs = []
    for he, lo, hi in era5_slots(a, b):
        tm = 0.5 * (lo + hi)
        dist = v_kn * (tm - a)
        s = dist if c["tag"] == "S" else NM_OPEN - dist
        ip = int(min(len(ROUTE) - 1, max(0, np.searchsorted(EDGES, s, side="right") - 1)))
        hs.append((hi - lo, ip, he))
    return v_kn, b - a, hs


def forecast_errors(n, lead0, level, rng):
    """AR(1) errors of speed (m/s) and direction (degrees) for n hours, with a standard deviation that grows linearly
    with the lead time (hours since the forecast was issued, at departure) up to 24 h."""
    f = FORECAST[level]
    eu, ed = np.zeros(n), np.zeros(n)
    zu, zd = rng.standard_normal(), rng.standard_normal()
    for i in range(n):
        if i:
            zu = RHO_AR * zu + math.sqrt(1 - RHO_AR ** 2) * rng.standard_normal()
            zd = RHO_AR * zd + math.sqrt(1 - RHO_AR ** 2) * rng.standard_normal()
        lead = min(lead0 + i, 24.0) / 24.0
        eu[i] = zu * (f["su0"] + (f["su1"] - f["su0"]) * lead)
        ed[i] = zd * (f["sd0"] + (f["sd1"] - f["sd0"]) * lead)
    return eu, ed


def evaluate(fields, cross, p, binned=None, forecast=None):
    """Benefit per unit (kW) of each state in each hour of each crossing (columns COLS).
    binned: groups of route_model.load_wind(); if given, the wind is moved to the centre of its bin (1 m/s, 15°) and
    the ratio u55/u40 to the mean of the bin, as in the histogram model (diagnostic).
    forecast: error level ('good' or 'poor'); if given, the evaluation uses the forecast wind (same hours)."""
    recs = []
    rng = np.random.default_rng(SEED)
    for ci, c in enumerate(cross):
        v_kn, t_open, hs = crossing_hours(c)
        v = v_kn * rm.KN
        if forecast is not None:
            eu, ed = forecast_errors(len(hs), 0.0, forecast, rng)
        rows = []
        for k, (w, ip, he) in enumerate(hs):
            pid, psi = ROUTE[ip][0], ROUTE[ip][4]
            fd = fields[pid]
            u, r, d100 = fd["u40"][he], fd["ratio"][he], fd["d100"][he]
            if not (np.isfinite(u) and np.isfinite(r) and np.isfinite(d100)):
                continue
            if forecast is not None:
                u, d100 = max(0.0, u + eu[k]), d100 + ed[k]
            g = float(eh.fold(d100 - psi))
            if binned is not None:
                iu = min(eh.NU - 1, int(math.floor(u)))
                s_leg = 0 if 4 <= c["month"] <= 11 else 1          # seasons of era5_route.json (version 4)
                u, r = iu + 0.5, float(binned[ROUTE[ip][1]]["ratio"][s_leg][iu])
                g = (min(eh.NG - 1, int(math.floor(g / 15.0))) + 0.5) * 15.0
            if c["tag"] == "N":
                g = 180.0 - g
            t = rm.bin_terms(float(u), float(r), g, v, p)
            rows.append((w, t["rot"], t["b_stop"], np.nan if t["ext_raw"] is None else t["ext_raw"], t["vawt_raw"],
                         t["vawt_ok"], np.nan if t["sav"] is None else t["sav"],
                         np.nan if t["sav3"] is None else t["sav3"], t["running"], t["cosb"], t["pen_r1"],
                         t["rot_run"], t["r0_wind"]))
        a = np.array(rows, float).reshape(-1, len(COLS))
        p_wind = float((a[:, 0] * a[:, 12]).sum() / a[:, 0].sum()) if len(a) else 0.0
        fuel, f_marg, load = rm.crossing_fuel(v_kn, t_open, p_wind, p)
        recs.append(dict(c=c, v_kn=v_kn, t_open=t_open, a=a, fuel=fuel, f_marg=f_marg))
    return recs


NEG = -1e18


def mode_matrix(a, topo):
    w, r1, bs, ext, vt, vok, sav, sav3 = (a[:, i] for i in range(8))
    run_v = a[:, 11]
    ext = np.where(np.isfinite(ext), ext, np.nan)
    m = {"R1": [r1],
         "R1 switched": [bs, run_v],
         "R1t": [bs, run_v, np.zeros_like(bs)],
         "T1": [r1 + vt],
         "T2": [r1 + vt, sav + vt],
         "T3": [r1, r1 + ext, sav, sav3],
         "T3x": [r1, r1 + ext],
         "T4": [r1, np.where(vok > 0, r1 + vt, np.nan)],
         "T5": [r1, sav]}[topo]
    B = np.column_stack(m)
    return np.where(np.isfinite(B), B, NEG), w, bs


def switch_cost(B, bs, t, i, j, tau):
    return tau * max(0.0, max(B[t, i], B[t, j]) - bs[t])


def path_value(B, w, bs, path, tau):
    """Energy (kWh per unit) of a sequence of modes, with departure and arrival in mode 0; and number of changes."""
    val, prev, n_sw = 0.0, 0, 0
    for t, m in enumerate(path):
        if m != prev:
            val -= switch_cost(B, bs, t, prev, m, tau)
            n_sw += 1
        val += w[t] * B[t, m]
        prev = m
    if prev != 0 and len(path):
        val -= switch_cost(B, bs, len(path) - 1, prev, 0, tau)
        n_sw += 1
    return val, n_sw


def _cost_matrix(B, bs, t, tau):
    """Cost of going from mode i to j in hour t (zero diagonal)."""
    mx = np.maximum(B[t][:, None], B[t][None, :])
    c = tau * np.maximum(0.0, mx - bs[t])
    np.fill_diagonal(c, 0.0)
    return c


def dp_window(B, w, bs, tau, m0, end_mode=None):
    """Optimal sequence in a window that starts in mode m0 (and ends in end_mode, if given)."""
    n, M = B.shape
    V = np.full(M, -np.inf)
    V[m0] = 0.0
    back = np.zeros((n, M), int)
    avail = B > NEG / 2
    for t in range(n):
        tot = V[:, None] - _cost_matrix(B, bs, t, tau)
        arg = np.argmax(tot, axis=0)
        best = tot[arg, np.arange(M)]
        nV = np.where(avail[t], best + w[t] * B[t], -np.inf)
        back[t] = arg
        V = nV
    if end_mode is not None:
        c_end = _cost_matrix(B, bs, n - 1, tau)[:, end_mode]
        final = np.where(np.isfinite(V), V - c_end, -np.inf)
    else:
        final = V
    j = int(np.argmax(final))
    path = [j]
    for t in range(n - 1, 0, -1):
        j = int(back[t, j])
        path.append(j)
    return path[::-1]


def dp(B, w, bs, tau):
    if B.shape[0] == 0:
        return 0.0, 0
    path = dp_window(B, w, bs, tau, 0, end_mode=0)
    return path_value(B, w, bs, path, tau)


def causal_path(B, w, bs, tau, L):
    """Causal control: in each hour only the wind of that hour is known and it is assumed to persist L hours (or what
    is left of the crossing); the mode changes to j if (B_j - B_current) * min(L, rest) exceeds the cost of the change.
    If the current mode is not available (structural limit, stop), it changes to the best available."""
    n, M = B.shape
    path, m, left = [], 0, float(w.sum())
    for t in range(n):
        if B[t, m] <= NEG / 2:
            m = int(np.argmax(B[t]))
        else:
            best_j, best_g = m, 0.0
            for j in range(M):
                if j == m or B[t, j] <= NEG / 2:
                    continue
                g = (B[t, j] - B[t, m]) * min(L, left) - switch_cost(B, bs, t, m, j, tau)
                if g > best_g:
                    best_j, best_g = j, g
            m = best_j
        path.append(m)
        left -= w[t]
    return path


def mpc_path(B, Bf, w, bs, bsf, tau, H):
    """Model predictive control: in hour t B[t] is known; the hours t+1 ... t+H-1 are taken from the forecast Bf. The
    window is solved starting in the current mode, the first step is applied and the clock advances one hour."""
    n, M = B.shape
    path, m = [], 0
    for t in range(n):
        hi = min(n, t + H)
        Bw = np.vstack([B[t:t + 1], Bf[t + 1:hi]])
        bw = np.concatenate([bs[t:t + 1], bsf[t + 1:hi]])
        end = 0 if hi == n else None
        win = dp_window(Bw, w[t:hi], bw, tau, m, end_mode=end)
        # the first step of the window starts from mode m: dp_window sets the mode of hour t
        m = win[0] if B[t, win[0]] > NEG / 2 else int(np.argmax(B[t]))
        path.append(m)
    return path


def strategies(recs, topo, taus=TAUS, Ls=LS, recs_fc=None, Hs=HS, mpc_taus=(0.25,)):
    """kWh per unit and crossing, and changes per crossing, for each strategy (model predictive control only with
    mpc_taus, for computing cost)."""
    out = {"inst": [], **{f"greedy {t}": [] for t in taus}, **{f"optimal {t}": [] for t in taus},
           **{f"switches greedy {t}": [] for t in taus}, **{f"switches optimal {t}": [] for t in taus},
           **{f"causal L{L:g} {t}": [] for t in taus for L in Ls},
           **{f"switches causal L{L:g} {t}": [] for t in taus for L in Ls}}
    if recs_fc:
        for lev in recs_fc:
            out.update({f"mpc {lev} H{H} {t}": [] for t in mpc_taus for H in Hs})
            out.update({f"switches mpc {lev} H{H} {t}": [] for t in mpc_taus for H in Hs})
    if topo == "T5":
        out.update({f"{r} {t}": [] for r in ("claim1", "intended") for t in (0.0,) + taus})
    for k, rec in enumerate(recs):
        a = rec["a"]
        B, w, bs = mode_matrix(a, topo)
        out["inst"].append(float((w * B.max(axis=1)).sum()))
        greedy = list(np.argmax(B, axis=1))
        for t in taus:
            gv, gs = path_value(B, w, bs, greedy, t)
            out[f"greedy {t}"].append(gv)
            out[f"switches greedy {t}"].append(gs)
            dv, ds = dp(B, w, bs, t) if B.shape[1] > 1 else (gv, gs)
            out[f"optimal {t}"].append(dv)
            out[f"switches optimal {t}"].append(ds)
            for L in Ls:
                cv, cs = path_value(B, w, bs, causal_path(B, w, bs, t, L), t) if B.shape[1] > 1 else (gv, gs)
                out[f"causal L{L:g} {t}"].append(cv)
                out[f"switches causal L{L:g} {t}"].append(cs)
            if recs_fc and t in mpc_taus:
                for lev, rf in recs_fc.items():
                    af = rf[k]["a"]
                    if af.shape != a.shape or B.shape[1] == 1:
                        for H in Hs:
                            out[f"mpc {lev} H{H} {t}"].append(gv)
                            out[f"switches mpc {lev} H{H} {t}"].append(gs)
                        continue
                    Bf, _, bsf = mode_matrix(af, topo)
                    for H in Hs:
                        mv, ms = path_value(B, w, bs, mpc_path(B, Bf, w, bs, bsf, t, H), t)
                        out[f"mpc {lev} H{H} {t}"].append(mv)
                        out[f"switches mpc {lev} H{H} {t}"].append(ms)
        if topo == "T5":
            run, cosb, sav_ok = a[:, 8] > 0, a[:, 9], B[:, 1] > NEG / 2
            rules = {"claim1": (~run) & sav_ok, "intended": (~run) & (cosb < 0) & sav_ok}
            for r, mask in rules.items():
                for t in (0.0,) + taus:
                    out[f"{r} {t}"].append(path_value(B, w, bs, list(mask.astype(int)), t)[0])
    return out


def annualise(recs, values, by=None, seasons=None, n_units=None):
    """t/yr of fuel (n units) from kWh per unit and crossing, with the marginal consumption of each crossing: mean by
    type of crossing times the weeks of each season (or of each year if by='year', only seasons of the base case).
    by = 0, 1 or 2: a single season."""
    vals = np.asarray(values, float) * np.array([r["f_marg"] for r in recs])
    n_units = n_units or rm.N_UNITS
    seasons = rm.BASE_SEASONS if seasons is None else seasons
    names = sorted({r["c"]["name"] for r in recs})
    if by == "year":
        res = {}
        for y in sorted({r["c"]["year"] for r in recs}):
            tot = 0.0
            for s in seasons:
                for nme in names:
                    sel = [k for k, r in enumerate(recs) if r["c"]["name"] == nme and r["c"]["year"] == y
                           and r["c"]["season"] == s]
                    if sel:
                        tot += vals[sel].mean() * rm.DAYS[s] / 7.0
            res[y] = n_units * tot * rm.SFOC_ME / 1e6
        return res
    tot = 0.0
    for s in (seasons if by is None else (by,)):
        for nme in names:
            sel = [k for k, r in enumerate(recs) if r["c"]["name"] == nme and r["c"]["season"] == s]
            if sel:
                tot += vals[sel].mean() * rm.DAYS[s] / 7.0
    return n_units * tot * rm.SFOC_ME / 1e6


def r0_fuel(recs, p, seasons=None):
    """Fuel of the ship without rotors (t/yr) in the given seasons, with the model of version 4."""
    seasons = rm.BASE_SEASONS if seasons is None else seasons
    names = sorted({r["c"]["name"] for r in recs})
    sea = 0.0
    for s in seasons:
        for nme in names:
            sel = [r["fuel"] for r in recs if r["c"]["name"] == nme and r["c"]["season"] == s]
            if sel:
                sea += float(np.mean(sel)) * rm.DAYS[s] / 7.0
    berth = p.hotel_berth * sum(rm.PORT_HOURS_WEEK.values()) * sum(rm.DAYS[s] for s in seasons) / 7.0 * rm.SFOC_AE / 1e6
    return sea + berth


def berth_generation(fields, p, seasons=None):
    """Energy at berth per unit in the actual hours of the timetable (seasons of the base case); t/yr for the units of
    the ship."""
    seasons = rm.BASE_SEASONS if seasons is None else seasons
    stays = [st for st in schedule(BERTHS) if st["season"] in seasons]
    weeks = sum(rm.DAYS[s] for s in seasons) / 7.0
    curves = {"vawt": (None, None, "u55"), "sav": (p.ds * rm.A_F, p.cp_s, "u40"),
              "sav3": (p.ds * (rm.A_F + rm.A_X), p.cp_s, "u40")}
    res = {k: {} for k in curves}
    for st in stays:
        fd = fields[st["name"]]
        e = {k: 0.0 for k in curves}
        for he, lo, hi in era5_slots(st["t1"], st["t2"]):
            for k, (area, cp, z) in curves.items():
                u = fd[z][he]
                if np.isfinite(u):
                    pw = rm.vawt_power_curve(np.array([u]), p)[0] if k == "vawt" else \
                        rm.power_curve(area, cp, np.array([u]))[0]
                    e[k] += (hi - lo) * float(pw)
        for k in curves:
            res[k].setdefault((st["name"], st["tag"]), []).append(e[k])
    out = {}
    for k in curves:
        kwh_year = sum(np.mean(v) * weeks for v in res[k].values())
        out[k] = dict(kW_mean_unit=kwh_year / (sum(rm.PORT_HOURS_WEEK.values()) * weeks),
                      t_yr=p.n_units * kwh_year * rm.SFOC_AE / 1e6)
    return out


def stats(d):
    v = np.array(list(d.values()), float)
    return dict(mean=float(v.mean()), sd=float(v.std(ddof=1)), min=float(v.min()), max=float(v.max()),
                cv=float(v.std(ddof=1) / v.mean()) if v.mean() else float("nan"))


def run(path="../data/era5_hourly_route.csv.gz", p=None, taus=TAUS, quiet=False, forecast=True):
    p = p or rm.P()
    data, meta = eh.load_hourly(path)
    fields = wind_fields(data)
    cross = schedule(TIMETABLE)
    recs = evaluate(fields, cross, p)
    recs_fc = {lev: evaluate(fields, cross, p, forecast=lev) for lev in FORECAST} if forecast else None
    base = [r for r in recs if r["c"]["season"] in rm.BASE_SEASONS]
    res = {"crossings": len(recs), "crossings_base_case": len(base),
           "hours_sea": float(sum(r["a"][:, 0].sum() for r in base)),
           "hours_missing": int(sum(len(crossing_hours(r["c"])[2]) - len(r["a"]) for r in recs)),
           "R0_t_yr": r0_fuel(recs, p), "R0_t_yr_whole_year": r0_fuel(recs, p, seasons=rm.ALL_SEASONS),
           "R0_t_Apr_Nov": r0_fuel(recs, p, seasons=rm.NOV_SEASONS),
           "tau_h": list(taus), "L_h": list(LS), "H_h": list(HS), "tau_R1": list(TAU_R1), "tau_tilt": list(TAU_TILT),
           "forecast": FORECAST}
    per = {topo: strategies(recs, topo, taus, recs_fc=recs_fc) for topo in ("R1", "T1", "T2", "T3", "T3x", "T4", "T5")}
    r1 = per["R1"]["inst"]
    res["R1"] = {"inst": annualise(recs, r1), "Apr-Oct": annualise(recs, r1, by=0), "Nov": annualise(recs, r1, by=2),
                 "Dec-Mar": annualise(recs, r1, by=1), "Apr-Nov": annualise(recs, r1, seasons=rm.NOV_SEASONS),
                 "whole year": annualise(recs, r1, seasons=rm.ALL_SEASONS), "per_year": annualise(recs, r1, by="year")}
    res["R1"]["per_year_stats"] = stats(res["R1"]["per_year"])
    pen = [float((r["a"][:, 0] * r["a"][:, 10]).sum()) for r in recs]
    res["R1"]["induced"] = annualise(recs, pen)
    # start and stop cost of R1, and tiltable rotor
    st_r1 = strategies(recs, "R1 switched", taus=TAU_R1, recs_fc=recs_fc, mpc_taus=(TAU_R1[-1],))
    res["R1"]["switched"] = {k: (float(np.mean(v)) if k.startswith("switches") else annualise(recs, v))
                             for k, v in st_r1.items()}
    st_t = strategies(recs, "R1t", taus=TAU_TILT, recs_fc=recs_fc, mpc_taus=(TAU_TILT[0],))
    res["R1t"] = {k: (float(np.mean(v)) if k.startswith("switches") else annualise(recs, v) - res["R1"]["inst"])
                  for k, v in st_t.items()}
    res["sea_increments"] = {}
    for topo in ("T1", "T2", "T3", "T3x", "T4", "T5"):
        d = {}
        for k, v in per[topo].items():
            if k.startswith("switches"):
                d[k] = float(np.mean(v))
            else:
                d[k] = annualise(recs, v) - res["R1"]["inst"]
        d["Apr-Oct inst"] = annualise(recs, per[topo]["inst"], by=0) - res["R1"]["Apr-Oct"]
        d["Apr-Nov inst"] = annualise(recs, per[topo]["inst"], seasons=rm.NOV_SEASONS) - res["R1"]["Apr-Nov"]
        d["Dec-Mar inst"] = annualise(recs, per[topo]["inst"], by=1) - res["R1"]["Dec-Mar"]
        for k in ["inst"] + [f"optimal {t}" for t in taus]:
            yr = {y: v - res["R1"]["per_year"][y] for y, v in annualise(recs, per[topo][k], by="year").items()}
            d[f"per_year {k}"] = yr
            d[f"per_year {k} stats"] = stats(yr)
        res["sea_increments"][topo] = d
    res["berth"] = berth_generation(fields, p)
    res["era5_meta"] = meta[:3]
    if not quiet:
        report(res)
    return res


def report(res):
    print(f"Simulated crossings: {res['crossings']} ({res['crossings_base_case']} from April to October); hours at "
          f"sea {res['hours_sea']:.0f}; hours without data {res['hours_missing']}; R0 {res['R0_t_yr']:.0f} t/yr "
          f"(whole year {res['R0_t_yr_whole_year']:.0f})")
    r1 = res["R1"]
    st = r1["per_year_stats"]
    print(f"R1 (four rotors): {r1['inst']:.1f} t/yr (Apr-Oct {r1['Apr-Oct']:.1f}, Apr-Nov {r1['Apr-Nov']:.1f}, "
          f"Dec-Mar {r1['Dec-Mar']:.1f}, "
          f"whole year {r1['whole year']:.1f}); induced resistance paid {r1['induced']:.1f} t/yr")
    print("  per year: " + ", ".join(f"{y} {v:.0f}" for y, v in r1["per_year"].items()) +
          f"  | mean {st['mean']:.1f}, sd {st['sd']:.1f} (CV {100 * st['cv']:.1f} %), {st['min']:.0f}-{st['max']:.0f}")
    c = r1["switched"]
    print("  R1 with start and stop: " + ", ".join(f"{k} {v:.1f}" for k, v in c.items()
                                                  if not k.startswith("switches") and ("L" not in k or "L3" in k)))
    print("R1t (increment over R1): " + ", ".join(f"{k} {v:.1f}" for k, v in res["R1t"].items()
                                                if not k.startswith("switches") and ("L" not in k or "L3" in k)))
    for topo, d in res["sea_increments"].items():
        print(f"{topo:4s} inst {d['inst']:7.1f} (Apr-Nov {d['Apr-Nov inst']:.1f}, Dec-Mar {d['Dec-Mar inst']:.1f})")
        for t in res["tau_h"]:
            line = (f"     tau {t}: optimal {d[f'optimal {t}']:6.1f} ({d[f'switches optimal {t}']:.1f} changes)"
                    f"  greedy {d[f'greedy {t}']:6.1f} ({d[f'switches greedy {t}']:.1f})")
            for L in (3.0,):
                line += f"  causal L{L:g} {d[f'causal L{L:g} {t}']:6.1f} ({d[f'switches causal L{L:g} {t}']:.1f})"
            for lev in res["forecast"]:
                for H in res["H_h"]:
                    key = f"mpc {lev} H{H} {t}"
                    if key in d:
                        line += f"  mpc {lev} H{H} {d[key]:6.1f}"
            print(line)
        s = d["per_year inst stats"]
        print(f"     per year (inst): " + ", ".join(f"{y} {v:.1f}" for y, v in d["per_year inst"].items()) +
              f" | mean {s['mean']:.1f}, sd {s['sd']:.1f}, {s['min']:.1f}-{s['max']:.1f}")
        if topo == "T5":
            print("     rules: " + ", ".join(f"{k} {d[k]:.1f}" for k in d if k.startswith(("claim1", "intended"))))
    print("Berth (actual hours, April-October): " + ", ".join(
        f"{k} {v['kW_mean_unit']:.2f} kW/unit, {v['t_yr']:.1f} t/yr" for k, v in res["berth"].items()))


if __name__ == "__main__":
    out = run(sys.argv[1] if len(sys.argv) > 1 else "../data/era5_hourly_route.csv.gz")
    json.dump(out, open("hourly_results.json", "w"), indent=1, ensure_ascii=False, default=float)
