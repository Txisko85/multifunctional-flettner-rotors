#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check of the model of the ship without rotors (route_model.py, case R0) against the public EU MRV data of the Galicia
(IMO 9856189) downloaded from EMSA THETIS-MRV on 18-10-2026 (mrv_galicia.json), and calibration of the sea margin.
Version 4 (22-10-2026). Version 5 (27-10-2026): whole-year wind from the histograms conditioned on the timetable
(era5_conditioned.py) and induced resistance from the sway-yaw balance (route_model.induced).

Scope of EU MRV for a UK-Spain ferry after Brexit: it includes all the crossings (each one departs from or arrives at
an EU port) and the berth in EU ports, but not the berth in UK ports. The years 2021-2023 are used, with an almost
exclusive UK-EU service (emissions between EU ports < 2 % of the total):
  1. navigation: fuel outside the berth, hours at sea, distance and mean speed;
  2. berth: the fuel at berth in the EU divided by all the hours off the sea gives a LOWER BOUND of the mean
     consumption at berth.
Check at equal speed (version 4; up to version 3 the hourly consumption of the model at its own mean speed, 17.4 kn,
was compared with that of EU MRV at 16.6 kn): the speeds of the four crossings of the 2026 timetable are scaled by a
common factor until the mean speed of the model (distance over hours at sea, with the port approaches) equals that of
EU MRV 2021-2023, with the whole-year wind. The sea margin mu, which multiplies the calm-water power by (1 + mu), is
calibrated so that the hourly consumption matches. It is calibrated for the power-speed curves of exponent 3 (base
case) and 4 (sensitivity); the result is copied into route_model.SEA_MARGIN, and the script warns if it does not
match.
Output: mrv_calibration.json.
"""
import json
import os
import route_model as rm

HERE = os.path.dirname(os.path.abspath(__file__))
YEARS_UK_EU = ("2021", "2022", "2023")
SEASONS_YEAR = rm.ALL_SEASONS          # version 5: April-October, December-March and November


def load_mrv():
    with open(os.path.join(HERE, "mrv_galicia.json")) as f:
        return json.load(f)


def per_year(rows):
    out = {}
    for y, r in rows.items():
        ef = r["co2"] / r["fuel"]
        f_berth = r["berth"] / ef
        h_year = 8784.0 if y in ("2020", "2024") else 8760.0
        h_off = h_year - r["sea_h"]
        dist = r["fuel"] / r["fpd"] * 1e3            # nm (consumption per distance = total consumption / distance)
        f_sea = r["fuel"] - f_berth
        out[y] = dict(ef=ef, f_berth=f_berth, h_off=h_off, dist=dist, v_mean=dist / r["sea_h"], f_sea=f_sea,
                      rate_sea=f_sea / r["sea_h"], fpd_sea=f_sea / dist * 1e3, intra_eu=r["between"] / r["co2"],
                      berth_lb_kgph=f_berth * 1e3 / h_off, berth_lb_kw=f_berth * 1e6 / h_off / rm.SFOC_AE)
    return out


def pooled_mrv(rows, py):
    ys = YEARS_UK_EU
    pooled = dict(f_sea=sum(py[y]["f_sea"] for y in ys), sea_h=sum(rows[y]["sea_h"] for y in ys),
                  dist=sum(py[y]["dist"] for y in ys))
    pooled["rate_sea"] = pooled["f_sea"] / pooled["sea_h"]
    pooled["fpd_sea"] = pooled["f_sea"] / pooled["dist"] * 1e3
    pooled["v_mean"] = pooled["dist"] / pooled["sea_h"]
    return pooled


_PW = {}


def model_at(p, factor, wind):
    """Consumption of R0 at sea with the open-sea speeds of the timetable multiplied by 'factor', with the whole-year
    wind: t/h, kg/nm and mean speed (distance over hours at sea, with the approaches)."""
    d, groups, nu, ng = wind
    fuel = hours = dist = 0.0
    days = rm.season_days(groups)
    for s in (SEASONS_YEAR if rm.is_cond(groups) else (0, 1)):
        weeks = days[s] / 7.0
        for name, direction, h in rm.crossings(p):
            v0, _ = rm.crossing_speed(h, p)
            v = v0 * factor
            t_open = p.nm_open / v
            key = (s, name, round(v, 9), p.r0_wind, p.axc, p.ayc, p.e_span, p.eta_d, p.sir_model, p.x_s, p.c_h,
                   rm.is_cond(groups))
            if key not in _PW:
                _PW[key] = rm.r0_wind_power(groups, nu, ng, s, direction, v, p, name=name)
            f, _, _ = rm.crossing_fuel(v, t_open, _PW[key], p)
            fuel += f * weeks
            hours += (t_open + rm.NM_APPROACH / rm.V_APPROACH) * weeks
            dist += (p.nm_open + rm.NM_APPROACH) * weeks
    return dict(rate=fuel / hours, fpd=fuel / dist * 1e3, v_mean=dist / hours)


def speed_factor(p, v_target, wind):
    lo, hi = 0.5, 1.5
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if model_at(rm.P(**{**vars(p), "r0_wind": False}), mid, wind)["v_mean"] < v_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def calibrate(n_prop, v_target, rate_target, wind):
    """Sea margin mu with which the hourly consumption of the model at the EU MRV speed equals the measured one."""
    p0 = rm.P(n_prop=n_prop, sea_margin=0.0)
    f = speed_factor(p0, v_target, wind)
    lo, hi = -0.3, 0.6
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if model_at(rm.P(n_prop=n_prop, sea_margin=mid), f, wind)["rate"] < rate_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), f


def main():
    mrv = load_mrv()
    rows = mrv["rows"]
    py = per_year(rows)
    pooled = pooled_mrv(rows, py)
    print("EU MRV, GALICIA (IMO 9856189)")
    print(f"{'year':5s} {'fuel t':>8s} {'CO2/fuel':>8s} {'intraEU':>7s} {'h sea':>7s} {'nm':>8s} {'kn':>6s} "
          f"{'t/h sea':>8s} {'kg/nm':>7s} {'brth EU':>8s} {'bound kg/h':>9s}")
    for y in sorted(rows):
        r, q = rows[y], py[y]
        print(f"{y:5s} {r['fuel']:8.0f} {q['ef']:8.3f} {100*q['intra_eu']:6.1f}% {r['sea_h']:7.0f} {q['dist']:8.0f} "
              f"{q['v_mean']:6.2f} {q['rate_sea']:8.3f} {q['fpd_sea']:7.1f} {q['f_berth']:8.0f} {q['berth_lb_kgph']:9.0f}")
    ef_all = sum(rows[y]["co2"] for y in rows if y != "2020") / sum(rows[y]["fuel"] for y in rows if y != "2020")
    print(f"\n2021-2023 pooled: {pooled['rate_sea']:.3f} t/h at sea, {pooled['fpd_sea']:.1f} kg/nm, "
          f"{pooled['v_mean']:.2f} kn, {pooled['sea_h']/3:.0f} h/yr at sea")
    print(f"CO2/fuel factor of the ship 2021-2025: {ef_all:.4f} t/t")

    wind = rm.load_wind()
    out = dict(mrv_2021_2023=pooled, ef_2021_2025=ef_all, per_year=py, cases={})
    v_t, r_t = pooled["v_mean"], pooled["rate_sea"]
    # comparison as in version 3: the model at its own mean speed (whole year, without wind or SFOC curve)
    p3 = rm.P(v3=True)
    own = model_at(p3, 1.0, wind)
    out["cases"]["version 3, at its own speed"] = own
    print(f"\nVersion 3 at its own speed ({own['v_mean']:.2f} kn): {own['rate']:.3f} t/h "
          f"({100 * (own['rate'] / r_t - 1):+.1f} %), {own['fpd']:.1f} kg/nm ({100 * (own['fpd'] / pooled['fpd_sea'] - 1):+.1f} %)")
    # at equal speed
    steps = [("cubic law, no wind, constant SFOC (version 3)", rm.P(v3=True, nm_open=rm.NM_OPEN_TOTAL)),
             ("+ wind resistance of the ship", rm.P(sfoc_curve=False, sea_margin=0.0)),
             ("+ SFOC with load", rm.P(sea_margin=0.0)),
             ("exponent 4, with wind and SFOC with load", rm.P(n_prop=4.0, sea_margin=0.0))]
    for label, p in steps:
        f = speed_factor(p, v_t, wind)
        m = model_at(p, f, wind)
        m["speed_factor"] = f
        out["cases"][label + ", at equal speed"] = m
        print(f"{label:58s} at {m['v_mean']:.2f} kn: {m['rate']:.3f} t/h ({100 * (m['rate'] / r_t - 1):+.1f} %), "
              f"{m['fpd']:.1f} kg/nm")
    out["sea_margin"] = {}
    for n in (3.0, 4.0):
        mu, f = calibrate(n, v_t, r_t, wind)
        out["sea_margin"][str(n)] = dict(mu=mu, speed_factor=f)
        chk = model_at(rm.P(n_prop=n, sea_margin=mu), f, wind)
        print(f"Calibrated sea margin, exponent {n:g}: mu = {mu:.4f} -> {chk['rate']:.3f} t/h at {chk['v_mean']:.2f} kn")
        if abs(rm.SEA_MARGIN.get(n, -1.0) - round(mu, 3)) > 5e-4:
            print(f"  WARNING: route_model.SEA_MARGIN[{n}] = {rm.SEA_MARGIN.get(n)}; the calibrated value is {mu:.3f}")
    lb = [py[y]["berth_lb_kgph"] for y in ("2021", "2022")]
    out["berth_bound_kgph"] = lb
    print(f"\nberth: model {rm.P().hotel_berth * rm.SFOC_AE / 1e3:.0f} kg/h (version 3: "
          f"{rm.P_HOTEL_BERTH_V3 * rm.SFOC_AE / 1e3:.0f} kg/h); MRV lower bound 2021-2022 "
          f"{min(lb):.0f}-{max(lb):.0f} kg/h")
    json.dump(out, open(os.path.join(HERE, "mrv_calibration.json"), "w"), indent=1, ensure_ascii=False)
    return out


if __name__ == "__main__":
    main()
