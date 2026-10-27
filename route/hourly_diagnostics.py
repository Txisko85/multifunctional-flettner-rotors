#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A1 - Diagnostic of the histogram model against the timetable (version 5, 27-10-2026; version 1: 19-10-2026).
Compares, with the same physical model (route_model.bin_terms):
  all hours        - histograms of all the hours of 2015-2024 at each point (era5_route.json, version 4; only
                     April-November and December-March);
  conditioned      - histograms of the hours of the timetable (era5_route_conditioned.json, version 5);
  hour by hour     - simulation of each crossing with the exact hourly wind (hourly_voyages.py).
The difference between the first two is the effect of the hours of the day sampled by the timetable; what remains
between the conditioned histograms and the hour-by-hour simulation is that of the bins (1 m/s x 15 degrees) and of the
mean marginal consumption.
Usage: python3 hourly_diagnostics.py [../data/era5_hourly_route.csv.gz]  ->  hourly_diagnostics.json
"""
import json
import sys
import numpy as np
import era5_hourly as eh
import hourly_voyages as rh
import route_model as rm

TOPOS = ("R1", "T1", "T3x", "T4", "R1t")


def hourly(recs, seasons):
    per = {t: rh.strategies(recs, t, taus=())["inst"] for t in ("R1", "T1", "T3x", "T4")}
    out = {"R1": rh.annualise(recs, per["R1"], seasons=seasons)}
    for t in ("T1", "T3x", "T4"):
        out[t] = rh.annualise(recs, per[t], seasons=seasons) - out["R1"]
    st = rh.strategies(recs, "R1t", taus=())["inst"]
    out["R1t"] = rh.annualise(recs, st, seasons=seasons) - out["R1"]
    return out


def hist(d, g, nu, ng, seasons):
    r = rm.annual(d, g, nu, ng, rm.P(), seasons=seasons)
    tot = {row[0]: row[1] for row in rm.table(r)}
    return {"R1": r["sea"]["R1"], "T1": tot["T1"], "T3x": tot["T3x"], "T4": tot["T4"], "R1t": r["sea"]["R1t"]}


def main(path):
    data, _ = eh.load_hourly(path)
    fields = rh.wind_fields(data)
    cross = rh.schedule(rh.TIMETABLE)
    p = rm.P()
    hod = np.zeros(24)
    for c in cross:
        for w, ip, he in rh.crossing_hours(c)[2]:
            hod[he % 24] += w
    hod /= hod.sum()
    recs = rh.evaluate(fields, cross, p)
    d4, g4, nu, ng = rm.load_wind("era5_route.json")
    d5, g5, _, _ = rm.load_wind()
    res = {"hours_UTC_fraction": [round(float(x), 4) for x in hod],
           "Apr-Nov": {"all hours": hist(d4, g4, nu, ng, (0,)), "conditioned": hist(d5, g5, nu, ng, rm.NOV_SEASONS),
                       "hour by hour": hourly(recs, rm.NOV_SEASONS)},
           "Apr-Oct": {"conditioned": hist(d5, g5, nu, ng, rm.BASE_SEASONS),
                       "hour by hour": hourly(recs, rm.BASE_SEASONS)}}
    print("Fraction of sea hours by UTC hour (%):", " ".join(f"{100 * x:.1f}" for x in hod))
    for per_, dd in res.items():
        if per_.startswith("hours"):
            continue
        print(per_)
        for lab, v in dd.items():
            print(f"  {lab:16s} " + "  ".join(f"{k} {v[k]:7.1f}" for k in TOPOS))
    json.dump(res, open("hourly_diagnostics.json", "w"), indent=1, ensure_ascii=False)
    return res


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "../data/era5_hourly_route.csv.gz")
