#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wind histograms conditioned on the timetable (version 5, 27-10-2026; simulated review v2, R3-2 and R1-4).

Up to version 4, the histogram model used all the hours of 2015-2024 at each route point (era5_route.json) and gave 6 %
more saving of R1 than the hour-by-hour simulation, because the timetable only passes each point at certain hours of
the day (hourly_diagnostics.py). Here only the hours in which the ship is at each point according to the 2026
timetable are aggregated (hourly_voyages.schedule and crossing_hours, with the same split of each hourly ERA5 value),
by type of crossing and month (that of the midpoint of
the crossing or stay), and the wind in port only in the actual berth hours. With them, the
histogram model reproduces the hour-by-hour simulation except for the effect of the bins (1 m/s x 15 degrees) and of
the mean marginal consumption, and it stays cheap for the sensitivity and uncertainty sweeps.

Format (era5_route_conditioned.json):
  crossings[name][month] = dict(h = open-sea hours (10 years), F = hours per bin (NU x NG, in rows of speed at 40 m
      and columns of the true wind angle relative to the heading towards Bilbao, folded to 0-180), rs, rc = sum of
      u55/u40 and number of hours per speed row, n = number of crossings);
  port_months[port][month] = dict(hr, ht = hours per speed row at 40 and 55 m, h = berth hours, n = calls).
Usage: python3 era5_conditioned.py [../data/era5_hourly_route.csv.gz]
"""
import json
import math
import os
import sys
import numpy as np
import era5_hourly as eh
import hourly_voyages as rh

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, "era5_route_conditioned.json")
NU, NG = eh.NU, eh.NG


def mid_month(c):
    """Month (UTC) of the midpoint of a crossing or berth stay: the month under which its hours are filed."""
    return (rh.T0 + rh.dt.timedelta(hours=0.5 * (c["t1"] + c["t2"]))).month


def aggregate(fields):
    trav, ports = {}, {}
    for c in rh.schedule(rh.TIMETABLE):
        _, _, hs = rh.crossing_hours(c)
        m = mid_month(c)
        e = trav.setdefault(c["name"], {}).setdefault(m, dict(h=0.0, F=np.zeros((NU, NG)), rs=np.zeros(NU),
                                                              rc=np.zeros(NU), n=0))
        e["n"] += 1
        for w, ip, he in hs:
            pid, psi = rh.ROUTE[ip][0], rh.ROUTE[ip][4]
            fd = fields[pid]
            u, r, d100 = fd["u40"][he], fd["ratio"][he], fd["d100"][he]
            if not (np.isfinite(u) and np.isfinite(r) and np.isfinite(d100)):
                continue
            iu = min(NU - 1, int(math.floor(u)))
            ig = min(NG - 1, int(math.floor(float(eh.fold(d100 - psi)) / 15.0)))
            e["F"][iu, ig] += w
            e["h"] += w
            if u > 0.1:
                e["rs"][iu] += w * r
                e["rc"][iu] += w
    for st in rh.schedule(rh.BERTHS):
        m = mid_month(st)
        fd = fields[st["name"]]
        e = ports.setdefault(st["name"], {}).setdefault(m, dict(hr=np.zeros(NU), ht=np.zeros(NU), h=0.0, n=0))
        e["n"] += 1
        for he, lo, hi in rh.era5_slots(st["t1"], st["t2"]):
            w = hi - lo
            for z, key in (("u40", "hr"), ("u55", "ht")):
                u = fd[z][he]
                if np.isfinite(u):
                    e[key][min(NU - 1, int(math.floor(u)))] += w
            e["h"] += w
    js = lambda e: {k: (np.round(v, 6).tolist() if isinstance(v, np.ndarray) else v) for k, v in e.items()}
    return dict(meta=dict(NU=NU, NG=NG, zr=eh.ZR, zt=eh.ZT, years="2015-2024", uBin=1, gammaBinDeg=15,
                          source="era5_conditioned.py: hours of the 2026 timetable (hourly_voyages.TIMETABLE, "
                                 "BERTHS) over hourly ERA5 2015-2024"),
                crossings={n: {str(m): js(e) for m, e in sorted(d.items())} for n, d in trav.items()},
                port_months={n: {str(m): js(e) for m, e in sorted(d.items())} for n, d in ports.items()})


def main(path):
    data, meta = eh.load_hourly(path)
    out = aggregate(rh.wind_fields(data))
    json.dump(out, open(OUTPUT, "w"), ensure_ascii=False)
    for n, d in out["crossings"].items():
        print(f"{n:28s} crossings per month: " + " ".join(f"{m}:{v['n']}" for m, v in d.items()))
    for n, d in out["port_months"].items():
        print(f"{n:12s} berth hours (10 years): {sum(v['h'] for v in d.values()):.0f}")
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "data", "era5_hourly_route.csv.gz"))
