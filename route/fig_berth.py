#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure of generation in port (A1, section 4.4): (a) electrical power curves per unit of the VAWT, of the rotor opened
as a Savonius and of the two cylinders as a Savonius (T3); (b) distribution of the wind at berth (ERA5 2015-2024) in
Portsmouth and Bilbao at 40 m and 55 m, with the mean power per unit.
Same assumptions as route_model.port_power: cut-in 3 m/s, rated 12 m/s, cut-out 20 m/s, eta_gen*eta_conv; from
21-10-2026 the curve of the VAWT comes from the DMST model (route_model.vawt_power_curve), with the tip speed limited
to 40 m/s. Version 4 (26-10-2026): only the seasons of the base case (rm.BASE_SEASONS). Version 5 (27-10-2026):
April-October and wind only in the berth hours of the timetable (conditioned histograms).
Writes berth_generation.pdf/png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import route_model as rm
from matplotlib import font_manager

# Typeface of the figures of the article (Inter); if not installed, the Matplotlib default.
if any(f.name == "Inter" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Inter"

HERE = os.path.dirname(os.path.abspath(__file__))


def power_curve(area, cp, u):
    p_r = 0.5 * rm.RHO * area * cp * 12.0 ** 3
    pw = np.where((u >= 3.0) & (u < 20.0), np.minimum(0.5 * rm.RHO * area * cp * u ** 3, p_r), 0.0)
    return pw * rm.ETA_GEN * rm.ETA_CONV / 1e3


def curve(key, area, cp, u, p):
    return rm.vawt_power_curve(u, p) if key == "vawt" else power_curve(area, cp, u)


def summary(d, nu, p, seasons=rm.BASE_SEASONS):
    pp = rm.port_power(d, nu, p)
    days = sum(rm.DAYS[s] for s in seasons)
    w_s = {s: rm.DAYS[s] / days for s in seasons}
    per_port = {prt: {k: sum(w_s[s] * pp[prt][s][k] for s in seasons) for k in ("vawt", "sav", "sav3")}
                for prt in pp}
    h = rm.PORT_HOURS_WEEK
    w_p = {prt: h[prt] / sum(h.values()) for prt in h}
    overall = {k: sum(w_p[prt] * per_port[prt][k] for prt in h) for k in ("vawt", "sav", "sav3")}
    hours = sum(h.values()) * days / 7.0
    return pp, per_port, overall, hours


def main(seasons=rm.BASE_SEASONS):
    d, groups, nu, ng = rm.load_wind()
    p = rm.P()
    pp, per_port, overall, hours = summary(d, nu, p, seasons)
    devices = (("vawt", "VAWT\n(T1, T2, T4), 55 m", rm.A_X, p.cp_t, "-", 0.8, "bottom"),
               ("sav", "Rotor as Savonius\n(T2, T5), 40 m", p.ds * rm.A_F, p.cp_s, "--", -0.8, "top"),
               ("sav3", "Rotor + extension as Savonius\n(T3), 40 m", p.ds * (rm.A_F + rm.A_X), p.cp_s, ":",
                0.8, "bottom"))
    u = np.linspace(0.0, 25.0, 501)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.48, 3.0), constrained_layout=True)
    for key, label, area, cp, ls, dy, va in devices:
        pc = curve(key, area, cp, u, p)
        ax1.plot(u, pc, ls, color="k", lw=1.3)
        p_r = curve(key, area, cp, np.array([15.0]), p)[0]
        ax1.text(19.7, p_r + dy, label, ha="right", va=va, fontsize=6.8)
    ax1.set_xlabel("Wind speed at device height (m/s)")
    ax1.set_ylabel("Electrical power per unit (kW)")
    ax1.set_xlim(0, 25)
    ax1.set_ylim(0, 1.22 * max(curve(k, a, c, np.array([12.0]), p)[0] for k, _, a, c, *_ in devices))
    ax1.set_title("(a)", loc="left", fontsize=9)
    colors = {"PORTSMOUTH": "#2a78d6", "BILBAO": "#eb6834"}
    names = {"PORTSMOUTH": "Portsmouth", "BILBAO": "Bilbao"}
    edges = np.arange(nu + 1)
    for prt in ("PORTSMOUTH", "BILBAO"):
        pdat = d["ports"][prt]
        n = sum(pdat["ns"][s] for s in seasons)
        for key, ls, z in (("hr", "-", "40 m"), ("ht", "--", "55 m")):
            f = sum(np.array(pdat[key][s], float) for s in seasons) / n
            ax2.step(edges[:-1], 100 * f, where="post", ls=ls, color=colors[prt], lw=1.1,
                     label=f"{names[prt]}, {z}")
    ax2.set_xlabel("Hourly wind speed at berth (m/s)")
    ax2.set_ylabel("Frequency (% of hours)")
    ax2.set_xlim(0, 25)
    ax2.set_ylim(0, None)
    txt = "\n".join(f"{names[p]}: {per_port[p]['vawt']:.1f} / {per_port[p]['sav']:.1f} / {per_port[p]['sav3']:.1f} kW"
                    for p in ("PORTSMOUTH", "BILBAO"))
    ax2.text(0.98, 0.62, "Mean power per unit\nVAWT / Savonius / T3\n" + txt, transform=ax2.transAxes,
             ha="right", va="top", fontsize=6.5)
    ax2.legend(frameon=False, fontsize=7, loc="upper right")
    ax2.set_title("(b)", loc="left", fontsize=9)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"berth_generation.{ext}"), dpi=300)
    print("Mean power per unit (kW), weighted by season:")
    for prt in per_port:
        print(f"  {prt:10s} " + "  ".join(f"{k} {v:.3f}" for k, v in per_port[prt].items()))
    print("  overall   " + "  ".join(f"{k} {v:.3f}" for k, v in overall.items()), f"(berth hours {hours:.1f} h/yr)")
    for k, v in overall.items():
        e = rm.N_UNITS * v * hours / 1e3
        print(f"  4 units {k}: {rm.N_UNITS*v:.2f} kW mean, {e:.1f} MWh/yr, "
              f"{e*rm.SFOC_AE/1e3:.2f} t/yr at {rm.SFOC_AE:.0f} g/kWh "
              f"(190-220 g/kWh: {e*0.190:.2f}-{e*0.220:.2f} t), covers {100*rm.N_UNITS*v/rm.P_HOTEL_BERTH:.2f} % of {rm.P_HOTEL_BERTH:.0f} kW")
    print("  rated power per unit (12 m/s): " + "  ".join(
        f"{k} {curve(k, a, c, np.array([12.0]), p)[0]:.1f} kW" for k, _, a, c, *_ in devices))


if __name__ == "__main__":
    main()
