#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LaTeX tables of the supplementary material (section S8) from uncertainty.json (version 5, 27-10-2026): total and
first-order Sobol indices, percentiles and ranking probabilities by scenario.
Usage: python3 sobol_tables_tex.py > sobol_tables.tex
"""
import json

d = json.load(open("uncertainty.json", encoding="utf-8"))
LAB = {"cd_mult": r"Rotor $C_D$ multiplier", "cd0": r"Stopped-rotor $C_{D,0}$", "cf": r"Spin friction coefficient $C_f$",
       "cp_s": r"Savonius $C_{p,S}$", "u_ct": r"Savonius $C_{t,S}$, position in its range",
       "f_cp_t": r"Turbine power factor $f_{cp}$", "f_ct_t": r"Turbine thrust factor $f_{ct}$",
       "f_cs_t": r"Turbine side-force factor $f_{cs}$", "cl_ext": r"Extension lift factor $f_X$",
       "f30": r"Lift factor $f_{30}$", "eta_d": r"Propulsive efficiency $\eta_D$", "sfoc_me": r"SFOC of the main engines",
       "c_h": r"Hull term of the induced resistance $c_H$", "cpa": r"Aft-body prismatic coefficient $C_{pa}$",
       "cwa": r"Aft-body waterplane coefficient $C_{wa}$", "x_s": r"Centre of the lateral windage $x_s$",
       "ayc": r"Lateral windage $A_LC_{Y,\max}$", "axc": r"Frontal windage $A_{\mathrm{fr}}C_{X,\max}$"}
OUT = [("R1", "R1"), ("pen R1", "R1 ind."), ("T3x", "T3x"), ("T4", "T4"), ("T1", "T1"), ("T5", "T5"),
       ("R1-30", "R1-30"), ("R1t", "R1t"), ("R1-D5", "R1-D5")]
ESC = [("A", "A"), ("A-VG", "A-VG"), ("A+", "A+"), ("AB", "AB"), ("B", "B"), ("B'", "B$'$")]


def num(x, nd=2):
    s = f"{max(x, 0.0) if nd == 2 else x:.{nd}f}"
    return s.replace("-", "$-$")


def sobol_table(key, caption, label):
    keys = [k for k in LAB if k in d["sobol"]["R1"]]
    cols = "R{0.85cm}" * len(OUT)
    print(r"\begin{table}[htbp]" "\n" r"\centering")
    print(rf"\caption{{{caption}}}" "\n" rf"\label{{{label}}}" "\n" r"\footnotesize" "\n" r"\setlength{\tabcolsep}{3pt}")
    print(rf"\begin{{tabularx}}{{\textwidth}}{{@{{}}L{cols}@{{}}}}" "\n" r"\toprule")
    print("Input & " + " & ".join(h for _, h in OUT) + r" \\" "\n" r"\midrule")
    for k in keys:
        cells = []
        for o, _ in OUT:
            v = d["sobol"][o][k][key]
            s = num(v)
            cells.append(rf"\textbf{{{s}}}" if v >= 0.1 else s)
        print(LAB[k] + " & " + " & ".join(cells) + r" \\")
    print(r"\bottomrule" "\n" r"\end{tabularx}" "\n" r"\end{table}" "\n")


def g(x):
    a = abs(x)
    s = f"{x:.0f}" if a >= 9.95 else f"{x:.1f}"
    return s.replace("-", "$-$")


sobol_table("ST", r"Total Sobol indices $S_{Ti}$ within scenario A ($N=1{,}024$, 18 inputs), for the annual saving of R1 "
            r"relative to R0, the induced resistance it pays, and the savings of the topologies and comparators relative to "
            r"R1. Indices of 0.1 or more in bold.", "tab:sobolT")
sobol_table("S1", r"First-order Sobol indices $S_i$ within scenario A, for the same outputs as Table~\ref{tab:sobolT}. "
            r"Negative estimates, all within their confidence intervals of zero, are shown as 0.00. Indices of 0.1 or more "
            r"in bold.", "tab:sobolS1")

esc = [e for e in ESC if e[0] in d["scenarios"]]
print(r"\begin{table}[htbp]" "\n" r"\centering")
print(r"\caption{Median and, in brackets, 5th and 95th percentiles of the annual savings (t/yr) by coefficient scenario: R1 relative to "
      r"R0, with the induced resistance it pays; topologies and comparators relative to R1, including berth.}")
print(r"\label{tab:pct}" "\n" r"\scriptsize" "\n" r"\setlength{\tabcolsep}{2pt}")
print(r"\begin{tabularx}{\textwidth}{@{}L" + "R{1.95cm}" * len(esc) + "@{}}" "\n" r"\toprule")
print("Output & " + " & ".join(h for _, h in esc) + r" \\" "\n" r"\midrule")
for o, h in [("R1", "R1"), ("pen R1", "R1, induced resistance"), ("T3x", "T3x"), ("T3", "T3"), ("T5", "T5"),
             ("T4", "T4"), ("T2", "T2"), ("T1", "T1"), ("R1-30", "R1-30"), ("R1t", "R1t"), ("R1-D5", "R1-D5")]:
    cells = []
    for e, _ in esc:
        r = d["scenarios"][e]["summary"][o]
        cells.append(f"{g(r['p50'])} [{g(r['p5'])}, {g(r['p95'])}]")
    print(h + " & " + " & ".join(cells) + r" \\")
cells = []
for e, _ in esc:
    r = d["scenarios"][e]["summary"]["R1 %"]
    cells.append(f"{r['p50']:.1f} [{r['p5']:.1f}, {r['p95']:.1f}]".replace("-", "$-$"))
print(r"R1, share of R0 fuel (\%) & " + " & ".join(cells) + r" \\")
print(r"\bottomrule" "\n" r"\end{tabularx}" "\n" r"\end{table}" "\n")

PR = [("ranking T3 >= T3x > T5 > T4 > T2 > T1", r"T3 $\geq$ T3x $>$ T5 $>$ T4 $>$ T2 $>$ T1"),
      ("T3x > T4", r"T3x $>$ T4"), ("T1 < 0", r"T1 $<$ 0"), ("T2 < 0", r"T2 $<$ 0"),
      ("T3x > R1t", r"T3x $>$ R1t"), ("T3x > R1-30", r"T3x $>$ R1-30"), ("T3x > R1-D5", r"T3x $>$ R1-D5"),
      ("T3 > R1t", r"T3 $>$ R1t"), ("T3 > R1-30", r"T3 $>$ R1-30"), ("R1-30 > R1t", r"R1-30 $>$ R1t")]
print(r"\begin{table}[htbp]" "\n" r"\centering")
print(r"\caption{Probability of each ranking by coefficient scenario, with its 95\% bootstrap confidence interval.}")
print(r"\label{tab:prob}" "\n" r"\scriptsize" "\n" r"\setlength{\tabcolsep}{2pt}")
print(r"\begin{tabularx}{\textwidth}{@{}L" + "R{1.95cm}" * len(esc) + "@{}}" "\n" r"\toprule")
print("Ranking & " + " & ".join(h for _, h in esc) + r" \\" "\n" r"\midrule")
for k, h in PR:
    cells = []
    for e, _ in esc:
        p = d["scenarios"][e]["probabilities"][k]
        cells.append(f"{p['p']:.2f} ({p['ci95'][0]:.2f}--{p['ci95'][1]:.2f})")
    print(h + " & " + " & ".join(cells) + r" \\")
print(r"\bottomrule" "\n" r"\end{tabularx}" "\n" r"\end{table}")
