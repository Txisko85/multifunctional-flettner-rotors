#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A1 - Diagnóstico del modelo de histogramas frente al horario (versión 5, 27-10-2026; versión 1: 19-10-2026).
Compara, con el mismo modelo físico (ruta_modelo.bin_terms):
  todas las horas  - histogramas de todas las horas de 2015-2024 de cada punto (era5_ruta.json, versión 4; solo
                     abril-noviembre y diciembre-marzo);
  condicionados    - histogramas de las horas del horario (era5_ruta_horario.json, versión 5);
  hora a hora      - simulación de cada travesía con el viento horario exacto (ruta_horaria.py).
La diferencia entre los dos primeros es el efecto de las horas del día que muestrea el horario; la que queda entre
los condicionados y la simulación hora a hora es la de las casillas (1 m/s x 15 grados) y la del consumo marginal medio.
Uso: python3 ruta_horaria_diag.py [../datos/era5_horario_ruta.csv.gz]  ->  diagnostico_horario.json
"""
import json
import sys
import numpy as np
import era5_horario as eh
import ruta_horaria as rh
import ruta_modelo as rm

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
    cross = rh.schedule(rh.HORARIO)
    p = rm.P()
    hod = np.zeros(24)
    for c in cross:
        for w, ip, he in rh.crossing_hours(c)[2]:
            hod[he % 24] += w
    hod /= hod.sum()
    recs = rh.evaluate(fields, cross, p)
    d4, g4, nu, ng = rm.load_wind("era5_ruta.json")
    d5, g5, _, _ = rm.load_wind()
    res = {"horas_UTC_fraccion": [round(float(x), 4) for x in hod],
           "abr-nov": {"todas las horas": hist(d4, g4, nu, ng, (0,)), "condicionados": hist(d5, g5, nu, ng, rm.NOV_SEASONS),
                       "hora a hora": hourly(recs, rm.NOV_SEASONS)},
           "abr-oct": {"condicionados": hist(d5, g5, nu, ng, rm.BASE_SEASONS),
                       "hora a hora": hourly(recs, rm.BASE_SEASONS)}}
    print("Fracción de horas de mar por hora UTC (%):", " ".join(f"{100 * x:.1f}" for x in hod))
    for per_, dd in res.items():
        if per_.startswith("horas"):
            continue
        print(per_)
        for lab, v in dd.items():
            print(f"  {lab:16s} " + "  ".join(f"{k} {v[k]:7.1f}" for k in TOPOS))
    json.dump(res, open("diagnostico_horario.json", "w"), indent=1, ensure_ascii=False)
    return res


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "era5_horario_ruta.csv.gz")
