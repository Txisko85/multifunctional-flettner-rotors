#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A1 - Diagnóstico del modelo horario frente al de histogramas (19-10-2026).
Separa la diferencia entre ambos en:
  muestreo  - el horario evalúa solo las horas en que el buque está en cada punto (hora del día, día de la semana);
  casillas  - el de histogramas evalúa cada hora en el centro de su casilla (1 m/s x 15°) y con la razón u55/u40
              media de la casilla; el horario usa los valores exactos.
Pasos: histogramas (todas las horas, casillas) -> horario con casillas (horas de travesía) -> horario exacto.
También comprueba que las horas UTC muestreadas se reparten por igual (pares e impares).
Uso: python3 ruta_horaria_diag.py [era5_horario_ruta.csv.gz]  ->  diagnostico_horario.json
"""
import json
import sys
import numpy as np
import era5_horario as eh
import ruta_horaria as rh
import ruta_modelo as rm


def sums(recs, topo_keys=("R1", "T1", "T3x", "T4")):
    out = {}
    per = {t: rh.strategies(recs, t, taus=())["inst"] for t in topo_keys}
    out["R1"] = rh.annualise(recs, per["R1"])
    out["R1 abr-nov"] = rh.annualise(recs, per["R1"], by=0)
    out["R1 dic-mar"] = rh.annualise(recs, per["R1"], by=1)
    for t in topo_keys[1:]:
        out[t] = rh.annualise(recs, per[t]) - out["R1"]
    return out


def main(path):
    data, _ = eh.load_hourly(path)
    fields = rh.wind_fields(data)
    cross = rh.schedule(rh.HORARIO)
    p = rm.P()
    # 1) horas UTC muestreadas
    hod = np.zeros(24)
    for c in cross:
        for w, ip, he in rh.crossing_hours(c)[2]:
            hod[he % 24] += w
    hod /= hod.sum()
    # 2) modelo de histogramas (referencia)
    d, groups, nu, ng = rm.load_wind()
    base = rm.annual(d, groups, nu, ng, p)
    base_s0 = rm.annual(d, groups, nu, ng, p, seasons=(0,))
    tot = {r[0]: r[1] for r in rm.table(base)}
    hist = {"R1": base["sea"]["R1"], "R1 abr-nov": base_s0["sea"]["R1"],
            "R1 dic-mar": base["sea"]["R1"] - base_s0["sea"]["R1"],
            "T1": tot["T1"], "T3x": tot["T3x"], "T4": tot["T4"]}
    # 3) horario con casillas y 4) horario exacto
    binned = sums(rh.evaluate(fields, cross, p, binned=groups))
    exact = sums(rh.evaluate(fields, cross, p))
    # 5) horario desplazado 0-23 h, con casillas (bin_terms en caché: las entradas son discretas)
    orig = rm.bin_terms
    cache = {}

    def cached(u, r, g, v, pp):
        key = (u, r, g, round(v, 9))
        if key not in cache:
            cache[key] = orig(u, r, g, v, pp)
        return cache[key]
    rm.bin_terms = cached
    shifts = {}
    for sh in range(24):
        shifts[sh] = sums(rh.evaluate(fields, rh.schedule(rh.HORARIO, shift=float(sh)), p, binned=groups))
    # 6) horario desplazado a lo largo de toda la semana (0-162 h cada 6 h): debe reproducir los histogramas
    week = {}
    for sh in range(0, 168, 6):
        week[sh] = sums(rh.evaluate(fields, rh.schedule(rh.HORARIO, shift=float(sh)), p, binned=groups))
    rm.bin_terms = orig
    mean_wk = {k: float(np.mean([s[k] for s in week.values()])) for k in hist}
    mean_sh = {k: float(np.mean([s[k] for s in shifts.values()])) for k in hist}
    res = {"horas_UTC_fraccion": [round(float(x), 4) for x in hod],
           "histogramas": hist, "horario_casillas": binned, "horario_exacto": exact,
           "efecto_muestreo": {k: binned[k] - hist[k] for k in hist},
           "efecto_casillas": {k: exact[k] - binned[k] for k in hist},
           "desplazamientos_casillas": shifts, "media_24_desplazamientos": mean_sh,
           "desplazamientos_semana": week, "media_semana": mean_wk}
    print("Fracción de horas de mar por hora UTC (%):", " ".join(f"{100 * x:.1f}" for x in hod))
    print(f"{'':12s} {'histogramas':>12s} {'media sem.':>11s} {'media 24 h':>11s} {'h. casillas':>12s} "
          f"{'h. exacto':>12s} {'muestreo':>9s} {'casillas':>9s}")
    for k in hist:
        print(f"{k:12s} {hist[k]:12.1f} {mean_wk[k]:11.1f} {mean_sh[k]:11.1f} {binned[k]:12.1f} {exact[k]:12.1f} "
              f"{binned[k] - hist[k]:+9.1f} {exact[k] - binned[k]:+9.1f}")
    print("R1 según el desplazamiento semanal (h): " + ", ".join(f"{sh} {v['R1']:.0f}" for sh, v in week.items()))
    print("R1 según el desplazamiento del horario (h): " +
          ", ".join(f"{sh} {v['R1']:.0f}" for sh, v in shifts.items()))
    print("R1 abr-nov según el desplazamiento: " + ", ".join(f"{sh} {v['R1 abr-nov']:.0f}" for sh, v in shifts.items()))
    print("T3x según el desplazamiento: " + ", ".join(f"{sh} {v['T3x']:.1f}" for sh, v in shifts.items()))
    json.dump(res, open("diagnostico_horario.json", "w"), indent=1, ensure_ascii=False)
    return res


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "era5_horario_ruta.csv.gz")
