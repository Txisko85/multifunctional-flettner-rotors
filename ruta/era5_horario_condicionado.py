#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Histogramas de viento condicionados al horario (versión 5, 27-10-2026; revisión simulada v2, R3-2 y R1-4).

Hasta la versión 4, el modelo de histogramas usaba todas las horas de 2015-2024 de cada punto de la ruta (era5_ruta.json)
y daba un 6 % más de ahorro de R1 que la simulación hora a hora, porque el horario solo pasa por cada punto a ciertas
horas del día (ruta_horaria_diag.py). Aquí se agregan solo las horas en que el buque está en cada punto según el
horario de 2026 (ruta_horaria.schedule y crossing_hours, con el mismo reparto de cada valor horario de ERA5), por tipo
de travesía y mes de salida, y el viento en puerto solo en las horas reales de atraque. Con ellos, el modelo de histogramas
reproduce la simulación hora a hora salvo el efecto de las casillas (1 m/s x 15 grados) y del consumo marginal
medio, y sigue siendo barato para los barridos de sensibilidad e incertidumbre.

Formato (era5_ruta_horario.json):
  travesias[nombre][mes] = dict(h = horas de mar abierto (10 años), F = horas por casilla (NU x NG, en filas de
      velocidad a 40 m y columnas del ángulo del viento verdadero respecto al rumbo hacia Bilbao, plegado a 0-180),
      rs, rc = suma de u55/u40 y número de horas por fila de velocidad, n = número de travesías);
  puertos[puerto][mes] = dict(hr, ht = horas por fila de velocidad a 40 y 55 m, h = horas de atraque, n = estancias).
Uso: python3 era5_horario_condicionado.py [../datos/era5_horario_ruta.csv.gz]
"""
import json
import math
import os
import sys
import numpy as np
import era5_horario as eh
import ruta_horaria as rh

HERE = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(HERE, "era5_ruta_horario.json")
NU, NG = eh.NU, eh.NG


def aggregate(fields):
    trav, ports = {}, {}
    for c in rh.schedule(rh.HORARIO):
        _, _, hs = rh.crossing_hours(c)
        m = c["month"]
        e = trav.setdefault(c["name"], {}).setdefault(m, dict(h=0.0, F=np.zeros((NU, NG)), rs=np.zeros(NU),
                                                              rc=np.zeros(NU), n=0))
        e["n"] += 1
        for w, ip, he in hs:
            pid, psi = rh.RUTA[ip][0], rh.RUTA[ip][4]
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
    for st in rh.schedule(rh.ATRAQUES):
        m = st["month"]
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
                          fuente="era5_horario_condicionado.py: horas del horario de 2026 (ruta_horaria.HORARIO, "
                                 "ATRAQUES) sobre ERA5 horario 2015-2024"),
                travesias={n: {str(m): js(e) for m, e in sorted(d.items())} for n, d in trav.items()},
                puertos={n: {str(m): js(e) for m, e in sorted(d.items())} for n, d in ports.items()})


def main(path):
    data, meta = eh.load_hourly(path)
    out = aggregate(rh.wind_fields(data))
    json.dump(out, open(SALIDA, "w"), ensure_ascii=False)
    for n, d in out["travesias"].items():
        print(f"{n:28s} travesías por mes: " + " ".join(f"{m}:{v['n']}" for m, v in d.items()))
    for n, d in out["puertos"].items():
        print(f"{n:12s} horas de atraque (10 años): {sum(v['h'] for v in d.values()):.0f}")
    return out


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "datos", "era5_horario_ruta.csv.gz"))
