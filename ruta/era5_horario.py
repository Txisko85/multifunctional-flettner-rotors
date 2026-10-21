#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lectura del ERA5 horario de la ruta (era5_horario_ruta.csv.gz, generado con descargar_era5_horario.py) y
agregación en histogramas idéntica a la del navegador (era5_agregacion_navegador.js) del 18-10-2026.

Uso: python3 era5_horario.py [era5_horario_ruta.csv.gz]
     Agrega el fichero, escribe era5_ruta_desde_horario.json y lo compara con era5_ruta.json.
"""
import csv
import gzip
import json
import math
import sys
import datetime as dt
import numpy as np

T0 = dt.datetime(2015, 1, 1, tzinfo=dt.timezone.utc)
N_HORAS = 87672
ZR, ZT, NU, NG = 40.0, 55.0, 30, 12
# id, grupo, latitud, longitud, rumbo hacia Bilbao (grados), millas que representa (como en el navegador)
PUNTOS = [
    ("C1", "C", 50.362, -1.850, 249.5, 45.36), ("C2", "C", 50.088, -2.950, 249.5, 45.36),
    ("C3", "C", 49.625, -4.062, 229.3, 58.84), ("C4", "C", 48.975, -5.188, 229.3, 58.84),
    ("B1", "B", 48.217, -5.527, 159.4, 55.23), ("B2", "B", 47.350, -5.082, 159.4, 55.23),
    ("B3", "B", 46.483, -4.638, 159.4, 55.23), ("B4", "B", 45.617, -4.192, 159.4, 55.23),
    ("B5", "B", 44.750, -3.748, 159.4, 55.23), ("B6", "B", 43.883, -3.303, 159.4, 55.23),
    ("PORTSMOUTH", "P", 50.790, -1.110, None, 0.0), ("BILBAO", "P", 43.360, -3.030, None, 0.0),
]


def load_hourly(path="era5_horario_ruta.csv.gz"):
    """Devuelve {punto: dict(u10, u100, d100)} con vectores de N_HORAS (NaN si falta) y los metadatos '#'."""
    data = {p[0]: {k: np.full(N_HORAS, np.nan) for k in ("u10", "u100", "d100")} for p in PUNTOS}
    meta = []
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                meta.append(line.strip())
                continue
            if line.startswith("punto,"):
                continue
            pid, h, a, b, c = line.rstrip("\n").split(",")
            h = int(h)
            if pid not in data or not 0 <= h < N_HORAS:
                continue
            dd = data[pid]
            if a:
                dd["u10"][h] = float(a)
            if b:
                dd["u100"][h] = float(b)
            if c:
                dd["d100"][h] = float(c)
    return data, meta


def shear(u10, u100):
    """Exponente horario de la ley potencial entre 10 y 100 m, como en el navegador (0,11 con viento flojo)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        al = np.where((u10 > 0.5) & (u100 > 0.5), np.log10(u100 / u10), 0.11)
    return np.clip(al, 0.0, 0.6)


def heights(dd):
    al = shear(dd["u10"], dd["u100"])
    return dd["u100"] * (ZR / 100.0) ** al, dd["u100"] * (ZT / 100.0) ** al, al


def fold(diff):
    g = np.mod(diff, 360.0)
    return np.where(g > 180.0, 360.0 - g, g)


def months():
    return np.array([(T0 + dt.timedelta(hours=int(h))).month for h in range(N_HORAS)])


def aggregate(data):
    """Histogramas por grupo (Canal, golfo) y estación, y viento en atraque, con el formato de era5_ruta.json."""
    mo = months()
    season = np.where((mo >= 4) & (mo <= 11), 0, 1)
    raw, pts = {}, {}
    for pid, g, lat, lon, psi, w in PUNTOS:
        dd = data[pid]
        ok = ~(np.isnan(dd["u10"]) | np.isnan(dd["u100"]) | np.isnan(dd["d100"]))
        ur, ut, al = heights(dd)
        iu = np.minimum(NU - 1, np.floor(np.nan_to_num(ur)).astype(int))
        it = np.minimum(NU - 1, np.floor(np.nan_to_num(ut)).astype(int))
        r = dict(cnt=[np.zeros(NU * NG), np.zeros(NU * NG)], ns=[0, 0], rs=[np.zeros(NU), np.zeros(NU)],
                 rc=[np.zeros(NU), np.zeros(NU)], hr=[np.zeros(NU, int), np.zeros(NU, int)],
                 ht=[np.zeros(NU, int), np.zeros(NU, int)])
        for s in (0, 1):
            m = ok & (season == s)
            r["ns"][s] = int(m.sum())
            r["hr"][s] = np.bincount(iu[m], minlength=NU)
            r["ht"][s] = np.bincount(it[m], minlength=NU)
            mr = m & (ur > 0.1)
            r["rs"][s] = np.bincount(iu[mr], weights=(ut / ur)[mr], minlength=NU)
            r["rc"][s] = np.bincount(iu[mr], minlength=NU).astype(float)
            if psi is not None:
                ig = np.minimum(NG - 1, np.floor(fold(dd["d100"] - psi) / 15.0).astype(int))
                r["cnt"][s] = np.bincount((iu * NG + ig)[m], minlength=NU * NG).astype(float)
        raw[pid] = r
        pts[pid] = dict(n=N_HORAS, miss=int((~ok).sum()), u100=round(float(np.nanmean(dd["u100"][ok])), 3),
                        u10=round(float(np.nanmean(dd["u10"][ok])), 3), alpha=round(float(al[ok].mean()), 4), w=w)
    out = dict(meta=dict(zr=ZR, zt=ZT, NU=NU, NG=NG, years="2015-2024", seasons=["Apr-Nov", "Dec-Mar"],
                         gammaBinDeg=15, uBin=1, source="era5_horario.py"), pts=pts, groups={}, ports={})
    for g in ("C", "B"):
        ids = [p for p in PUNTOS if p[1] == g]
        wsum = sum(p[5] for p in ids)
        out["groups"][g] = dict(wnm=wsum, F=[[], []], ratio=[[], []])
        for s in (0, 1):
            F, RS, RC = np.zeros(NU * NG), np.zeros(NU), np.zeros(NU)
            for p in ids:
                rr, wn = raw[p[0]], p[5] / wsum
                F += wn * rr["cnt"][s] / rr["ns"][s]
                RS += wn * rr["rs"][s]
                RC += wn * rr["rc"][s]
            out["groups"][g]["F"][s] = [int(round(v * 1e6)) for v in F]
            out["groups"][g]["ratio"][s] = [round(float(a / b), 4) if b > 0 else None for a, b in zip(RS, RC)]
    for pid in ("PORTSMOUTH", "BILBAO"):
        rr = raw[pid]
        out["ports"][pid] = dict(hr=[rr["hr"][s].tolist() for s in (0, 1)], ht=[rr["ht"][s].tolist() for s in (0, 1)],
                                 ns=rr["ns"])
    return out


def compare(a, b):
    """Diferencias entre dos agregaciones (máxima diferencia de F en ppm, de ratio y de los puertos)."""
    rep = {}
    for g in ("C", "B"):
        for s in (0, 1):
            fa, fb = np.array(a["groups"][g]["F"][s]), np.array(b["groups"][g]["F"][s])
            ra = np.array([np.nan if v is None else v for v in a["groups"][g]["ratio"][s]], float)
            rb = np.array([np.nan if v is None else v for v in b["groups"][g]["ratio"][s]], float)
            rep[f"{g}{s} F max |dif| ppm"] = int(np.abs(fa - fb).max())
            rep[f"{g}{s} ratio max |dif|"] = float(np.nanmax(np.abs(ra - rb)))
    for pid in ("PORTSMOUTH", "BILBAO"):
        for k in ("hr", "ht"):
            rep[f"{pid} {k} max |dif| h"] = int(max(np.abs(np.array(a["ports"][pid][k][s]) -
                                                         np.array(b["ports"][pid][k][s])).max() for s in (0, 1)))
        rep[f"{pid} ns"] = (a["ports"][pid]["ns"], b["ports"][pid]["ns"])
    return rep


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "era5_horario_ruta.csv.gz"
    data, meta = load_hourly(path)
    for m in meta:
        print(m)
    agg = aggregate(data)
    json.dump(agg, open("era5_ruta_desde_horario.json", "w"), indent=1, ensure_ascii=False)
    ref = json.load(open("era5_ruta.json", encoding="utf-8"))
    print("Comparación con era5_ruta.json (agregación del navegador, 18-10-2026):")
    for k, v in compare(agg, ref).items():
        print(f"  {k}: {v}")
    print("Medias por punto (u10, u100, alfa): nuevo / navegador")
    for pid in agg["pts"]:
        a, b = agg["pts"][pid], ref["pts"][pid]
        print(f"  {pid:10s} {a['u10']:.3f}/{b['u10']:.3f}  {a['u100']:.3f}/{b['u100']:.3f}  {a['alpha']:.4f}/{b['alpha']:.4f}"
              f"  faltan {a['miss']}")
