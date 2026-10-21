#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Monte Carlo de la incertidumbre de coeficientes (A1), versión 3 (21-10-2026). Muestreo uniforme e independiente
en los 14 rangos de coeficientes.RANGES: desde la versión 2 incluyen la eficacia e de la resistencia inducida y el
producto A_L C_Y,max del viento sobre el buque, y desde la versión 3 la turbina usa el modelo DMST con los factores
de validación f_cp_t, f_ct_t y f_cs_t en lugar de Cp y Cp/Ct supuestos; 1000 muestras, semilla fija.
Para que los resultados que no dependen de la turbina coincidan con los de la versión 2, los 13 parámetros de entonces
se muestrean con la misma semilla y en el mismo orden (los dos de la turbina supuesta se descartan) y los tres factores
de la turbina, con una segunda semilla (SEED + 1).
Uso: python3 ruta_montecarlo.py [N] [v1|v2]   ("v1" y "v2" repiten los cálculos del 18-10-2026 y del 19-10-2026)
Salida: montecarlo_ruta.json (o montecarlo_ruta_v1.json, _v2.json): percentiles 5/50/95, probabilidades de orden
y correlaciones de Spearman.
"""
import json
import sys
import numpy as np
from multiprocessing import Pool
import coeficientes as cf
import ruta_modelo as rm

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
VER = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] in ("v1", "v2") else "v3"
V1 = VER == "v1"
SEED = 20261018
RANGES = dict({"v1": cf.RANGES_V1, "v2": cf.RANGES_V2}.get(VER, cf.RANGES))
KEYS = list(RANGES)

_D = None


def _init():
    global _D
    _D = rm.load_wind()


def run(sample):
    d, groups, nu, ng = _D
    p = rm.P(v1=V1, turb="supuesta" if VER != "v3" else "dmst", **sample)
    r = rm.annual(d, groups, nu, ng, p)
    tot = {row[0]: row[3] for row in rm.table(r)}
    out = {"R1": r["sea"]["R1"], "fuel_total": r["fuel_total"], "pen R1": r["sea"]["_pen R1"]}
    for key in ("T1", "T2", "T3", "T3x", "T4", "T5"):
        out[key] = tot[key]
    for key in ("T5 literal", "T2 literal", "T3 literal", "T5 popa", "T5 inventor"):
        out[key] = r["sea"][key]
    return out


def spearman(x, y):
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def pct(v):
    return {"p5": float(np.percentile(v, 5)), "p50": float(np.percentile(v, 50)), "p95": float(np.percentile(v, 95)),
            "media": float(v.mean())}


def main():
    rng = np.random.default_rng(SEED)
    if VER == "v3":
        rng_t = np.random.default_rng(SEED + 1)
        samples = []
        for _ in range(N):
            smp = {k: float(rng.uniform(*cf.RANGES_V2[k])) for k in cf.RANGES_V2}
            del smp["cp_t"], smp["cpct_t"]
            smp.update({k: float(rng_t.uniform(*RANGES[k])) for k in ("f_cp_t", "f_ct_t", "f_cs_t")})
            samples.append({k: smp[k] for k in KEYS})
    else:
        samples = [{k: float(rng.uniform(*RANGES[k])) for k in KEYS} for _ in range(N)]
    with Pool(2, initializer=_init) as pool:
        res = pool.map(run, samples, chunksize=10)
    arr = {k: np.array([r[k] for r in res]) for k in res[0]}
    summary = {k: pct(v) for k, v in arr.items()}
    prob = {
        "T3x > T4": float((arr["T3x"] > arr["T4"]).mean()),
        "T3 >= T3x": float((arr["T3"] >= arr["T3x"] - 1e-9).mean()),
        "T5 > T4": float((arr["T5"] > arr["T4"]).mean()),
        "T1 < 0": float((arr["T1"] < 0).mean()),
        "T2 < 0": float((arr["T2"] < 0).mean()),
        "T5 literal < 0 (mar)": float((arr["T5 literal"] < 0).mean()),
        "T5 inventor = 0 (mar)": float((np.abs(arr["T5 inventor"]) < 0.05).mean()),
        "orden T3 >= T3x > T5 > T4 > T2 > T1": float(((arr["T3"] >= arr["T3x"] - 1e-9) & (arr["T3x"] > arr["T5"])
                                                      & (arr["T5"] > arr["T4"]) & (arr["T4"] > arr["T2"])
                                                      & (arr["T2"] > arr["T1"])).mean()),
    }
    summary["T3x/R1"] = pct(arr["T3x"] / arr["R1"])
    summary["T3x/T4"] = pct(arr["T3x"] / arr["T4"])
    summary["R1 %"] = pct(100 * arr["R1"] / arr["fuel_total"])
    summary["pen R1 / (R1 + pen R1)"] = pct(arr["pen R1"] / (arr["R1"] + arr["pen R1"]))
    corr = {out: {k: spearman(np.array([s[k] for s in samples]), arr[out]) for k in KEYS}
            for out in ("R1", "T3x", "T5 literal", "T4", "T5", "pen R1")}
    pairs = [("T3", "T3x"), ("T3x", "T5"), ("T5", "T4"), ("T4", "T2"), ("T2", "T1"), ("T3x", "T4"), ("T5", "T2")]
    prob["pares"] = {f"{a} > {b}": float((arr[a] > arr[b] + (0 if a != "T3" else -1e-9)).mean()) for a, b in pairs}
    sfx = "" if VER == "v3" else "_" + VER
    np.savez_compressed(f"montecarlo_muestras{sfx}.npz",
                        **{f"in_{k}": np.array([s[k] for s in samples]) for k in KEYS},
                        **{f"out_{k}": v for k, v in arr.items()})
    name = f"montecarlo_ruta{sfx}.json"
    json.dump({"N": N, "semilla": SEED, "version": int(VER[1]), "rangos": RANGES, "resumen": summary,
               "probabilidades": prob, "spearman": corr}, open(name, "w"), indent=1, ensure_ascii=False)
    print(json.dumps({"resumen": {k: {kk: round(vv, 2) for kk, vv in v.items()} for k, v in summary.items()},
                      "prob": prob}, indent=1, ensure_ascii=False))
    for out, c in corr.items():
        print(out, {k: round(v, 2) for k, v in sorted(c.items(), key=lambda kv: -abs(kv[1]))})


if __name__ == "__main__":
    main()
