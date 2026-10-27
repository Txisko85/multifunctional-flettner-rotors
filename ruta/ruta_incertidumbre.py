#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Incertidumbre de la versión 4 (22-10-2026; revisión simulada R2-1 y R3-4).

1. Incertidumbre epistémica de los coeficientes del rotor: escenarios (coeficientes.ESCENARIOS_COEF), no una
   variable muestreada. El escenario de referencia es A (validacion_viking_grace.py).
2. Incertidumbre paramétrica dentro de cada escenario: muestreo uniforme e independiente de los 15 rangos de
   coeficientes.RANGES_V4_MC. C_t,S se muestrea entre su cota cinemática, KAPPA_S C_p,S, y 1,0 mediante la variable
   auxiliar u_ct (independiente), para que el diseño de Sobol tenga entradas independientes.
3. Índices de Sobol de primer orden (estimador de Saltelli et al. 2010) y totales (estimador de Jansen), con el
   diseño de matrices A, B y A_B^(i) (N (d + 2) evaluaciones) e intervalos de confianza del 95 % por bootstrap
   (1000 remuestreos). Las salidas se centran (se resta su media) antes de aplicar los estimadores: el de primer
   orden de Saltelli es muy ruidoso con salidas de media grande frente a su dispersión (26-10-2026; la versión del
   22-10-2026 no centraba y daba índices de primer orden negativos o mayores que 1). Semilla fija. Las mismas 2N
   evaluaciones de A y B dan los percentiles y las probabilidades de orden, también con intervalos bootstrap.
Uso: python3 ruta_incertidumbre.py [N] [N_escenarios]   (por defecto N = 256 para A y 256 muestras para cada uno de
     los escenarios A+, AB, B y B')
     python3 ruta_incertidumbre.py --muestras   (repite solo el análisis con incertidumbre_muestras.npz)
Salida: incertidumbre_ruta.json y incertidumbre_muestras.npz.
"""
import json
import sys
import numpy as np
from multiprocessing import Pool
import coeficientes as cf
import ruta_modelo as rm

_ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
N = int(_ARGS[0]) if len(_ARGS) > 0 else 256
N_ESC = int(_ARGS[1]) if len(_ARGS) > 1 else 256
SEED = 20261022
N_BOOT = 1000
ESC_OTROS = ("A+", "AB", "B", "B'")
RANGES = {k: v for k, v in cf.RANGES_V4_MC.items() if k != "ct_s"}
RANGES["u_ct"] = (0.0, 1.0)
KEYS = list(RANGES)
OUTS = ("R1", "T3x", "T3", "T5", "T4", "T2", "T1", "R1t", "R1-30", "fuel_total", "pen R1")
_D = None


def _init():
    global _D
    _D = rm.load_wind()


def to_params(x):
    s = dict(zip(KEYS, x))
    u = s.pop("u_ct")
    lo = cf.KAPPA_S * s["cp_s"]
    s["ct_s"] = lo + u * (cf.CT_S_MAX - lo)
    return s


def run(args):
    x, esc = args
    d, groups, nu, ng = _D
    p = rm.P(**cf.ESCENARIOS_COEF[esc], **to_params(x))
    r = rm.annual(d, groups, nu, ng, p)
    tot = {row[0]: row[3] for row in rm.table(r)}
    out = {"R1": r["sea"]["R1"], "fuel_total": r["fuel_total"], "pen R1": r["sea"]["_pen R1"],
           "R1t": r["sea"]["R1t"], "R1-30": r["sea"]["R1-30"]}
    for key in ("T1", "T2", "T3", "T3x", "T4", "T5"):
        out[key] = tot[key]
    return [out[k] for k in OUTS]


def scale(u):
    lo = np.array([RANGES[k][0] for k in KEYS])
    hi = np.array([RANGES[k][1] for k in KEYS])
    return lo + u * (hi - lo)


def sobol(fa, fb, fab):
    """Índices de primer orden (Saltelli et al. 2010) y totales (Jansen) para una salida, con las salidas
    centradas en la media de A y B."""
    y = np.concatenate([fa, fb])
    var = np.var(y, ddof=1)
    if var <= 0:
        return np.zeros(fab.shape[1]), np.zeros(fab.shape[1])
    m = y.mean()
    fa, fb, fab = fa - m, fb - m, fab - m
    s1 = np.mean(fb[:, None] * (fab - fa[:, None]), axis=0) / var
    st = 0.5 * np.mean((fa[:, None] - fab) ** 2, axis=0) / var
    return s1, st


ORDER = [("T3", "T3x", True), ("T3x", "T5", False), ("T5", "T4", False), ("T4", "T2", False), ("T2", "T1", False),
         ("T3x", "T4", False), ("T3x", "R1t", False), ("T3", "R1t", False), ("T3", "R1-30", False),
         ("T3x", "R1-30", False), ("R1-30", "R1t", False)]


def probabilities(arr):
    pr = {}
    for a, b, ge in ORDER:
        pr[f"{a} {'>=' if ge else '>'} {b}"] = (arr[a] >= arr[b] - 1e-9) if ge else (arr[a] > arr[b])
    pr["T1 < 0"] = arr["T1"] < 0
    pr["T2 < 0"] = arr["T2"] < 0
    pr["orden T3 >= T3x > T5 > T4 > T2 > T1"] = ((arr["T3"] >= arr["T3x"] - 1e-9) & (arr["T3x"] > arr["T5"])
                                                 & (arr["T5"] > arr["T4"]) & (arr["T4"] > arr["T2"])
                                                 & (arr["T2"] > arr["T1"]))
    return pr


def boot_ci(values, stat, rng, n=N_BOOT):
    m = len(values)
    reps = [stat(values[rng.integers(0, m, m)]) for _ in range(n)]
    return float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5))


def summarize(arr, rng):
    out = {}
    for k, v in arr.items():
        out[k] = {"p5": float(np.percentile(v, 5)), "p50": float(np.percentile(v, 50)),
                  "p95": float(np.percentile(v, 95)), "media": float(v.mean()),
                  "ic95_p50": boot_ci(v, np.median, rng)}
    out["R1 %"] = {q: float(np.percentile(100 * arr["R1"] / arr["fuel_total"], int(q[1:]))) for q in ("p5", "p50", "p95")}
    probs = {}
    for k, mask in probabilities(arr).items():
        mask = mask.astype(float)
        probs[k] = {"p": float(mask.mean()), "ic95": boot_ci(mask, np.mean, rng)}
    return out, probs


def main():
    d = len(KEYS)
    if "--muestras" in sys.argv:
        z = np.load("incertidumbre_muestras.npz")
        assert list(z["claves"]) == KEYS and list(z["salidas"]) == list(OUTS)
        xs, ya, y_esc = z["x"], z["y"], z["y_esc"]
        esc_x = dict(zip(ESC_OTROS, z["x_esc"]))
        n_esc = y_esc.shape[0] // len(ESC_OTROS)
        n = xs.shape[0] // (d + 2)
        res = np.concatenate([ya, y_esc])
    else:
        n, n_esc = N, N_ESC
        rng = np.random.default_rng(SEED)
        A = scale(rng.random((n, d)))
        B = scale(rng.random((n, d)))
        ABs = []
        for i in range(d):
            ab = A.copy()
            ab[:, i] = B[:, i]
            ABs.append(ab)
        design = [A, B] + ABs
        xs = np.concatenate(design)
        jobs = [(x, "A") for x in xs]
        esc_x = {e: scale(rng.random((n_esc, d))) for e in ESC_OTROS}
        for e, xe in esc_x.items():
            jobs += [(x, e) for x in xe]
        with Pool(2, initializer=_init) as pool:
            res = np.array(pool.map(run, jobs, chunksize=8))
        ya = res[:len(xs)]
        np.savez_compressed("incertidumbre_muestras.npz", x=xs, y=ya, x_esc=np.stack([esc_x[e] for e in ESC_OTROS]),
                            y_esc=res[len(xs):], escenarios=np.array(ESC_OTROS), claves=np.array(KEYS),
                            salidas=np.array(OUTS))
    analyse(xs, ya, res, n, n_esc, d)


def analyse(xs, ya, res, N, N_ESC, d):
    n_a = len(xs)
    fa, fb = ya[:N], ya[N:2 * N]
    fab = np.stack([ya[(2 + i) * N:(3 + i) * N] for i in range(d)], axis=1)      # N x d x salidas
    out = {"N": N, "N_escenarios": N_ESC, "semilla": SEED, "rangos": RANGES, "escenario_sobol": "A",
           "sobol": {}, "escenarios": {}}
    rng_b = np.random.default_rng(SEED + 1)
    for j, name in enumerate(OUTS):
        if name in ("fuel_total",):
            continue
        s1, st = sobol(fa[:, j], fb[:, j], fab[:, :, j])
        ci1, cit = [], []
        for i in range(d):
            idx = [rng_b.integers(0, N, N) for _ in range(N_BOOT)]
            b1 = [sobol(fa[ix, j], fb[ix, j], fab[ix][:, i:i + 1, j])[0][0] for ix in idx]
            bt = [sobol(fa[ix, j], fb[ix, j], fab[ix][:, i:i + 1, j])[1][0] for ix in idx]
            ci1.append((float(np.percentile(b1, 2.5)), float(np.percentile(b1, 97.5))))
            cit.append((float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))))
        out["sobol"][name] = {k: {"S1": float(s1[i]), "S1_ic95": ci1[i], "ST": float(st[i]), "ST_ic95": cit[i]}
                              for i, k in enumerate(KEYS)}
    arr_a = {k: np.concatenate([fa[:, j], fb[:, j]]) for j, k in enumerate(OUTS)}
    s, p = summarize(arr_a, rng_b)
    out["escenarios"]["A"] = {"resumen": s, "probabilidades": p, "n": 2 * N}
    pos = n_a
    for e in ESC_OTROS:
        ye = res[pos:pos + N_ESC]
        pos += N_ESC
        arr = {k: ye[:, j] for j, k in enumerate(OUTS)}
        s, p = summarize(arr, rng_b)
        out["escenarios"][e] = {"resumen": s, "probabilidades": p, "n": N_ESC}
    json.dump(out, open("incertidumbre_ruta.json", "w"), indent=1, ensure_ascii=False)
    for e, v in out["escenarios"].items():
        r = v["resumen"]
        print(f"Escenario {e} (n = {v['n']}): " + "; ".join(
            f"{k} {r[k]['p50']:.1f} ({r[k]['p5']:.1f}-{r[k]['p95']:.1f})" for k in ("R1", "T3x", "T3", "T5", "R1t", "R1-30")))
        print("   " + "; ".join(f"P({k}) = {pp['p']:.2f}" for k, pp in v["probabilidades"].items()))
    for name, sv in out["sobol"].items():
        top = sorted(sv.items(), key=lambda kv: -kv[1]["ST"])[:5]
        print(f"Sobol {name}: " + ", ".join(f"{k} S1 {v['S1']:.2f} ST {v['ST']:.2f}" for k, v in top))


if __name__ == "__main__":
    main()
