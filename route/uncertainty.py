#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Uncertainty of version 5 (27-10-2026; second simulated review: R2-1, R3-m1). Version 4: 22-10-2026.

Changes of version 5: N = 1024 (previously 256); C_t,S is sampled from the strict bound C_p,S/lambda; the induced
resistance comes from the sway-yaw balance, so e_H is no longer sampled and c_H, C_pa, C_wa and the centre of the side
force of the wind on the ship x_s are sampled instead (18 inputs); scenario A-VG (fitted to the Viking Grace) added to
the others; comparator R1-D5.

1. Epistemic uncertainty of the rotor coefficients: scenarios (coefficients.COEF_SCENARIOS), not a sampled variable.
   The reference scenario is A (viking_grace_check.py).
2. Parametric uncertainty within each scenario: uniform and independent sampling of the ranges of
   coefficients.RANGES_V5_MC. C_t,S is sampled between its kinematic bound, KAPPA_S C_p,S, and 1.0 through the
   auxiliary variable u_ct (independent), so that the Sobol design has independent inputs.
3. First-order Sobol indices (estimator of Saltelli et al. 2010) and total indices (Jansen estimator), with the design
   of matrices A, B and A_B^(i) (N (d + 2) evaluations) and 95 % bootstrap confidence intervals (1000 resamples). The
   outputs are centred (their mean is subtracted) before applying the estimators: the first-order estimator of
   Saltelli is very noisy with outputs whose mean is large compared with their spread (26-10-2026; the version of
   22-10-2026 did not centre and gave first-order indices that were negative or larger than 1). Fixed seed. The same
   2N evaluations of A and B give the percentiles and the ranking probabilities, also with bootstrap intervals.
Usage: python3 uncertainty.py [N] [N_scenarios]   (by default N = 1024 for A and 256 samples for each of the
       scenarios A+, AB, B, B' and A-VG)
       python3 uncertainty.py --samples   (repeats only the analysis with uncertainty_samples.npz)
Output: uncertainty.json and uncertainty_samples.npz.
"""
import json
import sys
import numpy as np
from multiprocessing import Pool
import coefficients as cf
import route_model as rm

_ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
N = int(_ARGS[0]) if len(_ARGS) > 0 else 1024
N_SCEN = int(_ARGS[1]) if len(_ARGS) > 1 else 256
SEED = 20261022
N_BOOT = 1000
OTHER_SCEN = ("A+", "AB", "B", "B'", "A-VG")
RANGES = {k: v for k, v in cf.RANGES_V5_MC.items() if k != "ct_s"}
RANGES["u_ct"] = (0.0, 1.0)
KEYS = list(RANGES)
OUTS = ("R1", "T3x", "T3", "T5", "T4", "T2", "T1", "R1t", "R1-30", "R1-D5", "fuel_total", "pen R1")
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
    p = rm.P(**cf.COEF_SCENARIOS[esc], **to_params(x))
    r = rm.annual(d, groups, nu, ng, p)
    tot = {row[0]: row[3] for row in rm.table(r)}
    out = {"R1": r["sea"]["R1"], "fuel_total": r["fuel_total"], "pen R1": r["sea"]["_pen R1"],
           "R1t": r["sea"]["R1t"], "R1-30": r["sea"]["R1-30"], "R1-D5": r["sea"]["R1-D5"]}
    for key in ("T1", "T2", "T3", "T3x", "T4", "T5"):
        out[key] = tot[key]
    return [out[k] for k in OUTS]


def scale(u):
    lo = np.array([RANGES[k][0] for k in KEYS])
    hi = np.array([RANGES[k][1] for k in KEYS])
    return lo + u * (hi - lo)


def sobol(fa, fb, fab):
    """First-order (Saltelli et al. 2010) and total (Jansen) indices for one output, with the outputs centred on the
    mean of A and B."""
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
         ("T3x", "R1-30", False), ("R1-30", "R1t", False), ("T3x", "R1-D5", False)]


def probabilities(arr):
    pr = {}
    for a, b, ge in ORDER:
        pr[f"{a} {'>=' if ge else '>'} {b}"] = (arr[a] >= arr[b] - 1e-9) if ge else (arr[a] > arr[b])
    pr["T1 < 0"] = arr["T1"] < 0
    pr["T2 < 0"] = arr["T2"] < 0
    pr["ranking T3 >= T3x > T5 > T4 > T2 > T1"] = ((arr["T3"] >= arr["T3x"] - 1e-9) & (arr["T3x"] > arr["T5"])
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
                  "p95": float(np.percentile(v, 95)), "mean": float(v.mean()),
                  "ci95_p50": boot_ci(v, np.median, rng)}
    out["R1 %"] = {q: float(np.percentile(100 * arr["R1"] / arr["fuel_total"], int(q[1:]))) for q in ("p5", "p50", "p95")}
    probs = {}
    for k, mask in probabilities(arr).items():
        mask = mask.astype(float)
        probs[k] = {"p": float(mask.mean()), "ci95": boot_ci(mask, np.mean, rng)}
    return out, probs


CKPT = "uncertainty_partial.npy"
BLOCK = 256


def _run_checkpointed(jobs):
    """Evaluates the jobs in blocks and saves the partial results in CKPT, so that the computation can be resumed if it
    is interrupted (27-10-2026). The order of the jobs is fixed by the seed, so resuming gives the same result."""
    import os
    done = np.load(CKPT) if os.path.exists(CKPT) else np.empty((0, len(OUTS)))
    with Pool(2, initializer=_init) as pool:
        while len(done) < len(jobs):
            blk = jobs[len(done):len(done) + BLOCK]
            done = np.concatenate([done, np.array(pool.map(run, blk, chunksize=8))])
            np.save(CKPT, done)
            print(f"{len(done)}/{len(jobs)}", flush=True)
    return done


def main():
    d = len(KEYS)
    if "--samples" in sys.argv:
        z = np.load("uncertainty_samples.npz")
        assert list(z["keys"]) == KEYS and list(z["outputs"]) == list(OUTS)
        xs, ya, y_esc = z["x"], z["y"], z["y_esc"]
        esc_x = dict(zip(OTHER_SCEN, z["x_esc"]))
        n_esc = y_esc.shape[0] // len(OTHER_SCEN)
        n = xs.shape[0] // (d + 2)
        res = np.concatenate([ya, y_esc])
    else:
        n, n_esc = N, N_SCEN
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
        esc_x = {e: scale(rng.random((n_esc, d))) for e in OTHER_SCEN}
        for e, xe in esc_x.items():
            jobs += [(x, e) for x in xe]
        res = _run_checkpointed(jobs)
        ya = res[:len(xs)]
        np.savez_compressed("uncertainty_samples.npz", x=xs, y=ya, x_esc=np.stack([esc_x[e] for e in OTHER_SCEN]),
                            y_esc=res[len(xs):], scenarios=np.array(OTHER_SCEN), keys=np.array(KEYS),
                            outputs=np.array(OUTS))
    analyse(xs, ya, res, n, n_esc, d)


def analyse(xs, ya, res, N, N_SCEN, d):
    n_a = len(xs)
    fa, fb = ya[:N], ya[N:2 * N]
    fab = np.stack([ya[(2 + i) * N:(3 + i) * N] for i in range(d)], axis=1)      # N x d x outputs
    out = {"N": N, "N_scenarios": N_SCEN, "seed": SEED, "ranges": RANGES, "sobol_scenario": "A",
           "sobol": {}, "scenarios": {}}
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
        out["sobol"][name] = {k: {"S1": float(s1[i]), "S1_ci95": ci1[i], "ST": float(st[i]), "ST_ci95": cit[i]}
                              for i, k in enumerate(KEYS)}
    arr_a = {k: np.concatenate([fa[:, j], fb[:, j]]) for j, k in enumerate(OUTS)}
    s, p = summarize(arr_a, rng_b)
    out["scenarios"]["A"] = {"summary": s, "probabilities": p, "n": 2 * N}
    pos = n_a
    for e in OTHER_SCEN:
        ye = res[pos:pos + N_SCEN]
        pos += N_SCEN
        arr = {k: ye[:, j] for j, k in enumerate(OUTS)}
        s, p = summarize(arr, rng_b)
        out["scenarios"][e] = {"summary": s, "probabilities": p, "n": N_SCEN}
    json.dump(out, open("uncertainty.json", "w"), indent=1, ensure_ascii=False)
    for e, v in out["scenarios"].items():
        r = v["summary"]
        print(f"Scenario {e} (n = {v['n']}): " + "; ".join(
            f"{k} {r[k]['p50']:.1f} ({r[k]['p5']:.1f}-{r[k]['p95']:.1f})" for k in ("R1", "T3x", "T3", "T5", "R1t", "R1-30",
                                                                                 "R1-D5")))
        print("   " + "; ".join(f"P({k}) = {pp['p']:.2f}" for k, pp in v["probabilities"].items()))
    for name, sv in out["sobol"].items():
        top = sorted(sv.items(), key=lambda kv: -kv[1]["ST"])[:5]
        print(f"Sobol {name}: " + ", ".join(f"{k} S1 {v['S1']:.2f} ST {v['ST']:.2f}" for k, v in top))


if __name__ == "__main__":
    main()
