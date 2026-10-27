#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Three-dimensional validation of the lifting-line model (version 5, 27-10-2026; simulated review v2, R2-3).

1. Dependence of the lift on the aspect ratio: De Marco et al. (2016, Fig. 10; URANS with the rotor on a symmetry
   plane, as in the tests of Badalamenti and Prince) give C_L for AR = 2, 4, 6 and 8 and endplates de/d = 1 (no
   plate), 2 and 3, at k = 1.5, 2 and 3 (coef_data/demarco2016_fig10.json). The model (deck as a perfect mirror,
   m_d = 1; plate with m_p = 0 without plate and m_p = 0.5 or 1 with plate) predicts C_L(AR)/C_L(AR 6) with the two
   section models. The ratio is compared, not the absolute value: the model gives relative factors (f_X, f30,
   eta_B), and the absolute level is set by the coefficient scenario.
2. Empirical factor of the 30 m rotor: with plate de/d = 2, f30 = [(AR'/6) r - 1]/((AR' - 6)/6), with r = C_L(AR')/
   C_L(AR 6) interpolated at AR' = 7.5 (30 m / 4 m).
3. Insensitivity of the factors to the level of C_l,2D: f_X and f30 with the C_l,2D of set B multiplied by 1.7 (the
   ratio A/B at k = 3) and with the shape of set A. The wing analogy is linear: its factors depend neither on the
   level nor on the shape; the spin model does not depend on the level, but does on the shape at k >= 2.5.
4. Estimate of f_X: empirical f30 times the ratio f_X/f30 of the model (nearly constant in all variants).
Output: lifting_line_validation.json.
"""
import json
import os
import numpy as np
import coefficients as cf
import lifting_line as ls

HERE = os.path.dirname(os.path.abspath(__file__))
DM = json.load(open(os.path.join(HERE, "coef_data", "demarco2016_fig10.json"), encoding="utf-8"))["data"]
D = 4.0
KS = (1.5, 2.0, 3.0)
ARS = (2, 4, 6, 8)


def cl3d(ar, k, model, m_p):
    h = ar * D
    r = ls.Rotor([(0.0, h, D)], z_plate=h if m_p > 0 else None, m_d=1.0, m_p=m_p)
    s = r.solve(k, model)
    return r.segment(s, 0, h)[0] / (0.5 * D * h)


def ratios_obs(de, k):
    return {ar: DM[de][f"AR{ar}"]["CL"][str(k)] / DM[de]["AR6"]["CL"][str(k)] for ar in ARS}


def validation():
    out = {}
    for model in ("wing", "spin"):
        for de, m_ps in (("de1", (0.0,)), ("de2", (0.5, 1.0)), ("de3", (0.5, 1.0))):
            for m_p in m_ps:
                rows, err = {}, []
                for k in KS:
                    cl = {ar: cl3d(ar, k, model, m_p) for ar in ARS}
                    mod = {ar: cl[ar] / cl[6] for ar in ARS}
                    obs = ratios_obs(de, k)
                    rows[str(k)] = dict(model={str(a): mod[a] for a in ARS}, URANS={str(a): obs[a] for a in ARS},
                                        eta_AR6=cl[6] / float(ls.cl2d(k)))
                    err += [mod[a] - obs[a] for a in ARS if a != 6]
                out[f"{model}, {de}, m_p {m_p:g}"] = dict(k=rows, rms=float(np.sqrt(np.mean(np.square(err)))),
                                                        bias=float(np.mean(err)))
    return out


def f30_empirical():
    out = {}
    for de in ("de2", "de3"):
        for k in KS:
            r6, r8 = 1.0, ratios_obs(de, k)[8]
            r75 = r6 + (r8 - r6) * 0.75
            out[f"{de}, k {k:g}"] = ((7.5 / 6.0) * r75 - 1.0) / (1.5 / 6.0)
    return out


def insensitivity():
    base_cl = ls.CL_2D.copy()
    k_a = np.array(cf.K_A_EXT)
    cl_a = np.array(cf.CL_A_EXT)
    variants = {"B": base_cl, "B x 1.7": 1.7 * base_cl, "shape of A": np.interp(ls.K_2D, k_a, cl_a)}
    out = {}
    for lab, arr in variants.items():
        ls.CL_2D = arr
        for model in ("wing", "spin"):
            for m_p in (0.5, 1.0):
                f = ls.factors(1.0, m_p, model)
                out[f"{lab}; {model}; m_d 1, m_p {m_p:g}"] = dict(f_X=f["f_X"], f30=f["f30"], eta0=f["eta0"], k=f["k"],
                                                                 f_X_f30=[a / b for a, b in zip(f["f_X"], f["f30"])])
    ls.CL_2D = base_cl
    return out


def main():
    res = dict(validation=validation(), f30_empirical=f30_empirical(), insensitivity=insensitivity())
    # estimated f_X: empirical f30 (without the outlier de/d = 2, k = 2, AR 8) times the ratio f_X/f30 of the model,
    # which hardly depends on the section model, the level or the shape of C_l,2D
    f30e = [v for k, v in res["f30_empirical"].items() if k != "de2, k 2"]
    ratio = [x for v in res["insensitivity"].values() for x in v["f_X_f30"]]
    res["f_X_estimate"] = dict(f30_empirical=(min(f30e), max(f30e)), f_X_f30=(min(ratio), max(ratio)),
                               f_X=(min(f30e) * min(ratio), max(f30e) * max(ratio)))
    json.dump(res, open(os.path.join(HERE, "lifting_line_validation.json"), "w"), indent=1, ensure_ascii=False)
    print("C_L(AR)/C_L(AR 6): model against URANS (De Marco et al. 2016)")
    for key, v in res["validation"].items():
        print(f"  {key:22s} RMS {v['rms']:.3f}  bias {v['bias']:+.3f}  " + "; ".join(
            f"k {k}: " + "/".join(f"{v['k'][k]['model'][a]:.2f}~{v['k'][k]['URANS'][a]:.2f}" for a in ("2", "4", "8"))
            for k in v["k"]))
    print("empirical f30 (De Marco, AR 7.5):", {k: round(x, 2) for k, x in res["f30_empirical"].items()})
    print("estimated f_X:", {k: tuple(round(x, 2) for x in v) for k, v in res["f_X_estimate"].items()})
    print("Factors with different C_l,2D:")
    for k, v in res["insensitivity"].items():
        print(f"  {k:36s} f_X {np.round(v['f_X'], 3)}  f30 {np.round(v['f30'], 3)}  eta0 {np.round(v['eta0'], 3)}  "
              f"f_X/f30 {np.round(v['f_X_f30'], 3)}")
    return res


if __name__ == "__main__":
    main()
