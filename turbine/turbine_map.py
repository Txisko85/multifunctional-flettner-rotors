#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Maps of the case-study turbine for the route model (A1, 21-10-2026).

Three-bladed NACA 0021 H-Darrieus turbine on top of the rotor: R 1.7 m and H 6 m (the dimensions of the VAWT of the
route model), chord 0.283 m for the same solidity as the rotor validated in the wind tunnel (N c / D = 0.25); two
struts per blade (chord 0.10 m, CD 0.02) and a central shaft of 0.30 m (CD 0.6), assumed.
Central DMST model of dmst.py (Gormont-Berg with K1 = 0 when |alpha| falls), without correction factors: the route
model applies the factors f_cp_t, f_ct_t and f_cs_t that result from the validation (dmst_validation.py).
Output: turbine_map.json with Cp, Ct and |Cs| on a grid of lambda (0.5-5.0) and wind U (3-20 m/s), and Ct of the
parked turbine.
"""
import json
import os
import time
import numpy as np
import dmst

HERE = os.path.dirname(os.path.abspath(__file__))
LAMS = [round(x, 2) for x in np.arange(0.5, 5.001, 0.25)]
US = [3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]


def main():
    tur = dmst.Turbine()
    t0 = time.time()
    cp = np.zeros((len(US), len(LAMS)))
    ct = np.zeros_like(cp)
    cs = np.zeros_like(cp)
    park = []
    for i, u in enumerate(US):
        for j, lam in enumerate(LAMS):
            r = dmst.solve(tur, u, lam)
            cp[i, j], ct[i, j], cs[i, j] = r["cp"], r["ct"], abs(r["cs"])
        park.append(dmst.parked(tur, u)["ct"])
        jm = int(np.argmax(cp[i]))
        print(f"U {u:4.1f} m/s: max Cp {cp[i, jm]:.3f} at λ {LAMS[jm]:.2f} (Ct {ct[i, jm]:.3f}, |Cs| {cs[i, jm]:.3f}); "
              f"parked Ct {park[-1]:.3f}  [{time.time() - t0:.0f} s]", flush=True)
    out = dict(
        description="H-Darrieus turbine of the case study (A1): Cp, Ct and |Cs| of the central DMST model (dmst.py), "
                    "without correction factors; reference area 2 R H. Cs has the direction of the Magnus lift of a "
                    "cylinder turning in the same sense as the turbine.",
        geometry=dict(n=tur.n, R=tur.R, H=tur.H, c=tur.c, tc=tur.tc, solidity=tur.solidity, n_struts=tur.n_struts,
                      c_strut=tur.c_strut, cd_strut=tur.cd_strut, d_shaft=tur.d_shaft, cd_shaft=tur.cd_shaft,
                      area=tur.area),
        model="DMST, Gormont-Berg (A_M 6, K1 1/0), Buhl, finite aspect ratio, NACA 0021 polars of Sheldahl and Klimas",
        lam=LAMS, U=US, cp=cp.round(5).tolist(), ct=ct.round(5).tolist(), cs=cs.round(5).tolist(),
        ct_parked=[round(x, 5) for x in park])
    json.dump(out, open(os.path.join(HERE, "turbine_map.json"), "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
