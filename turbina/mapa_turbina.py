#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapas de la turbina del caso de estudio para el modelo de ruta (A1, 21-10-2026).

Turbina H-Darrieus de tres palas NACA 0021 sobre el rotor: R 1,7 m y H 6 m (las dimensiones de la VAWT del modelo
de ruta), cuerda 0,283 m para la misma solidez que el rotor validado en túnel (N c / D = 0,25); dos brazos por pala
(cuerda 0,10 m, CD 0,02) y eje central de 0,30 m (CD 0,6), supuestos.
Modelo DMST central de dmst.py (Gormont-Berg con K1 = 0 al decrecer |alfa|), sin factores de corrección: el
modelo de ruta aplica los factores f_cp_t, f_ct_t y f_cs_t que resultan de la validación (validacion_dmst.py).
Salida: mapa_turbina.json con Cp, Ct y |Cs| en una malla de lambda (0,5-5,0) y viento U (3-20 m/s), y Ct de la
turbina parada.
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
        print(f"U {u:4.1f} m/s: Cp máx {cp[i, jm]:.3f} a λ {LAMS[jm]:.2f} (Ct {ct[i, jm]:.3f}, |Cs| {cs[i, jm]:.3f}); "
              f"parada Ct {park[-1]:.3f}  [{time.time() - t0:.0f} s]", flush=True)
    out = dict(
        descripcion="Turbina H-Darrieus del caso (A1): Cp, Ct y |Cs| del modelo DMST central (dmst.py), sin factores "
                    "de corrección; área de referencia 2 R H. Cs tiene el sentido de la sustentación Magnus de un "
                    "cilindro que gira en el mismo sentido que la turbina.",
        geometria=dict(n=tur.n, R=tur.R, H=tur.H, c=tur.c, tc=tur.tc, solidez=tur.solidity, n_brazos=tur.n_struts,
                       c_brazo=tur.c_strut, cd_brazo=tur.cd_strut, d_eje=tur.d_shaft, cd_eje=tur.cd_shaft,
                       area=tur.area),
        modelo="DMST, Gormont-Berg (A_M 6, K1 1/0), Buhl, alargamiento finito, polares NACA 0021 de Sheldahl y Klimas",
        lam=LAMS, U=US, cp=cp.round(5).tolist(), ct=ct.round(5).tolist(), cs=cs.round(5).tolist(),
        ct_parada=[round(x, 5) for x in park])
    json.dump(out, open(os.path.join(HERE, "mapa_turbina.json"), "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
