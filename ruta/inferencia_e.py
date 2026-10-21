#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inferencia de la eficacia e del modelo de resistencia inducida, dR = Fy^2 / (pi e q_w (2T)^2), a partir de los
resultados de Bickert, Gentaz, Delvoye y Babarit (MARINE 2025): dR/R0 a escala real con Fy = 2 R0.
El artículo no da R0; se estima con ITTC-78 sin resistencia del aire: C_T = (1+k) C_F + C_W + dC_F,
con dimensiones de referencia del KVLCC2 (SIMMAN) y del KCS (Tokio 2015). Incertidumbre de C_T: +-10 %.
"""
import math

NU, RHO_W = 1.19e-6, 1025.0
CASES = {  # nombre: (Lpp m, T m, S m2, V kn, 1+k, C_W, dC_F, dR/R0 con Fy = 2 R0)
    "KVLCC2, 15.5 kn, 3 DOF": (320.0, 20.8, 27194.0, 15.5, 1.20, 0.0001, 0.000125, 0.067),
    "KCS, 20 kn, 3+1 DOF": (230.0, 10.8, 9424.0, 20.0, 1.10, 0.0002, 0.000100, 0.100),
    "KCS, 20 kn, 4 DOF": (230.0, 10.8, 9424.0, 20.0, 1.10, 0.0002, 0.000100, 0.211),
}


def cf_ittc57(re):
    return 0.075 / (math.log10(re) - 2.0) ** 2


def infer(L, T, S, vkn, ff, cw, dcf, rel, fac=1.0):
    v = vkn * 0.514444
    q = 0.5 * RHO_W * v * v
    ct = fac * (ff * cf_ittc57(v * L / NU) + cw + dcf)
    r0 = ct * q * S
    return ct, r0, 4.0 * r0 / (rel * math.pi * q * (2.0 * T) ** 2)


if __name__ == "__main__":
    for name, args in CASES.items():
        ct, r0, e = infer(*args)
        lo, hi = infer(*args, fac=0.9)[2], infer(*args, fac=1.1)[2]
        print(f"{name:24s} C_T {ct:.5f}  R0 {r0 / 1e3:5.0f} kN  e = {e:.2f} ({lo:.2f}-{hi:.2f})")
