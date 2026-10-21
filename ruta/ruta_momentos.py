#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Indicadores estructurales de A1 (apartado 4.6), con la base de coeficientes (19-10-2026; antes, nivel 0).

Momento flector en la base de una unidad, en operación, con el mismo viento aparente en toda la altura:
- rotor girando a k = min(3, U_s,max/U_a): resultante q A_F (CL^2 + CD^2)^0,5 aplicada a H_F/2;
- prolongación (T3, T3x) a k_X = k D_X/D_F, con el factor f_X en la sustentación, aplicada a H_F + H_X/2;
- VAWT (T1, T2, T4), desde el 21-10-2026: resultante de su empuje y su fuerza lateral (modelo DMST de ../turbina)
  en el estado de mayor fuerza permitido a ese viento, aplicada en H_F + H_X/2 (hasta el 19-10-2026, solo el
  empuje T_t = q A_X Cp/(Cp/Ct) con coeficientes supuestos);
- rotor abierto en Savonius (T2, T3, T5): T_S = C_t,S q D_S H_F en H_F/2; su sustentación no se conoce.
Las fuerzas de rotor y dispositivo se suman en módulo (cota superior). Referencia: momento del rotor solo en el
límite de operación supuesto, U_a = 20 m/s (k = 1,8). Se calcula también el viento aparente máximo con el que la
prolongación cabe en esa envolvente. La supervivencia con el rotor parado en temporal no se evalúa.
"""
import json
import math
import numpy as np
import ruta_modelo as rm

UA_DESIGN = 20.0


def coeffs(k, p):
    cl, cd = rm.rotor_coeffs(np.array([k]), p)
    return float(cl[0]), float(cd[0])


def moments(u_a, p):
    """Momentos en la base (N m) por unidad: rotor solo, prolongación, VAWT y Savonius abierto."""
    q = 0.5 * rm.RHO * u_a ** 2
    k = min(rm.SR_NOM, rm.US_MAX / u_a)
    cl, cd = coeffs(k, p)
    m_rot = q * rm.A_F * math.hypot(cl, cd) * rm.H_F / 2
    clx, cdx = coeffs(k * rm.D_X / rm.D_F, p)
    m_ext = q * rm.A_X * math.hypot(p.cl_ext * clx, cdx) * (rm.H_F + rm.H_X / 2)
    if p.turb == "dmst":
        _, th, sf = rm.vawt_states(u_a, p)
        f_vawt = float(np.max(np.hypot(th, sf)))
    else:
        f_vawt = q * rm.A_X * p.cp_t / p.cpct_t
    m_vawt = f_vawt * (rm.H_F + rm.H_X / 2)
    m_sav = p.ct_s * q * p.ds * rm.A_F * rm.H_F / 2
    return dict(k=k, rot=m_rot, ext=m_ext, vawt=m_vawt, sav=m_sav)


def ext_limit(p):
    m_design = moments(UA_DESIGN, p)["rot"]
    for u_a in np.arange(1.0, UA_DESIGN, 0.01):
        m = moments(u_a, p)
        if m["rot"] + m["ext"] > m_design:
            return round(float(u_a) - 0.01, 2), m_design
    return UA_DESIGN, m_design


def report(p, label):
    m12 = moments(12.0, p)                 # k = 3, límite de la relación de giro
    m20 = moments(UA_DESIGN, p)
    lim, m_design = ext_limit(p)
    m_sav_cut = moments(rm.UA_CUT_S, p)["sav"]
    out = dict(caso=label,
               M_rotor_kNm_12=m12["rot"] / 1e3, M_rotor_kNm_20=m20["rot"] / 1e3,
               ext_pct_k3=100 * m12["ext"] / m12["rot"], vawt_pct_k3=100 * m12["vawt"] / m12["rot"],
               ext_pct_20=100 * m20["ext"] / m20["rot"], vawt_pct_20=100 * m20["vawt"] / m20["rot"],
               sav_pct_of_design=100 * m_sav_cut / m_design,
               ua_ext_max=lim)
    return out


if __name__ == "__main__":
    cases = [("central", rm.P()), ("conjunto A (w = 1)", rm.P(w=1.0)), ("conjunto B (w = 0)", rm.P(w=0.0)),
             ("f_X = 1,0", rm.P(cl_ext=1.0)), ("f_X = 0,6", rm.P(cl_ext=0.6)),
             ("Ct Savonius 1,0", rm.P(ct_s=1.0)), ("Ct Savonius 0,4", rm.P(ct_s=0.4)),
             ("nivel 0", rm.P(nivel0=True))]
    res = [report(p, lab) for lab, p in cases]
    for r in res:
        print(f"{r['caso']:20s} M rotor {r['M_rotor_kNm_12']:6.0f} kN m (12 m/s, k 3), {r['M_rotor_kNm_20']:6.0f} (20 m/s); "
              f"prolongación +{r['ext_pct_k3']:4.1f} % (k 3) +{r['ext_pct_20']:4.1f} % (20 m/s); "
              f"VAWT +{r['vawt_pct_k3']:4.1f} % / +{r['vawt_pct_20']:4.1f} %; Savonius a 20 m/s "
              f"{r['sav_pct_of_design']:5.1f} % del de diseño; prolongación hasta Ua {r['ua_ext_max']:.2f} m/s")
    json.dump(res, open("momentos.json", "w"), indent=1, ensure_ascii=False)
