#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Indicadores estructurales de A1 (apartado 4.6). Versión 4 (22-10-2026): límites de giro del parámetro p (k_max,
U_s,max), rotor de 30 m (R1-30), calado aéreo absoluto y momento escorante del buque.

Momento flector en la base de una unidad, en operación, con el mismo viento aparente en toda la altura:
- rotor girando a k = min(k_max, U_s,max/U_a): resultante q A_F (CL^2 + CD^2)^0,5 aplicada a H_F/2;
- prolongación (T3, T3x) a k_X = k D_X/D_F, con el factor f_X en la sustentación, aplicada a H_F + H_X/2;
- tramo superior del rotor de 30 m (R1-30), con el diámetro del rotor y el factor f30, en H_F + H_X/2;
- VAWT (T1, T2, T4): resultante de su empuje y su fuerza lateral (modelo DMST de ../turbina) en el estado de mayor
  fuerza permitido a ese viento, aplicada en H_F + H_X/2;
- rotor abierto en Savonius (T2, T3, T5): T_S = C_t,S q D_S H_F en H_F/2; su sustentación no se conoce.
Las fuerzas de rotor y dispositivo se suman en módulo (cota superior). Referencia: momento del rotor solo en el
límite de operación supuesto, U_a = 20 m/s. Se calcula también el viento aparente máximo con el que la prolongación
cabe en esa envolvente. La supervivencia con el rotor parado en temporal no se evalúa.
Calado aéreo: base de los rotores a Z_DECK = Z_ROT - H_F/2 = 28 m sobre la flotación (supuesto, el mismo que fija el
centro del rotor a 40 m): 52 m en lo alto del rotor y 58 m con la prolongación, la VAWT o el rotor de 30 m. El perfil
de puerto de Portsmouth no da límite de calado aéreo (Cruise Britain) y el atraque de Bilbao está en el puerto
exterior, en Santurtzi, aguas afuera de los puentes de la ría (Autoridad Portuaria de Bilbao).
Momento escorante: fuerza lateral de las cuatro unidades con el viento de través aparente a U_a = 20 m/s, aplicada a
la altura de su centro respecto al centro de resistencia lateral (a T/2 bajo la flotación); se da también la escora
estática con un desplazamiento y una altura metacéntrica supuestos (no hay datos públicos de estabilidad).
"""
import json
import math
import numpy as np
import ruta_modelo as rm

UA_DESIGN = 20.0
Z_DECK = rm.Z_ROT - rm.H_F / 2.0            # 28 m: base de los rotores sobre la flotación (supuesto)
DISPL_T = (20000.0, 24000.0)               # t, supuesto: L_pp ~ 200 m x 27,8 m x 6,4 m x C_B 0,55-0,65
GM_M = (1.5, 2.5)                          # m, altura metacéntrica supuesta


def coeffs(k, p):
    cl, cd = rm.rotor_coeffs(np.array([k]), p)
    return float(cl[0]), float(cd[0])


def k_op(u_a, p):
    return min(p.k_max, p.us_max / u_a)


def moments(u_a, p):
    """Momentos en la base (N m) por unidad y fuerza lateral del rotor (N): rotor solo, prolongación, VAWT,
    Savonius abierto y tramo superior del rotor de 30 m."""
    q = 0.5 * rm.RHO * u_a ** 2
    k = k_op(u_a, p)
    cl, cd = coeffs(k, p)
    m_rot = q * rm.A_F * math.hypot(cl, cd) * rm.H_F / 2
    clx, cdx = coeffs(k * rm.D_X / rm.D_F, p)
    m_ext = q * rm.A_X * math.hypot(p.cl_ext * clx, cdx) * (rm.H_F + rm.H_X / 2)
    m_30 = q * rm.A_UP * math.hypot(p.f30 * cl, cd) * (rm.H_F + rm.H_X / 2)
    if p.turb == "dmst":
        _, th, sf = rm.vawt_states(u_a, p)
        f_vawt = float(np.max(np.hypot(th, sf)))
    else:
        f_vawt = q * rm.A_X * p.cp_t / p.cpct_t
    m_vawt = f_vawt * (rm.H_F + rm.H_X / 2)
    m_sav = p.ct_s * q * p.ds * rm.A_F * rm.H_F / 2
    return dict(k=k, rot=m_rot, ext=m_ext, r30=m_30, vawt=m_vawt, sav=m_sav,
                f_rot=q * rm.A_F * math.hypot(cl, cd), f_ext=q * rm.A_X * math.hypot(p.cl_ext * clx, cdx),
                f_30=q * rm.A_UP * math.hypot(p.f30 * cl, cd))


def ext_limit(p):
    m_design = moments(UA_DESIGN, p)["rot"]
    for u_a in np.arange(1.0, UA_DESIGN, 0.01):
        m = moments(u_a, p)
        if m["rot"] + m["ext"] > m_design:
            return round(float(u_a) - 0.01, 2), m_design
    return UA_DESIGN, m_design


def heeling(p, u_a=UA_DESIGN):
    """Momento escorante (MN m) de las unidades a U_a de través y escora estática (grados) para el rango supuesto
    de desplazamiento y GM: rotor solo (R1) y con la prolongación o el rotor de 30 m."""
    m = moments(u_a, p)
    arm_rot = rm.Z_ROT + p.t_draught / 2.0
    arm_top = Z_DECK + rm.H_F + rm.H_X / 2.0 + p.t_draught / 2.0
    # versión 5: (a) cota con la resultante de cada unidad; (b) viento aparente de través, con la fuerza lateral
    # de los rotores igual a su resistencia (C_D q A)
    q = 0.5 * rm.RHO * u_a ** 2
    k = k_op(u_a, p)
    _, cd = coeffs(k, p)
    _, cdx = coeffs(k * rm.D_X / rm.D_F, p)
    lat = dict(rot=q * rm.A_F * cd, ext=q * rm.A_X * cdx, r30=q * rm.A_UP * cd)
    out = {}
    for lab, mom, mom_t in (("R1", m["f_rot"] * arm_rot, lat["rot"] * arm_rot),
                            ("T3x", m["f_rot"] * arm_rot + m["f_ext"] * arm_top, lat["rot"] * arm_rot + lat["ext"] * arm_top),
                            ("R1-30", m["f_rot"] * arm_rot + m["f_30"] * arm_top, lat["rot"] * arm_rot + lat["r30"] * arm_top)):
        mh, mt = p.n_units * mom, p.n_units * mom_t
        heel = [math.degrees(math.asin(min(1.0, mh / (dsp * 1e3 * 9.81 * gm)))) for dsp in DISPL_T for gm in GM_M]
        heel_t = [math.degrees(math.asin(min(1.0, mt / (dsp * 1e3 * 9.81 * gm)))) for dsp in DISPL_T for gm in GM_M]
        out[lab] = dict(M_MNm=mh / 1e6, escora_grados=(min(heel), max(heel)), M_traves_MNm=mt / 1e6,
                        escora_traves_grados=(min(heel_t), max(heel_t)))
    return out


def report(p, label):
    k_lim = p.us_max / p.k_max                # viento aparente a partir del cual manda U_s,max
    m_k = moments(k_lim, p)
    m20 = moments(UA_DESIGN, p)
    lim, m_design = ext_limit(p)
    m_sav_cut = moments(rm.UA_CUT_S, p)["sav"]
    return dict(caso=label, ua_kmax=k_lim,
                M_rotor_kNm_kmax=m_k["rot"] / 1e3, M_rotor_kNm_20=m20["rot"] / 1e3,
                ext_pct_kmax=100 * m_k["ext"] / m_k["rot"], vawt_pct_kmax=100 * m_k["vawt"] / m_k["rot"],
                r30_pct_kmax=100 * m_k["r30"] / m_k["rot"],
                ext_pct_20=100 * m20["ext"] / m20["rot"], vawt_pct_20=100 * m20["vawt"] / m20["rot"],
                r30_pct_20=100 * m20["r30"] / m20["rot"],
                sav_pct_of_design=100 * m_sav_cut / m_design, ua_ext_max=lim,
                escorante=heeling(p))


def air_draft():
    return {"R1": Z_DECK + rm.H_F, "T3x, T3, T1, T2, T4, R1-30": Z_DECK + rm.H_F + rm.H_X}


if __name__ == "__main__":
    import coeficientes as cf
    cases = [("A (referencia)", rm.P()), ("A+ (k <= 4, U_s <= 46 m/s)", rm.P(**cf.ESCENARIOS_COEF["A+"])),
             ("A-VG (ajustado al Viking Grace)", rm.P(**cf.ESCENARIOS_COEF["A-VG"])),
             ("AB (versión 3)", rm.P(**cf.ESCENARIOS_COEF["AB"])), ("B", rm.P(**cf.ESCENARIOS_COEF["B"])),
             ("f_X = 1,30", rm.P(cl_ext=1.30)), ("f_X = 0,90", rm.P(cl_ext=0.90)), ("f_X = 0,5 (placa abajo)", rm.P(cl_ext=0.5)),
             ("Ct Savonius 1,0", rm.P(ct_s=1.0)), ("Ct Savonius en la cota", rm.P(ct_s=cf.KAPPA_S * cf.CENTRAL["cp_s"]))]
    res = [report(p, lab) for lab, p in cases]
    for r in res:
        print(f"{r['caso']:26s} M rotor {r['M_rotor_kNm_kmax']:6.0f} kN m (Ua {r['ua_kmax']:.0f} m/s, k_max), "
              f"{r['M_rotor_kNm_20']:6.0f} (20 m/s); prolongación +{r['ext_pct_kmax']:4.1f} % / +{r['ext_pct_20']:4.1f} %; "
              f"30 m +{r['r30_pct_kmax']:4.1f} % / +{r['r30_pct_20']:4.1f} %; VAWT +{r['vawt_pct_kmax']:4.1f} % / "
              f"+{r['vawt_pct_20']:4.1f} %; Savonius {r['sav_pct_of_design']:5.1f} % del de diseño; prolongación hasta "
              f"Ua {r['ua_ext_max']:.2f} m/s")
        e = r["escorante"]
        print("    momento escorante a 20 m/s (cota con la resultante | de través): " + "; ".join(
            f"{k} {v['M_MNm']:.1f} MN m ({v['escora_grados'][0]:.1f}-{v['escora_grados'][1]:.1f} grados) | "
            f"{v['M_traves_MNm']:.1f} MN m ({v['escora_traves_grados'][0]:.1f}-{v['escora_traves_grados'][1]:.1f})"
            for k, v in e.items()))
    print("calado aéreo (m sobre la flotación):", air_draft())
    json.dump(dict(casos=res, calado_aereo=air_draft(), z_deck=Z_DECK, desplazamiento_t=DISPL_T, gm_m=GM_M),
              open("momentos.json", "w"), indent=1, ensure_ascii=False)
