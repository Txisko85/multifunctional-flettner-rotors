#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Análisis de la ruta, versión 5 (27-10-2026), tras la segunda revisión simulada: caso de referencia (temporada de
abril a octubre del horario de 2026, histogramas condicionados al horario, resistencia inducida del equilibrio de
deriva y guiñada, coeficientes del escenario A), escenario A-VG ajustado al Viking Grace, comparadores no patentados
(rotor abatible, de 30 m, de 5 m de diámetro, quinta unidad y abatible con prolongación), envolvente estructural común,
desglose por travesía, paso a paso desde la versión 4, escenarios, sensibilidades de un factor, umbral de interacción,
reglas de la conversión en Savonius y su punto de equilibrio en C_t,S (resultados_ruta.json).
Versión 4 (22-10-2026): servicio de abril a noviembre con los histogramas de todas las horas (P(v4=True) y
era5_ruta.json la reproducen)."""
import json
import coeficientes as cf
import ruta_modelo as rm
import ruta_momentos as mom
from ruta_modelo import P, annual, table, break_even, load_wind, mean_benefits, crossing_speed

d, groups, nu, ng = load_wind()
base = annual(d, groups, nu, ng)
KEYS_SEA = ("R1", "T3x", "T3", "T4", "T1", "T2", "T5", "R1t", "R1-30", "R1-D5", "R1t+X", "_pen R1", "_pen propio",
            "_pen cruzado",
            "_R1 empuje", "_R1 giro", "_R1 parado", "_R0 viento", "T5 literal", "T5 inventor", "T5 popa",
            "T2 literal", "T2 inventor", "T3 literal", "T3 inventor")


def summ(r, extra=None):
    t = {row[0]: row[3] for row in table(r)}
    out = {"R0": round(r["fuel_total"], 1), "R0_mar": round(r["fuel_sea"], 1), "R0_atraque": round(r["fuel_berth"], 1),
           "R1": round(r["sea"]["R1"], 2), "R1_pct": round(100 * r["sea"]["R1"] / r["fuel_total"], 3),
           "pen_R1": round(r["sea"]["_pen R1"], 2)}
    for k in ("T3x", "T3", "T4", "T1", "T2", "T5"):
        out[k] = round(t[k], 2)
    for k in ("R1t", "R1-30", "R1-D5", "R1t+X", "T5 literal", "T5 inventor", "T5 popa"):
        out[k] = round(r["sea"][k], 2)
    out["T5_mar"] = round(r["sea"]["T5"], 2)
    out["empuje"], out["giro"], out["parado"] = (round(r["sea"][k], 2) for k in ("_R1 empuje", "_R1 giro", "_R1 parado"))
    out["pen_propio"], out["pen_cruzado"] = round(r["sea"]["_pen propio"], 2), round(r["sea"]["_pen cruzado"], 2)
    if extra:
        out.update(extra)
    return out


# ---------------- caso de referencia ----------------
r5 = annual(d, groups, nu, ng, P(n_units=5))
res = {"version": 5,
       "base": {"fuel_total": base["fuel_total"], "fuel_sea": base["fuel_sea"], "fuel_berth": base["fuel_berth"],
                "fuel_me": base["fuel_me"], "load_mean": base["load_mean"], "seasons": base["seasons"],
                "sea": {k: base["sea"][k] for k in KEYS_SEA}, "port": base["port"], "rows": table(base),
                "port_kw": base["port_kw"], "hours_sea": base["hours_sea"], "hours_port": base["hours_port"],
                "R1_quinta_unidad": r5["sea"]["R1"] - base["sea"]["R1"], "R1_cinco_unidades": r5["sea"]["R1"],
                "parametros": {k: v for k, v in vars(P()).items() if isinstance(v, (int, float, bool)) or v is None}}}

# ---------------- desglose por travesía (kW medios por unidad, mar abierto) ----------------
per = []
for name, direction, hours in rm.crossings(P()):
    v, t_open = crossing_speed(hours, P())
    row = {"travesia": name, "horas": hours, "v_kn": round(v, 2), "h_mar": round(t_open, 2)}
    for s_, lab in ((0, "abr-oct"), (2, "nov"), (1, "dic-mar")):
        mb = mean_benefits(groups, nu, ng, s_, direction, v, P(), name=name)
        row[lab] = {k: round(mb[k], 2) for k in ("R1", "_R1 empuje", "_R1 giro", "_R1 parado", "_pen R1",
                                                    "_R1 marcha", "T3x", "T4", "T1", "T5", "R1t", "R1-30")}
        row[lab]["_R0 viento kW buque"] = round(mb["_R0 viento"], 1)
    per.append(row)
res["por_travesia_kW_unidad"] = per

# ---------------- umbral de interacción (T4 = T3x) ----------------
res["eps_star"] = break_even(d, groups, nu, ng, base)

# ---------------- paso a paso desde la versión 4 ----------------
d4, g4, _, _ = load_wind("era5_ruta.json")
steps = [
    ("versión 4 (26-10-2026): abril-noviembre, histogramas de todas las horas", (d4, g4), P(v4=True), (0,)),
    ("+ histogramas condicionados al horario (abril-noviembre)", (d, groups), P(v4=True), rm.NOV_SEASONS),
    ("+ temporada abril-octubre", (d, groups), P(v4=True), None),
    ("+ resistencia inducida del equilibrio de deriva y guiñada", (d, groups), P(v4=True, sir_model="mmg"), None),
    ("+ margen de mar recalibrado (0,077)", (d, groups),
     P(v4=True, sir_model="mmg", sea_margin=rm.SEA_MARGIN[3.0]), None),
    ("+ cota estricta del Savonius (C_t,S central 0,62)", (d, groups),
     P(v4=True, sir_model="mmg", sea_margin=rm.SEA_MARGIN[3.0], ct_s=cf.CENTRAL_V4["ct_s"]), None),
    ("+ f30 contrastado con De Marco et al. (central 1,025): versión 5", (d, groups), P(), None),
]
res["paso_a_paso"] = {lab: summ(annual(dd[0], dd[1], nu, ng, p, seasons=ss)) for lab, dd, p, ss in steps}

# ---------------- escenarios epistémicos ----------------
UA_LIM = mom.ext_limit(P())[0]
res["ua_ext_max_estructural"] = UA_LIM
scen = {f"coeficientes {k}": (P(**v), None) for k, v in cf.ESCENARIOS_COEF.items()}
scen.update({
    "abril-noviembre (con noviembre)": (P(), rm.NOV_SEASONS),
    "año completo (dos rotaciones todo el año)": (P(), rm.ALL_SEASONS),
    "noviembre solo": (P(), (2,)),
    "diciembre-marzo solo": (P(), (1,)),
    "exponente 4 de la curva potencia-velocidad": (P(n_prop=4.0), None),
    "cinco unidades": (P(n_units=5), None),
    "velocidad -2 kn": (P(dv=-2.0), None),
    "velocidad -4 kn": (P(dv=-4.0), None),
    "velocidad -6 kn": (P(dv=-6.0), None),
    "interferencia en tándem: -10 % de sustentación a sotavento (|beta| < 30°)": (P(int_frac=0.10), None),
    "interferencia en tándem: -20 % de sustentación a sotavento (|beta| < 30°)": (P(int_frac=0.20), None),
    "interferencia en tándem: -10 % a sotavento (|beta| < 60°)": (P(int_frac=0.10, int_beta=60.0), None),
    "interferencia en tándem: -20 % a sotavento (|beta| < 60°)": (P(int_frac=0.20, int_beta=60.0), None),
    "inducida: ala de bajo alargamiento, e = 0,6 (versión 4)": (P(sir_model="ala", e_span=0.6), None),
    "inducida: ala de bajo alargamiento, e = 0,4": (P(sir_model="ala", e_span=0.4), None),
    "inducida: unidades con el centroide a 0,2 L a popa": (P(x_u=-0.2), None),
    "inducida: unidades con el centroide a 0,2 L a proa": (P(x_u=0.2), None),
    "inducida: c_H = 0": (P(c_h=0.0), None),
    "inducida: c_H = 0,25 (KVLCC2, Yasukawa y Yoshimura 2015)": (P(c_h=0.25), None),
    "inducida: c_H = 1": (P(c_h=1.0), None),
    "inducida: viento sobre el buque centrado en x_s = 0,05 L": (P(x_s=0.05), None),
    "inducida: viento sobre el buque centrado en x_s = 0,10 L": (P(x_s=0.10), None),
    "f_X 0,5 (placa de extremo sin trasladar)": (P(cl_ext=0.5), None),
    "f_X 0,90": (P(cl_ext=0.90), None),
    "f_X 1,30": (P(cl_ext=1.30), None),
    "f30 0,90": (P(f30=0.90), None),
    "f30 1,15": (P(f30=1.15), None),
    f"límite estructural de la prolongación Ua <= {UA_LIM:.2f} m/s (versión 4)": (P(ua_ext_max=UA_LIM), None),
    "envolvente común: momento del rotor de referencia a 20 m/s": (P(m_cap="ref"), None),
    "envolvente común: momento del rotor de 30 m a 20 m/s (cimentación reforzada)": (P(m_cap="30m"), None),
    "A-VG, envolvente común del rotor de referencia": (P(m_cap="ref", **cf.ESCENARIOS_COEF["A-VG"]), None),
    "A-VG, envolvente común del rotor de 30 m": (P(m_cap="30m", **cf.ESCENARIOS_COEF["A-VG"]), None),
    "interacción eps = -0,01": (P(eps=-0.01), None),
    "interacción eps = +0,01": (P(eps=0.01), None),
    "fuerza lateral del Savonius C_y,S = +0,5": (P(cy_s=0.5), None),
    "fuerza lateral del Savonius C_y,S = -0,5": (P(cy_s=-0.5), None),
    "viento en la turbina x0,9 (frenado por el rotor)": (P(u_turb=0.9), None),
    "viento en la turbina x1,17 (línea sustentadora)": (P(u_turb=1.17), None),
    "corriente en la turbina girada -4°": (P(turb_dang=-4.0), None),
    "corriente en la turbina girada +19°": (P(turb_dang=19.0), None),
    "viento en la turbina x1,17 y girada +19°": (P(u_turb=1.17, turb_dang=19.0), None),
    "kappa 1,2 (grupos auxiliares)": (P(kappa=1.2), None),
    "sin resistencia del rotor parado": (P(stop_drag=False), None),
    "sin resistencia inducida": (P(sir=False), None),
    "relación de giro fija k = min(k_max, U_s,max/Ua)": (P(k_opt=False), None),
    "Savonius sin ensanchar": (P(ds=1.0), None),
    "turbina supuesta (versión 2: Cp 0,25, Cp/Ct 0,4)": (P(turb="supuesta"), None),
    "turbina girando en sentido contrario al rotor": (P(corot=False), None),
})
res["escenarios"] = {k: summ(annual(d, groups, nu, ng, p, seasons=ss)) for k, (p, ss) in scen.items()}

# ---------------- conversión en Savonius: punto de equilibrio en C_t,S (versión 5, R2-1) ----------------
# C_t,S con el que la regla de la reivindicación 1 (T5 literal) deja de perder en el mar, y C_t,S por debajo del cual
# abrir el rotor en el mar de forma óptima ahorra al menos 1 t/año, para C_p,S = 0,14, 0,19 y 0,24; y el radio
# efectivo de los cangilones que da ese C_t,S (C_t,S = C_p,S R_S / (lambda r_ef)).


def _bisect(fun, lo, hi, n=30):
    flo = fun(lo)
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if (fun(mid) > 0) == (flo > 0):
            lo, flo = mid, fun(mid)
        else:
            hi = mid
    return 0.5 * (lo + hi)


sav = {}
for cp in (0.14, 0.19, 0.24):
    lit = lambda ct: annual(d, groups, nu, ng, P(cp_s=cp, ct_s=ct))["sea"]["T5 literal"]
    opt = lambda ct: annual(d, groups, nu, ng, P(cp_s=cp, ct_s=ct))["sea"]["T5"] - 1.0
    ct_lit, ct_opt = _bisect(lit, 0.05, 1.2), _bisect(opt, 0.05, 1.2)
    curve = {f"{ct:.2f}": {k: round(annual(d, groups, nu, ng, P(cp_s=cp, ct_s=ct))["sea"][k], 2)
                          for k in ("T5", "T5 literal", "T2 literal", "T3 literal")}
             for ct in (cp / cf.LAMBDA_S, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0)}
    sav[f"{cp:.2f}"] = dict(ct_equilibrio_reivindicacion1=ct_lit, ct_equilibrio_optima=ct_opt,
                            cota_estricta=cp / cf.LAMBDA_S, estimacion_centroide=cf.KAPPA_S_CENTROIDE * cp,
                            radio_efectivo_reiv1=cp * cf.R_S / (cf.LAMBDA_S * ct_lit),
                            radio_efectivo_optima=cp * cf.R_S / (cf.LAMBDA_S * ct_opt), curva_mar=curve)
res["savonius_equilibrio"] = sav

# ---------------- sensibilidades de un factor (extremos de los rangos muestreados) ----------------
sens = {}
for k, (lo, hi) in cf.RANGES_V5_MC.items():
    for lab, val in (("min", lo), ("max", hi)):
        sens[f"{k} = {val:g}"] = summ(annual(d, groups, nu, ng, P(**{k: val})))
res["sensibilidad"] = sens
json.dump(res, open("resultados_ruta.json", "w"), indent=1, ensure_ascii=False)

# ---------------- salida ----------------
b = res["base"]
print(f"Caso de referencia (escenario A, abril-octubre): R0 {b['fuel_total']:.0f} t (mar {b['fuel_sea']:.0f}, "
      f"atraque {b['fuel_berth']:.0f}); carga media {b['load_mean']:.2f}")
print(f"R1 {b['sea']['R1']:.1f} t ({100 * b['sea']['R1'] / b['fuel_total']:.2f} %): empuje {b['sea']['_R1 empuje']:.1f}, "
      f"giro {b['sea']['_R1 giro']:.1f}, parado {b['sea']['_R1 parado']:.1f}, inducida {b['sea']['_pen R1']:.1f} "
      f"(propia {b['sea']['_pen propio']:.1f}, cruzada {b['sea']['_pen cruzado']:.1f}); viento sobre el buque "
      f"{b['sea']['_R0 viento']:.0f} t")
for r in table(base):
    print(f"  {r[0]:4s} mar {r[1]:7.1f}  puerto {r[2]:5.1f}  total {r[3]:7.1f}  ({r[4]:.2f} %)")
print(f"  comparadores: R1t {b['sea']['R1t']:.1f}, R1-30 {b['sea']['R1-30']:.1f}, R1-D5 {b['sea']['R1-D5']:.1f}, "
      f"R1t+X {b['sea']['R1t+X']:.1f}, quinta unidad {b['R1_quinta_unidad']:.1f} t")
print("  reglas:", {k: round(b['sea'][k], 1) for k in ("T5 literal", "T5 inventor", "T5 popa", "T2 literal",
                                                       "T2 inventor", "T3 literal", "T3 inventor")})
print("eps* =", round(res["eps_star"], 3))
print("\nPASO A PASO")
for k, v in res["paso_a_paso"].items():
    print(f"  {k[:70]:70s} R0 {v['R0']:7.0f}  R1 {v['R1']:6.1f} ({v['R1_pct']:.2f} %)  T3x {v['T3x']:5.1f}  T3 {v['T3']:5.1f}  "
          f"T5 {v['T5']:5.1f}  T4 {v['T4']:4.1f}  T1 {v['T1']:6.1f}  R1t {v['R1t']:6.1f}  R1-30 {v['R1-30']:6.1f}")
print("\nPOR TRAVESIA (kW por unidad, abril-octubre)")
for row in per:
    print(f"  {row['travesia']:28s} {row['v_kn']:5.2f} kn  {row['abr-oct']}")
print("\nSAVONIUS: equilibrio en C_t,S")
for cp, v in res["savonius_equilibrio"].items():
    print(f"  C_p,S {cp}: reivindicación 1 {v['ct_equilibrio_reivindicacion1']:.3f} (r_ef {v['radio_efectivo_reiv1']:.2f} m), "
          f"óptima {v['ct_equilibrio_optima']:.3f} (r_ef {v['radio_efectivo_optima']:.2f} m); cota {v['cota_estricta']:.3f}, "
          f"centroide {v['estimacion_centroide']:.3f}")
print("\nESCENARIOS")
for k, v in res["escenarios"].items():
    print(f"  {k[:58]:58s} R0 {v['R0']:7.0f} R1 {v['R1']:6.1f} ({v['R1_pct']:.2f} %) T3x {v['T3x']:5.1f} T3 {v['T3']:5.1f} "
          f"T5 {v['T5']:5.1f} T4 {v['T4']:5.1f} T1 {v['T1']:6.1f} T2 {v['T2']:6.1f} R1t {v['R1t']:6.1f} R1-30 {v['R1-30']:6.1f}")
print("\nSENSIBILIDAD DE UN FACTOR")
for k, v in sens.items():
    print(f"  {k:22s} R1 {v['R1']:6.1f} T3x {v['T3x']:5.1f} T3 {v['T3']:5.1f} T5 {v['T5']:5.1f} T4 {v['T4']:5.1f} "
          f"T1 {v['T1']:6.1f} R1t {v['R1t']:6.1f} R1-30 {v['R1-30']:6.1f}")
