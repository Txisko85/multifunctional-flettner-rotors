#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Análisis de la ruta, versión 4 (22-10-2026), tras la revisión simulada: caso de referencia (servicio documentado
de abril a noviembre, coeficientes del escenario A), comparadores no patentados, desglose por travesía, paso a paso
desde la versión 3, escenarios epistémicos, sensibilidades de un factor, umbral de interacción y reglas de la
conversión en Savonius (resultados_ruta.json)."""
import json
import coeficientes as cf
import ruta_modelo as rm
import ruta_momentos as mom
from ruta_modelo import P, annual, table, break_even, load_wind, mean_benefits, crossing_speed

d, groups, nu, ng = load_wind()
base = annual(d, groups, nu, ng)
KEYS_SEA = ("R1", "T3x", "T3", "T4", "T1", "T2", "T5", "R1t", "R1-30", "_pen R1", "_pen propio", "_pen cruzado",
            "_R1 empuje", "_R1 giro", "_R1 parado", "_R0 viento", "T5 literal", "T5 inventor", "T5 popa",
            "T2 literal", "T2 inventor", "T3 literal", "T3 inventor")


def summ(r, extra=None):
    t = {row[0]: row[3] for row in table(r)}
    out = {"R0": round(r["fuel_total"], 1), "R0_mar": round(r["fuel_sea"], 1), "R0_atraque": round(r["fuel_berth"], 1),
           "R1": round(r["sea"]["R1"], 2), "R1_pct": round(100 * r["sea"]["R1"] / r["fuel_total"], 3),
           "pen_R1": round(r["sea"]["_pen R1"], 2)}
    for k in ("T3x", "T3", "T4", "T1", "T2", "T5"):
        out[k] = round(t[k], 2)
    for k in ("R1t", "R1-30", "T5 literal", "T5 inventor", "T5 popa"):
        out[k] = round(r["sea"][k], 2)
    if extra:
        out.update(extra)
    return out


# ---------------- caso de referencia ----------------
r5 = annual(d, groups, nu, ng, P(n_units=5))
res = {"version": 4,
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
    for s, lab in ((0, "abr-nov"), (1, "dic-mar")):
        mb = mean_benefits(groups, nu, ng, s, direction, v, P())
        row[lab] = {k: round(mb[k], 2) for k in ("R1", "_R1 empuje", "_R1 giro", "_R1 parado", "_pen R1",
                                                    "_R1 marcha", "T3x", "T4", "T1", "T5", "R1t", "R1-30")}
        row[lab]["_R0 viento kW buque"] = round(mb["_R0 viento"], 1)
    per.append(row)
res["por_travesia_kW_unidad"] = per

# ---------------- umbral de interacción (T4 = T3x) ----------------
res["eps_star"] = break_even(d, groups, nu, ng, base)

# ---------------- paso a paso desde la versión 3 ----------------
V3_CT_S, V3_CL_EXT = cf.CENTRAL["ct_s"], cf.CENTRAL["cl_ext"]
steps = [
    ("versión 3 (21-10-2026): año completo, mezcla AB", P(v3=True), (0, 1)),
    ("+ servicio documentado (abril-noviembre)", P(v3=True), (0,)),
    ("+ calado 6,4 m y hotel en atraque 2 000 kW (fuentes primarias, EU MRV)",
     P(v3=True, t_draught=rm.T_DRAUGHT, hotel_berth=rm.P_HOTEL_BERTH), (0,)),
    ("+ viento sobre el buque en R0, SFOC según la carga y margen de mar calibrado a igual velocidad",
     P(w=0.5, ct_s=V3_CT_S, cl_ext=V3_CL_EXT), (0,)),
    ("+ cota del Savonius en el centroide del cangilón (C_t,S central 0,77)", P(w=0.5, cl_ext=V3_CL_EXT), (0,)),
    ("+ f_X de la línea sustentadora (central 1,10)", P(w=0.5), (0,)),
    ("+ escenario de coeficientes A (versión 4, referencia)", P(), (0,)),
]
res["paso_a_paso"] = {lab: summ(annual(d, groups, nu, ng, p, seasons=ss)) for lab, p, ss in steps}

# ---------------- escenarios epistémicos ----------------
UA_LIM = mom.ext_limit(P())[0]
res["ua_ext_max_estructural"] = UA_LIM
scen = {f"coeficientes {k}": (P(**v), None) for k, v in cf.ESCENARIOS_COEF.items()}
scen.update({
    "año completo (dos rotaciones todo el año)": (P(), (0, 1)),
    "diciembre-marzo solo": (P(), (1,)),
    "exponente 4 de la curva potencia-velocidad": (P(n_prop=4.0), None),
    "cinco unidades": (P(n_units=5), None),
    "velocidad -2 kn": (P(dv=-2.0), None),
    "velocidad -4 kn": (P(dv=-4.0), None),
    "velocidad -6 kn": (P(dv=-6.0), None),
    "interferencia en tándem: -10 % de sustentación a sotavento (|beta| < 30°)": (P(int_frac=0.10), None),
    "interferencia en tándem: -20 % de sustentación a sotavento (|beta| < 30°)": (P(int_frac=0.20), None),
    "e = 1,3 (equilibrio de deriva y guiñada, unidades centradas)": (P(e_span=1.3), None),
    "e = 0,72 (equilibrio de deriva y guiñada, unidades a popa)": (P(e_span=0.72), None),
    "f_X 0,5 (placa de extremo sin trasladar)": (P(cl_ext=0.5), None),
    "f_X 0,95": (P(cl_ext=0.95), None),
    "f_X 1,25": (P(cl_ext=1.25), None),
    "f30 0,95": (P(f30=0.95), None),
    "f30 1,15": (P(f30=1.15), None),
    f"límite estructural de la prolongación Ua <= {UA_LIM:.2f} m/s": (P(ua_ext_max=UA_LIM), None),
    "interacción eps = -0,01": (P(eps=-0.01), None),
    "interacción eps = +0,01": (P(eps=0.01), None),
    "fuerza lateral del Savonius C_y,S = +0,5": (P(cy_s=0.5), None),
    "fuerza lateral del Savonius C_y,S = -0,5": (P(cy_s=-0.5), None),
    "viento en la turbina x0,9 (frenado por el rotor)": (P(u_turb=0.9), None),
    "viento en la turbina x1,1 (aceleración sobre la placa)": (P(u_turb=1.1), None),
    "kappa 1,2 (grupos auxiliares)": (P(kappa=1.2), None),
    "sin resistencia del rotor parado": (P(stop_drag=False), None),
    "sin resistencia inducida": (P(sir=False), None),
    "relación de giro fija k = min(k_max, U_s,max/Ua)": (P(k_opt=False), None),
    "Savonius sin ensanchar": (P(ds=1.0), None),
    "turbina supuesta (versión 2: Cp 0,25, Cp/Ct 0,4)": (P(turb="supuesta"), None),
    "turbina girando en sentido contrario al rotor": (P(corot=False), None),
})
res["escenarios"] = {k: summ(annual(d, groups, nu, ng, p, seasons=ss)) for k, (p, ss) in scen.items()}

# ---------------- sensibilidades de un factor (extremos de los rangos muestreados) ----------------
sens = {}
for k, (lo, hi) in cf.RANGES_V4_MC.items():
    for lab, val in (("min", lo), ("max", hi)):
        sens[f"{k} = {val:g}"] = summ(annual(d, groups, nu, ng, P(**{k: val})))
res["sensibilidad"] = sens
json.dump(res, open("resultados_ruta.json", "w"), indent=1, ensure_ascii=False)

# ---------------- salida ----------------
b = res["base"]
print(f"Caso de referencia (escenario A, abril-noviembre): R0 {b['fuel_total']:.0f} t (mar {b['fuel_sea']:.0f}, "
      f"atraque {b['fuel_berth']:.0f}); carga media {b['load_mean']:.2f}")
print(f"R1 {b['sea']['R1']:.1f} t ({100 * b['sea']['R1'] / b['fuel_total']:.2f} %): empuje {b['sea']['_R1 empuje']:.1f}, "
      f"giro {b['sea']['_R1 giro']:.1f}, parado {b['sea']['_R1 parado']:.1f}, inducida {b['sea']['_pen R1']:.1f} "
      f"(propia {b['sea']['_pen propio']:.1f}, cruzada {b['sea']['_pen cruzado']:.1f}); viento sobre el buque "
      f"{b['sea']['_R0 viento']:.0f} t")
for r in table(base):
    print(f"  {r[0]:4s} mar {r[1]:7.1f}  puerto {r[2]:5.1f}  total {r[3]:7.1f}  ({r[4]:.2f} %)")
print(f"  comparadores: R1t {b['sea']['R1t']:.1f}, R1-30 {b['sea']['R1-30']:.1f}, quinta unidad "
      f"{b['R1_quinta_unidad']:.1f} t")
print("  reglas:", {k: round(b['sea'][k], 1) for k in ("T5 literal", "T5 inventor", "T5 popa", "T2 literal",
                                                       "T2 inventor", "T3 literal", "T3 inventor")})
print("eps* =", round(res["eps_star"], 3))
print("\nPASO A PASO")
for k, v in res["paso_a_paso"].items():
    print(f"  {k[:70]:70s} R0 {v['R0']:7.0f}  R1 {v['R1']:6.1f} ({v['R1_pct']:.2f} %)  T3x {v['T3x']:5.1f}  T3 {v['T3']:5.1f}  "
          f"T5 {v['T5']:5.1f}  T4 {v['T4']:4.1f}  T1 {v['T1']:6.1f}  R1t {v['R1t']:6.1f}  R1-30 {v['R1-30']:6.1f}")
print("\nPOR TRAVESIA (kW por unidad, abril-noviembre)")
for row in per:
    print(f"  {row['travesia']:28s} {row['v_kn']:5.2f} kn  {row['abr-nov']}")
print("\nESCENARIOS")
for k, v in res["escenarios"].items():
    print(f"  {k[:58]:58s} R0 {v['R0']:7.0f} R1 {v['R1']:6.1f} ({v['R1_pct']:.2f} %) T3x {v['T3x']:5.1f} T3 {v['T3']:5.1f} "
          f"T5 {v['T5']:5.1f} T4 {v['T4']:5.1f} T1 {v['T1']:6.1f} T2 {v['T2']:6.1f} R1t {v['R1t']:6.1f} R1-30 {v['R1-30']:6.1f}")
print("\nSENSIBILIDAD DE UN FACTOR")
for k, v in sens.items():
    print(f"  {k:22s} R1 {v['R1']:6.1f} T3x {v['T3x']:5.1f} T3 {v['T3']:5.1f} T5 {v['T5']:5.1f} T4 {v['T4']:5.1f} "
          f"T1 {v['T1']:6.1f} R1t {v['R1t']:6.1f} R1-30 {v['R1-30']:6.1f}")
