#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Análisis de la ruta (versión 3, 21-10-2026: turbina con modelo DMST): caso central, desglose por travesía y
estación, umbral de interacción, efecto de las correcciones respecto a la versión 1 y sensibilidades de un factor
(resultados_ruta.json)."""
import json
import ruta_modelo as rm
import ruta_momentos as mom
from ruta_modelo import P, annual, table, break_even, load_wind, mean_benefits, crossing_speed

d, groups, nu, ng = load_wind()
base = annual(d, groups, nu, ng)
res = {"version": 3,
       "base": {"fuel_total": base["fuel_total"], "fuel_me": base["fuel_me"], "R1": base["sea"]["R1"],
                "pen_R1": base["sea"]["_pen R1"], "rows": table(base),
                "reglas_mar": {k: base["sea"][k] for k in ("T5 literal", "T5 inventor", "T5 popa", "T2 literal",
                                                            "T2 inventor", "T3 literal", "T3 inventor")},
                "port_kw": base["port_kw"], "hours_sea": base["hours_sea"], "hours_port": base["hours_port"],
                "parametros": {k: v for k, v in vars(P()).items() if isinstance(v, (int, float, bool)) or v is None}}}

# Desglose por travesía (kW medios por unidad)
per = []
for name, direction, hours in rm.crossings(P()):
    v, t_open = crossing_speed(hours, P())
    row = {"travesia": name, "horas": hours, "v_kn": round(v, 2), "h_mar": round(t_open, 2)}
    for s, lab in ((0, "abr-nov"), (1, "dic-mar")):
        mb = mean_benefits(groups, nu, ng, s, direction, v, P())
        row[lab] = {k: round(mb[k], 2) for k in ("R1", "_pen R1", "_R1 marcha", "T3x", "T4", "T1", "T5")}
    per.append(row)
res["por_travesia_kW_unidad"] = per

# Solo temporada de servicio 2026 (abril-noviembre)
season = annual(d, groups, nu, ng, seasons=(0,))
res["solo_abr_nov"] = {"fuel_total": season["fuel_total"], "R1": season["sea"]["R1"], "rows": table(season)}

# Umbral de interacción
res["eps_star"] = break_even(d, groups, nu, ng, base)


def summ(r):
    t = {row[0]: row[3] for row in table(r)}
    return {"R0": round(r["fuel_total"], 0), "R1": round(r["sea"]["R1"], 1), "pen_R1": round(r["sea"]["_pen R1"], 1),
            "T3x": round(t["T3x"], 1), "T3": round(t["T3"], 1), "T4": round(t["T4"], 1), "T1": round(t["T1"], 1),
            "T2": round(t["T2"], 1), "T5": round(t["T5"], 1),
            "T5_literal_mar": round(r["sea"]["T5 literal"], 1), "T5_inventor_mar": round(r["sea"]["T5 inventor"], 1),
            "T5_popa_mar": round(r["sea"]["T5 popa"], 1)}


# Paso a paso desde la versión 1
res["paso_a_paso"] = {k: summ(annual(d, groups, nu, ng, p)) for k, p in {
    "versión 1 (18-10-2026)": P(v1=True),
    "+ horario corregido": P(sir=False, turb="supuesta"),
    "+ resistencia inducida de los rotores": P(windage=False, turb="supuesta"),
    "+ término cruzado con el viento sobre el buque (versión 2)": P(turb="supuesta"),
    "+ turbina con modelo DMST (versión 3)": P(),
}.items()}

UA_LIM = mom.ext_limit(P())[0]
res["ua_ext_max_estructural"] = UA_LIM
cases = {
    "central": P(),
    "conjunto A (De Marco 2016, w = 1)": P(w=1.0),
    "conjunto B (Bordogna 2019 + inducida, w = 0)": P(w=0.0),
    "relacion de giro fija k = min(3, 36/Ua)": P(k_opt=False),
    "CD rotor x0,8": P(cd_mult=0.8),
    "CD rotor x1,2": P(cd_mult=1.2),
    "CD rotor parado 0,4": P(cd0=0.4),
    "CD rotor parado 0,7": P(cd0=0.7),
    "Cf giro 0,002": P(cf=0.002),
    "Cf giro 0,007": P(cf=0.007),
    "Ct Savonius 0,4": P(ct_s=0.4),
    "Ct Savonius 1,0": P(ct_s=1.0),
    "Cp Savonius 0,14": P(cp_s=0.14),
    "Cp Savonius 0,24": P(cp_s=0.24),
    "CL prolongacion x0,6": P(cl_ext=0.6),
    "CL prolongacion x1,0": P(cl_ext=1.0),
    "turbina Cp x0,7": P(f_cp_t=0.7),
    "turbina Cp x0,9": P(f_cp_t=0.9),
    "turbina empuje x1,0": P(f_ct_t=1.0),
    "turbina empuje x1,1": P(f_ct_t=1.1),
    "turbina fuerza lateral x0,6": P(f_cs_t=0.6),
    "turbina fuerza lateral x1,0": P(f_cs_t=1.0),
    "turbina sin fuerza lateral": P(f_cs_t=0.0),
    "turbina girando en sentido contrario al rotor": P(corot=False),
    "velocidad periferica de la turbina 30 m/s": P(v_tip_t=30.0),
    "velocidad periferica de la turbina 50 m/s": P(v_tip_t=50.0),
    "turbina supuesta (version 2: Cp 0,25, Cp/Ct 0,4)": P(turb="supuesta"),
    f"limite estructural Ua <= {UA_LIM:.2f}": P(ua_ext_max=UA_LIM),
    "kappa 1,2 (grupos auxiliares)": P(kappa=1.2),
    "Savonius sin ensanchar": P(ds=1.0),
    "velocidad -2 kn": P(dv=-2.0),
    "sin resistencia del rotor parado": P(stop_drag=False),
    "eficacia e 0,4": P(e_span=0.4),
    "eficacia e 0,8": P(e_span=0.8),
    "A_L C_Y 3000 m2": P(ayc=3000.0),
    "A_L C_Y 6000 m2": P(ayc=6000.0),
    "sin viento sobre el buque": P(windage=False),
    "sin resistencia inducida": P(sir=False),
    "nivel 0 (calculo anterior)": P(nivel0=True, sir=False),
    "version 1 (18-10-2026)": P(v1=True),
}
res["sensibilidad"] = {k: summ(annual(d, groups, nu, ng, p)) for k, p in cases.items()}
json.dump(res, open("resultados_ruta.json", "w"), indent=1, ensure_ascii=False)

print(f"R0 {base['fuel_total']:.0f} t/año; R1 {base['sea']['R1']:.1f} t/año "
      f"({100 * base['sea']['R1'] / base['fuel_total']:.2f} %); inducida pagada por R1 {base['sea']['_pen R1']:.1f}")
for r in table(base):
    print(f"  {r[0]:4s} mar {r[1]:7.1f}  puerto {r[2]:5.1f}  total {r[3]:7.1f}  ({r[4]:.2f} %)")
print("reglas", {k: round(v, 1) for k, v in res["base"]["reglas_mar"].items()})
print("eps* =", round(res["eps_star"], 3))
print("solo abr-nov: R1", round(season["sea"]["R1"], 1), {r[0]: round(r[3], 1) for r in table(season)})
print("\nPASO A PASO")
for k, v in res["paso_a_paso"].items():
    print(f"  {k:58s} {v}")
print("\nPOR TRAVESIA (kW por unidad)")
for row in per:
    print(f"  {row['travesia']:28s} {row['horas']:5.2f} h {row['v_kn']:5.2f} kn  verano {row['abr-nov']}  invierno {row['dic-mar']}")
print("\nSENSIBILIDAD")
for k, v in res["sensibilidad"].items():
    print(f"  {k:46s} R1 {v['R1']:6.1f} pen {v['pen_R1']:5.1f}  T3x {v['T3x']:6.1f}  T4 {v['T4']:4.1f}  T1 {v['T1']:6.1f}  "
          f"T5 {v['T5']:5.1f}  T5lit {v['T5_literal_mar']:6.1f}")
