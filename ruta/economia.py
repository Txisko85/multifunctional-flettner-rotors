#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Valor económico y regulatorio de los ahorros (versión 4, 22-10-2026; revisión simulada R3-1, R3-2 y R3-m5).

Para cada topología y comparador, el ahorro anual del modelo de ruta (ruta_modelo.annual, caso de referencia) se
convierte en euros y en la inversión máxima que justificaría (valor actual a 20 años y 8 %). No se conocen costes de
inversión ni de mantenimiento de ninguna de las opciones: el resultado es un umbral, no una rentabilidad.

Precios (medias recientes; rangos en la sensibilidad):
- Combustible: índices G20 de Ship & Bunker, medias de 2025: HSFO 461, VLSFO 535 y MGO 754 USD/t. En el mar se usa
  HSFO (el Galicia lleva un depurador de circuito cerrado por motor principal); en atraque, MGO (supuesto).
- Cambio: 1,12 USD/EUR (BCE, referencia del 21-10-2026: 1,1234).
- EU ETS: 73,5 EUR/t, media de enero-noviembre de 2025 (rango 60-85; Qin, ClearBlue Markets, 13-11-2025).
- UK ETS: futuro ICE de diciembre de 2026, unos 70 USD/t de media en 2026 (Ship & Bunker, 2026).
- Electricidad de la conexión a red en atraque (OPS): 150-250 EUR/MWh (supuesto).
Regulación:
- EU ETS (Directiva (UE) 2023/959): 50 % de las emisiones de los viajes entre un puerto de la UE y uno de fuera, y
  100 % en atraque en puertos de la UE; 100 % de las emisiones verificadas desde las de 2026.
- UK ETS (desde el 1-07-2026; DNV, 05-06-2026): viajes entre puertos del Reino Unido y emisiones en sus puertos;
  los viajes internacionales no. Aquí: atraque en Portsmouth al 100 %, travesías al 0 %.
- FuelEU Maritime (Reglamento (UE) 2023/1805): 50 % de la energía de los viajes con un puerto de fuera de la UE y
  100 % en atraque en puertos de la UE (art. 2); límite de intensidad 91,16 gCO2e/MJ con reducciones del 2 % (2025),
  6 % (2030) y 14,5 % (2035) (art. 4); factores por defecto del anexo II (HFO: 0,0405 MJ/g, WtT 13,5 g/MJ, 3,114 g
  CO2/g, 0,00005 g CH4/g, 0,00018 g N2O/g; MGO: 0,0427, 14,4, 3,206) y GWP a 100 años de la Directiva (UE)
  2018/2001 (CH4 25, N2O 298); penalización de 2 400 EUR por t de VLSFO equivalente (41 000 MJ) del déficit
  (anexo IV). El ahorro de energía reduce el déficit; aquí se valora a esa tasa, que es una cota superior del coste de
  cumplimiento (se puede cumplir con biocombustibles o con un pool a menor coste).
  Factor de recompensa del viento f_wind (anexo I): 0,99, 0,97 y 0,95 para P_wind/P_prop = 0,05, 0,1 y >= 0,15, con
  P_wind = f_eff P_eff (MEPC.1/Circ.896) y P_prop = P_ME (75 % del MCR). Aquí se estima P_wind con la potencia
  equivalente al empuje media de la ruta, como orden de magnitud: el valor oficial necesita la matriz de vientos de
  la Circ.896.
  Atraque: desde 2030, los buques de pasaje en puertos del art. 9 del Reglamento (UE) 2023/1804 deben usar la OPS
  para toda su demanda eléctrica, salvo que la cubran entera con tecnologías de cero emisiones del anexo III, entre
  ellas la generación a bordo con viento (art. 6.5.b). La producción en puerto de las turbinas cubre menos del 1 %
  de esa demanda: no exime.
- EEXI y CII: la reducción del CII alcanzado es la del CO2 anual; la del EEXI, P_wind/(P_ME + P_AE), con P_AE =
  0,025 MCR + 250 kW, como orden de magnitud.
Salida: economia.json.
"""
import json
import os
import ruta_modelo as rm

HERE = os.path.dirname(os.path.abspath(__file__))
USD_EUR = 1.12
FUEL_USD = {"HSFO": 461.0, "VLSFO": 535.0, "MGO": 754.0}
EUA, EUA_RANGE = 73.5, (60.0, 85.0)
UKA = 70.0 / USD_EUR
ELEC = (150.0, 250.0)                       # EUR/MWh, supuesto
FUELS = {"HFO": dict(lcv=0.0405, wtt=13.5, cf=3.114, ch4=0.00005, n2o=0.00018),
         "MGO": dict(lcv=0.0427, wtt=14.4, cf=3.206, ch4=0.00005, n2o=0.00018)}
GWP_CH4, GWP_N2O = 25.0, 298.0
GHG_REF = 91.16
LIMIT_CUT = {"2025-2029": 0.02, "2030-2034": 0.06, "2035-2039": 0.145}
PENALTY, VLSFO_MJ = 2400.0, 41000.0
COV = dict(ets_sea=0.5, ets_berth_eu=1.0, uk_ets_berth=1.0, fueleu_sea=0.5, fueleu_berth_eu=1.0)
YEARS, RATE = 20, 0.08
P_ME = 0.75 * rm.P_MCR
P_AE_EEXI = 0.025 * rm.P_MCR + 250.0
SFOC_AE = rm.SFOC_AE
KEYS = ("R1", "T3", "T3x", "T5", "T4", "T2", "T1", "R1t", "R1-30", "quinta unidad")


def ghg_intensity(f):
    """Intensidad de GEI del pozo a la estela (gCO2e/MJ) con los factores por defecto del anexo II."""
    return f["wtt"] + (f["cf"] + GWP_CH4 * f["ch4"] + GWP_N2O * f["n2o"]) / f["lcv"]


def fueleu_value_per_t(fuel, period):
    """Valor (EUR) de no consumir 1 t de combustible dentro del ámbito de FuelEU, a la tasa de la penalización."""
    f = FUELS[fuel]
    i_act = ghg_intensity(f)
    target = GHG_REF * (1.0 - LIMIT_CUT[period])
    e_mj = f["lcv"] * 1e6
    return max(i_act - target, 0.0) * e_mj / (i_act * VLSFO_MJ) * PENALTY


def annuity():
    return (1.0 - (1.0 + RATE) ** -YEARS) / RATE


def port_split(res, n_units):
    """Producción en puerto por topología y puerto: MWh/año y t/año de combustible de los auxiliares evitado."""
    out = {}
    for key, field in (("T1", ("vawt",)), ("T4", ("vawt",)), ("T2", ("vawt", "sav")), ("T5", ("sav",)), ("T3", ("sav3",))):
        out[key] = {}
        for prt, h in rm.PORT_HOURS_WEEK.items():
            mwh = 0.0
            for s in res["seasons"]:
                weeks = rm.DAYS[s] / 7.0
                kw = sum(res["port_kw"][prt][s][f] for f in field)
                mwh += n_units * kw * h * weeks / 1e3
            out[key][prt] = dict(mwh=mwh, fuel_t=mwh * SFOC_AE / 1e3)
    return out


def values(res, res5, period="2025-2029", eua=EUA, uka=UKA, fuel_sea="HSFO", fuel_berth="MGO", ops=False,
           elec=ELEC[0]):
    """Valor anual (EUR) por topología y comparador: combustible, EU ETS, UK ETS y FuelEU, en el mar y en atraque."""
    cf_sea = rm.CF_FUEL                                   # t CO2/t (EU MRV del Galicia)
    v_sea_t = (FUEL_USD[fuel_sea] / USD_EUR + cf_sea * COV["ets_sea"] * eua
               + COV["fueleu_sea"] * fueleu_value_per_t("HFO", period))
    cf_b = FUELS["MGO"]["cf"] if fuel_berth == "MGO" else FUELS["HFO"]["cf"]
    fb = "MGO" if fuel_berth == "MGO" else "HFO"
    v_bilbao_t = FUEL_USD[fuel_berth] / USD_EUR + cf_b * COV["ets_berth_eu"] * eua \
        + COV["fueleu_berth_eu"] * fueleu_value_per_t(fb, period)
    v_ports_t = FUEL_USD[fuel_berth] / USD_EUR + cf_b * COV["uk_ets_berth"] * uka
    sea = dict(res["sea"])
    sea["quinta unidad"] = res5["sea"]["R1"] - res["sea"]["R1"]
    ps = port_split(res, P_UNITS)
    out = {}
    for k in KEYS:
        d_sea = sea[k]
        port_val, port_t, port_mwh = 0.0, 0.0, 0.0
        if k in ps:
            for prt, pv in ps[k].items():
                port_t += pv["fuel_t"]
                port_mwh += pv["mwh"]
                if ops:
                    port_val += pv["mwh"] * elec
                else:
                    port_val += pv["fuel_t"] * (v_bilbao_t if prt == "BILBAO" else v_ports_t)
        total = d_sea * v_sea_t + port_val
        out[k] = dict(mar_t=d_sea, puerto_t=port_t, puerto_mwh=port_mwh, eur_mar=d_sea * v_sea_t,
                      eur_puerto=port_val, eur_total=total, inversion_max=total * annuity(),
                      co2_t=d_sea * cf_sea + (0.0 if ops else port_t * cf_b))
    return out, dict(eur_t_mar=v_sea_t, eur_t_bilbao=v_bilbao_t, eur_t_portsmouth=v_ports_t)


P_UNITS = rm.N_UNITS


def fuel_eu_wind_ratio(groups, nu, ng, p, n_units):
    """P_wind/P_prop con la potencia equivalente al empuje media de la ruta (orden de magnitud)."""
    acc, hrs = 0.0, 0.0
    for s in rm.BASE_SEASONS:
        for name, direction, h in rm.crossings(p):
            v_kn, t_open = rm.crossing_speed(h, p)
            mb = rm.mean_benefits(groups, nu, ng, s, direction, v_kn, p)
            acc += (mb["_R1 empuje"] - mb["_R1 giro"]) * t_open
            hrs += t_open
    p_wind = n_units * acc / hrs
    return p_wind, p_wind / P_ME


def main():
    d, groups, nu, ng = rm.load_wind(os.path.join(HERE, "era5_ruta.json"))
    res = rm.annual(d, groups, nu, ng)
    res5 = rm.annual(d, groups, nu, ng, rm.P(n_units=5))
    out = dict(precios=dict(fuel_usd=FUEL_USD, usd_eur=USD_EUR, eua=EUA, eua_rango=EUA_RANGE, uka=UKA, elec=ELEC),
               intensidad_ghg={k: ghg_intensity(v) for k, v in FUELS.items()},
               fueleu_eur_t={per: {f: fueleu_value_per_t(f, per) for f in FUELS} for per in LIMIT_CUT},
               anualidad=annuity(), casos={})
    cases = {
        "2026, sin OPS": dict(period="2025-2029"),
        "2026, sin OPS, EUA 60": dict(period="2025-2029", eua=60.0),
        "2026, sin OPS, EUA 85": dict(period="2025-2029", eua=85.0),
        "2026, sin OPS, VLSFO en el mar": dict(period="2025-2029", fuel_sea="VLSFO"),
        "2030, con OPS en los dos puertos (150 EUR/MWh)": dict(period="2030-2034", ops=True, elec=ELEC[0]),
        "2030, con OPS en los dos puertos (250 EUR/MWh)": dict(period="2030-2034", ops=True, elec=ELEC[1]),
        "2035, con OPS": dict(period="2035-2039", ops=True, elec=ELEC[0]),
    }
    for lab, kw in cases.items():
        v, unit = values(res, res5, **kw)
        out["casos"][lab] = dict(valores=v, por_tonelada=unit)
    p_wind, ratio = fuel_eu_wind_ratio(groups, nu, ng, rm.P(), rm.N_UNITS)
    p_wind5, ratio5 = fuel_eu_wind_ratio(groups, nu, ng, rm.P(), 5)
    out["fueleu_viento"] = dict(p_wind_kw=p_wind, ratio=ratio, f_wind=1.0 if ratio < 0.05 else 0.99,
                                p_wind_kw_5=p_wind5, ratio_5=ratio5)
    out["eexi_reduccion_aprox"] = p_wind / (P_ME + P_AE_EEXI)
    co2_total = res["fuel_sea"] * rm.CF_FUEL + res["fuel_berth"] * FUELS["MGO"]["cf"]
    out["cii_reduccion"] = {k: (res["sea"][k] * rm.CF_FUEL) / co2_total for k in ("R1", "T3x", "R1t", "R1-30")}
    json.dump(out, open(os.path.join(HERE, "economia.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"Intensidad GEI (gCO2e/MJ): " + ", ".join(f"{k} {v:.2f}" for k, v in out["intensidad_ghg"].items()))
    print("FuelEU, valor de 1 t ahorrada en su ámbito a la tasa de la penalización (EUR): " +
          "; ".join(f"{per}: HFO {v['HFO']:.0f}, MGO {v['MGO']:.0f}" for per, v in out["fueleu_eur_t"].items()))
    print(f"Anualidad (20 años, 8 %): {annuity():.3f}")
    for lab, c in out["casos"].items():
        u = c["por_tonelada"]
        print(f"\n{lab}: EUR por t ahorrada en el mar {u['eur_t_mar']:.0f}; en atraque Bilbao {u['eur_t_bilbao']:.0f}, "
              f"Portsmouth {u['eur_t_portsmouth']:.0f}")
        for k, v in c["valores"].items():
            print(f"  {k:14s} mar {v['mar_t']:6.1f} t  puerto {v['puerto_t']:4.1f} t ({v['puerto_mwh']:5.1f} MWh)  "
                  f"{v['eur_total'] / 1e3:7.1f} kEUR/año  inversión máx. {v['inversion_max'] / 1e6:5.2f} MEUR "
                  f"(CO2 {v['co2_t']:6.0f} t)")
    fw = out["fueleu_viento"]
    print(f"\nFuelEU: P_wind ~ {fw['p_wind_kw']:.0f} kW (4 unidades), P_wind/P_prop ~ {fw['ratio']:.3f} -> f_wind "
          f"{fw['f_wind']}; con cinco unidades {fw['ratio_5']:.3f}")
    print(f"EEXI: reducción aproximada {100 * out['eexi_reduccion_aprox']:.1f} %")
    print("CII: reducción del CO2 anual " + ", ".join(f"{k} {100 * v:.2f} %" for k, v in out["cii_reduccion"].items()))
    return out


if __name__ == "__main__":
    main()
