#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Valor económico y regulatorio de los ahorros (versión 5, 27-10-2026; segunda revisión simulada R3-1, R3-3 y menores;
versión 4: 22-10-2026).

Para cada topología y comparador, el ahorro anual del modelo de ruta (ruta_modelo.annual, temporada de abril a
octubre) se convierte en euros y en la inversión máxima que justificaría (valor actual a 20 años y 8 %), en los
escenarios de coeficientes A (referencia), A+ y A-VG (ajustado al Viking Grace). El resultado principal es ese umbral;
como orden de magnitud del coste se dan los costes publicados de un rotor comercial (EMSA 2023, tabla 7: precio de
0,56-1,05 MEUR por unidad, instalación en un buque existente de 0,14-0,26 MEUR, formación 10 kEUR y mantenimiento de
12-22 kEUR/año, en EUR de 2021). No hay costes publicados de ninguna de las topologías de la patente.

Precios (rangos en la sensibilidad, +-30 % en combustible y derechos de emisión):
- Combustible: índices G20 de Ship & Bunker, medias de 2025: HSFO 461, VLSFO 535 y MGO 754 USD/t. En el mar se usa
  HSFO (el Galicia lleva un depurador de circuito cerrado por motor principal); en atraque, MGO (supuesto).
- Cambio: 1,12 USD/EUR (BCE, referencia del 21-10-2026: 1,1234).
- EU ETS: 73,43 EUR/t, media de las subastas de 2025 (ICAP, ficha del EU ETS; secundario 74,35).
- UK ETS: 48,05 GBP/t (63,27 USD/t), media de las subastas de 2025 (ICAP, ficha del UK ETS).
- Electricidad de la conexión a red en atraque (OPS): 150-250 EUR/MWh (supuesto).
Regulación:
- EU ETS (Directiva (UE) 2023/959): 50 % de las emisiones de los viajes entre un puerto de la UE y uno de fuera, y
  100 % en atraque en puertos de la UE; 100 % de las emisiones verificadas desde las de 2026.
- UK ETS (Orden SI 2026/392, en vigor desde el 1-07-2026, anexo 1, párrafo 7): viajes que empiezan y acaban en
  puertos del Reino Unido y «in-port activity», que incluye el atraque en un puerto del Reino Unido; buques de 5 000 GT
  o más. Aquí: atraque en Portsmouth al 100 %, travesías al 0 %.
- FuelEU Maritime (Reglamento (UE) 2023/1805): 50 % de la energía de los viajes con un puerto de fuera de la UE y
  100 % en atraque en puertos de la UE (art. 2); límite de intensidad 91,16 gCO2e/MJ con reducciones del 2 % (2025),
  6 % (2030) y 14,5 % (2035) (art. 4); factores por defecto del anexo II (HFO: 0,0405 MJ/g, WtT 13,5 g/MJ, 3,114 g
  CO2/g, 0,00005 g CH4/g, 0,00018 g N2O/g; MGO: 0,0427, 14,4, 3,206) y GWP a 100 años de la Directiva (UE)
  2018/2001 (CH4 25, N2O 298); penalización de 2 400 EUR por t de VLSFO equivalente (41 000 MJ) del déficit
  (anexo IV). El ahorro de energía reduce el déficit; aquí se valora a esa tasa, que es una cota superior del coste de
  cumplimiento, y al precio de mercado de los excedentes cedidos en un pool (OceanScore Pool Price Index, marzo de
  2026: unos 200 EUR/t CO2e; ofertas de 175-216).
  Factor de recompensa del viento f_wind (anexo I): con P_Wind de la MEPC.1/Circ.896 (fueleu_viento.py). Un escalón
  de 0,01 vale 0,01 E/41 000 x 2 400 EUR a la tasa de la penalización y 0,01 I E x precio en un pool (E, energía en el
  ámbito de FuelEU en MJ; I, intensidad real).
  Atraque: desde 2030, los buques de pasaje en puertos del art. 9 del Reglamento (UE) 2023/1804 deben usar la OPS
  para toda su demanda eléctrica, salvo que la cubran entera con tecnologías de cero emisiones del anexo III, entre
  ellas la generación a bordo con viento (art. 6.5.b). La producción en puerto de las turbinas cubre menos del 1 %
  de esa demanda: no exime.
- EEXI y CII: la reducción del EEXI con la Circ.896 es del orden de P_Wind/(P_ME + P_AE), con P_AE = 0,025 MCR +
  250 kW; la del CII alcanzado es la del CO2 del año, y el CII es anual: con solo la temporada modelada, la reducción
  relativa al CO2 de la temporada es una cota superior si el buque navega el resto del año.
Salida: economia.json.
"""
import json
import os
import ruta_modelo as rm

HERE = os.path.dirname(os.path.abspath(__file__))
USD_EUR = 1.12
FUEL_USD = {"HSFO": 461.0, "VLSFO": 535.0, "MGO": 754.0}
EUA, EUA_RANGE = 73.43, (73.43 * 0.7, 73.43 * 1.3)
UKA = 63.27 / USD_EUR
POOL_EUR_T = 200.0                          # EUR/t CO2e, OceanScore Pool Price Index (marzo de 2026)
CAPEX_ROTOR = (560e3 + 140e3 + 10e3, 1050e3 + 263e3 + 10e3)    # EUR por rotor (EMSA 2023, tabla 7; 2021)
OPEX_ROTOR = (12e3, 22e3)                   # EUR/año por rotor (mantenimiento y reparación)
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
KEYS = ("R1", "T3", "T3x", "T5", "T4", "T2", "T1", "R1t", "R1-30", "R1-D5", "R1t+X", "quinta unidad")


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
           elec=ELEC[0], fuel_mult=1.0):
    """Valor anual (EUR) por topología y comparador: combustible, EU ETS, UK ETS y FuelEU, en el mar y en atraque."""
    cf_sea = rm.CF_FUEL                                   # t CO2/t (EU MRV del Galicia)
    v_sea_t = (fuel_mult * FUEL_USD[fuel_sea] / USD_EUR + cf_sea * COV["ets_sea"] * eua
               + COV["fueleu_sea"] * fueleu_value_per_t("HFO", period))
    cf_b = FUELS["MGO"]["cf"] if fuel_berth == "MGO" else FUELS["HFO"]["cf"]
    fb = "MGO" if fuel_berth == "MGO" else "HFO"
    v_bilbao_t = fuel_mult * FUEL_USD[fuel_berth] / USD_EUR + cf_b * COV["ets_berth_eu"] * eua \
        + COV["fueleu_berth_eu"] * fueleu_value_per_t(fb, period)
    v_ports_t = fuel_mult * FUEL_USD[fuel_berth] / USD_EUR + cf_b * COV["uk_ets_berth"] * uka
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


def energy_scope(res):
    """Energía (MJ/año) en el ámbito de FuelEU del servicio modelado e intensidad media real (gCO2e/MJ): 50 % de la de
    las travesías (HFO) y 100 % de la del atraque en Bilbao (MGO)."""
    h_b = rm.PORT_HOURS_WEEK["BILBAO"] / sum(rm.PORT_HOURS_WEEK.values())
    e_sea = COV["fueleu_sea"] * res["fuel_sea"] * FUELS["HFO"]["lcv"] * 1e6
    e_b = COV["fueleu_berth_eu"] * res["fuel_berth"] * h_b * FUELS["MGO"]["lcv"] * 1e6
    i_mean = (e_sea * ghg_intensity(FUELS["HFO"]) + e_b * ghg_intensity(FUELS["MGO"])) / (e_sea + e_b)
    return e_sea + e_b, i_mean


def wind_step_value(e_mj, i_mean, df=0.01):
    """Valor anual (EUR) de mejorar f_wind en df: a la tasa de la penalización y al precio del pool."""
    return dict(penalizacion=df * e_mj / VLSFO_MJ * PENALTY, pool=df * i_mean * e_mj / 1e6 * POOL_EUR_T,
                t_co2e=df * i_mean * e_mj / 1e6)


def capex_pv(n_units):
    """Coste actual (EUR) de n rotores comerciales con su mantenimiento a 20 años y 8 % (EMSA 2023)."""
    return tuple(n_units * (c + o * annuity()) for c, o in zip(CAPEX_ROTOR, OPEX_ROTOR))


ESC_VALOR = ("A", "A+", "A-VG", "A-VG 207", "A-VG 282")


def main():
    import coeficientes as cf
    fw = json.load(open(os.path.join(HERE, "fueleu_viento.json"), encoding="utf-8"))["resultados"]
    d, groups, nu, ng = rm.load_wind()
    out = dict(precios=dict(fuel_usd=FUEL_USD, usd_eur=USD_EUR, eua=EUA, eua_rango=EUA_RANGE, uka=UKA, elec=ELEC,
                            pool_eur_t=POOL_EUR_T),
               intensidad_ghg={k: ghg_intensity(v) for k, v in FUELS.items()},
               fueleu_eur_t={per: {f: fueleu_value_per_t(f, per) for f in FUELS} for per in LIMIT_CUT},
               anualidad=annuity(), capex_rotor=CAPEX_ROTOR, opex_rotor=OPEX_ROTOR,
               coste_actual_4_rotores=capex_pv(4), coste_actual_5_rotores=capex_pv(5), escenarios={})
    for esc in ESC_VALOR:
        p = rm.P(**cf.ESCENARIOS_COEF[esc])
        res = rm.annual(d, groups, nu, ng, p)
        res5 = rm.annual(d, groups, nu, ng, rm.P(n_units=5, **cf.ESCENARIOS_COEF[esc]))
        cases = {"2026, sin OPS": dict(period="2025-2029")}
        if esc == "A":
            cases.update({
                "2026, sin OPS, EUA -30 %": dict(period="2025-2029", eua=EUA_RANGE[0]),
                "2026, sin OPS, EUA +30 %": dict(period="2025-2029", eua=EUA_RANGE[1]),
                "2026, sin OPS, combustible -30 %": dict(period="2025-2029", fuel_mult=0.7),
                "2026, sin OPS, combustible +30 %": dict(period="2025-2029", fuel_mult=1.3),
                "2026, sin OPS, VLSFO en el mar": dict(period="2025-2029", fuel_sea="VLSFO"),
                "2030, con OPS en los dos puertos (150 EUR/MWh)": dict(period="2030-2034", ops=True, elec=ELEC[0]),
                "2030, con OPS en los dos puertos (250 EUR/MWh)": dict(period="2030-2034", ops=True, elec=ELEC[1]),
                "2035, con OPS": dict(period="2035-2039", ops=True, elec=ELEC[0]),
            })
        e_mj, i_mean = energy_scope(res)
        step = wind_step_value(e_mj, i_mean)
        fwe = fw[esc]
        f_r0 = 1.0
        wind = {}
        for k, key in (("R1", "R1"), ("quinta unidad", "R1x5"), ("T3x", "T3x"), ("T3", "T3x"), ("R1-30", "R1-30"),
                       ("R1t", "R1t")):
            f = fwe[key]["f_wind"]
            base = f_r0 if k == "R1" else fwe["R1"]["f_wind"]
            wind[k] = dict(f_wind=f, ratio=fwe[key]["ratio"], df=base - f,
                           eur_penalizacion=(base - f) / 0.01 * step["penalizacion"],
                           eur_pool=(base - f) / 0.01 * step["pool"])
        e = dict(energia_fueleu_MJ=e_mj, intensidad_media=i_mean, escalon_001=step, f_wind=wind, casos={},
                 R0=res["fuel_total"], R0_mar=res["fuel_sea"], R0_atraque=res["fuel_berth"])
        for lab, kw in cases.items():
            v, unit = values(res, res5, **kw)
            e["casos"][lab] = dict(valores=v, por_tonelada=unit)
        co2_total = res["fuel_sea"] * rm.CF_FUEL + res["fuel_berth"] * FUELS["MGO"]["cf"]
        e["cii_reduccion_temporada"] = {k: (res["sea"][k] * rm.CF_FUEL) / co2_total
                                        for k in ("R1", "T3x", "R1t", "R1-30")}
        e["eexi_reduccion_aprox"] = fwe["R1"]["p_wind_kw"] / (P_ME + P_AE_EEXI)
        out["escenarios"][esc] = e
    # CII con el año completo (dos rotaciones todo el año), escenario A
    ra = rm.annual(d, groups, nu, ng, rm.P(), seasons=rm.ALL_SEASONS)
    co2a = ra["fuel_sea"] * rm.CF_FUEL + ra["fuel_berth"] * FUELS["MGO"]["cf"]
    out["cii_reduccion_anio_completo_A"] = (ra["sea"]["R1"] * rm.CF_FUEL) / co2a
    json.dump(out, open(os.path.join(HERE, "economia.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"Intensidad GEI (gCO2e/MJ): " + ", ".join(f"{k} {v:.2f}" for k, v in out["intensidad_ghg"].items()))
    print(f"Anualidad (20 años, 8 %): {annuity():.3f}; coste actual de 4 rotores comerciales "
          f"{out['coste_actual_4_rotores'][0] / 1e6:.2f}-{out['coste_actual_4_rotores'][1] / 1e6:.2f} MEUR")
    for esc, e in out["escenarios"].items():
        st = e["escalon_001"]
        print(f"\n== {esc}: energía FuelEU {e['energia_fueleu_MJ'] / 1e6:.0f} TJ; escalón de 0,01 = {st['t_co2e']:.0f} t CO2e "
              f"= {st['penalizacion'] / 1e3:.0f} kEUR (penalización) o {st['pool'] / 1e3:.0f} kEUR (pool)")
        print("   f_wind: " + "; ".join(f"{k} {v['f_wind']:.2f} ({v['ratio']:.3f}, df {v['df']:.2f}, "
                                        f"{v['eur_penalizacion'] / 1e3:.0f}/{v['eur_pool'] / 1e3:.0f} kEUR)"
                                        for k, v in e["f_wind"].items()))
        for lab, c in e["casos"].items():
            u = c["por_tonelada"]
            print(f"  {lab}: EUR/t mar {u['eur_t_mar']:.0f}; Bilbao {u['eur_t_bilbao']:.0f}; Portsmouth {u['eur_t_portsmouth']:.0f}")
            for k, v in c["valores"].items():
                print(f"    {k:14s} mar {v['mar_t']:6.1f} t  puerto {v['puerto_t']:4.1f} t  {v['eur_total'] / 1e3:7.1f} kEUR/año"
                      f"  inversión máx. {v['inversion_max'] / 1e6:5.2f} MEUR (CO2 {v['co2_t']:5.0f} t)")
        print(f"  CII (temporada): " + ", ".join(f"{k} {100 * v:.2f} %" for k, v in e["cii_reduccion_temporada"].items())
              + f"; EEXI ~{100 * e['eexi_reduccion_aprox']:.1f} %")
    print(f"CII año completo, R1 (A): {100 * out['cii_reduccion_anio_completo_A']:.2f} %")
    return out


if __name__ == "__main__":
    main()
