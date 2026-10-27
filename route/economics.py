#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Economic and regulatory value of the savings (version 5, 27-10-2026; second simulated review R3-1, R3-3 and minor
comments; version 4: 22-10-2026).

For each topology and comparator, the annual saving of the route model (route_model.annual, April-October season) is
converted into euros and into the largest investment it would justify (present value over 20 years at 8 %), in the
coefficient scenarios A (reference), A+ and A-VG (fitted to the Viking Grace). The main result is that threshold; as
an order of magnitude of the cost, the published costs of a commercial rotor are given (EMSA 2023, table 7: price of
0.56-1.05 MEUR per unit, installation on an existing ship 0.14-0.26 MEUR, training 10 kEUR and maintenance of
12-22 kEUR/yr, in EUR of 2021). There are no published costs of any of the topologies of the patent.

Prices (ranges in the sensitivity, +-30 % on fuel and emission allowances):
- Fuel: G20 indices of Ship & Bunker, means of 2025: HSFO 461, VLSFO 535 and MGO 754 USD/t. At sea HSFO is used (the
  Galicia has a closed-loop scrubber per main engine); at berth, MGO (assumed).
- Exchange rate: 1.12 USD/EUR (ECB, reference rate of 21-10-2026: 1.1234).
- EU ETS: 73.43 EUR/t, mean of the 2025 auctions (ICAP, EU ETS factsheet; secondary market 74.35).
- UK ETS: 48.05 GBP/t (63.27 USD/t), mean of the 2025 auctions (ICAP, UK ETS factsheet).
- Electricity of the onshore power supply at berth (OPS): 150-250 EUR/MWh (assumed).
Regulation:
- EU ETS (Directive (EU) 2023/959): 50 % of the emissions of voyages between an EU port and a non-EU port, and 100 %
  at berth in EU ports; 100 % of the verified emissions from those of 2026.
- UK ETS (Order SI 2026/392, in force from 1-07-2026, schedule 1, paragraph 7): voyages that start and end in UK ports
  and «in-port activity», which includes the berth in a UK port; ships of 5,000 GT or more. Here: berth in Portsmouth
  at 100 %, crossings at 0 %.
- FuelEU Maritime (Regulation (EU) 2023/1805): 50 % of the energy of voyages with a non-EU port and 100 % at berth in
  EU ports (art. 2); intensity limit 91.16 gCO2e/MJ with reductions of 2 % (2025), 6 % (2030) and 14.5 % (2035)
  (art. 4); default factors of annex II (HFO: 0.0405 MJ/g, WtT 13.5 g/MJ, 3.114 g CO2/g, 0.00005 g CH4/g, 0.00018 g
  N2O/g; MGO: 0.0427, 14.4, 3.206) and 100-year GWP of Directive (EU) 2018/2001 (CH4 25, N2O 298); penalty of
  2,400 EUR per t of VLSFO equivalent (41,000 MJ) of the deficit (annex IV). The energy saving reduces the deficit;
  here it is valued at that rate, which is an upper bound of the compliance cost, and at the market price of surpluses
  transferred in a pool (OceanScore Pool Price Index, March 2026: about 200 EUR/t CO2e; offers of 175-216).
  Wind reward factor f_wind (annex I): with P_Wind from MEPC.1/Circ.896 (fueleu_wind.py). A step of 0.01 is worth
  0.01 E/41,000 x 2,400 EUR at the penalty rate and 0.01 I E x price in a pool (E, energy within the scope of FuelEU
  in MJ; I, actual intensity).
  Berth: from 2030, passenger ships in ports of art. 9 of Regulation (EU) 2023/1804 must use OPS for all their
  electrical demand, unless they cover it entirely with the zero-emission technologies of annex III, among them
  on-board wind generation (art. 6.5.b). The production of the turbines in port covers less than 1 % of that demand:
  it does not exempt.
- EEXI and CII: the reduction of the EEXI with Circ.896 is of the order of P_Wind/(P_ME + P_AE), with P_AE =
  0.025 MCR + 250 kW; that of the attained CII is that of the CO2 of the year, and the CII is annual: with only the
  modelled season, the reduction relative to the CO2 of the season is an upper bound if the ship sails the rest of
  the year.
Output: economics.json.
"""
import json
import os
import route_model as rm

HERE = os.path.dirname(os.path.abspath(__file__))
USD_EUR = 1.12
FUEL_USD = {"HSFO": 461.0, "VLSFO": 535.0, "MGO": 754.0}
EUA, EUA_RANGE = 73.43, (73.43 * 0.7, 73.43 * 1.3)
UKA = 63.27 / USD_EUR
POOL_EUR_T = 200.0                          # EUR/t CO2e, OceanScore Pool Price Index (March 2026)
CAPEX_ROTOR = (560e3 + 140e3 + 10e3, 1050e3 + 263e3 + 10e3)    # EUR per rotor (EMSA 2023, table 7; 2021)
OPEX_ROTOR = (12e3, 22e3)                   # EUR/yr per rotor (maintenance and repair)
ELEC = (150.0, 250.0)                       # EUR/MWh, assumed
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
KEYS = ("R1", "T3", "T3x", "T5", "T4", "T2", "T1", "R1t", "R1-30", "R1-D5", "R1t+X", "fifth unit")


def ghg_intensity(f):
    """Well-to-wake GHG intensity (gCO2e/MJ) with the default factors of annex II."""
    return f["wtt"] + (f["cf"] + GWP_CH4 * f["ch4"] + GWP_N2O * f["n2o"]) / f["lcv"]


def fueleu_value_per_t(fuel, period):
    """Value (EUR) of not consuming 1 t of fuel within the scope of FuelEU, at the penalty rate."""
    f = FUELS[fuel]
    i_act = ghg_intensity(f)
    target = GHG_REF * (1.0 - LIMIT_CUT[period])
    e_mj = f["lcv"] * 1e6
    return max(i_act - target, 0.0) * e_mj / (i_act * VLSFO_MJ) * PENALTY


def annuity():
    return (1.0 - (1.0 + RATE) ** -YEARS) / RATE


def port_split(res, n_units):
    """Production in port by topology and port: MWh/yr and t/yr of auxiliary fuel avoided."""
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
    """Annual value (EUR) by topology and comparator: fuel, EU ETS, UK ETS and FuelEU, at sea and at berth."""
    cf_sea = rm.CF_FUEL                                   # t CO2/t (EU MRV of the Galicia)
    v_sea_t = (fuel_mult * FUEL_USD[fuel_sea] / USD_EUR + cf_sea * COV["ets_sea"] * eua
               + COV["fueleu_sea"] * fueleu_value_per_t("HFO", period))
    cf_b = FUELS["MGO"]["cf"] if fuel_berth == "MGO" else FUELS["HFO"]["cf"]
    fb = "MGO" if fuel_berth == "MGO" else "HFO"
    v_bilbao_t = fuel_mult * FUEL_USD[fuel_berth] / USD_EUR + cf_b * COV["ets_berth_eu"] * eua \
        + COV["fueleu_berth_eu"] * fueleu_value_per_t(fb, period)
    v_ports_t = fuel_mult * FUEL_USD[fuel_berth] / USD_EUR + cf_b * COV["uk_ets_berth"] * uka
    sea = dict(res["sea"])
    sea["fifth unit"] = res5["sea"]["R1"] - res["sea"]["R1"]
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
        out[k] = dict(sea_t=d_sea, port_t=port_t, port_mwh=port_mwh, eur_sea=d_sea * v_sea_t,
                      eur_port=port_val, eur_total=total, max_investment=total * annuity(),
                      co2_t=d_sea * cf_sea + (0.0 if ops else port_t * cf_b))
    return out, dict(eur_t_sea=v_sea_t, eur_t_bilbao=v_bilbao_t, eur_t_portsmouth=v_ports_t)


P_UNITS = rm.N_UNITS


def energy_scope(res):
    """Energy (MJ/yr) within the scope of FuelEU of the modelled service and actual mean intensity (gCO2e/MJ): 50 % of
    that of the crossings (HFO) and 100 % of that of the berth in Bilbao (MGO)."""
    h_b = rm.PORT_HOURS_WEEK["BILBAO"] / sum(rm.PORT_HOURS_WEEK.values())
    e_sea = COV["fueleu_sea"] * res["fuel_sea"] * FUELS["HFO"]["lcv"] * 1e6
    e_b = COV["fueleu_berth_eu"] * res["fuel_berth"] * h_b * FUELS["MGO"]["lcv"] * 1e6
    i_mean = (e_sea * ghg_intensity(FUELS["HFO"]) + e_b * ghg_intensity(FUELS["MGO"])) / (e_sea + e_b)
    return e_sea + e_b, i_mean


def wind_step_value(e_mj, i_mean, df=0.01):
    """Annual value (EUR) of improving f_wind by df: at the penalty rate and at the pool price."""
    return dict(penalty=df * e_mj / VLSFO_MJ * PENALTY, pool=df * i_mean * e_mj / 1e6 * POOL_EUR_T,
                t_co2e=df * i_mean * e_mj / 1e6)


def capex_pv(n_units):
    """Present cost (EUR) of n commercial rotors with their maintenance over 20 years at 8 % (EMSA 2023)."""
    return tuple(n_units * (c + o * annuity()) for c, o in zip(CAPEX_ROTOR, OPEX_ROTOR))


VALUE_SCEN = ("A", "A+", "A-VG", "A-VG 207", "A-VG 282")


def main():
    import coefficients as cf
    fw = json.load(open(os.path.join(HERE, "fueleu_wind.json"), encoding="utf-8"))["results"]
    d, groups, nu, ng = rm.load_wind()
    out = dict(prices=dict(fuel_usd=FUEL_USD, usd_eur=USD_EUR, eua=EUA, eua_range=EUA_RANGE, uka=UKA, elec=ELEC,
                           pool_eur_t=POOL_EUR_T),
               ghg_intensity={k: ghg_intensity(v) for k, v in FUELS.items()},
               fueleu_eur_t={per: {f: fueleu_value_per_t(f, per) for f in FUELS} for per in LIMIT_CUT},
               annuity=annuity(), capex_rotor=CAPEX_ROTOR, opex_rotor=OPEX_ROTOR,
               present_cost_4_rotors=capex_pv(4), present_cost_5_rotors=capex_pv(5), scenarios={})
    for esc in VALUE_SCEN:
        p = rm.P(**cf.COEF_SCENARIOS[esc])
        res = rm.annual(d, groups, nu, ng, p)
        res5 = rm.annual(d, groups, nu, ng, rm.P(n_units=5, **cf.COEF_SCENARIOS[esc]))
        cases = {"2026, no OPS": dict(period="2025-2029")}
        if esc == "A":
            cases.update({
                "2026, no OPS, EUA -30 %": dict(period="2025-2029", eua=EUA_RANGE[0]),
                "2026, no OPS, EUA +30 %": dict(period="2025-2029", eua=EUA_RANGE[1]),
                "2026, no OPS, fuel -30 %": dict(period="2025-2029", fuel_mult=0.7),
                "2026, no OPS, fuel +30 %": dict(period="2025-2029", fuel_mult=1.3),
                "2026, no OPS, VLSFO at sea": dict(period="2025-2029", fuel_sea="VLSFO"),
                "2030, with OPS in both ports (150 EUR/MWh)": dict(period="2030-2034", ops=True, elec=ELEC[0]),
                "2030, with OPS in both ports (250 EUR/MWh)": dict(period="2030-2034", ops=True, elec=ELEC[1]),
                "2035, with OPS": dict(period="2035-2039", ops=True, elec=ELEC[0]),
            })
        e_mj, i_mean = energy_scope(res)
        step = wind_step_value(e_mj, i_mean)
        fwe = fw[esc]
        f_r0 = 1.0
        wind = {}
        for k, key in (("R1", "R1"), ("fifth unit", "R1x5"), ("T3x", "T3x"), ("T3", "T3x"), ("R1-30", "R1-30"),
                       ("R1t", "R1t")):
            f = fwe[key]["f_wind"]
            base = f_r0 if k == "R1" else fwe["R1"]["f_wind"]
            wind[k] = dict(f_wind=f, ratio=fwe[key]["ratio"], df=base - f,
                           eur_penalty=(base - f) / 0.01 * step["penalty"],
                           eur_pool=(base - f) / 0.01 * step["pool"])
        e = dict(energy_fueleu_MJ=e_mj, mean_intensity=i_mean, step_001=step, f_wind=wind, cases={},
                 R0=res["fuel_total"], R0_sea=res["fuel_sea"], R0_berth=res["fuel_berth"])
        for lab, kw in cases.items():
            v, unit = values(res, res5, **kw)
            e["cases"][lab] = dict(values=v, per_tonne=unit)
        co2_total = res["fuel_sea"] * rm.CF_FUEL + res["fuel_berth"] * FUELS["MGO"]["cf"]
        e["cii_reduction_season"] = {k: (res["sea"][k] * rm.CF_FUEL) / co2_total
                                     for k in ("R1", "T3x", "R1t", "R1-30")}
        e["eexi_reduction_approx"] = fwe["R1"]["p_wind_kw"] / (P_ME + P_AE_EEXI)
        out["scenarios"][esc] = e
    # CII with the whole year (two rotations all year), scenario A
    ra = rm.annual(d, groups, nu, ng, rm.P(), seasons=rm.ALL_SEASONS)
    co2a = ra["fuel_sea"] * rm.CF_FUEL + ra["fuel_berth"] * FUELS["MGO"]["cf"]
    out["cii_reduction_whole_year_A"] = (ra["sea"]["R1"] * rm.CF_FUEL) / co2a
    json.dump(out, open(os.path.join(HERE, "economics.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"GHG intensity (gCO2e/MJ): " + ", ".join(f"{k} {v:.2f}" for k, v in out["ghg_intensity"].items()))
    print(f"Annuity (20 years, 8 %): {annuity():.3f}; present cost of 4 commercial rotors "
          f"{out['present_cost_4_rotors'][0] / 1e6:.2f}-{out['present_cost_4_rotors'][1] / 1e6:.2f} MEUR")
    for esc, e in out["scenarios"].items():
        st = e["step_001"]
        print(f"\n== {esc}: FuelEU energy {e['energy_fueleu_MJ'] / 1e6:.0f} TJ; step of 0.01 = {st['t_co2e']:.0f} t CO2e "
              f"= {st['penalty'] / 1e3:.0f} kEUR (penalty) or {st['pool'] / 1e3:.0f} kEUR (pool)")
        print("   f_wind: " + "; ".join(f"{k} {v['f_wind']:.2f} ({v['ratio']:.3f}, df {v['df']:.2f}, "
                                        f"{v['eur_penalty'] / 1e3:.0f}/{v['eur_pool'] / 1e3:.0f} kEUR)"
                                        for k, v in e["f_wind"].items()))
        for lab, c in e["cases"].items():
            u = c["per_tonne"]
            print(f"  {lab}: EUR/t sea {u['eur_t_sea']:.0f}; Bilbao {u['eur_t_bilbao']:.0f}; Portsmouth {u['eur_t_portsmouth']:.0f}")
            for k, v in c["values"].items():
                print(f"    {k:14s} sea {v['sea_t']:6.1f} t  port {v['port_t']:4.1f} t  {v['eur_total'] / 1e3:7.1f} kEUR/yr"
                      f"  max. investment {v['max_investment'] / 1e6:5.2f} MEUR (CO2 {v['co2_t']:5.0f} t)")
        print(f"  CII (season): " + ", ".join(f"{k} {100 * v:.2f} %" for k, v in e["cii_reduction_season"].items())
              + f"; EEXI ~{100 * e['eexi_reduction_approx']:.1f} %")
    print(f"CII whole year, R1 (A): {100 * out['cii_reduction_whole_year_A']:.2f} %")
    return out


if __name__ == "__main__":
    main()
