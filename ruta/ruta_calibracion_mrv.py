#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contraste del modelo de buque (ruta_modelo.py, caso R0) con los datos públicos EU MRV del Galicia
(IMO 9856189) descargados de EMSA THETIS-MRV el 18-10-2026 (mrv_galicia.json). Versión 2 (19-10-2026): usa el
horario corregido de ruta_modelo (P().tz_fix) y muestra también la versión 1.

Alcance de EU MRV para un ferry Reino Unido-España después del Brexit: incluye todas las travesías
(cada una sale de o llega a un puerto de la UE) y el atraque en puertos de la UE, pero no el atraque
en puertos del Reino Unido. Por eso la comparación se hace por partes:
  1. navegación: combustible fuera del atraque, horas en el mar, distancia y velocidad media;
  2. atraque: el combustible en atraque en la UE dividido por todas las horas fuera del mar da una
     COTA INFERIOR del consumo medio en atraque (las horas reales de atraque en la UE son menos).
Se usan los años 2021-2023, con servicio Reino Unido-UE casi exclusivo (emisiones entre puertos de
la UE < 2 % del total). 2024 (sobre todo travesías intra-UE) y 2025 (mixto) se muestran aparte.
"""
import json
import os
import ruta_modelo as rm

HERE = os.path.dirname(os.path.abspath(__file__))
EF_MGO = 3.206    # t CO2/t, gasóleo marino (Anexo I del Reglamento (UE) 2015/757)
YEARS_UK_EU = ("2021", "2022", "2023")


def load_mrv():
    with open(os.path.join(HERE, "mrv_galicia.json")) as f:
        return json.load(f)


def per_year(rows):
    out = {}
    for y, r in rows.items():
        ef = r["co2"] / r["fuel"]
        f_berth = r["berth"] / ef
        h_year = 8784.0 if y in ("2020", "2024") else 8760.0
        h_off = h_year - r["sea_h"]
        dist = r["fuel"] / r["fpd"] * 1e3            # mn (consumo por distancia = consumo total / distancia)
        f_sea = r["fuel"] - f_berth
        out[y] = dict(ef=ef, f_berth=f_berth, h_off=h_off, dist=dist, v_mean=dist / r["sea_h"], f_sea=f_sea,
                      rate_sea=f_sea / r["sea_h"], fpd_sea=f_sea / dist * 1e3, intra_eu=r["between"] / r["co2"],
                      berth_lb_kgph=f_berth * 1e3 / h_off, berth_lb_kw=f_berth * 1e6 / h_off / rm.SFOC_AE)
    return out


def model_r0(p):
    d, groups, nu, ng = rm.load_wind()
    res = rm.annual(d, groups, nu, ng, p)
    weeks = sum(rm.DAYS.values()) / 7.0
    n_cross = len(rm.crossings(p)) * weeks
    dist = n_cross * (540.0 + rm.NM_APPROACH)
    f_hotel_sea = rm.P_HOTEL_SEA * res["hours_sea"] * p.sfoc_me / rm.ETA_PTO / 1e6
    f_sea = res["fuel_me"] + f_hotel_sea
    h_bilbao = rm.PORT_HOURS_WEEK["BILBAO"] * weeks
    f_berth_eu = rm.P_HOTEL_BERTH * h_bilbao * rm.SFOC_AE / 1e6
    return dict(fuel_total=res["fuel_total"], fuel_me=res["fuel_me"], f_sea=f_sea, hours_sea=res["hours_sea"],
                hours_port=res["hours_port"], dist=dist, v_mean=dist / res["hours_sea"],
                rate_sea=f_sea / res["hours_sea"], fpd_sea=f_sea / dist * 1e3,
                f_mrv_scope=f_sea + f_berth_eu, berth_kgph=rm.P_HOTEL_BERTH * rm.SFOC_AE / 1e3)


def main():
    mrv = load_mrv()
    rows = mrv["rows"]
    py = per_year(rows)
    print("EU MRV, GALICIA (IMO 9856189)")
    print(f"{'año':5s} {'comb t':>8s} {'CO2/comb':>8s} {'intraUE':>7s} {'h mar':>7s} {'mn':>8s} {'kn':>6s} "
          f"{'t/h mar':>8s} {'kg/mn':>7s} {'atr UE t':>8s} {'cota kg/h':>9s}")
    for y in sorted(rows):
        r, q = rows[y], py[y]
        print(f"{y:5s} {r['fuel']:8.0f} {q['ef']:8.3f} {100*q['intra_eu']:6.1f}% {r['sea_h']:7.0f} {q['dist']:8.0f} "
              f"{q['v_mean']:6.2f} {q['rate_sea']:8.3f} {q['fpd_sea']:7.1f} {q['f_berth']:8.0f} {q['berth_lb_kgph']:9.0f}")
    ys = YEARS_UK_EU
    pooled = dict(f_sea=sum(py[y]["f_sea"] for y in ys), sea_h=sum(rows[y]["sea_h"] for y in ys),
                  dist=sum(py[y]["dist"] for y in ys))
    pooled["rate_sea"] = pooled["f_sea"] / pooled["sea_h"]
    pooled["fpd_sea"] = pooled["f_sea"] / pooled["dist"] * 1e3
    pooled["v_mean"] = pooled["dist"] / pooled["sea_h"]
    ef_all = sum(rows[y]["co2"] for y in rows if y != "2020") / sum(rows[y]["fuel"] for y in rows if y != "2020")
    print(f"\n2021-2023 agregados: {pooled['rate_sea']:.3f} t/h en el mar, {pooled['fpd_sea']:.1f} kg/mn, "
          f"{pooled['v_mean']:.2f} kn, {pooled['sea_h']/3:.0f} h/año en el mar")
    print(f"Factor CO2/combustible del buque 2021-2025: {ef_all:.4f} t/t")
    out = {}
    for label, p in (("versión 1 (horario sin corregir)", rm.P(v1=True)), ("versión 2", rm.P())):
        m = model_r0(p)
        out[label] = m
        print(f"\nMODELO R0, {label}")
        print(f"  combustible total {m['fuel_total']:.0f} t/año (propulsión {m['fuel_me']:.0f}); en el mar {m['f_sea']:.0f} t/año")
        print(f"  {m['hours_sea']:.0f} h en el mar, {m['hours_port']:.0f} h atracado, {m['dist']:.0f} mn/año, {m['v_mean']:.2f} kn")
        print(f"  {m['rate_sea']:.3f} t/h en el mar ({100*(m['rate_sea']/pooled['rate_sea']-1):+.1f} % frente a MRV 2021-2023); "
              f"{m['fpd_sea']:.1f} kg/mn ({100*(m['fpd_sea']/pooled['fpd_sea']-1):+.1f} %)")
        print(f"  alcance MRV equivalente (mar + atraque en Bilbao): {m['f_mrv_scope']:.0f} t/año; "
              f"MRV 2021 {rows['2021']['fuel']:.0f}, 2022 {rows['2022']['fuel']:.0f}")
    lb = [py[y]["berth_lb_kgph"] for y in ("2021", "2022")]
    print(f"\natraque: modelo {rm.P_HOTEL_BERTH * rm.SFOC_AE / 1e3:.0f} kg/h; cota inferior MRV 2021-2022 {min(lb):.0f}-{max(lb):.0f} kg/h")
    return dict(per_year=py, pooled=pooled, ef_all=ef_all, model=out)


if __name__ == "__main__":
    main()
