#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contraste del modelo del buque sin rotores (ruta_modelo.py, caso R0) con los datos públicos EU MRV del Galicia
(IMO 9856189) descargados de EMSA THETIS-MRV el 18-10-2026 (mrv_galicia.json), y calibración del margen de mar.
Versión 4 (22-10-2026).

Alcance de EU MRV para un ferry Reino Unido-España después del Brexit: incluye todas las travesías (cada una sale de o
llega a un puerto de la UE) y el atraque en puertos de la UE, pero no el atraque en puertos del Reino Unido. Se usan
los años 2021-2023, con servicio Reino Unido-UE casi exclusivo (emisiones entre puertos de la UE < 2 % del total):
  1. navegación: combustible fuera del atraque, horas en el mar, distancia y velocidad media;
  2. atraque: el combustible en atraque en la UE dividido por todas las horas fuera del mar da una COTA INFERIOR del
     consumo medio en atraque.
Contraste a igual velocidad (versión 4; hasta la versión 3 se comparaba el consumo horario del modelo a su propia
velocidad media, 17,4 kn, con el de EU MRV a 16,6 kn): las velocidades de las cuatro travesías del horario de 2026 se
escalan con un mismo factor hasta que la velocidad media del modelo (distancia entre horas en el mar, con las
aproximaciones a puerto) iguala la de EU MRV 2021-2023, con el viento de todo el año. El margen de mar mu, que
multiplica la potencia en aguas tranquilas por (1 + mu), se calibra para que el consumo horario coincida. Se calibra
para las curvas potencia-velocidad de exponente 3 (caso base) y 4 (sensibilidad); el resultado se copia en
ruta_modelo.SEA_MARGIN, y el script avisa si no coincide.
Salida: calibracion_mrv.json.
"""
import json
import os
import ruta_modelo as rm

HERE = os.path.dirname(os.path.abspath(__file__))
YEARS_UK_EU = ("2021", "2022", "2023")
SEASONS_YEAR = (0, 1)


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


def pooled_mrv(rows, py):
    ys = YEARS_UK_EU
    pooled = dict(f_sea=sum(py[y]["f_sea"] for y in ys), sea_h=sum(rows[y]["sea_h"] for y in ys),
                  dist=sum(py[y]["dist"] for y in ys))
    pooled["rate_sea"] = pooled["f_sea"] / pooled["sea_h"]
    pooled["fpd_sea"] = pooled["f_sea"] / pooled["dist"] * 1e3
    pooled["v_mean"] = pooled["dist"] / pooled["sea_h"]
    return pooled


_PW = {}


def model_at(p, factor, wind):
    """Consumo de R0 en el mar con las velocidades de mar abierto del horario multiplicadas por 'factor', con el
    viento de todo el año: t/h, kg/mn y velocidad media (distancia entre horas en el mar, con las aproximaciones)."""
    d, groups, nu, ng = wind
    fuel = hours = dist = 0.0
    for s in SEASONS_YEAR:
        weeks = rm.DAYS[s] / 7.0
        for name, direction, h in rm.crossings(p):
            v0, _ = rm.crossing_speed(h, p)
            v = v0 * factor
            t_open = p.nm_open / v
            key = (s, direction, round(v, 9), p.r0_wind, p.axc, p.ayc, p.e_span, p.eta_d)
            if key not in _PW:
                _PW[key] = rm.r0_wind_power(groups, nu, ng, s, direction, v, p)
            f, _, _ = rm.crossing_fuel(v, t_open, _PW[key], p)
            fuel += f * weeks
            hours += (t_open + rm.NM_APPROACH / rm.V_APPROACH) * weeks
            dist += (p.nm_open + rm.NM_APPROACH) * weeks
    return dict(rate=fuel / hours, fpd=fuel / dist * 1e3, v_mean=dist / hours)


def speed_factor(p, v_target, wind):
    lo, hi = 0.5, 1.5
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if model_at(rm.P(**{**vars(p), "r0_wind": False}), mid, wind)["v_mean"] < v_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def calibrate(n_prop, v_target, rate_target, wind):
    """Margen de mar mu con el que el consumo horario del modelo a la velocidad de EU MRV iguala al medido."""
    p0 = rm.P(n_prop=n_prop, sea_margin=0.0)
    f = speed_factor(p0, v_target, wind)
    lo, hi = -0.3, 0.6
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if model_at(rm.P(n_prop=n_prop, sea_margin=mid), f, wind)["rate"] < rate_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), f


def main():
    mrv = load_mrv()
    rows = mrv["rows"]
    py = per_year(rows)
    pooled = pooled_mrv(rows, py)
    print("EU MRV, GALICIA (IMO 9856189)")
    print(f"{'año':5s} {'comb t':>8s} {'CO2/comb':>8s} {'intraUE':>7s} {'h mar':>7s} {'mn':>8s} {'kn':>6s} "
          f"{'t/h mar':>8s} {'kg/mn':>7s} {'atr UE t':>8s} {'cota kg/h':>9s}")
    for y in sorted(rows):
        r, q = rows[y], py[y]
        print(f"{y:5s} {r['fuel']:8.0f} {q['ef']:8.3f} {100*q['intra_eu']:6.1f}% {r['sea_h']:7.0f} {q['dist']:8.0f} "
              f"{q['v_mean']:6.2f} {q['rate_sea']:8.3f} {q['fpd_sea']:7.1f} {q['f_berth']:8.0f} {q['berth_lb_kgph']:9.0f}")
    ef_all = sum(rows[y]["co2"] for y in rows if y != "2020") / sum(rows[y]["fuel"] for y in rows if y != "2020")
    print(f"\n2021-2023 agregados: {pooled['rate_sea']:.3f} t/h en el mar, {pooled['fpd_sea']:.1f} kg/mn, "
          f"{pooled['v_mean']:.2f} kn, {pooled['sea_h']/3:.0f} h/año en el mar")
    print(f"Factor CO2/combustible del buque 2021-2025: {ef_all:.4f} t/t")

    wind = rm.load_wind()
    out = dict(mrv_2021_2023=pooled, ef_2021_2025=ef_all, por_anio=py, casos={})
    v_t, r_t = pooled["v_mean"], pooled["rate_sea"]
    # comparación como en la versión 3: el modelo a su propia velocidad media (todo el año, sin viento ni curva SFOC)
    p3 = rm.P(v3=True)
    own = model_at(p3, 1.0, wind)
    out["casos"]["versión 3, a su propia velocidad"] = own
    print(f"\nVersión 3 a su propia velocidad ({own['v_mean']:.2f} kn): {own['rate']:.3f} t/h "
          f"({100 * (own['rate'] / r_t - 1):+.1f} %), {own['fpd']:.1f} kg/mn ({100 * (own['fpd'] / pooled['fpd_sea'] - 1):+.1f} %)")
    # a igual velocidad
    steps = [("ley cúbica, sin viento, SFOC constante (versión 3)", rm.P(v3=True, nm_open=rm.NM_OPEN_TOTAL)),
             ("+ resistencia del viento sobre el buque", rm.P(sfoc_curve=False, sea_margin=0.0)),
             ("+ SFOC según la carga", rm.P(sea_margin=0.0)),
             ("exponente 4, con viento y SFOC según la carga", rm.P(n_prop=4.0, sea_margin=0.0))]
    for label, p in steps:
        f = speed_factor(p, v_t, wind)
        m = model_at(p, f, wind)
        m["factor_velocidad"] = f
        out["casos"][label + ", a igual velocidad"] = m
        print(f"{label:58s} a {m['v_mean']:.2f} kn: {m['rate']:.3f} t/h ({100 * (m['rate'] / r_t - 1):+.1f} %), "
              f"{m['fpd']:.1f} kg/mn")
    out["margen_de_mar"] = {}
    for n in (3.0, 4.0):
        mu, f = calibrate(n, v_t, r_t, wind)
        out["margen_de_mar"][str(n)] = dict(mu=mu, factor_velocidad=f)
        chk = model_at(rm.P(n_prop=n, sea_margin=mu), f, wind)
        print(f"Margen de mar calibrado, exponente {n:g}: mu = {mu:.4f} -> {chk['rate']:.3f} t/h a {chk['v_mean']:.2f} kn")
        if abs(rm.SEA_MARGIN.get(n, -1.0) - round(mu, 3)) > 5e-4:
            print(f"  AVISO: ruta_modelo.SEA_MARGIN[{n}] = {rm.SEA_MARGIN.get(n)}; el calibrado es {mu:.3f}")
    lb = [py[y]["berth_lb_kgph"] for y in ("2021", "2022")]
    out["atraque_cota_kgph"] = lb
    print(f"\natraque: modelo {rm.P().hotel_berth * rm.SFOC_AE / 1e3:.0f} kg/h (versión 3: "
          f"{rm.P_HOTEL_BERTH_V3 * rm.SFOC_AE / 1e3:.0f} kg/h); cota inferior MRV 2021-2022 "
          f"{min(lb):.0f}-{max(lb):.0f} kg/h")
    json.dump(out, open(os.path.join(HERE, "calibracion_mrv.json"), "w"), indent=1, ensure_ascii=False)
    return out


if __name__ == "__main__":
    main()
