#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Caso de validación del rotor convencional (R1) a escala real: Viking Grace (Viking Line), con un rotor Norsepower
de 24 m x 4 m (Lloyd's Register, nota de prensa del 30-04-2018), el mismo tamaño que el de este artículo, en la ruta
Turku-Mariehamn-Estocolmo. Versión 4 (22-10-2026), en respuesta a la revisión simulada (R1-1, R2-6).

Dato de contraste (Norsepower, nota del 16-05-2019): con un análisis de regresión de NAPA, el desglose de la potencia
propulsora con datos de a bordo de ABB y galgas extensométricas en el rotor, el ahorro verificado es de 231-315 t de
GNL al año, «equalling an average propulsion power between 207kW and 282kW», a una velocidad de servicio de 21 kn.
La nota no dice si la cifra es neta del consumo eléctrico del rotor; por cómo se obtuvo (potencia propulsora y
empuje medido), aquí se trata como bruta de ese consumo. Norsepower (presentación de 2018) da para sus rotores un
consumo eléctrico medio de 15-35 kW y una velocidad de giro máxima de 250 rpm.

Método: el mismo código del rotor que en Portsmouth-Bilbao (ruta_modelo.bin_terms con r1_only = True): coeficientes,
relación de giro óptima limitada por k_max y por la velocidad periférica máxima, potencia de giro y rotor parado.
Viento: ERA5 horario 2015-2024 en 10 puntos de la ruta (descargar_era5_viking_grace.py; Open-Meteo, models=era5),
extrapolado a 40 m con el exponente horario entre 10 y 100 m (como en la ruta del Galicia), en histogramas de 1 m/s y
15 grados por punto, tramo y sentido. Cada punto representa la mitad de los tramos que lo rodean; los puertos de Turku
y Estocolmo se asignan al punto más próximo. Las coordenadas aproximan la ruta, no el track real.
Buque (Wärtsilä, ficha del Viking Grace): propulsión eléctrica con cuatro motores 8L50DF de doble combustible
(30 400 kW), 21,8 kn de servicio, 6,80 m de calado de proyecto. El rotor se alimenta de la red del buque: kappa = 1.
Magnitudes (kW, medias de todas las horas de navegación, en ambos sentidos):
  empuje      F_x V / eta_D con el rotor en marcha (potencia propulsora equivalente al empuje, bruta);
  empuje_ef   F_x V (potencia efectiva del empuje);
  giro        consumo eléctrico del rotor; parado: resistencia del rotor parado / eta_D;
  inducida    resistencia inducida por la fuerza lateral del rotor (solo el término propio: la del viento sobre el
              Viking Grace no se conoce), con el calado de 6,80 m;
  neto        empuje - giro - inducida - parado (lo que el modelo de ruta cuenta como ahorro de R1).
Salida: validacion_viking_grace.json.
"""
import datetime as dt
import gzip
import json
import math
import os
import numpy as np
import ruta_modelo as rm
import era5_horario as eh
import coeficientes as cf

HERE = os.path.dirname(os.path.abspath(__file__))
DATOS = os.path.join(HERE, "..", "datos", "era5_horario_viking_grace.csv.gz")
N_HORAS = 87672
T0 = dt.datetime(2015, 1, 1, tzinfo=dt.timezone.utc)
NU, NG = 30, 12
NM_KM = 1.852

# puntos ERA5 (id, latitud, longitud, celda de tierra según la elevación de Open-Meteo) y puertos
PUNTOS = {"VG1": (60.380, 22.000), "VG2": (60.240, 21.620), "VG3": (60.150, 21.150), "VG4": (60.050, 20.650),
          "VG5": (60.050, 20.200), "MARIEHAMN": (60.095, 19.935), "VG6": (59.950, 19.500), "VG7": (59.800, 19.150),
          "VG8": (59.620, 18.850), "VG9": (59.420, 18.400)}
PUERTOS = {"TURKU": (60.434, 22.222), "ESTOCOLMO": (59.316, 18.090)}
RUTA = ["TURKU", "VG1", "VG2", "VG3", "VG4", "VG5", "MARIEHAMN", "VG6", "VG7", "VG8", "VG9", "ESTOCOLMO"]

T_VG = 6.80                 # m, calado de proyecto (Wärtsilä)
V_SERVICE = 21.0            # kn, velocidad de servicio de la nota de Norsepower (21,8 kn según Wärtsilä)
REF_KW = (207.0, 282.0)     # Norsepower (2019): potencia propulsora media verificada
REF_T = (231.0, 315.0)      # t de GNL al año
PERIODO_ENSAYO = (dt.datetime(2018, 5, 1, tzinfo=dt.timezone.utc), dt.datetime(2019, 5, 1, tzinfo=dt.timezone.utc))


def load_vg(path=DATOS):
    data = {p: {k: np.full(N_HORAS, np.nan) for k in ("u10", "u100", "d100")} for p in PUNTOS}
    meta = []
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                meta.append(line.strip())
                continue
            if line.startswith("punto,"):
                continue
            pid, h, a, b, c = line.rstrip("\n").split(",")
            h = int(h)
            if pid in data and 0 <= h < N_HORAS:
                dd = data[pid]
                if a:
                    dd["u10"][h] = float(a)
                if b:
                    dd["u100"][h] = float(b)
                if c:
                    dd["d100"][h] = float(c)
    return data, meta


def land_cells(meta):
    """Puntos cuya celda ERA5 tiene elevación > 0 en la cabecera (celda con tierra)."""
    out = set()
    for m in meta:
        if m.startswith("# VG") or m.startswith("# MARIEHAMN"):
            pid = m[2:].split(":")[0]
            elev = float(m.split("elevación")[1].split(";")[0])
            if elev > 0.0:
                out.add(pid)
    return out


def coord(n):
    return PUNTOS.get(n) or PUERTOS[n]


def dist_bearing(a, b):
    """Distancia (mn) y rumbo inicial (grados) de a a b (ortodrómica)."""
    la1, lo1, la2, lo2 = map(math.radians, (a[0], a[1], b[0], b[1]))
    dlo = lo2 - lo1
    h = math.sin((la2 - la1) / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin(dlo / 2) ** 2
    d_km = 2 * 6371.0 * math.asin(math.sqrt(h))
    y = math.sin(dlo) * math.cos(la2)
    x = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(dlo)
    return d_km / NM_KM, (math.degrees(math.atan2(y, x)) + 360.0) % 360.0


def subtramos():
    """(punto ERA5, rumbo hacia Estocolmo, millas): cada tramo se parte en dos mitades, cada una con su punto."""
    out = []
    for a, b in zip(RUTA[:-1], RUTA[1:]):
        d, psi = dist_bearing(coord(a), coord(b))
        first = a if a in PUNTOS else b
        second = b if b in PUNTOS else a
        out += [(first, psi, d / 2.0), (second, psi, d / 2.0)]
    return out


def histograms(data, mask_hours=None, excluir=()):
    """Histogramas normalizados (NU x NG) por subtramo y sentido, con el peso de sus millas."""
    hist = []
    sub = [s for s in subtramos() if s[0] not in excluir]
    w_tot = sum(s[2] for s in sub)
    for pid, psi, nm in sub:
        dd = data[pid]
        ok = ~(np.isnan(dd["u10"]) | np.isnan(dd["u100"]) | np.isnan(dd["d100"]))
        if mask_hours is not None:
            ok &= mask_hours
        ur, _, _ = eh.heights(dd)
        iu = np.minimum(NU - 1, np.floor(np.nan_to_num(ur)).astype(int))
        for sentido, rumbo in (("O", psi), ("E", (psi + 180.0) % 360.0)):
            ig = np.minimum(NG - 1, np.floor(eh.fold(dd["d100"] - rumbo) / 15.0).astype(int))
            cnt = np.bincount((iu * NG + ig)[ok], minlength=NU * NG).reshape(NU, NG).astype(float)
            hist.append((pid, sentido, 0.5 * nm / w_tot, cnt / cnt.sum()))
    return hist, w_tot


def evaluate(hist, v_kn, p):
    """Medias (kW) de las magnitudes del rotor en todas las horas de navegación."""
    v = v_kn * rm.KN
    acc = dict(empuje=0.0, empuje_ef=0.0, giro=0.0, parado=0.0, inducida=0.0, neto=0.0, marcha=0.0, k_medio=0.0)
    cache = {}
    for pid, sentido, w, F in hist:
        for iu in range(NU):
            for ig in range(NG):
                f = F[iu, ig]
                if f == 0.0:
                    continue
                key = (iu, ig)
                if key not in cache:
                    t = rm.bin_terms(iu + 0.5, 1.0, (ig + 0.5) * 15.0, v, p)
                    run = t["running"]
                    cache[key] = dict(empuje=t["rot_fx"], empuje_ef=t["rot_fx"] * p.eta_d, giro=t["spin_r1"],
                                      parado=(-t["b_stop"] if not run else 0.0),
                                      inducida=(t["pen_r1"] if run else 0.0), neto=t["rot"],
                                      marcha=1.0 if run else 0.0, k_medio=t["k"])
                for k2, val in cache[key].items():
                    acc[k2] += w * f * val
    acc["k_medio"] = acc["k_medio"] / acc["marcha"] if acc["marcha"] > 0 else 0.0
    return acc


def p_vg(**kw):
    base = dict(n_units=1, r1_only=True, t_draught=T_VG, kappa=1.0, windage=False, r0_wind=False)
    base.update(kw)
    return rm.P(**base)


ESCENARIOS = [  # escenarios epistémicos de coeficientes (coeficientes.ESCENARIOS_COEF)
    ("A: De Marco et al. (2016), k <= 3, U_s <= 36 m/s", "A"),
    ("A+: A prolongado, k <= 4, U_s <= 46 m/s", "A+"),
    ("AB: mezcla de A y B a partes iguales (versión 3)", "AB"),
    ("B: Bordogna et al. (2019), AR_eff 12", "B"),
    ("B': B con la línea sustentadora (eta 0,83)", "B'"),
    ("B'': B sin reflexión en la cubierta (eta 0,725, AR_eff 6)", "B''"),
]


def main():
    data, meta = load_vg()
    land = land_cells(meta)
    sub = subtramos()
    tot = sum(s[2] for s in sub)
    print("Viking Grace: ruta Turku-Mariehamn-Estocolmo aproximada")
    print(f"  {tot:.1f} mn; celdas ERA5 con tierra: {sorted(land)}")
    hours = np.arange(N_HORAS)
    t_ms = np.array([(T0 + dt.timedelta(hours=int(h))) for h in hours])
    ensayo = np.array([PERIODO_ENSAYO[0] <= t < PERIODO_ENSAYO[1] for t in t_ms])
    h_all, _ = histograms(data)
    h_test, _ = histograms(data, mask_hours=ensayo)
    h_sea, _ = histograms(data, excluir=land)
    u40 = {pid: float(np.nanmean(eh.heights(data[pid])[0])) for pid in PUNTOS}
    print("  viento medio a 40 m (m/s): " + ", ".join(f"{k} {v:.1f}" for k, v in u40.items()))
    out = dict(fuente="Norsepower (2019); LR (2018); Wärtsilä", ref_kw=REF_KW, ref_t=REF_T, millas=tot,
               celdas_tierra=sorted(land), viento_40m=u40, escenarios={}, sensibilidad={})
    print(f"\nReferencia: {REF_KW[0]:.0f}-{REF_KW[1]:.0f} kW de potencia propulsora media (Norsepower, 2019)")
    hdr = f"{'escenario':58s} {'empuje':>7s} {'ef.':>6s} {'giro':>6s} {'parado':>6s} {'induc':>6s} {'neto':>6s} {'marcha':>6s} {'k':>5s}"
    print(f"\n21 kn, ERA5 2015-2024 (kW por rotor)\n{hdr}")
    for label, esc in ESCENARIOS:
        r = evaluate(h_all, V_SERVICE, p_vg(**cf.ESCENARIOS_COEF[esc]))
        out["escenarios"][label] = r
        print(f"{label:58s} {r['empuje']:7.1f} {r['empuje_ef']:6.1f} {r['giro']:6.1f} {r['parado']:6.1f} "
              f"{r['inducida']:6.1f} {r['neto']:6.1f} {100 * r['marcha']:5.0f}% {r['k_medio']:5.2f}")
    print(f"\nSensibilidad (escenario A)\n{hdr}")
    sens = [("velocidad 21,8 kn", h_all, 21.8, {}), ("velocidad 18 kn", h_all, 18.0, {}),
            ("velocidad 16 kn", h_all, 16.0, {}), ("periodo del ensayo (05-2018 a 04-2019)", h_test, V_SERVICE, {}),
            ("sin celdas con tierra", h_sea, V_SERVICE, {}), ("eta_D 0,65", h_all, V_SERVICE, dict(eta_d=0.65)),
            ("eta_D 0,75", h_all, V_SERVICE, dict(eta_d=0.75)), ("Cf 0,002 (giro)", h_all, V_SERVICE, dict(cf=0.002)),
            ("Cf 0,007 (giro)", h_all, V_SERVICE, dict(cf=0.007)),
            ("A+, sin celdas con tierra", h_sea, V_SERVICE, dict(k_max=4.0, us_max=46.0)),
            ("A+, sin celdas con tierra, viento x1,1 en el rotor", h_sea, V_SERVICE, dict(k_max=4.0, us_max=46.0, f_u=1.1))]
    for label, hh, v, kw in sens:
        r = evaluate(hh, v, p_vg(**kw))
        out["sensibilidad"][label] = r
        print(f"{label:58s} {r['empuje']:7.1f} {r['empuje_ef']:6.1f} {r['giro']:6.1f} {r['parado']:6.1f} "
              f"{r['inducida']:6.1f} {r['neto']:6.1f} {100 * r['marcha']:5.0f}% {r['k_medio']:5.2f}")
    # el mismo cálculo para el Galicia (por rotor, abril-noviembre, travesías del horario), para comparar
    d, groups, nu, ng = rm.load_wind(os.path.join(HERE, "era5_ruta.json"))
    gal = {}
    for label, esc in ESCENARIOS:
        p = rm.P(**cf.ESCENARIOS_COEF[esc])
        acc, hrs = {}, 0.0
        for name, direction, h in rm.crossings(p):
            v_kn, t_open = rm.crossing_speed(h, p)
            mb = rm.mean_benefits(groups, nu, ng, 0, direction, v_kn, p)
            for k2 in ("R1", "_R1 empuje", "_R1 giro", "_R1 parado", "_pen R1", "_R1 marcha"):
                acc[k2] = acc.get(k2, 0.0) + mb[k2] * t_open
            hrs += t_open
        gal[label] = {k2: v2 / hrs for k2, v2 in acc.items()}
    out["galicia_abr_nov"] = gal
    print("\nGalicia, abril-noviembre, mar abierto (kW por rotor, mismo código)")
    for label, r in gal.items():
        print(f"{label:58s} empuje {r['_R1 empuje']:6.1f}  giro {r['_R1 giro']:5.1f}  parado {r['_R1 parado']:5.1f}  "
              f"inducida {r['_pen R1']:5.1f}  neto {r['R1']:6.1f}  marcha {100 * r['_R1 marcha']:3.0f}%")
    json.dump(out, open(os.path.join(HERE, "validacion_viking_grace.json"), "w", encoding="utf-8"), indent=1,
              ensure_ascii=False)
    return out


if __name__ == "__main__":
    main()
