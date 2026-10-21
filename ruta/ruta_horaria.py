#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A1 - Modelo horario de travesías (19-10-2026). Simula cada travesía del horario de 2026 sobre el viento ERA5
horario de 2015-2024 (era5_horario_ruta.csv.gz): posición, rumbo y velocidad del buque en cada hora, viento
verdadero y aparente, y beneficio de cada modo de cada topología con ruta_modelo.bin_terms (la misma física que el
modelo de histogramas, incluida la resistencia inducida).

Estrategias, por travesía:
  inst    el mejor modo en cada hora, sin coste de conmutación (equivale al modelo de histogramas);
  voraz   el mejor modo en cada hora, pagando cada cambio;
  optima  secuencia óptima con coste de conmutación por programación dinámica, conociendo el viento de toda la
          travesía (cota superior de cualquier control real);
  reglas  reglas de la patente para el Savonius (lectura literal y prevista), pagando cada cambio.
Coste de conmutación: cada cambio de modo dura tau horas, durante las que la unidad rinde como rotor parado; se
cobra la diferencia con el mejor de los dos modos (cota prudente). El dispositivo sale recogido y llega recogido.
La generación en puerto se calcula en las horas reales de atraque del horario.

Uso: python3 ruta_horaria.py [era5_horario_ruta.csv.gz]   ->  resultados_horarios.json
"""
import json
import math
import sys
import datetime as dt
from zoneinfo import ZoneInfo
import numpy as np
import ruta_modelo as rm
import era5_horario as eh

LON, MAD = ZoneInfo("Europe/London"), ZoneInfo("Europe/Madrid")
T0, N_HORAS = eh.T0, eh.N_HORAS
RUTA = [p for p in eh.PUNTOS if p[1] in ("C", "B")]          # sentido sur: C1 ... C4, B1 ... B6
NM_OPEN = 540.0
EDGES = np.concatenate([[0.0], np.cumsum([p[5] for p in RUTA])]) * NM_OPEN / sum(p[5] for p in RUTA)
TAUS = (0.1, 0.25, 0.5)
WEEKS = sum(rm.DAYS.values()) / 7.0


def gc_nm(a, b):
    p1, p2, dl = math.radians(a[0]), math.radians(b[0]), math.radians(b[1] - a[1])
    d = 2 * math.asin(math.sqrt(math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2))
    return d * 3440.065


# Aproximaciones: Portsmouth-Nab-St Catherine y aproximación-Bilbao (ruta_geometria.py), escaladas a 32 mn a 10 kn
_NM_P = gc_nm((50.790, -1.110), (50.670, -0.950)) + gc_nm((50.670, -0.950), (50.500, -1.300))
_NM_B = gc_nm((43.450, -3.080), (43.360, -3.030))
T_APP = {"P": _NM_P * rm.NM_APPROACH / (_NM_P + _NM_B) / rm.V_APPROACH,
         "B": _NM_B * rm.NM_APPROACH / (_NM_P + _NM_B) / rm.V_APPROACH}

# Horario de 2026 en horas locales: (nombre, sentido, salida (día ISO, h, min, zona), llegada (días después, h, min, zona))
HORARIO = [("Portsmouth-Bilbao, jueves", "S", (4, 19, 0, LON), (2, 8, 0, MAD)),
           ("Portsmouth-Bilbao, domingo", "S", (7, 22, 0, LON), (2, 8, 0, MAD)),
           ("Bilbao-Portsmouth, martes", "N", (2, 13, 0, MAD), (1, 20, 45, LON)),
           ("Bilbao-Portsmouth, sabado", "N", (6, 13, 0, MAD), (1, 17, 30, LON))]
ATRAQUES = [("BILBAO", "sábado", (6, 8, 0, MAD), (0, 13, 0, MAD)), ("BILBAO", "martes", (2, 8, 0, MAD), (0, 13, 0, MAD)),
            ("PORTSMOUTH", "domingo", (7, 17, 30, LON), (0, 22, 0, LON)),
            ("PORTSMOUTH", "miércoles-jueves", (3, 20, 45, LON), (1, 19, 0, LON))]


def _hours(d, hh, mm, tz):
    return (dt.datetime(d.year, d.month, d.day, hh, mm, tzinfo=tz) - T0).total_seconds() / 3600.0


def schedule(table, shift=0.0):
    """Intervalos (horas desde 2015-01-01 UTC) de cada travesía o atraque que cae entero dentro de los datos.
    shift: desplazamiento (h) de todo el horario, para el diagnóstico de la hora del día."""
    out, d = [], dt.date(2014, 12, 20)
    while d <= dt.date(2025, 1, 10):
        for name, tag, (wd, hh, mm, tz), (nd, ah, am, atz) in table:
            if d.isoweekday() == wd:
                t1 = _hours(d, hh, mm, tz) + shift
                t2 = _hours(d + dt.timedelta(days=nd), ah, am, atz) + shift
                if t1 >= 0 and t2 <= N_HORAS - 1:
                    out.append(dict(name=name, tag=tag, t1=t1, t2=t2, year=d.year,
                                    season=0 if 4 <= d.month <= 11 else 1))
        d += dt.timedelta(days=1)
    return out


def wind_fields(data):
    f = {}
    for pid, dd in data.items():
        u40, u55, _ = eh.heights(dd)
        bad = np.isnan(dd["u10"]) | np.isnan(dd["u100"]) | np.isnan(dd["d100"])   # como en la agregación
        u40, u55 = np.where(bad, np.nan, u40), np.where(bad, np.nan, u55)
        f[pid] = dict(u40=u40, ratio=np.where(u40 > 0, u55 / np.where(u40 > 0, u40, 1.0), 1.0), u55=u55, d100=dd["d100"])
    return f


def era5_slots(t1, t2):
    """Reparto del intervalo [t1, t2] (horas desde T0) entre valores horarios de ERA5: cada valor es instantáneo,
    a la hora en punto, y representa la hora centrada en él. Devuelve (hora ERA5, inicio, fin) de cada tramo.
    (Corrige la versión del 19-10-2026 por la mañana, que asignaba cada hora con round() y, por el redondeo al par,
    solo usaba horas UTC pares.)"""
    out = []
    for h in range(int(math.floor(t1 + 0.5)), int(math.floor(t2 + 0.5)) + 1):
        lo, hi = max(t1, h - 0.5), min(t2, h + 0.5)
        if hi > lo:
            out.append((min(N_HORAS - 1, max(0, h)), lo, hi))
    return out


def crossing_hours(c):
    """(peso h, índice de punto, hora ERA5) de cada tramo de mar abierto, y velocidad (kn)."""
    o, de = ("P", "B") if c["tag"] == "S" else ("B", "P")
    a, b = c["t1"] + T_APP[o], c["t2"] - T_APP[de]
    v_kn = NM_OPEN / (b - a)
    hs = []
    for he, lo, hi in era5_slots(a, b):
        tm = 0.5 * (lo + hi)
        dist = v_kn * (tm - a)
        s = dist if c["tag"] == "S" else NM_OPEN - dist
        ip = int(min(len(RUTA) - 1, max(0, np.searchsorted(EDGES, s, side="right") - 1)))
        hs.append((hi - lo, ip, he))
    return v_kn, b - a, hs


def evaluate(fields, cross, p, binned=None):
    """Beneficio por unidad (kW) de cada estado en cada hora de cada travesía.
    binned: grupos de ruta_modelo.load_wind(); si se da, el viento se lleva al centro de su casilla (1 m/s, 15°) y
    la razón u55/u40 a la media de la casilla, como en el modelo de histogramas (diagnóstico)."""
    recs = []
    for ci, c in enumerate(cross):
        v_kn, t_open, hs = crossing_hours(c)
        v = v_kn * rm.KN
        rows = []
        for w, ip, he in hs:
            pid, psi = RUTA[ip][0], RUTA[ip][4]
            fd = fields[pid]
            u, r, d100 = fd["u40"][he], fd["ratio"][he], fd["d100"][he]
            if not (np.isfinite(u) and np.isfinite(r) and np.isfinite(d100)):
                continue
            g = float(eh.fold(d100 - psi))
            if binned is not None:
                iu = min(eh.NU - 1, int(math.floor(u)))
                u, r = iu + 0.5, float(binned[RUTA[ip][1]]["ratio"][c["season"]][iu])
                g = (min(eh.NG - 1, int(math.floor(g / 15.0))) + 0.5) * 15.0
            if c["tag"] == "N":
                g = 180.0 - g
            t = rm.bin_terms(float(u), float(r), g, v, p)
            rows.append((w, t["rot"], t["b_stop"], np.nan if t["ext_raw"] is None else t["ext_raw"], t["vawt_raw"],
                         t["vawt_ok"], np.nan if t["sav"] is None else t["sav"],
                         np.nan if t["sav3"] is None else t["sav3"], t["running"], t["cosb"], t["pen_r1"]))
        recs.append(dict(c=c, v_kn=v_kn, t_open=t_open, a=np.array(rows, float).reshape(-1, 11)))
    return recs


NEG = -1e18


def mode_matrix(a, topo):
    w, r1, bs, ext, vt, vok, sav, sav3 = (a[:, i] for i in range(8))
    ext = np.where(np.isfinite(ext), ext, np.nan)
    m = {"R1": [r1],
         "T1": [r1 + vt],
         "T2": [r1 + vt, sav + vt],
         "T3": [r1, r1 + ext, sav, sav3],
         "T3x": [r1, r1 + ext],
         "T4": [r1, np.where(vok > 0, r1 + vt, np.nan)],
         "T5": [r1, sav]}[topo]
    B = np.column_stack(m)
    return np.where(np.isfinite(B), B, NEG), w, bs


def switch_cost(B, bs, t, i, j, tau):
    return tau * max(0.0, max(B[t, i], B[t, j]) - bs[t])


def path_value(B, w, bs, path, tau):
    """Energía (kWh por unidad) de una secuencia de modos, con salida y llegada en el modo 0; y número de cambios."""
    val, prev, n_sw = 0.0, 0, 0
    for t, m in enumerate(path):
        if m != prev:
            val -= switch_cost(B, bs, t, prev, m, tau)
            n_sw += 1
        val += w[t] * B[t, m]
        prev = m
    if prev != 0 and len(path):
        val -= switch_cost(B, bs, len(path) - 1, prev, 0, tau)
        n_sw += 1
    return val, n_sw


def dp(B, w, bs, tau):
    n, M = B.shape
    if n == 0:
        return 0.0, 0
    V = np.full(M, -np.inf)
    V[0] = 0.0
    back = np.zeros((n, M), int)
    for t in range(n):
        nV = np.full(M, -np.inf)
        for j in range(M):
            if B[t, j] <= NEG / 2:
                continue
            best, arg = -np.inf, 0
            for i in range(M):
                if not np.isfinite(V[i]):
                    continue
                val = V[i] - (0.0 if i == j else switch_cost(B, bs, t, i, j, tau))
                if val > best:
                    best, arg = val, i
            nV[j] = best + w[t] * B[t, j]
            back[t, j] = arg
        V = nV
    final = np.array([V[j] - (0.0 if j == 0 else switch_cost(B, bs, n - 1, j, 0, tau)) if np.isfinite(V[j]) else -np.inf
                      for j in range(M)])
    j = int(np.argmax(final))
    path = [j]
    for t in range(n - 1, 0, -1):
        j = int(back[t, j])
        path.append(j)
    path = path[::-1]
    val, n_sw = path_value(B, w, bs, path, tau)
    return val, n_sw


def causal_path(B, w, bs, tau, L):
    """Control causal: en cada hora solo se conoce el viento de esa hora y se supone que persiste L horas (o lo que
    quede de travesía); se cambia al modo j si (B_j - B_actual) * min(L, resto) supera el coste del cambio. Si el modo
    actual no está disponible (límite estructural, parada), se cambia al mejor disponible."""
    n, M = B.shape
    path, m, left = [], 0, float(w.sum())
    for t in range(n):
        if B[t, m] <= NEG / 2:
            m = int(np.argmax(B[t]))
        else:
            best_j, best_g = m, 0.0
            for j in range(M):
                if j == m or B[t, j] <= NEG / 2:
                    continue
                g = (B[t, j] - B[t, m]) * min(L, left) - switch_cost(B, bs, t, m, j, tau)
                if g > best_g:
                    best_j, best_g = j, g
            m = best_j
        path.append(m)
        left -= w[t]
    return path


LS = (1.0, 2.0, 3.0, 6.0, 12.0)


def strategies(recs, topo, taus=TAUS, Ls=LS):
    """kWh por unidad y travesía, y cambios por travesía, para cada estrategia."""
    out = {"inst": [], **{f"voraz {t}": [] for t in taus}, **{f"optima {t}": [] for t in taus},
           **{f"cambios voraz {t}": [] for t in taus}, **{f"cambios optima {t}": [] for t in taus},
           **{f"causal L{L:g} {t}": [] for t in taus for L in Ls},
           **{f"cambios causal L{L:g} {t}": [] for t in taus for L in Ls}}
    if topo == "T5":
        out.update({f"{r} {t}": [] for r in ("literal", "prevista") for t in (0.0,) + taus})
    for rec in recs:
        a = rec["a"]
        B, w, bs = mode_matrix(a, topo)
        out["inst"].append(float((w * B.max(axis=1)).sum()))
        greedy = list(np.argmax(B, axis=1))
        for t in taus:
            gv, gs = path_value(B, w, bs, greedy, t)
            out[f"voraz {t}"].append(gv)
            out[f"cambios voraz {t}"].append(gs)
            dv, ds = dp(B, w, bs, t) if B.shape[1] > 1 else (gv, gs)
            out[f"optima {t}"].append(dv)
            out[f"cambios optima {t}"].append(ds)
            for L in Ls:
                cv, cs = path_value(B, w, bs, causal_path(B, w, bs, t, L), t) if B.shape[1] > 1 else (gv, gs)
                out[f"causal L{L:g} {t}"].append(cv)
                out[f"cambios causal L{L:g} {t}"].append(cs)
        if topo == "T5":
            run, cosb, sav_ok = a[:, 8] > 0, a[:, 9], B[:, 1] > NEG / 2
            rules = {"literal": (~run) & sav_ok, "prevista": (~run) & (cosb < 0) & sav_ok}
            for r, mask in rules.items():
                for t in (0.0,) + taus:
                    out[f"{r} {t}"].append(path_value(B, w, bs, list(mask.astype(int)), t)[0])
    return out


def annualise(recs, values, by=None):
    """t/año de combustible (cuatro unidades) a partir de kWh por unidad y travesía: media por tipo de travesía
    por semanas del año (o de la estación, o de cada año si by='year')."""
    vals = np.asarray(values, float)
    names = sorted({r["c"]["name"] for r in recs})
    if by == "year":
        res = {}
        for y in sorted({r["c"]["year"] for r in recs}):
            tot = 0.0
            for nme in names:
                sel = [k for k, r in enumerate(recs) if r["c"]["name"] == nme and r["c"]["year"] == y]
                if sel:
                    tot += vals[sel].mean() * WEEKS
            res[y] = rm.N_UNITS * tot * rm.SFOC_ME / 1e6
        return res
    tot = 0.0
    for s in ((0, 1) if by is None else (by,)):
        for nme in names:
            sel = [k for k, r in enumerate(recs) if r["c"]["name"] == nme and r["c"]["season"] == s]
            if sel:
                tot += vals[sel].mean() * rm.DAYS[s] / 7.0
    return rm.N_UNITS * tot * rm.SFOC_ME / 1e6


def r0_fuel(recs):
    """Combustible del buque sin rotores (t/año), con las mismas hipótesis que ruta_modelo.annual."""
    per = []
    for r in recs:
        v = r["v_kn"]
        p_b0 = rm.LOAD_DESIGN * rm.P_MCR * (v / rm.V_DESIGN) ** rm.N_PROP
        p_app = rm.LOAD_DESIGN * rm.P_MCR * (rm.V_APPROACH / rm.V_DESIGN) ** rm.N_PROP
        hours = r["t_open"] + rm.NM_APPROACH / rm.V_APPROACH
        per.append(((p_b0 * r["t_open"] + p_app * rm.NM_APPROACH / rm.V_APPROACH) * rm.SFOC_ME
                    + rm.P_HOTEL_SEA * hours * rm.SFOC_ME / rm.ETA_PTO) / 1e6)
    sea = annualise(recs, per) / rm.N_UNITS / rm.SFOC_ME * 1e6     # annualise multiplica por N y SFOC
    berth = rm.P_HOTEL_BERTH * sum(rm.PORT_HOURS_WEEK.values()) * WEEKS * rm.SFOC_AE / 1e6
    return sea + berth


def berth_generation(fields, p):
    """Energía en atraque por unidad (kWh por estancia) en las horas reales del horario; t/año para cuatro unidades."""
    stays = schedule(ATRAQUES)
    curves = {"vawt": (None, None, "u55"), "sav": (p.ds * rm.A_F, p.cp_s, "u40"),
              "sav3": (p.ds * (rm.A_F + rm.A_X), p.cp_s, "u40")}
    res = {k: {} for k in curves}
    for st in stays:
        fd = fields[st["name"]]
        e = {k: 0.0 for k in curves}
        for he, lo, hi in era5_slots(st["t1"], st["t2"]):
            for k, (area, cp, z) in curves.items():
                u = fd[z][he]
                if np.isfinite(u):
                    pw = rm.vawt_power_curve(np.array([u]), p)[0] if k == "vawt" else \
                        rm.power_curve(area, cp, np.array([u]))[0]
                    e[k] += (hi - lo) * float(pw)
        for k in curves:
            res[k].setdefault((st["name"], st["tag"]), []).append(e[k])
    out = {}
    for k in curves:
        kwh_year = sum(np.mean(v) * WEEKS for v in res[k].values())
        out[k] = dict(kW_medio_unidad=kwh_year / (sum(rm.PORT_HOURS_WEEK.values()) * WEEKS),
                      t_anio=rm.N_UNITS * kwh_year * rm.SFOC_AE / 1e6)
    return out


def stats(d):
    v = np.array(list(d.values()), float)
    return dict(media=float(v.mean()), dt=float(v.std(ddof=1)), min=float(v.min()), max=float(v.max()),
                cv=float(v.std(ddof=1) / v.mean()) if v.mean() else float("nan"))


def run(path="era5_horario_ruta.csv.gz", p=None, taus=TAUS, quiet=False):
    p = p or rm.P()
    data, meta = eh.load_hourly(path)
    fields = wind_fields(data)
    cross = schedule(HORARIO)
    recs = evaluate(fields, cross, p)
    res = {"travesias": len(recs), "horas_mar": float(sum(r["a"][:, 0].sum() for r in recs)),
           "horas_sin_dato": int(sum(len(crossing_hours(r["c"])[2]) - len(r["a"]) for r in recs)),
           "R0_t_anio": r0_fuel(recs), "tau_h": list(taus), "L_h": list(LS)}
    per = {topo: strategies(recs, topo, taus) for topo in ("R1", "T1", "T2", "T3", "T3x", "T4", "T5")}
    r1 = per["R1"]["inst"]
    res["R1"] = {"inst": annualise(recs, r1), "abr-nov": annualise(recs, r1, by=0), "dic-mar": annualise(recs, r1, by=1),
                 "por_anio": annualise(recs, r1, by="year")}
    res["R1"]["por_anio_estad"] = stats(res["R1"]["por_anio"])
    pen = [float((r["a"][:, 0] * r["a"][:, 10]).sum()) for r in recs]
    res["R1"]["inducida"] = annualise(recs, pen)
    res["incrementos_mar"] = {}
    for topo in ("T1", "T2", "T3", "T3x", "T4", "T5"):
        d = {}
        for k, v in per[topo].items():
            if k.startswith("cambios"):
                d[k] = float(np.mean(v))
            else:
                d[k] = annualise(recs, v) - res["R1"]["inst"]
        d["abr-nov inst"] = annualise(recs, per[topo]["inst"], by=0) - res["R1"]["abr-nov"]
        d["dic-mar inst"] = annualise(recs, per[topo]["inst"], by=1) - res["R1"]["dic-mar"]
        for k in ["inst"] + [f"optima {t}" for t in taus]:
            yr = {y: v - res["R1"]["por_anio"][y] for y, v in annualise(recs, per[topo][k], by="year").items()}
            d[f"por_anio {k}"] = yr
            d[f"por_anio {k} estad"] = stats(yr)
        res["incrementos_mar"][topo] = d
    res["atraque"] = berth_generation(fields, p)
    res["meta_era5"] = meta[:3]
    if not quiet:
        report(res)
    return res


def report(res):
    print(f"Travesías simuladas: {res['travesias']}; horas en el mar {res['horas_mar']:.0f}; "
          f"horas sin dato {res['horas_sin_dato']}; R0 {res['R0_t_anio']:.0f} t/año")
    r1 = res["R1"]
    st = r1["por_anio_estad"]
    print(f"R1 (cuatro rotores): {r1['inst']:.1f} t/año (abr-nov {r1['abr-nov']:.1f}, dic-mar {r1['dic-mar']:.1f}); "
          f"inducida pagada {r1['inducida']:.1f} t/año")
    print("  por año: " + ", ".join(f"{y} {v:.0f}" for y, v in r1["por_anio"].items()) +
          f"  | media {st['media']:.1f}, dt {st['dt']:.1f} (CV {100 * st['cv']:.1f} %), {st['min']:.0f}-{st['max']:.0f}")
    for topo, d in res["incrementos_mar"].items():
        print(f"{topo:4s} inst {d['inst']:7.1f} (abr-nov {d['abr-nov inst']:.1f}, dic-mar {d['dic-mar inst']:.1f})")
        for t in res["tau_h"]:
            line = (f"     tau {t}: óptima {d[f'optima {t}']:6.1f} ({d[f'cambios optima {t}']:.1f} cambios)"
                    f"  voraz {d[f'voraz {t}']:6.1f} ({d[f'cambios voraz {t}']:.1f})")
            for L in res["L_h"]:
                line += f"  causal L{L:g} {d[f'causal L{L:g} {t}']:6.1f} ({d[f'cambios causal L{L:g} {t}']:.1f})"
            print(line)
        s = d["por_anio inst estad"]
        print(f"     por año (inst): " + ", ".join(f"{y} {v:.1f}" for y, v in d["por_anio inst"].items()) +
              f" | media {s['media']:.1f}, dt {s['dt']:.1f}, {s['min']:.1f}-{s['max']:.1f}")
        if topo == "T5":
            print("     reglas: " + ", ".join(f"{k} {d[k]:.1f}" for k in d if k.startswith(("literal", "prevista"))))
    print("Atraque (horas reales): " + ", ".join(f"{k} {v['kW_medio_unidad']:.2f} kW/unidad, {v['t_anio']:.1f} t/año"
                                             for k, v in res["atraque"].items()))


if __name__ == "__main__":
    out = run(sys.argv[1] if len(sys.argv) > 1 else "era5_horario_ruta.csv.gz")
    json.dump(out, open("resultados_horarios.json", "w"), indent=1, ensure_ascii=False, default=float)
