#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A1 - Modelo horario de travesías. Versión 4 (22-10-2026; la versión 3 está en el historial del repositorio).
Simula cada travesía del horario de 2026 sobre el viento ERA5 horario de 2015-2024 (era5_horario_ruta.csv.gz):
posición, rumbo y velocidad del buque en cada hora, viento verdadero y aparente, y beneficio de cada modo de cada
topología con ruta_modelo.bin_terms (la misma física que el modelo de histogramas).

Cambios de la versión 4 (revisión simulada R3-5 y R3-m4):
- Caso base: servicio documentado de abril a noviembre (ruta_modelo.BASE_SEASONS); diciembre-marzo como escenario.
- R0 con el modelo de la versión 4 (resistencia del viento sobre el buque, curva SFOC, margen de mar) y ahorros
  convertidos con el consumo marginal de cada travesía, como en ruta_modelo.annual.
- Coste de arranque y parada también para el rotor convencional: R1 tiene dos modos (parado, en marcha) y cada cambio
  dura tau_R1 (0,05-0,1 h de aceleración o frenado).
- Rotor abatible (R1t, comparador comercial): modos parado, en marcha y abatido (sin resistencia ni fuerza lateral);
  cada abatimiento o izado dura tau_t (0,25-0,5 h, supuesto).
- Control con previsión imperfecta (control predictivo con horizonte H): en cada hora se conoce el viento actual y,
  para las horas siguientes, una previsión emitida a la salida con errores AR(1) de velocidad y dirección que crecen
  con el plazo; se resuelve la secuencia óptima en el horizonte, se aplica el primer paso y se repite. El resultado se
  evalúa con el viento real. Dos niveles de error (supuestos, del orden de los de un modelo numérico a corto plazo en
  el mar): bueno, 1,0 -> 2,0 m/s y 10 -> 20 grados en 24 h; malo, 2,0 -> 3,5 m/s y 20 -> 40 grados.

Estrategias, por travesía:
  inst     el mejor modo en cada hora, sin coste de conmutación (equivale al modelo de histogramas);
  voraz    el mejor modo en cada hora, pagando cada cambio;
  optima   secuencia óptima con coste de conmutación por programación dinámica, conociendo el viento de toda la
           travesía (cota superior de cualquier control real);
  causal   en cada hora se supone que el viento actual persiste L horas (persistencia);
  mpc      control predictivo con la previsión imperfecta (horizonte H);
  reglas   reglas de la patente para el Savonius (lectura literal y prevista), pagando cada cambio.
Coste de conmutación: cada cambio de modo dura tau horas, durante las que la unidad rinde como rotor parado; se
cobra la diferencia con el mejor de los dos modos (cota prudente). Las topologías de la patente parten del mejor
estado de R1 en cada hora (sin coste propio de R1); el coste de arranque y parada de R1 se da aparte.
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
NM_OPEN = rm.NM_OPEN_TOTAL
EDGES = np.concatenate([[0.0], np.cumsum([p[5] for p in RUTA])]) * NM_OPEN / sum(p[5] for p in RUTA)
TAUS = (0.1, 0.25, 0.5)
TAU_R1 = (0.05, 0.1)
TAU_TILT = (0.25, 0.5)
LS = (1.0, 2.0, 3.0, 6.0, 12.0)
HS = (3, 6, 12)
FORECAST = {"buena": dict(su0=1.0, su1=2.0, sd0=10.0, sd1=20.0), "mala": dict(su0=2.0, su1=3.5, sd0=20.0, sd1=40.0)}
RHO_AR = 0.9
SEED = 20261022
COLS = ("w", "rot", "b_stop", "ext", "vawt", "vawt_ok", "sav", "sav3", "running", "cosb", "pen_r1", "rot_run",
        "r0_wind")


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
    a la hora en punto, y representa la hora centrada en él. Devuelve (hora ERA5, inicio, fin) de cada tramo."""
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


def forecast_errors(n, lead0, level, rng):
    """Errores AR(1) de velocidad (m/s) y dirección (grados) para n horas, con desviación típica que crece
    linealmente con el plazo (horas desde la emisión de la previsión, a la salida) hasta 24 h."""
    f = FORECAST[level]
    eu, ed = np.zeros(n), np.zeros(n)
    zu, zd = rng.standard_normal(), rng.standard_normal()
    for i in range(n):
        if i:
            zu = RHO_AR * zu + math.sqrt(1 - RHO_AR ** 2) * rng.standard_normal()
            zd = RHO_AR * zd + math.sqrt(1 - RHO_AR ** 2) * rng.standard_normal()
        lead = min(lead0 + i, 24.0) / 24.0
        eu[i] = zu * (f["su0"] + (f["su1"] - f["su0"]) * lead)
        ed[i] = zd * (f["sd0"] + (f["sd1"] - f["sd0"]) * lead)
    return eu, ed


def evaluate(fields, cross, p, binned=None, forecast=None):
    """Beneficio por unidad (kW) de cada estado en cada hora de cada travesía (columnas COLS).
    binned: grupos de ruta_modelo.load_wind(); si se da, el viento se lleva al centro de su casilla (1 m/s, 15°) y
    la razón u55/u40 a la media de la casilla, como en el modelo de histogramas (diagnóstico).
    forecast: nivel de error ('buena' o 'mala'); si se da, se evalúa con el viento previsto (mismas horas)."""
    recs = []
    rng = np.random.default_rng(SEED)
    for ci, c in enumerate(cross):
        v_kn, t_open, hs = crossing_hours(c)
        v = v_kn * rm.KN
        if forecast is not None:
            eu, ed = forecast_errors(len(hs), 0.0, forecast, rng)
        rows = []
        for k, (w, ip, he) in enumerate(hs):
            pid, psi = RUTA[ip][0], RUTA[ip][4]
            fd = fields[pid]
            u, r, d100 = fd["u40"][he], fd["ratio"][he], fd["d100"][he]
            if not (np.isfinite(u) and np.isfinite(r) and np.isfinite(d100)):
                continue
            if forecast is not None:
                u, d100 = max(0.0, u + eu[k]), d100 + ed[k]
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
                         np.nan if t["sav3"] is None else t["sav3"], t["running"], t["cosb"], t["pen_r1"],
                         t["rot_run"], t["r0_wind"]))
        a = np.array(rows, float).reshape(-1, len(COLS))
        p_wind = float((a[:, 0] * a[:, 12]).sum() / a[:, 0].sum()) if len(a) else 0.0
        fuel, f_marg, load = rm.crossing_fuel(v_kn, t_open, p_wind, p)
        recs.append(dict(c=c, v_kn=v_kn, t_open=t_open, a=a, fuel=fuel, f_marg=f_marg))
    return recs


NEG = -1e18


def mode_matrix(a, topo):
    w, r1, bs, ext, vt, vok, sav, sav3 = (a[:, i] for i in range(8))
    run_v = a[:, 11]
    ext = np.where(np.isfinite(ext), ext, np.nan)
    m = {"R1": [r1],
         "R1 conmutado": [bs, run_v],
         "R1t": [bs, run_v, np.zeros_like(bs)],
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


def _cost_matrix(B, bs, t, tau):
    """Coste de pasar del modo i al j en la hora t (diagonal nula)."""
    mx = np.maximum(B[t][:, None], B[t][None, :])
    c = tau * np.maximum(0.0, mx - bs[t])
    np.fill_diagonal(c, 0.0)
    return c


def dp_window(B, w, bs, tau, m0, end_mode=None):
    """Secuencia óptima en una ventana que empieza en el modo m0 (y acaba en end_mode, si se da)."""
    n, M = B.shape
    V = np.full(M, -np.inf)
    V[m0] = 0.0
    back = np.zeros((n, M), int)
    avail = B > NEG / 2
    for t in range(n):
        tot = V[:, None] - _cost_matrix(B, bs, t, tau)
        arg = np.argmax(tot, axis=0)
        best = tot[arg, np.arange(M)]
        nV = np.where(avail[t], best + w[t] * B[t], -np.inf)
        back[t] = arg
        V = nV
    if end_mode is not None:
        c_end = _cost_matrix(B, bs, n - 1, tau)[:, end_mode]
        final = np.where(np.isfinite(V), V - c_end, -np.inf)
    else:
        final = V
    j = int(np.argmax(final))
    path = [j]
    for t in range(n - 1, 0, -1):
        j = int(back[t, j])
        path.append(j)
    return path[::-1]


def dp(B, w, bs, tau):
    if B.shape[0] == 0:
        return 0.0, 0
    path = dp_window(B, w, bs, tau, 0, end_mode=0)
    return path_value(B, w, bs, path, tau)


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


def mpc_path(B, Bf, w, bs, bsf, tau, H):
    """Control predictivo: en la hora t se conoce B[t]; las horas t+1 ... t+H-1 se toman de la previsión Bf. Se
    resuelve la ventana empezando en el modo actual, se aplica el primer paso y se avanza una hora."""
    n, M = B.shape
    path, m = [], 0
    for t in range(n):
        hi = min(n, t + H)
        Bw = np.vstack([B[t:t + 1], Bf[t + 1:hi]])
        bw = np.concatenate([bs[t:t + 1], bsf[t + 1:hi]])
        end = 0 if hi == n else None
        win = dp_window(Bw, w[t:hi], bw, tau, m, end_mode=end)
        # el primer paso de la ventana parte del modo m: dp_window fija el modo de la hora t
        m = win[0] if B[t, win[0]] > NEG / 2 else int(np.argmax(B[t]))
        path.append(m)
    return path


def strategies(recs, topo, taus=TAUS, Ls=LS, recs_fc=None, Hs=HS, mpc_taus=(0.25,)):
    """kWh por unidad y travesía, y cambios por travesía, para cada estrategia (el control predictivo solo con
    mpc_taus, por coste de cálculo)."""
    out = {"inst": [], **{f"voraz {t}": [] for t in taus}, **{f"optima {t}": [] for t in taus},
           **{f"cambios voraz {t}": [] for t in taus}, **{f"cambios optima {t}": [] for t in taus},
           **{f"causal L{L:g} {t}": [] for t in taus for L in Ls},
           **{f"cambios causal L{L:g} {t}": [] for t in taus for L in Ls}}
    if recs_fc:
        for lev in recs_fc:
            out.update({f"mpc {lev} H{H} {t}": [] for t in mpc_taus for H in Hs})
            out.update({f"cambios mpc {lev} H{H} {t}": [] for t in mpc_taus for H in Hs})
    if topo == "T5":
        out.update({f"{r} {t}": [] for r in ("literal", "prevista") for t in (0.0,) + taus})
    for k, rec in enumerate(recs):
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
            if recs_fc and t in mpc_taus:
                for lev, rf in recs_fc.items():
                    af = rf[k]["a"]
                    if af.shape != a.shape or B.shape[1] == 1:
                        for H in Hs:
                            out[f"mpc {lev} H{H} {t}"].append(gv)
                            out[f"cambios mpc {lev} H{H} {t}"].append(gs)
                        continue
                    Bf, _, bsf = mode_matrix(af, topo)
                    for H in Hs:
                        mv, ms = path_value(B, w, bs, mpc_path(B, Bf, w, bs, bsf, t, H), t)
                        out[f"mpc {lev} H{H} {t}"].append(mv)
                        out[f"cambios mpc {lev} H{H} {t}"].append(ms)
        if topo == "T5":
            run, cosb, sav_ok = a[:, 8] > 0, a[:, 9], B[:, 1] > NEG / 2
            rules = {"literal": (~run) & sav_ok, "prevista": (~run) & (cosb < 0) & sav_ok}
            for r, mask in rules.items():
                for t in (0.0,) + taus:
                    out[f"{r} {t}"].append(path_value(B, w, bs, list(mask.astype(int)), t)[0])
    return out


def annualise(recs, values, by=None, seasons=None, n_units=None):
    """t/año de combustible (n unidades) a partir de kWh por unidad y travesía, con el consumo marginal de cada
    travesía: media por tipo de travesía por semanas de cada estación (o de cada año si by='year', solo estaciones
    del caso base). by = 0 o 1: una sola estación."""
    vals = np.asarray(values, float) * np.array([r["f_marg"] for r in recs])
    n_units = n_units or rm.N_UNITS
    seasons = rm.BASE_SEASONS if seasons is None else seasons
    names = sorted({r["c"]["name"] for r in recs})
    if by == "year":
        res = {}
        for y in sorted({r["c"]["year"] for r in recs}):
            tot = 0.0
            for s in seasons:
                for nme in names:
                    sel = [k for k, r in enumerate(recs) if r["c"]["name"] == nme and r["c"]["year"] == y
                           and r["c"]["season"] == s]
                    if sel:
                        tot += vals[sel].mean() * rm.DAYS[s] / 7.0
            res[y] = n_units * tot * rm.SFOC_ME / 1e6
        return res
    tot = 0.0
    for s in (seasons if by is None else (by,)):
        for nme in names:
            sel = [k for k, r in enumerate(recs) if r["c"]["name"] == nme and r["c"]["season"] == s]
            if sel:
                tot += vals[sel].mean() * rm.DAYS[s] / 7.0
    return n_units * tot * rm.SFOC_ME / 1e6


def r0_fuel(recs, p, seasons=None):
    """Combustible del buque sin rotores (t/año) en las estaciones dadas, con el modelo de la versión 4."""
    seasons = rm.BASE_SEASONS if seasons is None else seasons
    names = sorted({r["c"]["name"] for r in recs})
    sea = 0.0
    for s in seasons:
        for nme in names:
            sel = [r["fuel"] for r in recs if r["c"]["name"] == nme and r["c"]["season"] == s]
            if sel:
                sea += float(np.mean(sel)) * rm.DAYS[s] / 7.0
    berth = p.hotel_berth * sum(rm.PORT_HOURS_WEEK.values()) * sum(rm.DAYS[s] for s in seasons) / 7.0 * rm.SFOC_AE / 1e6
    return sea + berth


def berth_generation(fields, p, seasons=None):
    """Energía en atraque por unidad en las horas reales del horario (estaciones del caso base); t/año para las
    unidades del buque."""
    seasons = rm.BASE_SEASONS if seasons is None else seasons
    stays = [st for st in schedule(ATRAQUES) if st["season"] in seasons]
    weeks = sum(rm.DAYS[s] for s in seasons) / 7.0
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
        kwh_year = sum(np.mean(v) * weeks for v in res[k].values())
        out[k] = dict(kW_medio_unidad=kwh_year / (sum(rm.PORT_HOURS_WEEK.values()) * weeks),
                      t_anio=p.n_units * kwh_year * rm.SFOC_AE / 1e6)
    return out


def stats(d):
    v = np.array(list(d.values()), float)
    return dict(media=float(v.mean()), dt=float(v.std(ddof=1)), min=float(v.min()), max=float(v.max()),
                cv=float(v.std(ddof=1) / v.mean()) if v.mean() else float("nan"))


def run(path="era5_horario_ruta.csv.gz", p=None, taus=TAUS, quiet=False, forecast=True):
    p = p or rm.P()
    data, meta = eh.load_hourly(path)
    fields = wind_fields(data)
    cross = schedule(HORARIO)
    recs = evaluate(fields, cross, p)
    recs_fc = {lev: evaluate(fields, cross, p, forecast=lev) for lev in FORECAST} if forecast else None
    base = [r for r in recs if r["c"]["season"] in rm.BASE_SEASONS]
    res = {"travesias": len(recs), "travesias_caso_base": len(base),
           "horas_mar": float(sum(r["a"][:, 0].sum() for r in base)),
           "horas_sin_dato": int(sum(len(crossing_hours(r["c"])[2]) - len(r["a"]) for r in recs)),
           "R0_t_anio": r0_fuel(recs, p), "R0_t_anio_completo": r0_fuel(recs, p, seasons=(0, 1)),
           "tau_h": list(taus), "L_h": list(LS), "H_h": list(HS), "tau_R1": list(TAU_R1), "tau_abatible": list(TAU_TILT),
           "prevision": FORECAST}
    per = {topo: strategies(recs, topo, taus, recs_fc=recs_fc) for topo in ("R1", "T1", "T2", "T3", "T3x", "T4", "T5")}
    r1 = per["R1"]["inst"]
    res["R1"] = {"inst": annualise(recs, r1), "abr-nov": annualise(recs, r1, by=0), "dic-mar": annualise(recs, r1, by=1),
                 "año completo": annualise(recs, r1, seasons=(0, 1)), "por_anio": annualise(recs, r1, by="year")}
    res["R1"]["por_anio_estad"] = stats(res["R1"]["por_anio"])
    pen = [float((r["a"][:, 0] * r["a"][:, 10]).sum()) for r in recs]
    res["R1"]["inducida"] = annualise(recs, pen)
    # coste de arranque y parada de R1 y rotor abatible
    st_r1 = strategies(recs, "R1 conmutado", taus=TAU_R1, recs_fc=recs_fc, mpc_taus=(TAU_R1[-1],))
    res["R1"]["conmutado"] = {k: (float(np.mean(v)) if k.startswith("cambios") else annualise(recs, v))
                              for k, v in st_r1.items()}
    st_t = strategies(recs, "R1t", taus=TAU_TILT, recs_fc=recs_fc, mpc_taus=(TAU_TILT[0],))
    res["R1t"] = {k: (float(np.mean(v)) if k.startswith("cambios") else annualise(recs, v) - res["R1"]["inst"])
                  for k, v in st_t.items()}
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
    print(f"Travesías simuladas: {res['travesias']} ({res['travesias_caso_base']} de abril a noviembre); horas en el "
          f"mar {res['horas_mar']:.0f}; horas sin dato {res['horas_sin_dato']}; R0 {res['R0_t_anio']:.0f} t/año "
          f"(año completo {res['R0_t_anio_completo']:.0f})")
    r1 = res["R1"]
    st = r1["por_anio_estad"]
    print(f"R1 (cuatro rotores): {r1['inst']:.1f} t/año (abr-nov {r1['abr-nov']:.1f}, dic-mar {r1['dic-mar']:.1f}, "
          f"año completo {r1['año completo']:.1f}); inducida pagada {r1['inducida']:.1f} t/año")
    print("  por año: " + ", ".join(f"{y} {v:.0f}" for y, v in r1["por_anio"].items()) +
          f"  | media {st['media']:.1f}, dt {st['dt']:.1f} (CV {100 * st['cv']:.1f} %), {st['min']:.0f}-{st['max']:.0f}")
    c = r1["conmutado"]
    print("  R1 con arranque y parada: " + ", ".join(f"{k} {v:.1f}" for k, v in c.items()
                                                    if not k.startswith("cambios") and ("L" not in k or "L3" in k)))
    print("R1t (incremento sobre R1): " + ", ".join(f"{k} {v:.1f}" for k, v in res["R1t"].items()
                                                  if not k.startswith("cambios") and ("L" not in k or "L3" in k)))
    for topo, d in res["incrementos_mar"].items():
        print(f"{topo:4s} inst {d['inst']:7.1f} (abr-nov {d['abr-nov inst']:.1f}, dic-mar {d['dic-mar inst']:.1f})")
        for t in res["tau_h"]:
            line = (f"     tau {t}: óptima {d[f'optima {t}']:6.1f} ({d[f'cambios optima {t}']:.1f} cambios)"
                    f"  voraz {d[f'voraz {t}']:6.1f} ({d[f'cambios voraz {t}']:.1f})")
            for L in (3.0,):
                line += f"  causal L{L:g} {d[f'causal L{L:g} {t}']:6.1f} ({d[f'cambios causal L{L:g} {t}']:.1f})"
            for lev in res["prevision"]:
                for H in res["H_h"]:
                    key = f"mpc {lev} H{H} {t}"
                    if key in d:
                        line += f"  mpc {lev} H{H} {d[key]:6.1f}"
            print(line)
        s = d["por_anio inst estad"]
        print(f"     por año (inst): " + ", ".join(f"{y} {v:.1f}" for y, v in d["por_anio inst"].items()) +
              f" | media {s['media']:.1f}, dt {s['dt']:.1f}, {s['min']:.1f}-{s['max']:.1f}")
        if topo == "T5":
            print("     reglas: " + ", ".join(f"{k} {d[k]:.1f}" for k in d if k.startswith(("literal", "prevista"))))
    print("Atraque (horas reales, abril-noviembre): " + ", ".join(
        f"{k} {v['kW_medio_unidad']:.2f} kW/unidad, {v['t_anio']:.1f} t/año" for k, v in res["atraque"].items()))


if __name__ == "__main__":
    out = run(sys.argv[1] if len(sys.argv) > 1 else "era5_horario_ruta.csv.gz")
    json.dump(out, open("resultados_horarios.json", "w"), indent=1, ensure_ascii=False, default=float)
