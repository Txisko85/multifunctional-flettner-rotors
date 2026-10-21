#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A1 - Modelo de ruta (nivel 1): ro-pax Galicia (Brittany Ferries), Portsmouth-Bilbao. Versión 2, 19-10-2026.

Cambios respecto a la versión del 18-10-2026, que se reproduce exactamente con P(v1=True):
1. Horario. Las duraciones de la versión 1 restaban horas locales, y España va una hora por delante del Reino
   Unido todo el año. Travesías hacia Bilbao: 36,0 y 33,0 h (antes 37,0 y 34,0); hacia Portsmouth: 32,75 y
   29,5 h (antes 31,75 y 28,5). El total semanal en el mar (131,25 h) y el atraque (36,75 h) no cambian.
2. Resistencia inducida por la fuerza lateral (deriva y timón): dR = Fy_tot^2 / (pi e q_w (2T)^2), ala de bajo
   alargamiento con envergadura efectiva 2T. Fy_tot es la fuerza lateral total del buque: las cuatro unidades
   más el viento sobre casco y superestructura, F_s = q_s (A_L C_Y,max) sin(beta_s). Solo cuenta el incremento
   respecto al buque sin rotores, cuyo efecto ya está en la curva de potencia contrastada con EU MRV. La relación
   de giro y los despliegues se eligen con ella incluida. Convención: Fy > 0 hacia sotavento; con viento aparente
   de popa del través, la componente lateral de la sustentación del rotor apunta a barlovento y reduce Fy_tot.
3. El factor de interacción eps escala la fuerza aerodinámica del rotor (sus dos componentes) con la VAWT desplegada.

Datos públicos:
- Buque: 214,5 m de eslora, 27,8 m de manga, 6,7 m de calado, 2 x MaK 12M43C de 12 600 kW (25 200 kW),
  22 kn de proyecto, 2 generadores de cola de 2 600 kW (niferry.co.uk; Wikipedia, E-Flexer).
- Horario 2026 Portsmouth-Bilbao: dos rotaciones semanales (Brittany Ferries newsroom).
- Viento: ERA5 horario 2015-2024 a 10 y 100 m (Open-Meteo, models=era5), en 10 puntos de la ruta y en los
  dos puertos; extrapolado con ley potencial a 40 m (centro del rotor) y 55 m (VAWT o prolongación).
Coeficientes: coeficientes.py. Supuestos: rendimientos, curva potencia-velocidad (ley cúbica, contrastada con
EU MRV), SFOC, cargas de hotel, altura de la cubierta, coeficientes de la VAWT y área lateral del buque.
"""
import json
import math
import os
import numpy as np
import coeficientes as cf

RHO, NU_AIR, RHO_W = 1.225, 1.48e-5, 1025.0
KN = 0.514444

# ---------------- Rotor y dispositivos ----------------
D_F, H_F, R_E = 4.0, 24.0, 4.0           # R_E: radio de la placa de extremo (De = 2D)
D_X, H_X = 3.4, 6.0                      # VAWT o prolongación
A_F, A_X = D_F * H_F, D_X * H_X
CL3, CD0_ROT, SR_NOM, US_MAX, CF_SKIN = 8.0, 2.0, 3.0, 36.0, 0.0025   # nivel 0 (solo P(nivel0=True))
CP_T, CPCT = 0.30, 0.40                  # turbina de nivel 0 (P(nivel0=True))
CT_T = CP_T / CPCT
CP_T_V2, CPCT_V2 = 0.25, 0.40            # turbina supuesta de las versiones 1 y 2 (centro de sus rangos)
UA_CUT_T = 20.0
# Turbina con modelo DMST (21-10-2026): mapas de ../turbina/mapa_turbina.json, factores de validación f_cp_t,
# f_ct_t y f_cs_t (coeficientes.py), velocidad periférica máxima y vientos nominal y de arranque del rotor H de
# 12 kW de Uppsala (Apelfröjd et al. 2016: 40 m/s y 12 m/s; arranque 3 m/s, supuesto).
MAPA_TURBINA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "turbina", "mapa_turbina.json")
V_TIP_T, U_RATED_T, U_IN_T = 40.0, 12.0, 3.0
DS_FACTOR, CP_S, CT_S, CD_STOP, UA_CUT_S = 1.8, 0.15, 0.9, 0.6, 20.0
ETA_M, ETA_MOT, ETA_CONV, ETA_GEN = 0.97, 0.95, 0.97, 0.94
N_UNITS = 4

# ---------------- Buque (público) y supuestos ----------------
P_MCR, V_DESIGN = 25200.0, 22.0
LOAD_DESIGN, N_PROP, ETA_D = 0.85, 3.0, 0.70       # supuestos
SFOC_ME, SFOC_AE, ETA_PTO = 180.0, 200.0, 0.95       # supuestos
P_HOTEL_SEA, P_HOTEL_BERTH = 2000.0, 1500.0          # kW, supuestos
T_DRAUGHT = 6.7                                      # m, calado
Z_SHIP, Z_ROT, ALPHA_SHIP = 20.0, 40.0, 0.08         # centro del área lateral, centro del rotor, exponente medio
CF_FUEL = 3.15    # t CO2/t de combustible: cociente CO2/combustible del Galicia en EU MRV 2021-2025
                  # (3,1518; anual 3,13-3,20). Anexo I del Reglamento (UE) 2015/757: HFO 3,114, LFO 3,151, MGO 3,206.

# ---------------- Ruta y horario (públicos) ----------------
NM_OPEN = {"C": 208.4, "B": 331.4}                   # mn de mar abierto: Canal y golfo
NM_APPROACH, V_APPROACH = 32.0, 10.0                 # aproximaciones a puerto (supuesto 10 kn)
CROSSINGS = [("Portsmouth-Bilbao, jueves", "S", 36.0), ("Portsmouth-Bilbao, domingo", "S", 33.0),
             ("Bilbao-Portsmouth, martes", "N", 32.75), ("Bilbao-Portsmouth, sabado", "N", 29.5)]
CROSSINGS_V1 = [("Portsmouth-Bilbao, jueves", "S", 37.0), ("Portsmouth-Bilbao, domingo", "S", 34.0),
                ("Bilbao-Portsmouth, martes", "N", 31.75), ("Bilbao-Portsmouth, sabado", "N", 28.5)]
PORT_HOURS_WEEK = {"BILBAO": 10.0, "PORTSMOUTH": 26.75}
DAYS = {0: 244.0, 1: 121.25}                         # 0: abril-noviembre, 1: diciembre-marzo


def apparent(u_w, gamma_deg, v_s):
    g = math.radians(gamma_deg)
    ux, uy = u_w * math.cos(g) + v_s, u_w * math.sin(g)
    return math.hypot(ux, uy), math.atan2(uy, ux)


def spin_power(u_a, k, cf_skin=CF_SKIN):
    """Potencia de giro (W) del rotor: fricción del cilindro y disco de la placa de extremo; k escalar o vector."""
    u_s = np.asarray(k, float) * u_a
    p_cyl = 0.5 * RHO * cf_skin * u_s ** 3 * math.pi * D_F * H_F
    w = u_s / (D_F / 2)
    re_d = np.maximum(w * R_E ** 2 / NU_AIR, 1.0)
    p_disc = 0.5 * 0.146 * re_d ** -0.2 * RHO * w ** 3 * R_E ** 5
    return p_cyl + p_disc


K_GRID = np.round(np.arange(0.25, 3.0001, 0.25), 2)


class P:
    """Parámetros que se varían en la sensibilidad y en el Monte Carlo."""
    def __init__(self, **kw):
        self.v1 = False        # True: versión del 18-10-2026 (horario sin corregir, sin resistencia inducida)
        self.nivel0 = False    # True: coeficientes de nivel 0 y relación de giro fija
        self.cd = CD0_ROT      # solo nivel 0
        self.cl_ext = 1.0
        self.kappa = 1.0 / ETA_PTO
        self.ds = DS_FACTOR
        self.eps = 0.0
        self.ua_ext_max = None
        self.dv = 0.0          # cambio de velocidad del buque (kn)
        self.stop_drag = True  # resistencia del rotor parado en R1
        self.k_opt = True      # relación de giro óptima en cada estado
        self.tz_fix = True     # horario en duraciones reales
        self.sir = True        # resistencia inducida por la fuerza lateral
        self.windage = True    # fuerza lateral del viento sobre el buque (término cruzado)
        self.eta_d = ETA_D
        self.sfoc_me = SFOC_ME
        self.turb = "dmst"     # "dmst": modelo de la turbina (21-10-2026); "supuesta": Cp y Cp/Ct supuestos
        self.cp_t, self.cpct_t = CP_T_V2, CPCT_V2      # solo con turb = "supuesta"
        self.v_tip_t = V_TIP_T
        self.corot = True      # la turbina gira con el rotor cuando lo asiste (su fuerza lateral suma a la Magnus)
        self.__dict__.update(cf.CENTRAL)
        self.__dict__.update(kw)
        if self.v1:
            self.tz_fix, self.sir, self.turb = False, False, "supuesta"
            self.__dict__.update({k: v for k, v in kw.items() if k not in ("v1", "tz_fix", "sir")})
        if self.nivel0:
            self.k_opt, self.turb = False, "supuesta"
            self.cd0, self.cf, self.cp_s, self.ct_s, self.cp_t, self.cpct_t = CD_STOP, CF_SKIN, CP_S, CT_S, CP_T, CPCT
            self.cl_ext, self.eta_d, self.sfoc_me = 1.0, ETA_D, SFOC_ME
            self.__dict__.update({k: v for k, v in kw.items() if k != "nivel0"})


def crossings(p):
    return CROSSINGS if p.tz_fix else CROSSINGS_V1


def rotor_coeffs(k, p):
    """CL y CD (vectores) del rotor a las relaciones de giro k."""
    k = np.asarray(k, float)
    if p.nivel0:
        return CL3 * np.minimum(k, SR_NOM) / SR_NOM, np.full_like(k, p.cd)
    kc = np.clip(k, 0.0, 3.0)
    cl_a = np.interp(kc, cf.K_A, cf.CL_A)
    cd_a = np.interp(kc, cf.K_A, [p.cd0] + cf.CD_A[1:])
    cl_b = np.interp(kc, cf.K_B, cf.CL_B)
    cd_b = np.interp(kc, cf.K_B, [p.cd0] + cf.CD_B2D[1:]) + cl_b ** 2 / (math.pi * cf.AR_EFF_B)
    cl = p.w * cl_a + (1.0 - p.w) * cl_b
    cd = p.w * cd_a + (1.0 - p.w) * cd_b
    return cl, p.cd0 + (cd - p.cd0) * p.cd_mult


def sir_coef(v_s, p):
    """Coeficiente c (1/N) de la resistencia inducida, dR = c Fy^2, a la velocidad v_s (m/s)."""
    if not p.sir or v_s <= 0.0:
        return 0.0
    return 1.0 / (math.pi * p.e_span * 0.5 * RHO_W * v_s ** 2 * (2.0 * T_DRAUGHT) ** 2)


def ship_side_force(u_rot, gamma, v_s, p):
    """Fuerza lateral (N, hacia sotavento) del viento sobre casco y superestructura."""
    if not (p.sir and p.windage):
        return 0.0
    u_sh, b_sh = apparent(u_rot * (Z_SHIP / Z_ROT) ** ALPHA_SHIP, gamma, v_s)
    return 0.5 * RHO * u_sh ** 2 * p.ayc * math.sin(b_sh)


class _MapaTurbina:
    """Cp, Ct y |Cs| de la turbina (modelo DMST, sin factores) en función de lambda y del viento U."""

    def __init__(self, path=MAPA_TURBINA):
        d = json.load(open(path, encoding="utf-8"))
        self.lam, self.u = np.array(d["lam"], float), np.array(d["U"], float)
        self.cp, self.ct, self.cs = (np.array(d[k], float) for k in ("cp", "ct", "cs"))
        self.ct_park = np.array(d["ct_parada"], float)
        self._cache = {}

    LAM_FINE = np.round(np.arange(0.5, 5.0001, 0.05), 2)

    def row(self, u):
        """Cp, Ct y |Cs| en la malla fina de lambda (interpolación lineal del mapa) y Ct parada, al viento u."""
        key = round(float(u), 3)
        if key not in self._cache:
            uc = min(max(float(u), self.u[0]), self.u[-1])
            j = int(np.clip(np.searchsorted(self.u, uc) - 1, 0, len(self.u) - 2))
            t = (uc - self.u[j]) / (self.u[j + 1] - self.u[j])
            rows = [a[j] * (1 - t) + a[j + 1] * t for a in (self.cp, self.ct, self.cs)]
            self._cache[key] = tuple(np.interp(self.LAM_FINE, self.lam, r) for r in rows) + \
                (float(self.ct_park[j] * (1 - t) + self.ct_park[j + 1] * t),)
        return self._cache[key]


_MT = None


def mapa_turbina():
    global _MT
    if _MT is None:
        _MT = _MapaTurbina()
    return _MT


def rated_power_t(p):
    """Potencia nominal (W, eje) de la turbina: la máxima a U_RATED_T con la velocidad periférica permitida."""
    mt = mapa_turbina()
    cp, _, _, _ = mt.row(U_RATED_T)
    ok = mt.LAM_FINE <= p.v_tip_t / U_RATED_T + 1e-9
    return p.f_cp_t * float(cp[ok].max()) * 0.5 * RHO * U_RATED_T ** 3 * A_X


def vawt_states(u, p):
    """Estados posibles de la turbina con viento u (m/s) a su altura: potencia en el eje (W), empuje según el viento
    (N) y módulo de la fuerza lateral (N). El primero es la turbina parada; los demás, en marcha a cada lambda de la
    malla fina (pasos de 0,05) permitida por la velocidad periférica máxima, la potencia nominal y los vientos de
    arranque y parada."""
    mt = mapa_turbina()
    q = 0.5 * RHO * u * u
    cp, ct, cs, ct_park = mt.row(u)
    pw, th, sf = [0.0], [ct_park * q * A_X], [0.0]
    if U_IN_T <= u <= UA_CUT_T:
        lam = mt.LAM_FINE
        l_max = min(p.v_tip_t / u, lam[-1])
        p_r = rated_power_t(p)
        pm = p.f_cp_t * cp * q * u * A_X
        # candidatos: la malla hasta l_max, l_max y los cruces con la potencia nominal (control continuo de lambda)
        cand = list(lam[lam <= l_max + 1e-9]) + [l_max]
        for i in range(len(lam) - 1):
            if lam[i] < l_max and (pm[i] - p_r) * (pm[i + 1] - p_r) < 0.0:
                lc = lam[i] + (p_r - pm[i]) / (pm[i + 1] - pm[i]) * (lam[i + 1] - lam[i])
                if lc <= l_max:
                    cand.append(lc)
        cand = np.array(cand)
        pmc = p.f_cp_t * np.interp(cand, lam, cp) * q * u * A_X
        ok = (pmc > 0.0) & (pmc <= p_r * (1.0 + 1e-6))
        pw += list(pmc[ok])
        th += list(p.f_ct_t * np.interp(cand[ok], lam, ct) * q * A_X)
        sf += list(p.f_cs_t * np.interp(cand[ok], lam, cs) * q * A_X)
    return np.array(pw), np.array(th), np.array(sf)


def vawt_power_curve(u, p):
    """Potencia eléctrica (kW) de la turbina generando (en puerto), para un vector de vientos u (m/s)."""
    if p.turb != "dmst":
        return power_curve(A_X, p.cp_t, np.asarray(u, float))
    return np.array([float(vawt_states(float(x), p)[0].max()) for x in np.atleast_1d(u)]) * ETA_GEN * ETA_CONV / 1e3


def bin_terms(u_rot, ratio, gamma, v_s, p):
    """Términos por unidad (kW equivalentes del motor principal) en un estado de viento, con la resistencia
    inducida de la fuerza lateral total ya descontada."""
    u_a, b = apparent(u_rot, gamma, v_s)
    out = dict(rot=0.0, rot_fx=0.0, ext=0.0, vawt_opt=0.0, vawt_raw=0.0, sav=None, sav3=None,
               running=False, cosb=1.0, k=0.0, pen_r1=0.0, fy_r1=0.0, f_s=0.0,
               b_stop=0.0, ext_raw=None, vawt_ok=False)
    if u_a < 0.5:
        return out
    q = 0.5 * RHO * u_a ** 2
    sb, cb = math.sin(b), math.cos(b)
    out["cosb"] = cb
    eta_d = p.eta_d
    # resistencia inducida: penalización por unidad (kW) para una fuerza lateral por unidad fy (N)
    c_sir = sir_coef(v_s, p)
    f_s = ship_side_force(u_rot, gamma, v_s, p)
    out["f_s"] = f_s
    k_pen = c_sir * v_s / eta_d / 1e3 / N_UNITS

    def pen(fy_unit):
        tot = f_s + N_UNITS * np.asarray(fy_unit, float)
        return k_pen * (tot ** 2 - f_s ** 2)

    k_max = min(SR_NOM, US_MAX / u_a)
    if p.k_opt:
        ks = K_GRID[K_GRID <= k_max + 1e-9]
        if ks.size == 0 or ks[-1] < k_max - 1e-9:
            ks = np.append(ks, k_max)
    else:
        ks = np.array([k_max])
    cl, cd = rotor_coeffs(ks, p)
    fx = q * A_F * (cl * sb - cd * cb)
    fy = q * A_F * (cl * cb + cd * sb)
    p_spin = spin_power(u_a, ks, p.cf)
    br_all = (fx * v_s / eta_d - p.kappa * p_spin / (ETA_MOT * ETA_CONV)) / 1e3
    run_all = br_all - pen(fy)
    i = int(np.argmax(run_all))
    br = float(run_all[i])
    # rotor parado: resistencia del cilindro (cuenta si p.stop_drag)
    b_stop = -p.cd0 * q * A_F * cb * v_s / eta_d / 1e3
    fy_stop = p.cd0 * q * A_F * sb
    b_ref, fy_ref = (b_stop, fy_stop) if p.stop_drag else (0.0, 0.0)
    ref = b_ref - float(pen(fy_ref))
    running = br > ref
    out["running"] = running
    out["rot"] = max(br, ref)
    out["b_stop"] = ref                       # rotor parado (con su resistencia inducida), para el coste de conmutación
    fy_r1 = float(fy[i]) if running else fy_ref
    out["fy_r1"] = fy_r1
    out["pen_r1"] = float(pen(fy_r1))
    if running:
        out["k"] = float(ks[i])
        out["rot_fx"] = float(fx[i]) * v_s / eta_d / 1e3
    # extremo superior: viento a 55 m
    u_at, bt_ang = apparent(u_rot * ratio, gamma, v_s)
    qt = 0.5 * RHO * u_at ** 2
    sbt, cbt = math.sin(bt_ang), math.cos(bt_ang)
    # prolongación solidaria: misma velocidad angular que el rotor; k elegido para el conjunto
    if u_at > 0.5 and (p.ua_ext_max is None or u_a <= p.ua_ext_max) and (running or not p.nivel0):
        u_sx = ks * u_a * D_X / D_F
        cl_x, cd_x = rotor_coeffs(u_sx / u_at, p)
        fx_x = qt * A_X * (p.cl_ext * cl_x * sbt - cd_x * cbt)
        fy_x = qt * A_X * (p.cl_ext * cl_x * cbt + cd_x * sbt)
        p_x = 0.5 * RHO * p.cf * u_sx ** 3 * math.pi * D_X * H_X
        tot = br_all + (fx_x * v_s / eta_d - p.kappa * p_x / (ETA_MOT * ETA_CONV)) / 1e3 - pen(fy + fy_x)
        out["ext_raw"] = float(tot.max()) - out["rot"]          # incremento sin recortar (modelo horario)
        out["ext"] = max(out["ext_raw"], 0.0)
    # VAWT en el extremo
    if p.turb != "dmst" and u_at <= UA_CUT_T:          # versiones 1 y 2: Cp y Cp/Ct supuestos, sin turbina parada
        p_t = qt * u_at * A_X * p.cp_t
        t_t = qt * A_X * p.cp_t / p.cpct_t
        fy_t = t_t * sbt                       # fuerza lateral propia de la turbina (C_s) no modelada
        if running:   # asistencia mecánica hasta cubrir la potencia de giro; el exceso, por el generador
            ps = float(p_spin[i])
            mech = min(ETA_M * p_t, ps)
            saved_e = mech / (ETA_MOT * ETA_CONV) + max(ETA_M * p_t - ps, 0.0) * ETA_GEN * ETA_CONV
            inter = p.eps * float(fx[i]) * v_s / eta_d
            fy_inter = p.eps * float(fy[i])
        else:
            saved_e = p_t * ETA_GEN * ETA_CONV
            inter, fy_inter = 0.0, 0.0
        bt = (-t_t * cbt * v_s / eta_d + p.kappa * saved_e + inter) / 1e3 \
            - float(pen(fy_r1 + fy_t + fy_inter) - pen(fy_r1))
        out["vawt_raw"], out["vawt_opt"] = bt, max(bt, 0.0)
        out["vawt_ok"] = True                                   # por encima de la parada no se despliega
    elif p.turb == "dmst" and u_at > 0.5:
        # modelo DMST: la turbina elige el estado (parada o lambda) que más ahorra; su empuje sigue al viento aparente
        # y su fuerza lateral es perpendicular, en el sentido de la sustentación del rotor si gira con él
        pw, th, sf = vawt_states(u_at, p)
        if running:   # asistencia mecánica hasta cubrir la potencia de giro; el exceso, por el generador
            ps = float(p_spin[i])
            saved_e = np.minimum(ETA_M * pw, ps) / (ETA_MOT * ETA_CONV) \
                + np.maximum(ETA_M * pw - ps, 0.0) * ETA_GEN * ETA_CONV
            inter = p.eps * float(fx[i]) * v_s / eta_d
            fy_inter = p.eps * float(fy[i])
            signs = (1.0,) if p.corot else (-1.0,)
        else:         # rotor parado: la turbina arranca en el sentido que más conviene
            saved_e = pw * ETA_GEN * ETA_CONV
            inter, fy_inter = 0.0, 0.0
            signs = (1.0, -1.0)
        best = -np.inf
        for sg in signs:
            fx_t = sg * sf * sbt - th * cbt
            fy_t = sg * sf * cbt + th * sbt
            bt = (fx_t * v_s / eta_d + p.kappa * saved_e + inter) / 1e3 - (pen(fy_r1 + fy_t + fy_inter) - pen(fy_r1))
            best = max(best, float(np.max(bt)))
        out["vawt_raw"], out["vawt_opt"] = best, max(best, 0.0)
        out["vawt_ok"] = u_at <= UA_CUT_T                       # por encima de la parada no se despliega (T4)
    # rotor abierto en Savonius (a la altura del rotor): beneficio absoluto, con el mismo criterio que R1
    # (si no se cuenta la resistencia del rotor parado en R1, se mide respecto al rotor parado con resistencia)
    if u_a <= UA_CUT_S:
        for key, area_c in (("sav", A_F), ("sav3", A_F + A_X)):
            a_s = p.ds * area_c
            p_e = ETA_GEN * ETA_CONV * p.cp_s * q * u_a * a_s
            s_abs = (-p.ct_s * a_s * q * cb * v_s / eta_d + p.kappa * p_e) / 1e3
            fy_s = p.ct_s * a_s * q * sb       # sustentación del Savonius girando no modelada
            out[key] = (s_abs - float(pen(fy_s))) - (b_stop - float(pen(fy_stop))) + ref
    return out


def topologies(t):
    """Beneficio incremental por unidad sobre R1 (kW) navegando, por topología; t = bin_terms.
    Reglas de la conversión en Savonius:
      literal  : se abre siempre que el rotor no compensa (texto de la patente, incluye viento de proa);
      inventor : se abre solo si el rotor no compensa y el viento aparente empuja (de popa del través);
      popa     : se abre siempre que el viento aparente empuja, aunque el rotor compense."""
    r1, be, run = t["rot"], t["ext"], t["running"]
    s, s3 = t["sav"], t["sav3"]
    push = t["cosb"] < 0.0
    sav_opt = max(0.0, s - r1) if s is not None else 0.0
    sav_lit = (s - r1) if (s is not None and not run) else 0.0
    sav_inv = (s - r1) if (s is not None and not run and push) else 0.0
    sav_popa = (s - r1) if (s is not None and push) else 0.0
    s_best = max(v for v in (s, s3) if v is not None) if s is not None else None
    sav3_opt = max(0.0, s_best - (r1 + be)) if s_best is not None else 0.0
    sav3_lit = (s3 - r1) if (s3 is not None and not run) else 0.0
    return {
        "R1": r1,
        "T1": t["vawt_raw"],
        "T2": t["vawt_raw"] + sav_opt,
        "T2 literal": t["vawt_raw"] + sav_lit,
        "T2 inventor": t["vawt_raw"] + sav_inv,
        "T3": be + sav3_opt,
        "T3 literal": be + sav3_lit,
        "T3 inventor": be + sav_inv,
        "T3x": be,
        "T4": t["vawt_opt"],
        "T5": sav_opt,
        "T5 literal": sav_lit,
        "T5 inventor": sav_inv,
        "T5 popa": sav_popa,
        "_pen R1": t["pen_r1"],          # diagnóstico: resistencia inducida que paga R1 (kW por unidad)
        "_R1 marcha": 1.0 if run else 0.0,
    }


def load_wind(path="era5_ruta.json"):
    d = json.load(open(path, encoding="utf-8"))
    nu, ng = d["meta"]["NU"], d["meta"]["NG"]

    def ffill(vals):
        out, last = [], vals[0]
        for v in vals:
            last = v if v is not None else last
            out.append(last)
        return np.array(out, float)
    groups = {}
    for g, gd in d["groups"].items():
        groups[g] = dict(F=[np.array(gd["F"][s], float).reshape(nu, ng) / 1e6 for s in (0, 1)],
                         ratio=[ffill(gd["ratio"][s]) for s in (0, 1)])
    return d, groups, nu, ng


def crossing_speed(hours, p):
    t_open = hours - NM_APPROACH / V_APPROACH
    return 540.0 / t_open + p.dv, t_open


def mean_benefits(groups, nu, ng, season, direction, v_kn, p):
    """Media temporal por unidad (kW) de cada topología en una travesía."""
    v = v_kn * KN
    acc, w_tot = {}, sum(NM_OPEN.values())
    for g, gd in groups.items():
        F, ratio = gd["F"][season], gd["ratio"][season]
        wg = NM_OPEN[g] / w_tot
        for iu in range(nu):
            u = iu + 0.5
            for ig in range(ng):
                f = F[iu, ig]
                if f == 0:
                    continue
                gam = (ig + 0.5) * 15.0
                if direction == "N":
                    gam = 180.0 - gam
                tt = topologies(bin_terms(u, ratio[iu], gam, v, p))
                for key, val in tt.items():
                    acc[key] = acc.get(key, 0.0) + wg * f * val
    return acc


def power_curve(area, cp, u):
    """Potencia eléctrica (kW) de un generador eólico: arranque 3 m/s, nominal 12 m/s, parada 20 m/s."""
    p_r = 0.5 * RHO * area * cp * 12.0 ** 3
    pw = np.where((u >= 3.0) & (u < 20.0), np.minimum(0.5 * RHO * area * cp * u ** 3, p_r), 0.0)
    return pw * ETA_GEN * ETA_CONV / 1e3


def port_power(d, nu, p=None):
    """Potencia eléctrica media por unidad en puerto (kW): VAWT a 55 m; Savonius a 40 m."""
    p = p or P()
    res = {}
    u = np.arange(nu) + 0.5
    for port, pd in d["ports"].items():
        res[port] = []
        for s in (0, 1):
            hr, ht, n = np.array(pd["hr"][s], float), np.array(pd["ht"][s], float), pd["ns"][s]
            vawt = float((vawt_power_curve(u, p) * ht).sum() / n)
            sav = float((power_curve(p.ds * A_F, p.cp_s, u) * hr).sum() / n)
            sav3 = float((power_curve(p.ds * (A_F + A_X), p.cp_s, u) * hr).sum() / n)
            res[port].append(dict(vawt=vawt, sav=sav, sav3=sav3))
    return res


def annual(d, groups, nu, ng, p=None, seasons=(0, 1)):
    p = p or P()
    sea = {}
    fuel_me = 0.0
    hours_sea = 0.0
    for s in seasons:
        weeks = DAYS[s] / 7.0
        for name, direction, hours in crossings(p):
            v_kn, t_open = crossing_speed(hours, p)
            mb = mean_benefits(groups, nu, ng, s, direction, v_kn, p)
            for key, val in mb.items():
                sea[key] = sea.get(key, 0.0) + N_UNITS * val * t_open * weeks * p.sfoc_me / 1e6
            p_b0 = LOAD_DESIGN * P_MCR * (v_kn / V_DESIGN) ** N_PROP
            p_app = LOAD_DESIGN * P_MCR * (V_APPROACH / V_DESIGN) ** N_PROP
            fuel_me += (p_b0 * t_open + p_app * NM_APPROACH / V_APPROACH) * weeks * p.sfoc_me / 1e6
            hours_sea += hours * weeks
    pp = port_power(d, nu, p)
    port = {"T1": 0.0, "T2": 0.0, "T3": 0.0, "T4": 0.0, "T5": 0.0}
    hours_port = 0.0
    for s in seasons:
        weeks = DAYS[s] / 7.0
        for prt, h in PORT_HOURS_WEEK.items():
            e = N_UNITS * h * weeks * SFOC_AE / 1e6
            port["T1"] += e * pp[prt][s]["vawt"]
            port["T4"] += e * pp[prt][s]["vawt"]
            port["T2"] += e * (pp[prt][s]["vawt"] + pp[prt][s]["sav"])
            port["T5"] += e * pp[prt][s]["sav"]
            port["T3"] += e * pp[prt][s]["sav3"]
            hours_port += h * weeks
    fuel_hotel = (P_HOTEL_SEA * hours_sea * p.sfoc_me / ETA_PTO + P_HOTEL_BERTH * hours_port * SFOC_AE) / 1e6
    return dict(sea=sea, port=port, fuel_total=fuel_me + fuel_hotel, fuel_me=fuel_me,
                hours_sea=hours_sea, hours_port=hours_port, port_kw=pp)


def table(res):
    rows = []
    for key in ("T3", "T3x", "T2", "T5", "T4", "T1"):
        s = res["sea"][key]
        pt = res["port"].get(key.replace("x", ""), 0.0) if key != "T3x" else 0.0
        rows.append((key, s, pt, s + pt, 100 * (s + pt) / res["fuel_total"]))
    return rows


def break_even(d, groups, nu, ng, base, lo=0.0, hi=0.6, **kw):
    """Factor de interacción eps* con el que T4 iguala a T3x navegando (mismos parámetros kw)."""
    target = base["sea"]["T3x"]
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if annual(d, groups, nu, ng, P(eps=mid, **kw))["sea"]["T4"] < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    d, groups, nu, ng = load_wind()
    for label, p in (("versión 1 (18-10-2026)", P(v1=True)), ("versión 2", P())):
        base = annual(d, groups, nu, ng, p)
        print(f"== {label}")
        print(f"Combustible R0 (año completo): {base['fuel_total']:.0f} t; motor principal {base['fuel_me']:.0f} t")
        print(f"R1, 4 rotores: {base['sea']['R1']:.1f} t/año ({100 * base['sea']['R1'] / base['fuel_total']:.2f} %); "
              f"resistencia inducida pagada por R1 {base['sea']['_pen R1']:.1f} t/año")
        for r in table(base):
            print(f"  {r[0]:4s} mar {r[1]:8.1f}  puerto {r[2]:6.1f}  total {r[3]:8.1f}  ({r[4]:.2f} %)")
        for key in ("T5 literal", "T5 inventor", "T5 popa", "T2 literal", "T2 inventor", "T3 literal", "T3 inventor"):
            print(f"  {key:12s} mar {base['sea'][key]:8.1f}")
