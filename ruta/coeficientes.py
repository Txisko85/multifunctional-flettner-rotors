#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base de coeficientes del modelo de ruta (A1). Versión 18-10-2026, ampliada el 19-10-2026 con la
resistencia inducida por la fuerza lateral (e_span) y la fuerza lateral del viento sobre el buque (ayc).

Rotor girando, CL(k) y CD(k), con k = U_superficie / U_aparente:
  Conjunto A  De Marco et al. (2016), Fig. 10, URANS: AR = 6 y de/d = 2, la geometria del caso de estudio
              (digitalizado de los vectores del PDF). Coincide con los ensayos de Badalamenti (AR 5,1,
              de/d = 2; Seifert 2012, Figs. 21-24) y de Ackeret (Seifert 2012, Fig. 27). Re de modelo.
  Conjunto B  Bordogna et al. (2019), Figs. 4 y 7: ensayos a Re 3,6e5-1e6 de pared a pared (casi 2D),
              media de los puntos supercriticos, mas la resistencia inducida de un rotor de AR = 6 sobre
              cubierta (AR efectivo 12 por la reflexion en la cubierta): CD_i = CL^2 / (pi * 12).
  El modelo mezcla ambos con un peso w (w = 1: A; w = 0: B). El caso central usa w = 0,5.
Rotor parado: CD0 0,4-0,7 (Bordogna, k = 0, Re >= 3,6e5: 0,44-0,56; Ackeret, Re bajo: 0,84).
Potencia de giro: friccion del cilindro P = Cf * rho/2 * U_s^3 * pi * D * H, con Cf entre 0,002
  (Theodorsen-Regier, cilindro liso a Re ~1e7) y 0,007 (Bordogna, ensayo con perdidas de extremo y mecanicas),
  mas el disco de la placa de extremo (disco libre, CM = 0,146 Re^-0,2).
Savonius: Cp 0,14-0,24 en tunel con cangilones semicirculares (Akwa et al. 2012, Tabla 1).
  Empuje Ct 0,4-1,0: cota inferior cinematica de una maquina de arrastre, Ct >= Cp U / v_cangilon (~2,2 Cp a
  lambda ~0,8); Chauvin y Benghrib (1989) dan 0,14-0,20 a lambda 0,6-1,0 integrando seis tomas de presion a media
  altura; cota superior, la resistencia de un cuerpo romo de la misma area (supuesto).
VAWT (21-10-2026): modelo DMST de ../turbina (mapa_turbina.json) con tres factores que salen de su validación
  (validacion_dmst.py). f_cp_t 0,7-0,9: el Cp máximo medido en campo en el rotor H de 12 kW de Uppsala es 0,68-0,73
  veces el del modelo y el de túnel de Battisti et al. (2018), 0,92 veces. f_ct_t 1,0-1,1: el modelo da el 89 % del
  empuje medido en túnel cerca del máximo, que incluye el eje y los brazos del ensayo, no modelados. f_cs_t 0,6-1,0:
  el modelo sobrestima la fuerza lateral del túnel en un 35 % cerca del máximo (en un 0-80 % en todo el ensayo).
  Hasta el 19-10-2026 (versiones 1 y 2): Cp 0,20-0,30 y Cp/Ct 0,3-0,5, supuestos (RANGES_V2).
Prolongacion: CL x0,6-1,0 respecto al rotor (la placa de extremo sube con ella en el caso 1,0).
Buque: rendimiento propulsivo eta_D 0,65-0,75 y SFOC del motor principal 170-190 g/kWh.
Resistencia inducida por la fuerza lateral (19-10-2026): dR = Fy^2 / (pi e q_w (2T)^2), ala de bajo alargamiento con
  envergadura efectiva 2T (reflexion en la superficie libre). Bickert et al. (MARINE 2025) dan dR/R0 con Fy = 2 R0 a
  escala real: KVLCC2 a 15,5 kn 6,7 % (3 GDL), KCS a 20 kn 10,0 % (3+1 GDL) y 21,1 % (4 GDL, escora de 15 grados).
  Con R0 estimado por ITTC-78 (+-10 %), equivalen a e = 0,51-0,62, 0,43-0,52 y 0,20-0,25 (inferencia_e.py).
  Rango 0,4-0,8 (central 0,6): la escora pesa poco con las fuerzas laterales de este buque.
Fuerza lateral del viento sobre casco y superestructura: F_s = q_s (A_L C_Y,max) sin(beta_s), con A_L C_Y,max
  3 000-6 000 m2 (supuesto: A_L ~ 5 000 m2 para 214,5 m de eslora y ~25 m de altura media sobre la flotacion,
  C_Y,max 0,6-1,2). Hasta la version 3 solo contaba por el termino cruzado de la resistencia inducida.
Version 4 (22-10-2026), tras la revision simulada:
  - Resistencia del viento sobre el buque sin rotores (R0), explicita: longitudinal X_s = q_s (A_T C_X,max) cos(beta_s),
    con A_T C_X,max 300-600 m2 (supuesto: area frontal ~ 27,8 m de manga x ~ 27 m de altura, C_X,max 0,4-0,8, como los
    ferris de Blendermann 1994), e inducida por su propia fuerza lateral.
  - Empuje del Savonius: la cota cinematica se evalua en el centroide de la proyeccion del cangilon, r_c = (d - e)/2
    = 1,6 m: C_t,S >= C_p,S R/(lambda r_c) = 2,81 C_p,S a lambda = 0,8 (antes ~2,2 C_p,S con r ~ 2 m). El Monte Carlo
    muestrea C_t,S entre esa cota y 1,0; el valor central es el punto medio con C_p,S central (0,77).
  - Conjuntos de coeficientes del rotor como escenarios epistemicos (ESCENARIOS_COEF), no como una variable aleatoria:
    el peso w deja de muestrearse. El escenario de referencia es A: con el mismo codigo, en la ruta del Viking Grace
    (validacion_viking_grace.py), A da 141 kW de potencia propulsora equivalente al empuje del rotor de 24 m x 4 m,
    frente a 207-282 kW verificados (Norsepower 2019); la mezcla de A y B de la version 3 da 86 kW; B, 37 kW; B', 22 kW.
    B y B' quedan como cotas inferiores incompatibles con el dato a escala real.
"""
import math
import numpy as np

K_A = [0.0, 1.0, 1.5, 2.0, 3.0]
CL_A = [0.0, 1.807, 3.298, 5.451, 9.729]
CD_A = [None, 0.891, 0.869, 1.465, 3.977]          # k = 0: CD0 del rotor parado

K_B = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
CL_B = [0.0, 1.12, 2.20, 3.05, 3.85, 5.04, 5.58]
CD_B2D = [0.51, 0.77, 1.03, 1.33, 1.61, 1.87, 2.16]
AR_EFF_B = 12.0

# Versión 4: prolongación más allá de k = 3, solo para la sensibilidad con k_max = 4 (y la validación con el Viking
# Grace). Conjunto B: puntos de Bordogna et al. (2019) a Re = 3,6e5 (CD, media de las dos series de la Fig. 7).
# Conjunto A: De Marco et al. (2016) no pasa de k = 3; se prolonga con la forma relativa de los ensayos de Ackeret con
# placas de 120 mm (Seifert 2012, Fig. 27), que se satura: CL x1,08 / 1,10 / 1,12 y CD x1,23 / 1,34 / 1,45 a k = 3,5 /
# 4,0 / 4,5 respecto a k = 3.
K_B_EXT = K_B + [3.65, 4.13, 4.60]
CL_B_EXT = CL_B + [6.50, 7.31, 8.46]
CD_B2D_EXT = CD_B2D + [2.37, 2.50, 2.68]
K_A_EXT = K_A + [3.5, 4.0, 4.5]
CL_A_EXT = CL_A + [round(9.729 * r, 3) for r in (1.0822, 1.1032, 1.1244)]
CD_A_EXT = CD_A + [round(3.977 * r, 3) for r in (1.2337, 1.3405, 1.4472)]

# Corrección tridimensional del conjunto B (versión 4, linea_sustentadora.py): factor de sustentación eta_b y
# alargamiento efectivo de la resistencia inducida. Por defecto, sin corrección (1,0 y 12, como hasta la versión 3);
# escenario B' con la analogía de ala y la cubierta como plano de simetría (eta 0,83, AR_eff 12) y B'' sin reflexión
# en la cubierta (eta 0,725, AR_eff 6).
ETA_B_SCEN = {"B": (1.0, 12.0), "B'": (0.83, 12.0), "B''": (0.725, 6.0)}

# Valor central = punto medio de cada rango (el Monte Carlo muestrea los rangos de forma uniforme).
# Cota cinematica del Savonius (version 4): C_t,S >= KAPPA_S C_p,S, con la velocidad del cangilon en el centroide
# de su proyeccion: r_c = (d - e)/2 = 1,6 m, R = D_S/2 = 3,6 m, lambda = 0,8.
R_S, RC_S, LAMBDA_S = 3.6, 1.6, 0.8
# Versión 5 (revisión simulada v2, R2-1): la cota estricta usa la velocidad máxima de los cangilones, la del borde
# exterior (omega R_S): C_t,S >= C_p,S / lambda = 1,25 C_p,S. Evaluada en el centroide de la proyección del cangilón
# (versión 4) daba 2,81 C_p,S, que es una estimación con un radio efectivo supuesto, no una cota.
KAPPA_S = 1.0 / LAMBDA_S
KAPPA_S_CENTROIDE = R_S / (LAMBDA_S * RC_S)          # 2,8125: estimación de la versión 4
CT_S_MAX = 1.0

RANGES = dict(w=(0.0, 1.0), cd_mult=(0.8, 1.2), cd0=(0.4, 0.7), cf=(0.002, 0.007), cp_s=(0.14, 0.24),
              ct_s=(0.4, 1.0), f_cp_t=(0.7, 0.9), f_ct_t=(1.0, 1.1), f_cs_t=(0.6, 1.0), cl_ext=(0.6, 1.0),
              eta_d=(0.65, 0.75), sfoc_me=(170.0, 190.0), e_span=(0.4, 0.8), ayc=(3000.0, 6000.0))
CENTRAL = {k: 0.5 * (a + b) for k, (a, b) in RANGES.items()}
# version 4: resistencia frontal del viento sobre el buque y C_t,S acotado por la cota cinematica en el centroide
# y factores de la prolongación y del rotor de 30 m de la línea sustentadora (linea_sustentadora.py): con la placa en
# lo alto de la prolongación, f_X = 0,95-1,25 (1,0 si la circulación la fija el giro; hasta 1,25-1,3 con la analogía
# de ala, porque el rotor principal deja de tener extremo); f30 = 0,95-1,15.
# Versión 5 (revisión simulada v2, R2-3): f_X y f30 contrastados con la dependencia de C_L con el alargamiento de
# De Marco et al. (2016, URANS con placa de/d = 2 y 3): f30 empírico 0,89-1,13 y f_X = f30 x (1,01-1,18), el cociente
# del modelo de línea sustentadora con los dos modelos de sección (linea_sustentadora_validacion.py): f_X 0,90-1,30.
RANGES_V4 = dict(RANGES, ct_s=(KAPPA_S * RANGES["cp_s"][0], CT_S_MAX), axc=(300.0, 600.0), cl_ext=(0.90, 1.30),
                 f30=(0.90, 1.15))
CENTRAL_V4 = dict(CENTRAL, ct_s=0.5 * (KAPPA_S * CENTRAL["cp_s"] + CT_S_MAX), axc=450.0, cl_ext=1.10, f30=1.025,
                  w=1.0)
RANGES_V4_ORIG = dict(RANGES, ct_s=(KAPPA_S_CENTROIDE * RANGES["cp_s"][0], CT_S_MAX), axc=(300.0, 600.0),
                      cl_ext=(0.95, 1.25), f30=(0.95, 1.15))           # versión 4 (26-10-2026)
CENTRAL_V4_ORIG = dict(CENTRAL, ct_s=0.5 * (KAPPA_S_CENTROIDE * CENTRAL["cp_s"] + CT_S_MAX), axc=450.0, cl_ext=1.10,
                       f30=1.05, w=1.0)
# Escenarios epistemicos de los coeficientes del rotor (version 4). A: referencia (el mas cercano al dato a escala
# real); A+: A prolongado hasta k = 4 con U_s <= 46 m/s; AB: mezcla de la version 3; B: Bordogna con AR_eff 12; B':
# B corregido con la linea sustentadora (analogia de ala, cubierta como plano de simetria); B'': sin reflexion.
# Versión 5 (revisión simulada v2, R1-1 y R2-2): escenario A-VG, A con C_L multiplicado para que el empuje del modelo
# en la ruta del Viking Grace, con las celdas de ERA5 de mar, iguale 207-282 kW (centro 244,5 kW);
# validacion_viking_grace.py comprueba estos valores.
CL_MULT_VG = {"207": 1.141, "244.5": 1.25, "282": 1.359}
ESCENARIOS_COEF = {"A": dict(w=1.0), "A+": dict(w=1.0, k_max=4.0, us_max=46.0), "AB": dict(w=0.5),
                   "B": dict(w=0.0), "B'": dict(w=0.0, eta_b=ETA_B_SCEN["B'"][0], ar_b=ETA_B_SCEN["B'"][1]),
                   "B''": dict(w=0.0, eta_b=ETA_B_SCEN["B''"][0], ar_b=ETA_B_SCEN["B''"][1]),
                   "A-VG": dict(w=1.0, cl_mult=CL_MULT_VG["244.5"]),
                   "A-VG 207": dict(w=1.0, cl_mult=CL_MULT_VG["207"]),
                   "A-VG 282": dict(w=1.0, cl_mult=CL_MULT_VG["282"])}
# rangos que se muestrean en la version 4 (escenario fijado): todos menos w
RANGES_V4_MC = {k: v for k, v in RANGES_V4.items() if k != "w"}
# version 5: la resistencia inducida sale del equilibrio de deriva y guiñada (maniobra.py); en lugar de e_span se
# muestrean el término del casco c_H (0-1; el KVLCC2 de Yasukawa y Yoshimura 2015, X'_vv = -0,040 e Y'_v = -0,315,
# equivale a c_H = 2 x 0,040 / 0,315 = 0,25), los coeficientes de forma de popa C_pa y C_wa y el centro de la fuerza
# lateral del viento sobre el buque x_s (fracción de Lpp a proa de la maestra).
RANGES_V5_MC = {k: v for k, v in RANGES_V4.items() if k not in ("w", "e_span")}
RANGES_V5_MC.update(c_h=(0.0, 1.0), cpa=(0.65, 0.75), cwa=(0.80, 0.90), x_s=(0.05, 0.10))


def ct_s_range(cp_s):
    """Rango de C_t,S compatible con la cota cinematica para un C_p,S dado (version 4)."""
    return KAPPA_S * cp_s, CT_S_MAX
# rangos de las versiones anteriores, con la turbina supuesta (para reproducirlas)
RANGES_V2 = dict(w=(0.0, 1.0), cd_mult=(0.8, 1.2), cd0=(0.4, 0.7), cf=(0.002, 0.007), cp_s=(0.14, 0.24),
                 ct_s=(0.4, 1.0), cp_t=(0.20, 0.30), cpct_t=(0.3, 0.5), cl_ext=(0.6, 1.0),
                 eta_d=(0.65, 0.75), sfoc_me=(170.0, 190.0), e_span=(0.4, 0.8), ayc=(3000.0, 6000.0))   # 19-10-2026
RANGES_V1 = {k: v for k, v in RANGES_V2.items() if k not in ("e_span", "ayc")}                       # 18-10-2026


def _interp(k, ks, vs):
    return float(np.interp(k, ks, vs))


def rotor_cl_cd(k, w=0.5, cd0=0.55, cd_mult=1.0):
    """CL y CD del rotor girando a la relacion de giro k (0 <= k <= 3)."""
    k = min(max(k, 0.0), 3.0)
    cd_a = [cd0] + CD_A[1:]
    cl_a, cd_a_k = _interp(k, K_A, CL_A), _interp(k, K_A, cd_a)
    cl_b = _interp(k, K_B, CL_B)
    cd_b = _interp(k, K_B, [cd0] + CD_B2D[1:]) + cl_b ** 2 / (math.pi * AR_EFF_B)
    cl = w * cl_a + (1.0 - w) * cl_b
    cd = w * cd_a_k + (1.0 - w) * cd_b
    return cl, cd0 + (cd - cd0) * cd_mult if k > 0 else cd0


if __name__ == "__main__":
    print(" k    CL_A  CD_A   CL_B  CD_B   central CL  CD   L/D")
    for k in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        a = rotor_cl_cd(k, w=1.0); b = rotor_cl_cd(k, w=0.0); c = rotor_cl_cd(k)
        print(f"{k:3.1f}  {a[0]:5.2f} {a[1]:5.2f}  {b[0]:5.2f} {b[1]:5.2f}   {c[0]:5.2f} {c[1]:5.2f} {c[0]/c[1]:5.2f}")
