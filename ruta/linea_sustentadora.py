#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Línea sustentadora de rotores Flettner sobre cubierta (A1, versión 4, 22-10-2026), para lo que el modelo de ruta
no resuelve: el factor de sustentación de la prolongación (f_X) y del rotor de 30 m (f30), la corrección
tridimensional del conjunto B, la velocidad inducida por el rotor sobre la turbina y el orden de magnitud de la
interacción eps.

Modelo: línea sustentadora discreta (herraduras de Prandtl sobre el eje del rotor, con estelas rectas aguas abajo),
con imágenes en la cubierta (factor m_d: 1 = plano de simetría perfecto, AR efectivo 2H/D; 0 = sin reflexión,
AR efectivo H/D) y en la placa de extremo (factor m_p). Dos modelos de sección:
  "ala":  analogía con un ala de pendiente 2 pi, C_l = 2 pi (alpha_g - alpha_i) con alpha_g = C_l,2D(k)/(2 pi): la
          corrección clásica de alargamiento finito. Es el modelo de los factores: maneja los extremos libres.
  "giro": la circulación la fija la relación de giro, Gamma = 1/2 U_ef D C_l(k_ef), con k_ef = omega R / U_ef: un
          cilindro es axisimétrico y la velocidad inducida solo gira la corriente y cambia su módulo, así que la
          envergadura finita añade resistencia inducida pero apenas quita sustentación. No reproduce la pérdida en
          un extremo libre (la circulación no se anula allí), así que solo se usa con cubierta y placa como
          espejos perfectos, como cota superior (sin pérdida).
C_l,2D(k) es el del conjunto B (Bordogna et al. 2019, ensayo de pared a pared, casi bidimensional; puntos hasta
k = 4,6). Corriente uniforme: la cortadura del viento ya está en el modelo de ruta (viento a 40 y 55 m).
Uso: python3 linea_sustentadora.py   ->  linea_sustentadora.json
"""
import json
import math
import os
import numpy as np
import coeficientes as cf

HERE = os.path.dirname(os.path.abspath(__file__))
D_F, H_F, D_X, H_X, D_E = 4.0, 24.0, 3.4, 6.0, 8.0
Z_PLATE_GAP = 0.35            # separación entre la placa (en lo alto del rotor) y la turbina
L_WAKE = 1.0e5                # longitud de las estelas (m), "infinitas"

# C_l,2D(k) del conjunto B, prolongado hasta k = 4,6 con los puntos de Re = 3,6e5 de Bordogna et al. (2019)
K_2D = np.array(cf.K_B + [3.65, 4.13, 4.60])
CL_2D = np.array(cf.CL_B + [6.50, 7.31, 8.46])


def cl2d(k):
    return np.interp(k, K_2D, CL_2D)


def _segment(p, a, b, gam):
    """Velocidad inducida (vector) en los puntos p (n x 3) por un segmento de vórtice de a a b con circulación gam."""
    r0 = b - a
    r1 = p - a
    r2 = p - b
    cr = np.cross(r1, r2)
    cr2 = np.einsum("ij,ij->i", cr, cr)
    n1 = np.linalg.norm(r1, axis=1)
    n2 = np.linalg.norm(r2, axis=1)
    ok = (cr2 > 1e-10) & (n1 > 1e-9) & (n2 > 1e-9)
    k = np.zeros(len(p))
    k[ok] = gam / (4 * math.pi * cr2[ok]) * (np.einsum("ij,j->i", r1[ok], r0) / n1[ok] - np.einsum("ij,j->i", r2[ok], r0) / n2[ok])
    return cr * k[:, None]


def horseshoe(p, z0, z1, gam=1.0):
    """Herradura con el tramo ligado en el eje (x = 0) de z0 a z1 y estelas hacia +x."""
    A, B = np.array([0.0, 0.0, z0]), np.array([0.0, 0.0, z1])
    Af, Bf = np.array([L_WAKE, 0.0, z0]), np.array([L_WAKE, 0.0, z1])
    return _segment(p, Af, A, gam) + _segment(p, A, B, gam) + _segment(p, B, Bf, gam)


class Rotor:
    """Rotor escalonado: tramos (z_inicio, z_fin, diámetro); placa en z_plate (o ninguna)."""

    def __init__(self, tramos, z_plate=None, m_d=1.0, m_p=0.5, n=80):
        self.tramos = tramos
        z_top = tramos[-1][1]
        # paneles con espaciado coseno en cada tramo
        zs, ds = [], []
        for (za, zb, d) in tramos:
            m = max(8, int(round(n * (zb - za) / z_top)))
            e = za + (zb - za) * 0.5 * (1 - np.cos(np.linspace(0, math.pi, m + 1)))
            zs.append(e)
            ds += [d] * m
        edges = np.unique(np.concatenate(zs))
        self.e0, self.e1 = edges[:-1], edges[1:]
        self.zc = 0.5 * (self.e0 + self.e1)
        self.dz = self.e1 - self.e0
        self.D = np.array([next(d for (za, zb, d) in tramos if za - 1e-9 <= z <= zb + 1e-9) for z in self.zc])
        self.z_plate, self.m_d, self.m_p = z_plate, m_d, m_p
        self.A = self._influence(np.column_stack([np.zeros_like(self.zc), np.zeros_like(self.zc), self.zc]))

    def _influence(self, pts):
        """Matriz (puntos x paneles) de la velocidad lateral (y) inducida por circulación unitaria, con imágenes."""
        A = np.zeros((len(pts), len(self.zc)))
        for j, (z0, z1) in enumerate(zip(self.e0, self.e1)):
            v = horseshoe(pts, z0, z1)
            if self.m_d:
                v += self.m_d * horseshoe(pts, -z1, -z0)
            if self.z_plate is not None and self.m_p:
                zp = self.z_plate
                v += self.m_p * horseshoe(pts, 2 * zp - z1, 2 * zp - z0)
            A[:, j] = v[:, 1]
        return A

    def solve(self, k_ref, modelo="giro", u=1.0, v_ext=None, it=200):
        """Circulación (m2/s) en cada panel con viento u (escalar o por panel) y velocidad angular tal que el tramo
        inferior va a k_ref con viento unitario; v_ext: velocidad lateral externa (p. ej., la inducida por la
        turbina). Convención: con Gamma > 0 a lo largo de +z y viento hacia +x, la sustentación apunta a -y, así
        que una velocidad inducida v > 0 se opone a ella (alpha_i = v/u). Devuelve sustentación y resistencia
        inducida por unidad de longitud (divididas por rho)."""
        u = np.broadcast_to(np.asarray(u, float), self.zc.shape).copy()
        v_ext = np.zeros_like(self.zc) if v_ext is None else np.asarray(v_ext, float)
        omega = k_ref / (self.tramos[0][2] / 2)
        k = omega * (self.D / 2) / u
        if modelo == "ala":
            alpha_g = cl2d(k) / (2 * math.pi)
            # Gamma = pi u D (alpha_g - alpha_i), alpha_i = (A Gamma + v_ext) / u
            M = np.eye(len(k)) + math.pi * np.diag(self.D) @ self.A
            gam = np.linalg.solve(M, math.pi * self.D * (u * alpha_g - v_ext))
        else:
            gam = 0.5 * u * self.D * cl2d(k)
            for _ in range(it):
                v = self.A @ gam + v_ext
                u_ef = np.sqrt(u ** 2 + v ** 2)
                new = 0.5 * u_ef * self.D * cl2d(omega * (self.D / 2) / u_ef)
                if np.max(np.abs(new - gam)) < 1e-12:
                    gam = new
                    break
                gam = 0.5 * gam + 0.5 * new
        v = self.A @ gam + v_ext
        # Kutta-Joukowski con la corriente local (u, v): F/rho = (u x + v y) x Gamma z = Gamma (-u y + v x)
        lift = u * gam                               # / rho: sustentación por unidad de longitud (hacia -y)
        drag_i = v * gam                             # / rho: resistencia inducida por unidad de longitud
        drag_p = 0.5 * u ** 2 * self.D * np.interp(k, cf.K_B, cf.CD_B2D)   # / rho: resistencia de perfil
        return dict(gam=gam, v=v, k=k, u=u, lift=lift, drag_i=drag_i, drag_p=drag_p, omega=omega)

    def fx(self, sol, beta_deg):
        """Fuerza longitudinal (/rho) de todo el rotor con viento aparente a beta del rumbo."""
        b = math.radians(beta_deg)
        return float(((sol["lift"] * math.sin(b) - (sol["drag_i"] + sol["drag_p"]) * math.cos(b)) * self.dz).sum())

    def tramo(self, sol, za, zb):
        sel = (self.zc >= za - 1e-9) & (self.zc <= zb + 1e-9)
        return float((sol["lift"][sel] * self.dz[sel]).sum()), float((sol["drag_i"][sel] * self.dz[sel]).sum())

    def lift_2d(self, sol, za, zb, u=1.0):
        sel = (self.zc >= za - 1e-9) & (self.zc <= zb + 1e-9)
        return float((0.5 * u ** 2 * self.D[sel] * cl2d(sol["k"][sel]) * self.dz[sel]).sum())


def configs(m_d, m_p):
    zt = H_F + Z_PLATE_GAP
    return {
        "C0": Rotor([(0.0, H_F, D_F)], z_plate=H_F, m_d=m_d, m_p=m_p),
        "C2": Rotor([(0.0, H_F, D_F), (H_F, H_F + H_X, D_X)], z_plate=H_F + H_X, m_d=m_d, m_p=m_p),
        "C2 placa abajo": Rotor([(0.0, H_F, D_F), (H_F, H_F + H_X, D_X)], z_plate=None, m_d=m_d, m_p=0.0),
        "R1-30": Rotor([(0.0, H_F + H_X, D_F)], z_plate=H_F + H_X, m_d=m_d, m_p=m_p),
    }


def factors(m_d, m_p, modelo, ks=(1.5, 2.0, 2.5, 3.0)):
    """f_X, f30 y eficiencia de C0 (sustentación 3D / 2D) con el criterio del modelo de ruta (ec. 5): el rotor
    principal conserva la sustentación de C0 y el factor recoge todo lo demás."""
    cfg = configs(m_d, m_p)
    out = {"eta0": [], "f_X": [], "f_X placa abajo": [], "f30": [], "cdi0": [], "cdi_C2_extra": []}
    for kr in ks:
        s0 = cfg["C0"].solve(kr, modelo)
        l0, d0 = cfg["C0"].tramo(s0, 0, H_F)
        l0_2d = cfg["C0"].lift_2d(s0, 0, H_F)
        eta0 = l0 / l0_2d
        out["eta0"].append(eta0)
        out["cdi0"].append(d0 / (0.5 * D_F * H_F))
        for key, name, top_d in (("f_X", "C2", D_X), ("f_X placa abajo", "C2 placa abajo", D_X), ("f30", "R1-30", D_F)):
            r = cfg[name]
            s = r.solve(kr, modelo)
            lt, dt = r.tramo(s, 0, H_F + H_X)
            l_ext_2d = r.lift_2d(s, H_F, H_F + H_X)
            out[key].append((lt - l0) / (eta0 * l_ext_2d))
            if key == "f_X":
                # resistencia inducida de C2 frente a la suma que supone la ec. (5): la de C0 más la del tramo
                # superior con el mismo coeficiente inducido que C0 a su sustentación
                cl_ext = (lt - l0) / (0.5 * D_X * H_X)
                cdi_ext_ec5 = out["cdi0"][-1] * (cl_ext / (l0 / (0.5 * D_F * H_F))) ** 2
                out["cdi_C2_extra"].append((dt - d0 - cdi_ext_ec5 * 0.5 * D_X * H_X) / (0.5 * D_F * H_F))
    return {k: [float(x) for x in v] for k, v in out.items()} | {"k": list(ks)}


def _field(r, s, pts):
    """Velocidad (vx, vy) inducida en pts por el sistema de vórtices del rotor resuelto, con sus imágenes."""
    vx = np.zeros(len(pts))
    vy = np.zeros(len(pts))
    for j, (z0, z1) in enumerate(zip(r.e0, r.e1)):
        v = horseshoe(pts, z0, z1, s["gam"][j])
        if r.m_d:
            v += r.m_d * horseshoe(pts, -z1, -z0, s["gam"][j])
        if r.z_plate is not None and r.m_p:
            zp = r.z_plate
            v += r.m_p * horseshoe(pts, 2 * zp - z1, 2 * zp - z0, s["gam"][j])
        vx += v[:, 0]
        vy += v[:, 1]
    return vx, vy


def turbine_inflow(m_d, m_p, modelo, k_ref=2.5):
    """Corriente inducida por el rotor (con su placa) en la turbina, de 24,35 a 30,35 m sobre la cubierta y radio
    1,7 m, con viento unitario: ángulo medio (grados; positivo hacia el lado de la sustentación del rotor) y módulo
    medio relativo, en el eje y en cuatro puntos del borde del disco."""
    r = configs(m_d, m_p)["C0"]
    s = r.solve(k_ref, modelo)
    zs = np.linspace(H_F + Z_PLATE_GAP, H_F + Z_PLATE_GAP + H_X, 25)
    out = {}
    for name, (x, y) in {"eje": (0.0, 0.0), "barlovento": (-1.7, 0.0), "sotavento": (1.7, 0.0),
                         "lado sustentación": (0.0, -1.7), "lado opuesto": (0.0, 1.7)}.items():
        pts = np.column_stack([np.full_like(zs, x), np.full_like(zs, y), zs])
        vx, vy = _field(r, s, pts)
        out[name] = dict(angulo_medio_grados=float(np.degrees(np.mean(np.arctan2(-vy, 1.0 + vx)))),
                         velocidad_media_rel=float(np.mean(np.hypot(1.0 + vx, vy))))
    return out


def turbine_on_rotor(m_d, m_p, modelo, k_ref=2.5, cs=0.19, ct=1.0, betas=(60.0, 90.0, 120.0)):
    """Orden de magnitud de eps, el cambio relativo del empuje del rotor con la turbina desplegada sobre la placa:
    (a) la fuerza lateral de la turbina (C_s, del lado de la sustentación del rotor) como una herradura de 24,35 a
        30,35 m, con circulación por unidad de velocidad Gamma_T = C_s D_X / 2, induce velocidad lateral en el rotor;
    (b) el empuje de la turbina (C_t) como disco actuador de radio D_X/2 frena el viento en lo alto del rotor: en el
        eje del disco la inducción axial es a (1 - x/(R^2 + x^2)^0,5) a la distancia x del disco, con
        a = (1 - (1 - C_t)^0,5)/2 (teoría de cantidad de movimiento).
    Devuelve eps para (a), (b) y ambos, con viento aparente a 60, 90 y 120 grados."""
    r = configs(m_d, m_p)["C0"]
    s0 = r.solve(k_ref, modelo)
    pts = np.column_stack([np.zeros_like(r.zc), np.zeros_like(r.zc), r.zc])
    zt0, zt1 = H_F + Z_PLATE_GAP, H_F + Z_PLATE_GAP + H_X
    v_t = horseshoe(pts, zt0, zt1, 0.5 * cs * D_X)[:, 1]
    a = 0.5 * (1 - math.sqrt(max(0.0, 1 - min(ct, 0.96))))
    x = (zt0 + zt1) / 2 - r.zc                     # distancia al centro del disco, medida a lo largo de su eje
    rt = D_X / 2
    u_b = 1.0 - a * (1 - x / np.sqrt(rt ** 2 + x ** 2))
    sa = r.solve(k_ref, modelo, v_ext=v_t)
    sb = r.solve(k_ref, modelo, u=u_b)
    sab = r.solve(k_ref, modelo, u=u_b, v_ext=v_t)
    out = dict(induccion_axial=a, gamma_T_rel=float(0.5 * cs * D_X / np.max(s0["gam"])))
    for b in betas:
        f0 = r.fx(s0, b)
        out[f"beta {b:g}"] = dict(eps_fuerza_lateral=r.fx(sa, b) / f0 - 1, eps_frenado=r.fx(sb, b) / f0 - 1,
                                  eps_total=r.fx(sab, b) / f0 - 1)
    return out


def main():
    res = {"modelos": {}}
    for modelo in ("giro", "ala"):
        res["modelos"][modelo] = {}
        for m_d in (0.0, 0.5, 1.0):
            for m_p in (0.0, 0.5, 1.0):
                key = f"m_d {m_d:g}, m_p {m_p:g}"
                res["modelos"][modelo][key] = dict(factores=factors(m_d, m_p, modelo),
                                                   turbina=turbine_inflow(m_d, m_p, modelo),
                                                   eps=turbine_on_rotor(m_d, m_p, modelo))
    json.dump(res, open(os.path.join(HERE, "linea_sustentadora.json"), "w"), indent=1, ensure_ascii=False)
    for modelo, dd in res["modelos"].items():
        print(f"== modelo de sección '{modelo}'")
        for key, v in dd.items():
            f = v["factores"]
            print(f"  {key}: eta0 {np.round(f['eta0'], 3)} f_X {np.round(f['f_X'], 3)} "
                  f"f_X placa abajo {np.round(f['f_X placa abajo'], 3)} f30 {np.round(f['f30'], 3)} "
                  f"CDi0 {np.round(f['cdi0'], 3)}")
            t = v["turbina"]
            e = v["eps"]
            print(f"     turbina: ángulo {t['eje']['angulo_medio_grados']:.2f}° (eje), velocidad {t['eje']['velocidad_media_rel']:.3f}; "
                  + "; ".join(f"{b}: eps {e[b]['eps_total']:+.4f} (lateral {e[b]['eps_fuerza_lateral']:+.4f}, "
                              f"frenado {e[b]['eps_frenado']:+.4f})" for b in e if b.startswith("beta")))


if __name__ == "__main__":
    main()
