#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lifting line of Flettner rotors on deck (A1, version 4, 22-10-2026), for what the route model does not resolve: the
lift factor of the extension (f_X) and of the 30 m rotor (f30), the three-dimensional correction of set B, the
velocity induced by the rotor at the turbine and the order of magnitude of the interaction eps.

Model: discrete lifting line (Prandtl horseshoes on the rotor axis, with straight wakes downstream), with images in
the deck (factor m_d: 1 = perfect symmetry plane, effective AR 2H/D; 0 = no reflection, effective AR H/D) and in the
endplate (factor m_p). Two section models:
  "wing": analogy with a wing of slope 2 pi, C_l = 2 pi (alpha_g - alpha_i) with alpha_g = C_l,2D(k)/(2 pi): the
          classical finite aspect-ratio correction. It is the model of the factors: it handles free ends.
  "spin": the circulation is set by the spin ratio, Gamma = 1/2 U_ef D C_l(k_ef), with k_ef = omega R / U_ef: a
          cylinder is axisymmetric and the induced velocity only turns the flow and changes its magnitude, so the
          finite span adds induced drag but hardly removes lift. It does not reproduce the loss at a free end (the
          circulation does not vanish there), so it is used only with deck and plate as perfect mirrors, as an upper
          bound (no loss).
C_l,2D(k) is that of set B (Bordogna et al. 2019, wall-to-wall test, nearly two-dimensional; points up to k = 4.6).
Uniform flow: the wind shear is already in the route model (wind at 40 and 55 m).
Usage: python3 lifting_line.py   ->  lifting_line.json
"""
import json
import math
import os
import numpy as np
import coefficients as cf

HERE = os.path.dirname(os.path.abspath(__file__))
D_F, H_F, D_X, H_X, D_E = 4.0, 24.0, 3.4, 6.0, 8.0
Z_PLATE_GAP = 0.35            # gap between the plate (at the top of the rotor) and the turbine
L_WAKE = 1.0e5                # length of the wakes (m), "infinite"

# C_l,2D(k) of set B, extended to k = 4.6 with the points at Re = 3.6e5 of Bordogna et al. (2019)
K_2D = np.array(cf.K_B + [3.65, 4.13, 4.60])
CL_2D = np.array(cf.CL_B + [6.50, 7.31, 8.46])


def cl2d(k):
    return np.interp(k, K_2D, CL_2D)


def _segment(p, a, b, gam):
    """Induced velocity (vector) at the points p (n x 3) by a vortex segment from a to b with circulation gam."""
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
    """Horseshoe with the bound segment on the axis (x = 0) from z0 to z1 and wakes towards +x."""
    A, B = np.array([0.0, 0.0, z0]), np.array([0.0, 0.0, z1])
    Af, Bf = np.array([L_WAKE, 0.0, z0]), np.array([L_WAKE, 0.0, z1])
    return _segment(p, Af, A, gam) + _segment(p, A, B, gam) + _segment(p, B, Bf, gam)


class Rotor:
    """Stepped rotor: segments (z_start, z_end, diameter); plate at z_plate (or none)."""

    def __init__(self, segments, z_plate=None, m_d=1.0, m_p=0.5, n=80):
        self.segments = segments
        z_top = segments[-1][1]
        # panels with cosine spacing in each segment
        zs, ds = [], []
        for (za, zb, d) in segments:
            m = max(8, int(round(n * (zb - za) / z_top)))
            e = za + (zb - za) * 0.5 * (1 - np.cos(np.linspace(0, math.pi, m + 1)))
            zs.append(e)
            ds += [d] * m
        edges = np.unique(np.concatenate(zs))
        self.e0, self.e1 = edges[:-1], edges[1:]
        self.zc = 0.5 * (self.e0 + self.e1)
        self.dz = self.e1 - self.e0
        self.D = np.array([next(d for (za, zb, d) in segments if za - 1e-9 <= z <= zb + 1e-9) for z in self.zc])
        self.z_plate, self.m_d, self.m_p = z_plate, m_d, m_p
        self.A = self._influence(np.column_stack([np.zeros_like(self.zc), np.zeros_like(self.zc), self.zc]))

    def _influence(self, pts):
        """Matrix (points x panels) of the lateral (y) velocity induced by unit circulation, with images."""
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

    def solve(self, k_ref, model="spin", u=1.0, v_ext=None, it=200):
        """Circulation (m2/s) on each panel with wind u (scalar or per panel) and an angular velocity such that the
        lower segment runs at k_ref with unit wind; v_ext: external lateral velocity (e.g., that induced by the
        turbine). Convention: with Gamma > 0 along +z and wind towards +x, the lift points to -y, so an induced
        velocity v > 0 opposes it (alpha_i = v/u). Returns lift and induced drag per unit length (divided by
        rho)."""
        u = np.broadcast_to(np.asarray(u, float), self.zc.shape).copy()
        v_ext = np.zeros_like(self.zc) if v_ext is None else np.asarray(v_ext, float)
        omega = k_ref / (self.segments[0][2] / 2)
        k = omega * (self.D / 2) / u
        if model == "wing":
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
        # Kutta-Joukowski with the local flow (u, v): F/rho = (u x + v y) x Gamma z = Gamma (-u y + v x)
        lift = u * gam                               # / rho: lift per unit length (towards -y)
        drag_i = v * gam                             # / rho: induced drag per unit length
        drag_p = 0.5 * u ** 2 * self.D * np.interp(k, cf.K_B, cf.CD_B2D)   # / rho: profile drag
        return dict(gam=gam, v=v, k=k, u=u, lift=lift, drag_i=drag_i, drag_p=drag_p, omega=omega)

    def fx(self, sol, beta_deg):
        """Longitudinal force (/rho) of the whole rotor with the apparent wind at beta from the heading."""
        b = math.radians(beta_deg)
        return float(((sol["lift"] * math.sin(b) - (sol["drag_i"] + sol["drag_p"]) * math.cos(b)) * self.dz).sum())

    def segment(self, sol, za, zb):
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
        "C2 plate below": Rotor([(0.0, H_F, D_F), (H_F, H_F + H_X, D_X)], z_plate=None, m_d=m_d, m_p=0.0),
        "R1-30": Rotor([(0.0, H_F + H_X, D_F)], z_plate=H_F + H_X, m_d=m_d, m_p=m_p),
    }


def factors(m_d, m_p, model, ks=(1.5, 2.0, 2.5, 3.0)):
    """f_X, f30 and efficiency of C0 (3D / 2D lift) with the criterion of the route model (eq. 5): the main rotor keeps
    the lift of C0 and the factor collects everything else."""
    cfg = configs(m_d, m_p)
    out = {"eta0": [], "f_X": [], "f_X plate below": [], "f30": [], "cdi0": [], "cdi_C2_extra": []}
    for kr in ks:
        s0 = cfg["C0"].solve(kr, model)
        l0, d0 = cfg["C0"].segment(s0, 0, H_F)
        l0_2d = cfg["C0"].lift_2d(s0, 0, H_F)
        eta0 = l0 / l0_2d
        out["eta0"].append(eta0)
        out["cdi0"].append(d0 / (0.5 * D_F * H_F))
        for key, name, top_d in (("f_X", "C2", D_X), ("f_X plate below", "C2 plate below", D_X), ("f30", "R1-30", D_F)):
            r = cfg[name]
            s = r.solve(kr, model)
            lt, dt = r.segment(s, 0, H_F + H_X)
            l_ext_2d = r.lift_2d(s, H_F, H_F + H_X)
            out[key].append((lt - l0) / (eta0 * l_ext_2d))
            if key == "f_X":
                # induced drag of C2 against the sum assumed by eq. (5): that of C0 plus that of the upper segment
                # with the same induced coefficient as C0 at its lift
                cl_ext = (lt - l0) / (0.5 * D_X * H_X)
                cdi_ext_ec5 = out["cdi0"][-1] * (cl_ext / (l0 / (0.5 * D_F * H_F))) ** 2
                out["cdi_C2_extra"].append((dt - d0 - cdi_ext_ec5 * 0.5 * D_X * H_X) / (0.5 * D_F * H_F))
    return {k: [float(x) for x in v] for k, v in out.items()} | {"k": list(ks)}


def _field(r, s, pts):
    """Velocity (vx, vy) induced at pts by the vortex system of the solved rotor, with its images."""
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


def turbine_inflow(m_d, m_p, model, k_ref=2.5):
    """Flow induced by the rotor (with its plate) at the turbine, from 24.35 to 30.35 m above the deck and radius
    1.7 m, with unit wind: mean angle (degrees; positive towards the lift side of the rotor) and mean relative
    magnitude, on the axis and at four points of the edge of the disc."""
    r = configs(m_d, m_p)["C0"]
    s = r.solve(k_ref, model)
    zs = np.linspace(H_F + Z_PLATE_GAP, H_F + Z_PLATE_GAP + H_X, 25)
    out = {}
    for name, (x, y) in {"axis": (0.0, 0.0), "windward": (-1.7, 0.0), "leeward": (1.7, 0.0),
                         "lift side": (0.0, -1.7), "opposite side": (0.0, 1.7)}.items():
        pts = np.column_stack([np.full_like(zs, x), np.full_like(zs, y), zs])
        vx, vy = _field(r, s, pts)
        out[name] = dict(mean_angle_deg=float(np.degrees(np.mean(np.arctan2(-vy, 1.0 + vx)))),
                         mean_speed_rel=float(np.mean(np.hypot(1.0 + vx, vy))))
    return out


def turbine_on_rotor(m_d, m_p, model, k_ref=2.5, cs=0.19, ct=1.0, betas=(60.0, 90.0, 120.0)):
    """Order of magnitude of eps, the relative change of the thrust of the rotor with the turbine deployed above the
    plate:
    (a) the side force of the turbine (C_s, on the lift side of the rotor) as a horseshoe from 24.35 to 30.35 m, with
        circulation per unit speed Gamma_T = C_s D_X / 2, induces lateral velocity at the rotor;
    (b) the thrust of the turbine (C_t) as an actuator disc of radius D_X/2 slows the wind at the top of the rotor: on
        the axis of the disc the axial induction is a (1 - x/(R^2 + x^2)^0.5) at the distance x from the disc, with
        a = (1 - (1 - C_t)^0.5)/2 (momentum theory).
    Returns eps for (a), (b) and both, with the apparent wind at 60, 90 and 120 degrees."""
    r = configs(m_d, m_p)["C0"]
    s0 = r.solve(k_ref, model)
    pts = np.column_stack([np.zeros_like(r.zc), np.zeros_like(r.zc), r.zc])
    zt0, zt1 = H_F + Z_PLATE_GAP, H_F + Z_PLATE_GAP + H_X
    v_t = horseshoe(pts, zt0, zt1, 0.5 * cs * D_X)[:, 1]
    a = 0.5 * (1 - math.sqrt(max(0.0, 1 - min(ct, 0.96))))
    x = (zt0 + zt1) / 2 - r.zc                     # distance to the centre of the disc, measured along its axis
    rt = D_X / 2
    u_b = 1.0 - a * (1 - x / np.sqrt(rt ** 2 + x ** 2))
    sa = r.solve(k_ref, model, v_ext=v_t)
    sb = r.solve(k_ref, model, u=u_b)
    sab = r.solve(k_ref, model, u=u_b, v_ext=v_t)
    out = dict(axial_induction=a, gamma_T_rel=float(0.5 * cs * D_X / np.max(s0["gam"])))
    for b in betas:
        f0 = r.fx(s0, b)
        out[f"beta {b:g}"] = dict(eps_side_force=r.fx(sa, b) / f0 - 1, eps_braking=r.fx(sb, b) / f0 - 1,
                                  eps_total=r.fx(sab, b) / f0 - 1)
    return out


def main():
    res = {"models": {}}
    for model in ("spin", "wing"):
        res["models"][model] = {}
        for m_d in (0.0, 0.5, 1.0):
            for m_p in (0.0, 0.5, 1.0):
                key = f"m_d {m_d:g}, m_p {m_p:g}"
                res["models"][model][key] = dict(factors=factors(m_d, m_p, model),
                                                 turbine=turbine_inflow(m_d, m_p, model),
                                                 eps=turbine_on_rotor(m_d, m_p, model))
    json.dump(res, open(os.path.join(HERE, "lifting_line.json"), "w"), indent=1, ensure_ascii=False)
    for model, dd in res["models"].items():
        print(f"== section model '{model}'")
        for key, v in dd.items():
            f = v["factors"]
            print(f"  {key}: eta0 {np.round(f['eta0'], 3)} f_X {np.round(f['f_X'], 3)} "
                  f"f_X plate below {np.round(f['f_X plate below'], 3)} f30 {np.round(f['f30'], 3)} "
                  f"CDi0 {np.round(f['cdi0'], 3)}")
            t = v["turbine"]
            e = v["eps"]
            print(f"     turbine: angle {t['axis']['mean_angle_deg']:.2f}° (axis), speed {t['axis']['mean_speed_rel']:.3f}; "
                  + "; ".join(f"{b}: eps {e[b]['eps_total']:+.4f} (side force {e[b]['eps_side_force']:+.4f}, "
                              f"braking {e[b]['eps_braking']:+.4f})" for b in e if b.startswith("beta")))


if __name__ == "__main__":
    main()
