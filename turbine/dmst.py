#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H-Darrieus turbine: double-multiple streamtube model (DMST). A1, 21-10-2026.

Model of Paraschivoiu (2002) for a straight-bladed rotor, with:
  - static polars of the NACA 0021 from Sheldahl and Klimas (1981), in the version distributed with CACTUS (Sandia,
    BSD-3 licence; data/NACA_0021_SK.dat), interpolated in log(Re) and in the angle of attack;
  - dynamic stall of Gormont (1973) with the modification of Berg (1983), A_M = 6, in the form used by QBlade;
    K1 = 1 when |alpha| grows and K1 = 0 when it falls, as Shires (2013) recommends for H-rotors
    (option k1_dec = 0.5 for the original model);
  - high-thrust correction of Buhl (2005) for a > 0.4;
  - finite aspect-ratio correction of the blades (induced angle and Prandtl induced drag, only before stall);
  - parasitic losses of the struts (power) and drag of the central shaft (thrust).
Validation: dmst_validation.py (Battisti et al. 2018, wind tunnel; Kjellin et al. 2011, 12 kW turbine in the field).
Outputs are the power coefficient Cp, the thrust coefficient along the wind Ct and the side-force coefficient Cs, all
referred to the swept area A = 2 R H.

Conventions: wind along +x; rotor seen from above turning counterclockwise; azimuth psi of each blade measured from +x
counterclockwise (upwind blade with cos(psi) < 0). Cs > 0 is a force along +y, that is, to the left looking
downwind. For clockwise rotation, Cs changes sign.
With the same sense of rotation, the Magnus lift of a cylinder points along -y.

Usage: python dmst.py   (table of Cp, Ct and Cs of the case-study turbine)
"""
import math
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
POLAR_FILE = os.path.join(HERE, "data", "NACA_0021_SK.dat")


class Polar:
    """Static 2D polars (CL, CD) on a common grid of angles, interpolated in log10(Re)."""

    def __init__(self, path=POLAR_FILE, step=0.25):
        blocks, cur = [], None
        with open(path, encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if s.startswith("Reynolds Number:"):
                    cur = {"re": float(s.split(":")[1]), "a": [], "cl": [], "cd": []}
                    blocks.append(cur)
                    continue
                if cur is None or not s or not (s[0].isdigit() or s[0] == "-"):
                    continue
                parts = s.split()
                if len(parts) >= 3:
                    cur["a"].append(float(parts[0]))
                    cur["cl"].append(float(parts[1]))
                    cur["cd"].append(float(parts[2]))
        self.alpha = np.arange(-180.0, 180.0 + 1e-9, step)
        self.re = np.array([b["re"] for b in blocks])
        self.lre = np.log10(self.re)
        self.cl = np.array([np.interp(self.alpha, b["a"], b["cl"]) for b in blocks])
        self.cd = np.array([np.interp(self.alpha, b["a"], b["cd"]) for b in blocks])
        # angle of maximum lift (static stall) at each Re, positive side between 5 and 30 degrees
        m = (self.alpha >= 5.0) & (self.alpha <= 30.0)
        self.a_ss = np.array([self.alpha[m][np.argmax(c[m])] for c in self.cl])
        self.step = step

    def _w(self, re):
        lr = np.clip(np.log10(np.maximum(re, 1.0)), self.lre[0], self.lre[-1])
        j = np.clip(np.searchsorted(self.lre, lr) - 1, 0, len(self.lre) - 2)
        t = (lr - self.lre[j]) / (self.lre[j + 1] - self.lre[j])
        return j, t

    def coeffs(self, alpha_deg, re):
        """Static CL and CD at alpha (degrees, any value) and Re (vectors of the same shape)."""
        a = (np.asarray(alpha_deg, float) + 180.0) % 360.0 - 180.0
        j, t = self._w(np.asarray(re, float))
        x = (a + 180.0) / self.step
        i = np.clip(np.floor(x).astype(int), 0, len(self.alpha) - 2)
        u = x - i
        out = []
        for tab in (self.cl, self.cd):
            v0 = tab[j, i] * (1 - u) + tab[j, i + 1] * u
            v1 = tab[j + 1, i] * (1 - u) + tab[j + 1, i + 1] * u
            out.append(v0 * (1 - t) + v1 * t)
        return out[0], out[1]

    def stall_angle(self, re):
        j, t = self._w(np.asarray(re, float))
        return self.a_ss[j] * (1 - t) + self.a_ss[j + 1] * t


_POLAR = None


def polar():
    global _POLAR
    if _POLAR is None:
        _POLAR = Polar()
    return _POLAR


def ct_momentum(a):
    """Thrust coefficient of the streamtube as a function of the induction: 4a(1-a), with the Buhl correction if a > 0.4."""
    a = np.asarray(a, float)
    return np.where(a <= 0.4, 4.0 * a * (1.0 - a), 8.0 / 9.0 - 4.0 / 9.0 * a + 14.0 / 9.0 * a * a)


def a_from_ct(ct):
    """Inverse of ct_momentum (low-induction branch; Buhl above 0.96)."""
    ct = np.asarray(ct, float)
    low = 0.5 * (1.0 - np.sqrt(np.maximum(1.0 - np.minimum(ct, 0.96), 0.0)))
    high = (4.0 / 9.0 + np.sqrt(np.maximum(16.0 / 81.0 - 4.0 * 14.0 / 9.0 * (8.0 / 9.0 - np.minimum(ct, 2.0)), 0.0))) \
        / (2.0 * 14.0 / 9.0)
    return np.where(ct <= 0.96, low, high)


class Turbine:
    """Geometry: n blades of chord c, radius R and height H; t/c of the airfoil; struts and shaft."""

    def __init__(self, n=3, R=1.7, H=6.0, c=0.283, tc=0.21, n_struts=2, c_strut=0.10, cd_strut=0.02,
                 d_shaft=0.30, cd_shaft=0.6):
        self.n, self.R, self.H, self.c, self.tc = n, R, H, c, tc
        self.n_struts, self.c_strut, self.cd_strut = n_struts, c_strut, cd_strut
        self.d_shaft, self.cd_shaft = d_shaft, cd_shaft

    @property
    def area(self):
        return 2.0 * self.R * self.H

    @property
    def solidity(self):
        """Solidity N c / D."""
        return self.n * self.c / (2.0 * self.R)


def _blade(alpha_rad, re, adot, w, c, tc, ar, opts, pol):
    """CL and CD of the blade (3D, with dynamic stall) at alpha (rad), Re, d(alpha)/dt (rad/s) and relative speed w."""
    a_deg = np.degrees(alpha_rad)
    a_ss = pol.stall_angle(re)
    # finite aspect ratio before stall: effective 2D alpha and induced drag
    if opts.get("ar", True) and ar is not None:
        cl0, _ = pol.coeffs(a_deg, re)
        k_ar = np.clip((2.0 * a_ss - np.abs(a_deg)) / a_ss, 0.0, 1.0)      # 1 before stall, 0 at 2 alpha_ss
        a_eff = a_deg - k_ar * np.degrees(cl0 / (math.pi * ar))
    else:
        k_ar, a_eff = 0.0, a_deg
    cl_s, cd_s = pol.coeffs(a_eff, re)
    if opts.get("ds", True):
        tu = c / (2.0 * np.maximum(w, 1e-6))
        grow = a_deg * adot >= 0.0
        k1 = np.where(grow, 1.0, opts.get("k1_dec", 0.0))
        root = np.sqrt(np.abs(tu * adot)) * np.sign(adot)
        y_l = 1.4 - 6.0 * (0.06 - tc)
        y_m = 1.0 - 2.5 * (0.06 - tc)
        a_rl = np.degrees(np.radians(a_eff) - k1 * y_l * root)
        a_rm = np.degrees(np.radians(a_eff) - k1 * y_m * root)
        cl_r, _ = pol.coeffs(a_rl, re)
        _, cd_r = pol.coeffs(a_rm, re)
        safe = np.abs(a_rl) > 0.5
        cl_gb = np.where(safe, cl_r / np.where(safe, a_rl, 1.0) * a_eff, cl_s)
        a_end = opts.get("a_m", 6.0) * a_ss
        gam = np.clip((a_end - np.abs(a_eff)) / (a_end - a_ss), 0.0, 1.0)
        cl = cl_s + gam * (cl_gb - cl_s)
        cd = cd_s + gam * (cd_r - cd_s)
    else:
        cl, cd = cl_s, cd_s
    if ar is not None and opts.get("ar", True):
        cd = cd + k_ar * cl * cl / (math.pi * ar * 0.9)
    return cl, cd


def solve(tur, U, lam, rho=1.225, nu=1.48e-5, n_half=36, opts=None, iters=40, return_detail=False):
    """Cp, Ct and Cs of the turbine with wind U (m/s) and tip-speed ratio lam = omega R / U."""
    opts = dict(opts or {})
    pol = polar()
    R, c, n = tur.R, tur.c, tur.n
    omega = lam * U / R
    ar = tur.H / c if opts.get("ar", True) else None
    dpsi = math.pi / n_half
    psi_u = math.pi / 2 + (np.arange(n_half) + 0.5) * dpsi          # upwind: cos(psi) < 0
    psi_d = math.pi - psi_u                                          # downwind, same streamtube (same y)
    kfac_u = n * c / (2.0 * math.pi * R * np.abs(np.cos(psi_u)))
    kfac_d = n * c / (2.0 * math.pi * R * np.abs(np.cos(psi_d)))
    adot_u = np.zeros(n_half)
    adot_d = np.zeros(n_half)

    def half(psi, u_in, kfac, adot, a0):
        a = a0.copy()
        for _ in range(400):
            ul = u_in * (1.0 - a)
            wx = omega * R + ul * np.sin(psi)
            wy = ul * np.cos(psi)
            w = np.hypot(wx, wy)
            al = np.arctan2(wy, wx)
            cl, cd = _blade(al, w * c / nu, adot, w, c, tur.tc, ar, opts, pol)
            ct_b = cl * np.sin(al) - cd * np.cos(al)                  # tangential (driving)
            cn_b = cl * np.cos(al) + cd * np.sin(al)                  # normal (outwards)
            cx = -ct_b * np.sin(psi) + cn_b * np.cos(psi)             # on the blade, along +x
            ct_tube = kfac * (w / np.maximum(u_in, 1e-9)) ** 2 * cx
            a_new = a_from_ct(ct_tube)
            if np.max(np.abs(a_new - a)) < 1e-6:
                a = a_new
                break
            a = a + 0.3 * (a_new - a)
        ul = u_in * (1.0 - a)
        wx = omega * R + ul * np.sin(psi)
        wy = ul * np.cos(psi)
        w = np.hypot(wx, wy)
        al = np.arctan2(wy, wx)
        cl, cd = _blade(al, w * c / nu, adot, w, c, tur.tc, ar, opts, pol)
        ct_b = cl * np.sin(al) - cd * np.cos(al)
        cn_b = cl * np.cos(al) + cd * np.sin(al)
        return dict(a=a, w=w, al=al, ct=ct_b, cn=cn_b)

    def alpha_rate(psi, ul):
        # d(alpha)/dt = omega d(alpha)/d(psi) with the local speed of each streamtube (analytical derivative, no noise)
        wx = omega * R + ul * np.sin(psi)
        wy = ul * np.cos(psi)
        return omega * (-ul * omega * R * np.sin(psi) - ul * ul) / (wx * wx + wy * wy)

    cp_prev = None
    a_u0, a_d0 = np.full(n_half, 0.1), np.full(n_half, 0.1)
    for it in range(iters):
        up = half(psi_u, np.full(n_half, U), kfac_u, adot_u, a_u0)
        u_e = np.maximum(U * (1.0 - 2.0 * up["a"]), 0.05 * U)
        dn = half(psi_d, u_e, kfac_d, adot_d, a_d0)
        a_u0, a_d0 = up["a"], dn["a"]
        if not opts.get("ds", True):
            break
        new_u = alpha_rate(psi_u, U * (1.0 - up["a"]))
        new_d = alpha_rate(psi_d, u_e * (1.0 - dn["a"]))
        adot_u = new_u if it == 0 else 0.5 * (adot_u + new_u)
        adot_d = new_d if it == 0 else 0.5 * (adot_d + new_d)
        cp_now = float(np.sum(up["w"] ** 2 * up["ct"]) + np.sum(dn["w"] ** 2 * dn["ct"]))
        if cp_prev is not None and abs(cp_now - cp_prev) < 1e-5 * abs(cp_prev):
            break
        cp_prev = cp_now

    q_inf = 0.5 * rho * U * U
    tot = {"Q": 0.0, "Fx": 0.0, "Fy": 0.0}
    for psi, h in ((psi_u, up), (psi_d, dn)):
        qw = 0.5 * rho * h["w"] ** 2 * c
        ft, fn = qw * h["ct"], qw * h["cn"]
        fx = -ft * np.sin(psi) + fn * np.cos(psi)
        fy = ft * np.cos(psi) + fn * np.sin(psi)
        tot["Q"] += n / (2 * math.pi) * np.sum(ft * R) * dpsi
        tot["Fx"] += n / (2 * math.pi) * np.sum(fx) * dpsi
        tot["Fy"] += n / (2 * math.pi) * np.sum(fy) * dpsi
    H, A = tur.H, tur.area
    cp_blades = tot["Q"] * H * omega / (q_inf * U * A)
    cp_struts = tur.n * tur.n_struts * tur.cd_strut * tur.c_strut * lam ** 3 / (8.0 * H)
    cp = cp_blades - cp_struts
    ct = tot["Fx"] * H / (q_inf * A) + tur.cd_shaft * tur.d_shaft / (2.0 * R)
    cs = tot["Fy"] * H / (q_inf * A)
    out = dict(cp=float(cp), ct=float(ct), cs=float(cs), cp_blades=float(cp_blades), cp_struts=float(cp_struts),
               a_up_max=float(up["a"].max()), a_dn_max=float(dn["a"].max()))
    if return_detail:
        out.update(psi_u=psi_u, psi_d=psi_d, up=up, dn=dn)
    return out


def parked(tur, U, rho=1.225, nu=1.48e-5, n_pos=24):
    """Coefficients of the parked turbine (Ct and Cs), averaged over the angular position of the rotor."""
    pol = polar()
    re = U * tur.c / nu
    cts, css = [], []
    for k in range(n_pos):
        psi0 = 2 * math.pi * k / (n_pos * tur.n)
        fx = fy = 0.0
        for b in range(tur.n):
            psi = psi0 + 2 * math.pi * b / tur.n
            al = math.atan2(math.cos(psi), math.sin(psi))
            cl, cd = pol.coeffs(np.degrees(al), re)
            cl, cd = float(cl), float(cd)
            ct_b = cl * math.sin(al) - cd * math.cos(al)
            cn_b = cl * math.cos(al) + cd * math.sin(al)
            fx += -ct_b * math.sin(psi) + cn_b * math.cos(psi)
            fy += ct_b * math.cos(psi) + cn_b * math.sin(psi)
        cts.append(fx * tur.c * tur.H / tur.area)
        css.append(fy * tur.c * tur.H / tur.area)
    ct = float(np.mean(cts)) + tur.cd_shaft * tur.d_shaft / (2.0 * tur.R)
    return dict(ct=ct, cs=float(np.mean(css)))


if __name__ == "__main__":
    t = Turbine()
    print(f"Case-study turbine: R {t.R} m, H {t.H} m, {t.n} NACA 0021 blades, c {t.c} m, solidity Nc/D {t.solidity:.3f}")
    for U in (8.0, 14.0):
        print(f"U = {U} m/s")
        for lam in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0):
            r = solve(t, U, lam)
            print(f"  lambda {lam:3.1f}  Cp {r['cp']:6.3f}  Ct {r['ct']:6.3f}  Cs {r['cs']:6.3f}  "
                  f"a_max {r['a_up_max']:.2f}/{r['a_dn_max']:.2f}")
    print("Parked:", parked(t, 14.0))
