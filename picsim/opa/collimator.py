"""Fast-axis collimator designs for the stacked-chip relay, and their ray-traced maps for the wave model.

All lenses are cylindrical (power in y only), traced in the (y, z) plane with the facet at the origin. The design
focal length f is paraxial: a ray leaving the facet with s = k_y/k0 exits at height y_S ~ f s. The wave model
(picsim.opa.stack.RelayY) needs, per ray s: the normalised exit height g(s) = y_S/f (g = s for a lens obeying
the Abbe sine condition), its derivative, the residual wavefront W(s) in waves and the amplitude transmission.

Designs:
  ellipse_immersed  glued or printed onto the facet (source inside glass of index n); one elliptical exit surface
                    of eccentricity 1/n collimates exactly (Cartesian oval), but does not obey the sine condition
  plano_asphere     in air, flat side to the facet at working distance wd, exit asphere constructed for exact
                    on-axis collimation (constant optical path); no sine condition
  aplanat           in air, two aspheres designed together (Wassermann-Wolf construction in 2D) for exact
                    collimation and exact sine condition, y_S = f s
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .. import raytrace as rt


@dataclass
class LensMap:
    s: np.ndarray  # s >= 0 grid
    g: np.ndarray
    dg: np.ndarray
    W: np.ndarray  # waves
    t: np.ndarray  # amplitude transmission (Fresnel, uncoated, or 1 if AR coated)
    s_max: float


@dataclass
class Design:
    name: str
    surfaces: list
    indices: list  # source medium, ..., air
    f: float
    z_out: float

    @property
    def n_m(self):
        return self.indices[0]

    def fan(self, s, y0=0.0, z0=0.0):
        th = np.arcsin(np.clip(np.asarray(s) / self.n_m, -1, 1))
        return rt.trace(y0, z0, th, self.surfaces, self.indices, self.z_out)

    def lensmap(self, lam, s_max, n=1201, coated=False):
        s = np.linspace(0, s_max, n)
        r, _ = self.fan(s)
        ok = np.isfinite(r.y)
        if not ok.all():
            s_max = float(s[ok][-1])
            s = np.linspace(0, s_max, n)
            r, _ = self.fan(s)
        g = r.y / self.f
        W = (r.opl - r.opl[0]) / lam
        t = np.ones_like(s) if coated else np.sqrt(r.T)
        return LensMap(s, g, np.gradient(g, s), W, t, s_max)

    def geometry(self, s_max):
        """Lens dimensions, steepest surface slope, largest incidence angle and exit-ray tilt over |s| <= s_max."""
        s = np.linspace(0, s_max, 801)
        r, inc = self.fan(s)
        out = {"exit_height_um": float(np.nanmax(r.y)), "exit_tilt_mrad_max": float(np.nanmax(np.abs(r.dy)) * 1e3),
               "max_incidence_deg": [float(np.nanmax(a)) for a in inc],
               "T_uncoated_min": float(np.nanmin(r.T)), "z_vertices_um": [float(sf.z_v) for sf in self.surfaces]}
        slopes, heights = [], []
        r0, _ = rt.launch(0.0, 0.0, np.arcsin(np.clip(s / self.n_m, -1, 1)), self.n_m)
        rr = r0
        for j, sf in enumerate(self.surfaces):
            rr, _ = rt.refract(rr, sf, self.indices[j], self.indices[j + 1])
            heights.append(float(np.nanmax(np.abs(rr.y))))
            slopes.append(float(np.degrees(np.arctan(np.nanmax(np.abs(sf.dsag(rr.y)))))))
        out["ray_height_at_surfaces_um"] = heights
        out["max_slope_deg"] = slopes
        return out

    def sensitivity(self, s):
        """Optical-path change per unit sag error of each surface, d OPL / d z_j = n_in dz_in - n_out dz_out, and
        the ray height on that surface (for form-error Monte Carlo)."""
        th = np.arcsin(np.clip(np.asarray(s) / self.n_m, -1, 1))
        r, _ = rt.launch(0.0, 0.0, th, self.n_m)
        F, Y = [], []
        for j, sf in enumerate(self.surfaces):
            dz_in = r.dz
            r, _ = rt.refract(r, sf, self.indices[j], self.indices[j + 1])
            F.append(self.indices[j] * dz_in - self.indices[j + 1] * r.dz)
            Y.append(r.y.copy())
        return F, Y


def ellipse_immersed(n, f, gap=60.0):
    """Source inside glass (index n) at the far focus of an ellipse with e = 1/n; paraxial focal length f
    requires a vertex distance r0 = n f."""
    e = 1.0 / n
    r0 = n * f
    a = r0 / (1 + e)
    R = -a * (1 - e**2)
    surf = rt.Conic(z_v=r0, R=R, k=-e**2, y_max=a * np.sqrt(1 - e**2) * 0.9999)
    return Design("ellipse, immersed", [surf], [n, 1.0], f, r0 + gap)


def plano_asphere(n, f, wd, gap=60.0, s_max=0.99, npts=1500):
    """Flat entrance at z = wd, exit surface for exact on-axis collimation; paraxial f = wd + T/n."""
    T = (f - wd) * n
    C = (n - 1) * T
    th = np.arcsin(np.linspace(0, s_max, npts))
    tg = np.arcsin(np.sin(th) / n)
    L = (C + wd - wd / np.cos(th)) / (n - np.cos(tg))
    y = wd * np.tan(th) + L * np.sin(tg)
    z = wd + L * np.cos(tg)
    ok = (L > 0) & np.concatenate([[True], np.diff(y) > 0])
    k = np.argmax(~ok) if (~ok).any() else len(y)
    y, z = y[:k], z[:k]
    z_v = wd + T
    return Design(f"plano-asphere n={n:.3f}", [rt.Flat(wd, y_max=1e4), rt.Tab(z_v, y, z - z_v)], [1.0, n, 1.0], f,
                  z_v + gap)


def aplanat(n, f, r0, T, s_max=0.93, n_steps=4000, gap=60.0):
    """Two-surface aplanatic singlet in air (2D Wassermann-Wolf construction).

    Front surface P1 = r(theta) (sin, cos) from the source; inside the glass the ray has direction phi; the back
    surface point P2 = (f sin(theta), z2). Constant optical path r + n |P2 - P1| - z2 = (n - 1) T and Snell's law at
    the front surface, r' = -r n sin(phi - theta) / (n cos(theta - phi) - 1), determine r(theta); the back surface
    follows from the path condition (Fermat), so it refracts every ray parallel to the axis.
    """
    C = (n - 1) * T
    phic = np.arccos(1.0 / n)
    qmin = np.sqrt(n**2 - 1)

    def phi_of(th, r):
        h = f * np.sin(th)
        a = h - r * np.sin(th)
        rhs = C - r * (1 - np.cos(th))
        if abs(a) < 1e-14:
            return 0.0
        target = rhs / abs(a)
        if target < qmin:
            raise ValueError("no refraction solution (lens too strong)")
        g = lambda p: (n - np.cos(p)) / np.sin(p) - target
        p = brentq(g, 1e-14, phic)
        return np.sign(a) * p

    def drdth(th, r):
        p = phi_of(th, r)
        return -r * n * np.sin(p - th) / (n * np.cos(th - p) - 1)

    th_max = np.arcsin(s_max)
    ths = np.linspace(0, th_max, n_steps)
    rs = np.empty_like(ths)
    rs[0] = r0
    for i in range(1, len(ths)):  # RK4
        h = ths[i] - ths[i - 1]
        t0, r_ = ths[i - 1], rs[i - 1]
        try:
            k1 = drdth(t0, r_) if t0 > 0 else 0.0
            k2 = drdth(t0 + h / 2, r_ + h / 2 * k1)
            k3 = drdth(t0 + h / 2, r_ + h / 2 * k2)
            k4 = drdth(t0 + h, r_ + h * k3)
        except ValueError:
            ths, rs = ths[:i], rs[:i]
            break
        rs[i] = r_ + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    phis = np.array([phi_of(t, r) if t > 0 else 0.0 for t, r in zip(ths, rs)])
    P1y, P1z = rs * np.sin(ths), rs * np.cos(ths)
    h = f * np.sin(ths)
    L2 = np.where(ths > 0, (h - P1y) / np.where(np.sin(phis) == 0, 1, np.sin(phis)), T)
    L2[0] = T
    P2z = P1z + L2 * np.cos(phis)
    s1 = rt.Tab(r0, P1y, P1z - r0)
    s2 = rt.Tab(r0 + T, h, P2z - (r0 + T))
    d = Design(f"aplanat n={n:.3f}", [s1, s2], [1.0, n, 1.0], f, r0 + T + gap)
    d.construction = {"theta_max_deg": float(np.degrees(ths[-1])), "s_max": float(np.sin(ths[-1])),
                      "edge_thickness_um": float(np.min(P2z - P1z)), "front_height_um": float(P1y[-1]),
                      "back_height_um": float(h[-1])}
    return d


def weierstrass_aplanat(n, f, R, gap, T2, n2=None, **kw):
    """Immersed aplanatic (Weierstrass) cylinder plus an air-spaced aplanat.

    A glass block (index n) glued or printed onto the facet ends in a circular cylinder of radius R whose centre
    lies R/n beyond the facet. Its aplanatic points make every ray from the facet appear to diverge, without
    aberration and with the sine condition, from a virtual point n R behind the centre, with sin(u') = s / n^2.
    An aplanat of focal length f n^2 (designed for that virtual point) then collimates the reduced numerical
    aperture. The pair obeys the sine condition with focal length f.
    """
    n2 = n if n2 is None else n2
    z_c = R / n
    z_vertex = z_c + R
    z_virt = z_c - n * R
    circle = rt.Conic(z_v=z_vertex, R=-R, k=0.0, y_max=0.999 * R)
    d2 = aplanat(n2, f * n**2, (z_vertex + gap) - z_virt, T2, **kw)
    s1, s2 = d2.surfaces
    s1.z_v += z_virt
    s2.z_v += z_virt
    d = Design(f"Weierstrass + aplanat n={n:.3f}", [circle, s1, s2], [n, 1.0, n2, 1.0], f, s2.z_v + 60.0)
    d.construction = {**d2.construction, "s_max_system": d2.construction["s_max"] * n**2, "R_um": R,
                      "block_length_um": z_vertex, "virtual_source_um": z_virt}
    return d


def fac_design(n_glass, n2=1.8, f=175.0):
    """Adopted collimator: Weierstrass block (R = 80 um) + 10 um air gap + 250 um aplanat of index n2."""
    return weierstrass_aplanat(n_glass, f, 80.0, 10.0, 250.0, n2=n2)


def back_focal_plane(d: Design, x0=1.0):
    """z where a chief ray (launched parallel to the axis from height x0) crosses the axis."""
    r, _ = rt.trace(x0, 0.0, np.array([0.0]), d.surfaces, d.indices, d.z_out)
    return float(r.z[0] - r.y[0] * r.dz[0] / r.dy[0])
