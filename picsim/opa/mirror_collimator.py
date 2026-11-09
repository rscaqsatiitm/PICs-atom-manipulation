"""All-reflective fast-axis collimator: two cylindrical mirrors in one glass block glued to the facet.

Why reflective. A cylindrical surface (no power in x) conserves k_x, so light of a tweezer at x_a, with x direction
sine s_x, sees the y-z problem of the 2D design with every index n replaced by n_eff = sqrt(n^2 - s_x^2). A refracting
surface then changes power (docs/megascale.md Sec. 12). A mirror does not: the law of reflection in the y-z
projection is the 2D law, and in one medium the phase is k0 n_eff L_yz along the whole path. A 2D design that
collimates exactly therefore collimates exactly for every s_x. If it also obeys the sine condition in the medium,
y_S - y_c = sign n f sin(theta), it stays aplanatic at every s_x, with focal length f n / n_eff, and the rim of the
mirrors fixes theta, so the stripe edge at S does not move.

Construction (2D, Wassermann-Wolf with mirrors). The facet is at the origin and radiates into glass of index n
along +z, with ray angle theta. M1 = r(theta) (sin theta, cos theta). After M1 the ray has direction
phi (sin phi, cos phi) and meets M2 at height y_S(theta) = y_c + sign n f sin(theta); M2 sends it along +z. Constant
optical path to a plane z = const gives, in closed form,
    tan(phi / 2) = (K - r (1 - cos theta)) / (y_S - r sin theta),   K = y_c tan(phi0 / 2),
and the law of reflection at M1 gives r' = r cot((phi - theta) / 2). M2 = M1 + L (sin phi, cos phi) follows
(Fermat), and it reflects every ray parallel to z. Free parameters: r0 (distance of M1 on axis), y_c (centre of the
collimated beam), phi0 (direction of the axial ray after M1) and sign (mapping orientation).

The design is traced independently through spline fits of the two mirror curves, in 3D (rays with k_x), with the
vector law of reflection, optionally with s/p Jones coefficients of the coating. Lengths in um.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import CubicSpline

from .collimator import LensMap


# ---------------------------------------------------------------------------------------------- mirror curves
class Curve:
    """A mirror profile (y(u), z(u)) in the y-z plane (a cylinder along x), parameter u (construction angle)."""

    def __init__(self, u, y, z, u_lim=None):
        self.u = np.asarray(u, float)
        self.cy = CubicSpline(self.u, y)
        self.cz = CubicSpline(self.u, z)
        self.u_lim = (self.u[0], self.u[-1]) if u_lim is None else u_lim  # clear aperture in u

    def point(self, u):
        return self.cy(u), self.cz(u)

    def normal(self, u):
        ty, tz = self.cy(u, 1), self.cz(u, 1)
        nn = np.hypot(ty, tz)
        return -tz / nn, ty / nn  # (N_y, N_z); sign is irrelevant for reflection

    def curvature_radius(self, u):
        y1, z1, y2, z2 = self.cy(u, 1), self.cz(u, 1), self.cy(u, 2), self.cz(u, 2)
        k = np.abs(y1 * z2 - z1 * y2) / np.hypot(y1, z1) ** 3
        return 1.0 / np.maximum(k, 1e-300)

    def intersect(self, py, pz, dy, dz, n_grid=1201, n_bisect=60):
        """First crossing (smallest t > 0) of rays p + t d (projected on y-z) with the curve. Returns (u, t); NaN if
        the ray misses the curve or meets it outside the clear aperture u_lim."""
        ug = np.linspace(self.u[0], self.u[-1], n_grid)
        Y, Z = self.cy(ug)[:, None], self.cz(ug)[:, None]
        F = (Y - py[None]) * dz[None] - (Z - pz[None]) * dy[None]
        T = ((Y - py[None]) * dy[None] + (Z - pz[None]) * dz[None]) / (dy**2 + dz**2)[None]
        sc = (np.sign(F[:-1]) * np.sign(F[1:]) <= 0) & (T[:-1] > 1e-6) & (T[1:] > 1e-6)
        Tm = np.where(sc, 0.5 * (T[:-1] + T[1:]), np.inf)
        k = np.argmin(Tm, axis=0)
        j = np.arange(len(py))
        found = np.isfinite(Tm[k, j])
        lo, hi = ug[k], ug[k + 1]
        flo = F[k, j]
        for _ in range(n_bisect):
            mid = 0.5 * (lo + hi)
            fm = (self.cy(mid) - py) * dz - (self.cz(mid) - pz) * dy
            same = np.sign(fm) == np.sign(flo)
            lo, flo = np.where(same, mid, lo), np.where(same, fm, flo)
            hi = np.where(same, hi, mid)
        u = 0.5 * (lo + hi)
        t = ((self.cy(u) - py) * dy + (self.cz(u) - pz) * dz) / (dy**2 + dz**2)
        ok = found & (u >= self.u_lim[0] - 1e-12) & (u <= self.u_lim[1] + 1e-12)
        return np.where(ok, u, np.nan), np.where(ok, t, np.nan)


def fresnel_metal(n1, n2, cos_i):
    """s and p amplitude reflection coefficients from medium n1 onto n2 (complex). Convention: with s = d x N and
    p = s x d (before) and p' = s x d' (after), a perfect conductor gives r_s = -1, r_p = +1 (checked in the tests)."""
    sin2 = 1 - cos_i**2
    cos_t = np.sqrt(1 - (n1 / n2) ** 2 * sin2 + 0j)
    rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    rp = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    return rs, rp


def _reflect(d, Ny, Nz, E=None, coat=None):
    """Vector law of reflection for directions d (3, N) on a cylinder with normal (0, Ny, Nz); optional field E
    (3, N, complex) carried with the s/p coefficients coat(cos_i) -> (r_s, r_p); coat None: perfect conductor."""
    N = np.stack([np.zeros_like(Ny), Ny, Nz])
    dn = np.sum(d * N, axis=0)
    d2 = d - 2 * dn * N
    cos_i = np.abs(dn)
    if E is None:
        return d2, None, cos_i
    if coat is None:
        E2 = -E + 2 * np.sum(E * N, axis=0) * N
        return d2, E2, cos_i
    s = np.cross(d.T, N.T).T
    sn = np.linalg.norm(s, axis=0)
    # normal incidence: any s works; take x
    s = np.where(sn > 1e-12, s / np.where(sn > 1e-12, sn, 1), np.array([[1.0], [0.0], [0.0]]))
    p = np.cross(s.T, d.T).T
    p2 = np.cross(s.T, d2.T).T
    rs, rp = coat(cos_i)
    Es, Ep = np.sum(E * s, axis=0), np.sum(E * p, axis=0)
    return d2, rs * Es * s + rp * Ep * p2, cos_i


# ---------------------------------------------------------------------------------------------- the design
@dataclass
class ReflectiveDesign:
    name: str
    n: float
    f: float
    y_c: float
    sign: int
    theta_m: float  # rim of M1: largest ray angle in the glass for s_x = 0
    M1: Curve
    M2: Curve
    z_e: float  # flat exit face (normal to z)
    R: float = 1.0  # power reflectance per mirror (wave model); 1: ideal
    construction: dict = field(default_factory=dict)
    indices: list = None  # [source medium, between the mirrors, after the exit face]; m26.skewed() maps these

    def __post_init__(self):
        if self.indices is None:
            self.indices = [self.n, self.n, 1.0]

    @property
    def n_m(self):
        return self.indices[0]

    @property
    def z_out(self):
        """Exit-field plane (the exit face; collimated light crosses it without change)."""
        return self.z_e

    # ------------------------------------------------------------------ tracing
    def trace(self, sx, sy, y0=0.0, z0=0.0, n=None, pol=False, coat=None):
        """3D trace from a source point (x = 0, y0, z0) in the glass. sx, sy: direction sines in air units
        (k / k0); n: medium index (default: the design index). Returns a dict with exit height y and x on the
        plane z = z_e (in the glass), optical path (n x 3D length), exit direction, incidence cosines, M1
        parameter and, with pol, the field after both mirrors (TE launch: x projected transverse to the ray).
        coat: None (perfect conductor), a function cos_i -> (r_s, r_p), or a list of two (M1, M2)."""
        n = self.n if n is None else n
        sx, sy = np.broadcast_arrays(np.asarray(sx, float), np.asarray(sy, float))
        sx, sy = sx.ravel(), sy.ravel()
        d = np.stack([sx / n, sy / n, np.sqrt(np.clip(1 - (sx**2 + sy**2) / n**2, 0, None))])
        p = np.stack([np.zeros_like(sx), y0 + 0 * sx, z0 + 0 * sx])
        E = None
        if pol:
            ex = np.array([[1.0], [0.0], [0.0]]) * np.ones_like(sx)
            E = (ex - np.sum(ex * d, axis=0) * d).astype(complex)
            E /= np.linalg.norm(E, axis=0)
        L = np.zeros_like(sx)
        cos_inc, us = [], []
        coats = coat if isinstance(coat, (list, tuple)) else [coat, coat]
        for M, cj in zip((self.M1, self.M2), coats):
            dyz = np.hypot(d[1], d[2])
            u, t = M.intersect(p[1], p[2], d[1] / dyz, d[2] / dyz)
            Ny, Nz = M.normal(np.nan_to_num(u))
            step = t / dyz  # 3D length
            p = p + step * d
            L = L + step
            d, E, ci = _reflect(d, Ny, Nz, E, cj)
            cos_inc.append(ci)
            us.append(u)
        step = (self.z_e - p[2]) / d[2]
        p = p + step * d
        L = L + step
        bad = ~np.isfinite(us[0]) | ~np.isfinite(us[1])
        nan = np.where(bad, np.nan, 1.0)
        return {"y": p[1] * nan, "x": p[0] * nan, "opl": n * L * nan, "d": d * nan, "cos_i": [c * nan for c in cos_inc],
                "u1": us[0], "u": us, "E": E}

    def lensmap(self, lam, s_max, n=1201, coated=True):
        """Map for the wave model (picsim.opa.stack.RelayY), in the y-z problem with effective source index
        indices[0] (m26.skewed): g = sign (y - y_c) / f, W from the projected path k0 n_eff L_yz, t = R."""
        ne = self.indices[0]
        s_rim = ne * np.sin(self.theta_m)
        s = np.linspace(0.0, min(s_max, s_rim * (1 - 1e-9)), n)
        r = self.trace(0.0, s, n=ne)
        ok = np.isfinite(r["y"])
        if not ok.all():
            s = np.linspace(0.0, float(s[ok][-1]), n)
            r = self.trace(0.0, s, n=ne)
        g = self.sign * (r["y"] - self.y_c) / self.f
        W = (r["opl"] - r["opl"][0]) / lam
        t = np.full_like(s, self.R if coated else self.R)
        return LensMap(s, g, np.gradient(g, s), W, t, float(s[-1]))

    def back_focal_plane(self, dy0=0.5):
        """S: where a ray launched parallel to the axis from height dy0 at the facet crosses the line y = y_c after
        the system; returned in the air-equivalent frame (a plane inside the glass at z appears at
        z_e - (z_e - z) / n to the collimated beam in air)."""
        r = self.trace(np.array([0.0]), np.array([0.0]), y0=dy0)
        d = r["d"][:, 0]
        zB = self.z_e - (r["y"][0] - self.y_c) * d[2] / d[1]
        return float(self.z_e - (self.z_e - zB) / self.n)

    def chief_glass_length(self):
        """Projected glass path of the axial ray, facet -> M1 -> M2 -> exit face (equal for all rays in 2D)."""
        r = self.trace(np.array([0.0]), np.array([0.0]))
        return float(r["opl"][0] / self.n)

    def sigma_z(self, zb_air):
        """Physical z of the x slit plane, given the back focal distance zb_air of the x lens measured from the
        emitters in air: the unfolded glass path L moves the x image by L (1 - 1/n), and the fold places the exit
        face at z_e after a glass path L, so z = z_e + zb_air - L / n."""
        L = self.chief_glass_length()
        return self.z_e + zb_air - L / self.n

    # ------------------------------------------------------------------ geometry
    def light_region(self, n_rays=161):
        """Shapely polygon (coordinates (z, y)) of all light in the block: cone, intermediate beam, output beam."""
        from shapely.geometry import Polygon
        from shapely.ops import unary_union

        th = np.linspace(-self.theta_m, self.theta_m, n_rays)
        P1 = np.array(self.M1.point(th)).T
        u2, _ = self.M2.intersect(P1[:, 0], P1[:, 1], *self._dir_after_m1(th))
        P2 = np.array(self.M2.point(u2)).T
        polys = []
        for i in range(n_rays - 1):
            a, b, c, e = P1[i], P1[i + 1], P2[i + 1], P2[i]
            polys.append(Polygon([(0, 0), (a[1], a[0]), (b[1], b[0])]))
            polys.append(Polygon([(a[1], a[0]), (b[1], b[0]), (c[1], c[0]), (e[1], e[0])]).buffer(0))
            polys.append(Polygon([(e[1], e[0]), (c[1], c[0]), (self.z_e, c[0]), (self.z_e, e[0])]).buffer(0))
        return unary_union([p for p in polys if p.area > 0]).buffer(0.01).buffer(-0.01)

    def _dir_after_m1(self, th):
        P1y, P1z = self.M1.point(th)
        Ny, Nz = self.M1.normal(th)
        d = np.stack([0 * th, np.sin(th), np.cos(th)])
        d2, _, _ = _reflect(d, Ny, Nz)
        return d2[1], d2[2]

    def geometry(self, pitch=500.0, margin=10.0):
        """Block dimensions, footprint against the chip pitch, obscuration, incidence and curvature."""
        from shapely.affinity import translate
        from shapely.geometry import LineString, box

        reg = self.light_region()
        blk = reg.buffer(margin)
        overlap = sum(blk.intersection(translate(blk, yoff=k * pitch)).area for k in (1, 2))  # (z, y) coordinates
        th = np.linspace(-self.theta_m, self.theta_m, 401)
        m1 = LineString(np.array(self.M1.point(th))[::-1].T)
        th2 = th  # clear part of M2 (rays of the rim of M1)
        m2 = LineString(np.array(self.M2.point(th2))[::-1].T)
        inner = reg.buffer(-0.5)
        obsc = float(m1.intersection(inner).length + m2.intersection(inner).length)
        behind = float(reg.intersection(box(-1e4, -1e4, -0.5, 1e4)).area)
        s = np.linspace(-0.9, 0.9, 721) * np.sin(self.theta_m) * self.n
        r = self.trace(0.0, s)
        i1 = np.degrees(np.arccos(np.clip(r["cos_i"][0], 0, 1)))
        i2 = np.degrees(np.arccos(np.clip(r["cos_i"][1], 0, 1)))
        sens = np.sqrt(np.nanmean((2 * self.n * r["cos_i"][0]) ** 2 + (2 * self.n * r["cos_i"][1]) ** 2))
        zs = np.arange(0.0, self.z_e + 1, 2.0)
        prof = []
        for z in zs:
            sl = reg.intersection(box(z - 0.5, -1e4, z + 0.5, 1e4))
            prof.append(sl.area)  # ~ measure in y (1 um slab)
        y1, z1 = self.M1.point(th)
        y2, z2 = self.M2.point(th2)
        crit = np.degrees(np.arcsin(1 / self.n))
        return {"block_length_um": float(self.z_e), "region_y_extent_um": [float(reg.bounds[1]), float(reg.bounds[3])],
                "max_cross_section_um": float(max(prof)), "pitch_overlap_area_um2": float(overlap),
                "margin_um": margin, "obscured_length_um": obsc, "light_behind_facet_um2": behind,
                "M1_y_um": [float(y1.min()), float(y1.max())], "M1_z_um": [float(z1.min()), float(z1.max())],
                "M2_y_um": [float(y2.min()), float(y2.max())], "M2_z_um": [float(z2.min()), float(z2.max())],
                "incidence_M1_deg": [float(np.nanmin(i1)), float(np.nanmax(i1))],
                "incidence_M2_deg": [float(np.nanmin(i2)), float(np.nanmax(i2))],
                "tir_angle_deg": float(crit),
                "TIR_fraction_M1": float(np.nanmean(i1 > crit)), "TIR_fraction_M2": float(np.nanmean(i2 > crit)),
                "min_radius_of_curvature_um": [float(self.M1.curvature_radius(th).min()),
                                               float(self.M2.curvature_radius(th2).min())],
                "form_sensitivity_rms": float(sens),
                "exit_tilt_max_rad": float(np.nanmax(np.abs(r["d"][1] / r["d"][2]))),
                "glass_path_um": self.chief_glass_length()}

    def ok(self, pitch=500.0, margin=10.0):
        g = self.geometry(pitch, margin)
        return (g["pitch_overlap_area_um2"] < 1.0 and g["obscured_length_um"] < 0.5
                and g["light_behind_facet_um2"] < 1.0), g


def reflective_aplanat(n, f, r0, y_c, phi0_deg, sign=-1, s_max=0.9, ext=1.04, n_steps=3001, exit_gap=20.0, R=1.0):
    """Two-mirror reflective aplanat in glass of index n (2D Wassermann-Wolf construction, see the module doc).

    s_max: rim in air units for s_x = 0 (theta_m = asin(s_max / n)); the mirrors are built to ext * theta_m so that
    traced rays near the rim are not lost to the spline end, and the clear aperture of M1 is set to theta_m."""
    th_m = float(np.arcsin(s_max / n))
    phi0 = np.radians(phi0_deg)
    K = y_c * np.tan(phi0 / 2)

    def phi_of(t, r):
        ys = y_c + sign * n * f * np.sin(t)
        return 2 * np.arctan2(K - r * (1 - np.cos(t)), ys - r * np.sin(t))

    def drdt(t, r):
        return r / np.tan((phi_of(t, r) - t) / 2)

    halves = []
    for th_end in (ext * th_m, -ext * th_m):
        ths = np.linspace(0.0, th_end, n_steps)
        rs = np.empty_like(ths)
        rs[0] = r0
        for i in range(1, n_steps):
            h = ths[i] - ths[i - 1]
            t, r = ths[i - 1], rs[i - 1]
            k1 = drdt(t, r)
            k2 = drdt(t + h / 2, r + h / 2 * k1)
            k3 = drdt(t + h / 2, r + h / 2 * k2)
            k4 = drdt(t + h, r + h * k3)
            rs[i] = r + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        halves.append((ths, rs))
    th = np.concatenate([halves[1][0][::-1], halves[0][0][1:]])
    r = np.concatenate([halves[1][1][::-1], halves[0][1][1:]])
    if not np.all(np.isfinite(r)) or np.any(r <= 0):
        raise ValueError("construction failed")
    ph = phi_of(th, r)
    P1y, P1z = r * np.sin(th), r * np.cos(th)
    ys = y_c + sign * n * f * np.sin(th)
    L = (ys - P1y) / np.sin(ph)
    if np.any(L <= 0) or not np.all(np.isfinite(L)):
        raise ValueError("negative or undefined M1-M2 distance")
    P2z = P1z + L * np.cos(ph)
    M1 = Curve(th, P1y, P1z, u_lim=(-th_m, th_m))
    M2 = Curve(th, ys, P2z)
    z_e = float(max(P1z.max(), P2z.max()) + exit_gap)
    d = ReflectiveDesign(f"reflective aplanat n={n:.4f} r0={r0:g} y_c={y_c:g} phi0={phi0_deg:g} sign={sign:+d}",
                         n, f, y_c, sign, th_m, M1, M2, z_e, R=R)
    d.construction = {"n": n, "f_um": f, "r0_um": r0, "y_c_um": y_c, "phi0_deg": phi0_deg, "sign": sign,
                      "s_max": s_max, "theta_m_deg": float(np.degrees(th_m)), "F_glass_um": n * f,
                      "M1_M2_distance_um": [float(L.min()), float(L.max())], "exit_gap_um": exit_gap}
    return d


def skewed(d: ReflectiveDesign, sx):
    """Copy of d seen by light with x direction sine sx (effective indices; same convention as m26.skewed)."""
    dd = copy.copy(d)
    dd.indices = [float(np.sqrt(n**2 - sx**2)) for n in d.indices]
    return dd
