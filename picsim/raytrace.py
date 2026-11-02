"""2D sequential ray tracing in the (y, z) plane for cylindrical (1D) lenses. Lengths in um.

A surface is z = z_v + sag(y); rays travel towards +z. Refraction by the vector form of Snell's law, optical path
accumulated as n * length, s-polarised Fresnel power transmission (the chip TE field is along x, normal to the
plane of incidence). Rays are traced in parallel (numpy arrays); a ray lost by total internal reflection or by
missing a surface carries NaN.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass
class Conic:
    z_v: float
    R: float  # vertex radius of curvature (>0: centre of curvature at larger z)
    k: float = 0.0
    y_max: float = np.inf

    def sag(self, y):
        c = 1.0 / self.R
        arg = 1 - (1 + self.k) * c**2 * y**2
        return np.where(arg >= 0, c * y**2 / (1 + np.sqrt(np.clip(arg, 0, None))), np.nan)

    def dsag(self, y):
        c = 1.0 / self.R
        arg = 1 - (1 + self.k) * c**2 * y**2
        return np.where(arg > 0, c * y / np.sqrt(np.clip(arg, 1e-300, None)), np.nan)


@dataclass
class Flat:
    z_v: float
    y_max: float = np.inf

    def sag(self, y):
        return np.zeros_like(np.asarray(y, dtype=float))

    def dsag(self, y):
        return np.zeros_like(np.asarray(y, dtype=float))


@dataclass
class Tab:
    """Surface tabulated as sag(y) (symmetric in y), cubic spline in y."""
    z_v: float
    y: np.ndarray
    s: np.ndarray
    y_max: float = field(default=None)

    def __post_init__(self):
        yy = np.concatenate([-self.y[::-1][:-1], self.y]) if self.y[0] == 0 else np.concatenate([-self.y[::-1], self.y])
        ss = np.concatenate([self.s[::-1][:-1], self.s]) if self.y[0] == 0 else np.concatenate([self.s[::-1], self.s])
        self._cs = CubicSpline(yy, ss)
        self._d = self._cs.derivative()
        if self.y_max is None:
            self.y_max = float(self.y[-1])

    def sag(self, y):
        y = np.asarray(y, dtype=float)
        return np.where(np.abs(y) <= self.y_max, self._cs(y), np.nan)

    def dsag(self, y):
        y = np.asarray(y, dtype=float)
        return np.where(np.abs(y) <= self.y_max, self._d(y), np.nan)


@dataclass
class Rays:
    y: np.ndarray
    z: np.ndarray
    dy: np.ndarray
    dz: np.ndarray
    opl: np.ndarray
    T: np.ndarray  # power transmission (Fresnel, s-pol)


def launch(y0, z0, angles, n0=1.0):
    a = np.asarray(angles, dtype=float)
    one = np.ones_like(a)
    return Rays(y0 * one, z0 * one, np.sin(a), np.cos(a), 0.0 * one, one.copy()), n0


def _intersect(r: Rays, surf, n_grid=800, n_bisect=60):
    """First crossing of the ray with the surface: sign change of F(t) = z(t) - z_v - sag(y(t)) on a grid of t,
    refined by bisection (robust on steep aspheres)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t_y = np.where(r.dy != 0, (np.sign(r.dy) * min(surf.y_max, 1e5) - r.y) / r.dy, np.inf)
    t_z = (surf.z_v - r.z) / r.dz
    t_max = np.minimum(t_y, np.abs(t_z) * 3 + 2000.0)
    u = np.linspace(0, 1, n_grid)[:, None]
    T = u * t_max[None, :]
    Y = r.y[None, :] + T * r.dy[None, :]
    F = r.z[None, :] + T * r.dz[None, :] - surf.z_v - surf.sag(Y)
    F = np.where(np.isfinite(F), F, np.nan)
    pos = F >= 0
    idx = np.argmax(pos, axis=0)
    found = pos[idx, np.arange(F.shape[1])] & (idx > 0)
    j = np.arange(F.shape[1])
    lo = np.where(found, T[np.maximum(idx - 1, 0), j], np.nan)
    hi = np.where(found, T[idx, j], np.nan)
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        Fm = r.z + mid * r.dz - surf.z_v - surf.sag(r.y + mid * r.dy)
        go_hi = Fm >= 0
        hi = np.where(go_hi, mid, hi)
        lo = np.where(go_hi, lo, mid)
    return 0.5 * (lo + hi)


def refract(r: Rays, surf, n1, n2):
    """Propagate to the surface and refract from index n1 to n2. Returns new Rays and incidence angles."""
    t = _intersect(r, surf)
    y = r.y + t * r.dy
    z = r.z + t * r.dz
    m = surf.dsag(y)
    Nn = np.hypot(m, 1.0)
    Ny, Nz = -m / Nn, 1.0 / Nn
    cos_i = Ny * r.dy + Nz * r.dz
    sgn = np.sign(cos_i)
    Ny, Nz, cos_i = Ny * sgn, Nz * sgn, cos_i * sgn  # normal along the ray
    mu = n1 / n2
    sin2t = mu**2 * (1 - cos_i**2)
    tir = sin2t > 1
    cos_t = np.sqrt(np.clip(1 - sin2t, 0, None))
    dy = mu * r.dy + (cos_t - mu * cos_i) * Ny
    dz = mu * r.dz + (cos_t - mu * cos_i) * Nz
    ts = 2 * n1 * cos_i / (n1 * cos_i + n2 * cos_t)
    Tf = (n2 * cos_t) / (n1 * cos_i) * ts**2
    nan = tir | ~np.isfinite(t)
    out = Rays(np.where(nan, np.nan, y), np.where(nan, np.nan, z), np.where(nan, np.nan, dy),
               np.where(nan, np.nan, dz), r.opl + n1 * t, r.T * Tf)
    return out, np.degrees(np.arccos(np.clip(cos_i, -1, 1)))


def to_plane(r: Rays, z_p, n):
    t = (z_p - r.z) / r.dz
    return Rays(r.y + t * r.dy, z_p + 0 * r.z, r.dy, r.dz, r.opl + n * t, r.T)


def trace(y0, z0, angles, surfaces, indices, z_out):
    """Trace a fan from (y0, z0). indices: n before the first surface, ..., n after the last. Returns the rays
    on the plane z = z_out and the incidence angles (deg) at each surface."""
    r, _ = launch(y0, z0, angles, indices[0])
    inc = []
    for j, s in enumerate(surfaces):
        r, a = refract(r, s, indices[j], indices[j + 1])
        inc.append(a)
    return to_plane(r, z_out, indices[-1]), inc
