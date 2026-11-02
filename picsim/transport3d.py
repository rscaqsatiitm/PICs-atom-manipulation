"""3D classical transport of atoms between a static tweezer array and a moving tweezer, with handoff.

Potential U(r, t) = U_static(r) + U_move(r, t) - m g y. The static array (1075 nm) is a set of Gaussian beams of
waist w0 and Rayleigh range z_R at the site positions, all at a common (global) depth. The moving tweezer is a
tabulated, normalised 3D intensity g(x, y, z) (peak 1), scaled by a time-dependent depth and translated along
a prescribed trajectory. Both are optical dipole potentials of different wavelengths, so their intensities add
incoherently. Atoms are propagated with velocity Verlet; trajectories and depths are sampled and held at the
actuator update rate.

Lengths in m, times in s, energies in J inside the integrator.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.constants import g as G0, k as kB

from . import params as prm

M_CS = prm.CS_MASS.value


@dataclass
class StaticArray:
    sites: np.ndarray  # (n, 3) positions (m)
    w0: float  # m
    lam: float  # m
    depth: float  # J (global)

    @property
    def zR(self):
        return np.pi * self.w0**2 / self.lam

    def potential_force(self, r):
        U = np.zeros(len(r))
        F = np.zeros_like(r)
        zR, w0 = self.zR, self.w0
        for s in self.sites:
            d = r - s
            zeta = d[:, 2] / zR
            q = 1 + zeta**2
            rho2 = d[:, 0] ** 2 + d[:, 1] ** 2
            A = np.exp(-2 * rho2 / (w0**2 * q)) / q
            U -= self.depth * A
            F[:, 0] += self.depth * A * (-4 * d[:, 0] / (w0**2 * q))
            F[:, 1] += self.depth * A * (-4 * d[:, 1] / (w0**2 * q))
            F[:, 2] += self.depth * A * (2 * zeta / (zR * q)) * (2 * rho2 / (w0**2 * q) - 1)
        return U, F

    def site_radial_freq(self):
        return np.sqrt(4 * self.depth / (M_CS * self.w0**2))

    def site_axial_freq(self):
        return np.sqrt(2 * self.depth / (M_CS * self.zR**2))


@dataclass
class GridTrap:
    """Normalised intensity g on a regular grid (x, y, z axes in m), peak 1 at the origin of its frame."""
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    g: np.ndarray  # (nx, ny, nz)
    _grad: np.ndarray = field(default=None, repr=False)

    def __post_init__(self):
        gx, gy, gz = np.gradient(self.g, self.x, self.y, self.z)
        self._grad = np.stack([gx, gy, gz, self.g], axis=0)  # (4, nx, ny, nz)
        self._o = np.array([self.x[0], self.y[0], self.z[0]])
        self._d = np.array([self.x[1] - self.x[0], self.y[1] - self.y[0], self.z[1] - self.z[0]])
        self._n = np.array(self.g.shape)

    def eval(self, r):
        """Return g and grad g at positions r (n, 3) in the trap frame; zero outside the grid."""
        u = (r - self._o) / self._d
        i0 = np.floor(u).astype(int)
        inside = np.all((i0 >= 0) & (i0 < self._n - 1), axis=1)
        i0 = np.clip(i0, 0, self._n - 2)
        f = np.clip(u - i0, 0, 1)
        out = np.zeros((4, len(r)))
        for dx in (0, 1):
            wx = f[:, 0] if dx else 1 - f[:, 0]
            for dy in (0, 1):
                wy = f[:, 1] if dy else 1 - f[:, 1]
                for dz in (0, 1):
                    wz = f[:, 2] if dz else 1 - f[:, 2]
                    out += self._grad[:, i0[:, 0] + dx, i0[:, 1] + dy, i0[:, 2] + dz] * (wx * wy * wz)
        out[:, ~inside] = 0
        return out[3], out[:3].T


def thermal_start(n, T, static: StaticArray, site, rng):
    """Thermal sample (harmonic approximation) restricted to atoms bound in the real Gaussian site
    (a truncated Boltzmann distribution, as after loading and cooling). Returns positions, velocities and the
    fraction of harmonic samples that were unbound."""
    wr, wz = static.site_radial_freq(), static.site_axial_freq()
    sr, sz = np.sqrt(kB * T / (M_CS * wr**2)), np.sqrt(kB * T / (M_CS * wz**2))
    R, V, tried = [], [], 0
    while sum(len(x) for x in R) < n:
        r = np.column_stack([rng.normal(0, sr, n), rng.normal(0, sr, n), rng.normal(0, sz, n)]) + site
        v = rng.normal(0, np.sqrt(kB * T / M_CS), (n, 3))
        U, _ = static.potential_force(r)
        ok = 0.5 * M_CS * np.sum(v**2, axis=1) + U < 0
        R.append(r[ok]), V.append(v[ok])
        tried += n
    r, v = np.concatenate(R)[:n], np.concatenate(V)[:n]
    return r, v, 1 - sum(len(x) for x in R) / tried


def min_jerk_u(u):
    u = np.clip(u, 0, 1)
    return 10 * u**3 - 15 * u**4 + 6 * u**5


@dataclass
class Sequence:
    """Pick-up at src, optional lane offset in y, move in x, drop at dst; times in s, positions in m."""
    src: np.ndarray
    dst: np.ndarray
    depth_move: float  # J
    t_ramp: float
    t_move: float
    lane_dy: float = 0.0
    t_lane: float = 100e-6
    t_hold: float = 100e-6
    f_update: float | None = 200e3
    offset: np.ndarray = field(default_factory=lambda: np.zeros(3))  # moving-trap misalignment (m)

    def segments(self):
        lane = self.lane_dy != 0
        T = [self.t_ramp, self.t_lane if lane else 0.0, self.t_move, self.t_lane if lane else 0.0, self.t_ramp, self.t_hold]
        return np.cumsum([0.0] + T)

    @property
    def duration(self):
        return self.segments()[-1]

    def state(self, t):
        """Moving-trap centre and depth at time t (sample-and-hold at f_update)."""
        if self.f_update:
            t = np.floor(t * self.f_update) / self.f_update
        b = self.segments()
        pos = self.src.astype(float).copy()
        dep = 0.0
        if t < b[1]:
            dep = np.sin(0.5 * np.pi * t / self.t_ramp) ** 2
        elif t < b[4]:
            dep = 1.0
            if t < b[2]:
                pos = pos + np.array([0, self.lane_dy * min_jerk_u((t - b[1]) / max(b[2] - b[1], 1e-12)), 0])
            elif t < b[3]:
                pos = pos + np.array([0, self.lane_dy, 0])
                pos[0] += (self.dst[0] - self.src[0]) * min_jerk_u((t - b[2]) / self.t_move)
            else:
                pos = self.dst.astype(float).copy() + np.array([0, self.lane_dy * (1 - min_jerk_u((t - b[3]) / max(b[4] - b[3], 1e-12))), 0])
        elif t < b[5]:
            pos = self.dst.astype(float).copy()
            dep = np.cos(0.5 * np.pi * (t - b[4]) / self.t_ramp) ** 2
        else:
            pos = self.dst.astype(float).copy()
            dep = 0.0
        return pos + self.offset, dep * self.depth_move


def run(static: StaticArray, trap: GridTrap, seq: Sequence, T_atom=10e-6, n=1500, dt=0.25e-6, seed=0, gravity=True,
        atom_site=None):
    """Propagate atoms that start thermal in the source site (or in atom_site, for a bystander atom that should stay
    where it is while the moving tweezer passes)."""
    rng = np.random.default_rng(seed)
    start = seq.src if atom_site is None else np.asarray(atom_site, float)
    r, v, unbound0 = thermal_start(n, T_atom, static, start, rng)
    Fg = np.array([0.0, -M_CS * G0 if gravity else 0.0, 0.0])

    def force(r, t):
        Us, Fs = static.potential_force(r)
        c, dep = seq.state(t)
        if dep > 0:
            gm, dg = trap.eval(r - c)
            return Fs + dep * dg + Fg, Us - dep * gm
        return Fs + Fg, Us

    F, U = force(r, 0.0)
    E0 = 0.5 * M_CS * np.sum(v**2, axis=1) + U
    t = 0.0
    a = F / M_CS
    while t < seq.duration:
        v += 0.5 * dt * a
        r += dt * v
        t += dt
        F, U = force(r, t)
        a = F / M_CS
        v += 0.5 * dt * a
    Us, _ = static.potential_force(r)
    E1 = 0.5 * M_CS * np.sum(v**2, axis=1) + Us
    near = np.argmin(np.linalg.norm(r[:, None, :] - static.sites[None, :, :], axis=2), axis=1)
    dst_idx = int(np.argmin(np.linalg.norm(static.sites - seq.dst, axis=1)))
    src_idx = int(np.argmin(np.linalg.norm(static.sites - seq.src, axis=1)))
    bound = E1 < 0
    if atom_site is not None:  # bystander: success means staying in its own site
        dst_idx = int(np.argmin(np.linalg.norm(static.sites - start, axis=1)))
    delivered = bound & (near == dst_idx)
    left = bound & (near == src_idx)
    e_above0 = (E0 + static.depth) / kB
    e_above1 = (E1[delivered] + static.depth) / kB
    return {"initially_unbound_fraction_rejected": float(unbound0), "delivered": float(delivered.mean()), "left_behind": float(left.mean()),
            "misplaced": float((bound & (near != dst_idx) & (near != src_idx)).mean()), "lost": float((~bound).mean()),
            "E_initial_uK": float(e_above0.mean() * 1e6),
            "E_final_uK": float(e_above1.mean() * 1e6) if delivered.any() else None,
            "dE_uK": float((e_above1.mean() - e_above0.mean()) * 1e6) if delivered.any() else None}
