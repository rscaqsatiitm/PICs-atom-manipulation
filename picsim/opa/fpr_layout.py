"""Geometry of a confocal star coupler (free-propagation region, FPR) with tapered ports and arms.

Coordinates: x lateral, z along the device (um). The input arc has its vertex at (0, z0) and is centred on
O_out = (0, z0 + R); the output arc has its vertex at (0, z0 + R) and is centred on O_in = (0, z0). Port k ends
on the input arc at angle theta_k = k p_port / R and points at O_out; arm n starts on the output arc at
angle phi_n = n p_a / R and points away from O_in. Ports and arms taper linearly between w_narrow and w_wide
over L_taper next to the arcs. The slab fills the lens-shaped region between the arcs, limited to
|x| <= x_slab; outside it, and between waveguides, is cladding.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class FPRLayout:
    N_ports: int
    N_arms: int
    p_port: float
    p_a: float
    R: float
    w_wide: float = 1.7
    w_narrow: float = 0.6
    L_taper: float = 30.0
    L_port: float = 55.0  # total port length behind the arc (taper + narrow)
    L_arm: float = 70.0  # total arm length after the arc
    z0: float = 60.0
    x_slab: float | None = None

    def __post_init__(self):
        kp = np.arange(self.N_ports) - (self.N_ports - 1) / 2
        ka = np.arange(self.N_arms) - (self.N_arms - 1) / 2
        self.theta = kp * self.p_port / self.R
        self.phi = ka * self.p_a / self.R
        self.Px = self.R * np.sin(self.theta)
        self.Pz = self.z0 + self.R - self.R * np.cos(self.theta)
        self.Ax = self.R * np.sin(self.phi)
        self.Az = self.z0 + self.R * np.cos(self.phi)
        if self.x_slab is None:
            self.x_slab = max(self.N_ports * self.p_port, self.N_arms * self.p_a)  # covers port-array replicas

    # ---------------------------------------------------------------- per-z slices (for BPM)
    def _width(self, s):
        """Guide width vs distance from the arc s >= 0 (taper over L_taper)."""
        f = np.clip(s / self.L_taper, 0, 1)
        return self.w_wide + (self.w_narrow - self.w_wide) * f

    def ports_at(self, z):
        s = (self.Pz - z) / np.cos(self.theta)  # distance behind the arc along the port axis
        ok = (s >= 0) & (s <= self.L_port)
        xc = self.Px + s * np.sin(self.theta)
        hw = self._width(s) / (2 * np.cos(self.theta))
        return xc[ok], hw[ok]

    def arms_at(self, z):
        s = (z - self.Az) / np.cos(self.phi)
        ok = (s >= 0) & (s <= self.L_arm)
        xc = self.Ax + s * np.sin(self.phi)
        hw = self._width(s) / (2 * np.cos(self.phi))
        return xc[ok], hw[ok]

    def slab_at(self, z):
        d = z - self.z0
        if d < 0 or d > self.R:
            return None
        h = min(np.sqrt(max(self.R**2 - (self.R - d) ** 2, 0)), np.sqrt(max(self.R**2 - d**2, 0)), self.x_slab)
        return (-h, h) if h > 0 else None

    def intervals(self, z):
        a, b = [], []
        sl = self.slab_at(z)
        if sl is not None:
            a.append(sl[0]), b.append(sl[1])
        for xc, hw in (self.ports_at(z), self.arms_at(z)):
            a.extend(xc - hw), b.extend(xc + hw)
        a, b = np.array(a), np.array(b)
        if len(a) > 1:  # merge overlaps (guides touching the slab near the arcs)
            o = np.argsort(a)
            a, b = a[o], b[o]
            ma, mb = [a[0]], [b[0]]
            for lo, hi in zip(a[1:], b[1:]):
                if lo <= mb[-1]:
                    mb[-1] = max(mb[-1], hi)
                else:
                    ma.append(lo), mb.append(hi)
            a, b = np.array(ma), np.array(mb)
        return a, b

    def port_axis(self, k, z):
        """Centre x and direction angle (radians, from +z) of port k at plane z."""
        s = (self.Pz[k] - z) / np.cos(self.theta[k])
        return self.Px[k] + s * np.sin(self.theta[k]), -self.theta[k]

    def arm_axis(self, n, z):
        s = (z - self.Az[n]) / np.cos(self.phi[n])
        return self.Ax[n] + s * np.sin(self.phi[n]), self.phi[n]

    @property
    def z_launch(self):
        return float(self.Pz.max() - (self.L_port - 2.0) * np.cos(self.theta).min())

    @property
    def z_measure(self):
        return float(self.Az.min() + (self.L_arm - 5.0) * np.cos(self.phi).min())

    # ---------------------------------------------------------------- polygons (for FDTD)
    def polygons(self):
        """List of polygons [(x, z), ...] for the slab, ports and arms."""
        polys = []
        zz = np.linspace(self.z0, self.z0 + self.R, 400)
        left, right = [], []
        for z in zz:
            sl = self.slab_at(z)
            if sl is not None:
                left.append((sl[0], z)), right.append((sl[1], z))
        polys.append(right + left[::-1])
        for k in range(self.N_ports):
            t = self.theta[k]
            u = np.array([-np.sin(t), np.cos(t)])  # axis direction (towards the slab)
            nrm = np.array([np.cos(t), np.sin(t)])
            P = np.array([self.Px[k], self.Pz[k]])
            pts = []
            for s, w in ((self.L_port, self.w_narrow), (self.L_taper, self.w_narrow), (-1.0, self.w_wide)):
                c = P - s * u
                pts.append(c + nrm * w / 2)
            for s, w in ((-1.0, self.w_wide), (self.L_taper, self.w_narrow), (self.L_port, self.w_narrow)):
                c = P - s * u
                pts.append(c - nrm * w / 2)
            polys.append([tuple(p) for p in pts])
        for n in range(self.N_arms):
            f = self.phi[n]
            u = np.array([np.sin(f), np.cos(f)])
            nrm = np.array([np.cos(f), -np.sin(f)])
            A = np.array([self.Ax[n], self.Az[n]])
            pts = []
            for s, w in ((-1.0, self.w_wide), (self.L_taper, self.w_narrow), (self.L_arm, self.w_narrow)):
                c = A + s * u
                pts.append(c + nrm * w / 2)
            for s, w in ((self.L_arm, self.w_narrow), (self.L_taper, self.w_narrow), (-1.0, self.w_wide)):
                c = A + s * u
                pts.append(c - nrm * w / 2)
            polys.append([tuple(p) for p in pts])
        return polys
