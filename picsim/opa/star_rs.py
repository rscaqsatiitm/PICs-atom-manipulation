"""Star coupler by 2D Rayleigh-Sommerfeld (Huygens) integration with exact path lengths.

Arms sit on a circle of radius R about the origin O (the vertex of the port curve), at arc spacing p_a, and
point radially away from O. Ports sit on a circle of radius R_p through O, centred at (0, R_p), at arc spacing
p_port, and point at the vertex of the arm arc, C = (0, R). R_p = R is the confocal mounting and R_p = R/2 the
Rowland mounting. Each port radiates its lateral mode u(xi), sampled across its aperture, into the slab with the
2D Green function
    g(r) = sqrt(k_s / (2 pi r)) exp(-i (k_s r - pi/4)) (cos a_in + cos a_out) / 2,
and each arm takes the overlap of the slab field with its own mode across its aperture. This keeps exact path
lengths (aberrations), tilts and obliquity, which the paraxial model (star.py) drops; coupling between
neighbouring ports or arms inside the tapers is not included (see fpr_bpm.py).
Lengths in micrometres.
"""
from __future__ import annotations

import numpy as np

from .star import _mode
from ..mmi import Lateral


class StarCouplerRS:
    def __init__(self, lam, t, N_ports, N_arms, p_port, p_a, R, w_port, w_arm, R_p=None, dxi=0.05, half_ap=2.2):
        self.lam, self.N_ports, self.N_arms = lam, N_ports, N_arms
        self.R = R
        self.R_p = R if R_p is None else R_p
        lat = Lateral(lam, t)
        self.n_s = lat.n_s
        self.k_s = 2 * np.pi * self.n_s / lam
        self.xi_p, self.u_p = _mode(lat, w_port)
        self.xi_a, self.e_a = _mode(lat, w_arm)
        self.tau = np.arange(-half_ap, half_ap + dxi / 2, dxi)
        self.dxi = dxi
        self.up = np.interp(self.tau, self.xi_p, self.u_p)
        self.ea = np.interp(self.tau, self.xi_a, self.e_a)
        kp = np.arange(N_ports) - (N_ports - 1) / 2
        ka = np.arange(N_arms) - (N_arms - 1) / 2
        phi = ka * p_a / R
        self.A_c = np.stack([R * np.sin(phi), R * np.cos(phi)], 1)  # arm centres
        self.A_t = np.stack([np.cos(phi), -np.sin(phi)], 1)  # tangent across each arm aperture
        self.A_d = np.stack([np.sin(phi), np.cos(phi)], 1)  # arm axis (away from O)
        a = kp * p_port / self.R_p
        self.P_c = np.stack([self.R_p * np.sin(a), self.R_p - self.R_p * np.cos(a)], 1)
        d = np.array([0.0, R])[None, :] - self.P_c
        self.P_d = d / np.linalg.norm(d, axis=1, keepdims=True)  # port axis (towards the arm-arc vertex)
        self.P_t = np.stack([self.P_d[:, 1], -self.P_d[:, 0]], 1)

    def arm_amplitudes(self, k):
        """Amplitudes in all arms for unit power in port k (index 0 ... N_ports-1)."""
        src = self.P_c[k][None, :] + self.tau[:, None] * self.P_t[k][None, :]  # (ns, 2)
        obs = (self.A_c[:, None, :] + self.tau[None, :, None] * self.A_t[:, None, :]).reshape(-1, 2)  # (N*nt, 2)
        dv = obs[None, :, :] - src[:, None, :]  # (ns, no, 2)
        r = np.linalg.norm(dv, axis=2)
        cos_in = (dv @ self.P_d[k]) / r
        A_d_rep = np.repeat(self.A_d, len(self.tau), axis=0)
        cos_out = np.einsum("ijk,jk->ij", dv, A_d_rep) / r
        g = np.sqrt(self.k_s / (2 * np.pi * r)) * np.exp(-1j * (self.k_s * r - np.pi / 4)) * 0.5 * (cos_in + cos_out)
        E = (self.up * self.dxi) @ g  # field at the arm aperture points
        E = E.reshape(self.N_arms, len(self.tau))
        return E @ (self.ea * self.dxi)

    def matrix(self):
        return np.array([self.arm_amplitudes(k) for k in range(self.N_ports)])
