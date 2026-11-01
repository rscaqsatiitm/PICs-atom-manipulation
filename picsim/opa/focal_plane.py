"""Generalised sliding-window focal-plane OPA (route E) for design optimisation.

N_ports dense input ports (pitch p_port, width w_port at the FPR) feed N_arms arms (pitch p_a, width w_arm at
the FPR) through a confocal star coupler of radius R = N_arms p_a p_port / lambda_s, so that adjacent ports land
one resolution cell lambda/(2 NA) apart at the atoms. The arms end in N_arms end-fire emitters at pitch p_out
whose facet width w_emit sets the element pattern. With N_arms > N_ports (oversampling f = N_ports/N_arms) the
ports use only the central fraction f of the arms' and emitters' angular range, which reduces scan roll-off.
The transfer matrix A (ports x arms) comes from the paraxial star-coupler model unless a BPM matrix is given.
"""
from __future__ import annotations

import numpy as np

from .. import modes as M
from .array import OPA, channel_response, spot_metrics
from .star import StarCoupler
from .star_rs import StarCouplerRS

_EMITTER_CACHE = {}


def emitter_profile(lam, width, t=0.4):
    key = (round(lam, 6), round(width, 4), t)
    if key not in _EMITTER_CACHE:
        ms = M.solve(M.XSection(cores=[(0, width)], t=t, x_left=-5 - width / 2, x_right=5 + width / 2), lam,
                     num_modes=2)
        te = M.fundamental(ms)
        x = np.linspace(-4.5 - width / 2, 4.5 + width / 2, 1201)
        ex = np.abs(M.field_on_grid(te.mode, x, np.array([t / 2]), 0))[0]
        _EMITTER_CACHE[key] = (x, ex / ex.max())
    x, ex = _EMITTER_CACHE[key]
    return lambda xx: np.interp(xx, x, ex, left=0.0, right=0.0)


class FocalPlane:
    def __init__(self, lam=1.064435, NA=0.5, N_ports=128, N_arms=128, p_port=2.0, p_a=2.0, w_port=1.7,
                 w_arm=1.7, p_out=1.5, w_emit=0.6, t=0.4, A=None, dx_atoms=0.02, model="paraxial", rp_ratio=1.0):
        self.lam, self.NA = lam, NA
        self.N_ports, self.N_arms = N_ports, N_arms
        self.opa = OPA(N=N_arms, pitch=p_out, lam=lam, element=emitter_profile(lam, w_emit, t), NA=NA)
        n_s = StarCoupler(lam, t, 2, p_a, 100.0, w_port, w_arm).n_s
        self.R = N_arms * p_a * p_port / (lam / n_s)
        self.k = np.arange(N_ports) - (N_ports - 1) / 2
        self.x = np.arange(-1.02 * self.opa.x_range, 1.02 * self.opa.x_range, dx_atoms)
        self.G = channel_response(self.opa, self.x)
        if A is None and model == "paraxial":
            sc = StarCoupler(lam, t, N_arms, p_a, self.R, w_port, w_arm)
            A = np.array([sc.arm_amplitudes(kk * p_port) for kk in self.k])[::-1]
        elif A is None:
            # exact-path model; arm phases calibrated with the two central ports (static arm trims)
            A = StarCouplerRS(lam, t, N_ports, N_arms, p_port, p_a, self.R, w_port, w_arm,
                              R_p=rp_ratio * self.R).matrix()
            c0 = N_ports // 2
            A = A * np.exp(-1j * np.angle(A[c0 - 1] + A[c0]))[None, :]
            xs = [spot_metrics(self.x, np.abs(A[k] @ self.G) ** 2)["x_peak"] for k in (N_ports // 8, 7 * N_ports // 8)]
            A = A if xs[1] > xs[0] else A[::-1]
        self.A = A
        self.dx_port = lam / (2 * NA)
        self.e_pow = float(np.sum(self.opa.element(self.opa.x) ** 2) * self.opa.dx)

    def weights(self, u, sigma, amp=1.0, phase=0.0):
        w = np.exp(-((self.k - u) ** 2) / (2 * sigma**2))
        return amp * np.exp(1j * phase) * w / np.sqrt(np.sum(w**2))

    def intensity(self, port_amps):
        c = port_amps @ self.A
        return np.abs(c @ self.G) ** 2 / (self.opa.lam * self.opa.x_scale * self.e_pow)

    def spot(self, u, sigma):
        I = self.intensity(self.weights(u, sigma))
        s = spot_metrics(self.x, I)
        nb = max(np.interp(s["x_peak"] + 4.0, self.x, I), np.interp(s["x_peak"] - 4.0, self.x, I)) / s["peak"]
        return {**s, "nb4": float(nb), "delivered": float(np.trapezoid(I, self.x))}

    def field_scan(self, sigma, x_half=60.0, n=13):
        """Peak per unit port power across the field of use (|x| <= x_half)."""
        us = np.linspace(-x_half / self.dx_port, x_half / self.dx_port, n)
        return [self.spot(u, sigma) for u in us]
