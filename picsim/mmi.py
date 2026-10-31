"""Multimode-interference (MMI) couplers by effective-index reduction and eigenmode expansion (EME).

Vertical step: the SiN film (thickness t, oxide above and below) is reduced to its TE slab index n_s.
Lateral step: a 1D index profile n(x) = n_s in the cores, n_SiO2 outside. The chip TE mode has E along x,
normal to the lateral interfaces, hence the lateral problem uses the 1D TM operator.

EME: the access-guide mode is expanded at the junction onto the guided modes of the MMI section, propagated
over length L, and projected onto the output access-guide modes. Power in unguided MMI content is counted as
radiated (excess loss). Back-reflection is neglected (abrupt junctions between guides of similar index).
Lengths in micrometres.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import materials as mat
from . import slab


@dataclass
class Lateral:
    lam: float
    t: float
    dt: float = 0.0  # film thickness error (um)
    dx: float = 0.01
    margin: float = 3.0

    def __post_init__(self):
        self.n_cl = float(mat.n_sio2(self.lam))
        self.n_s = slab.slab_neff(self.t + self.dt, self.lam, float(mat.n_sin(self.lam)), self.n_cl, npts=3001)

    def grid(self, half_width):
        X = half_width + self.margin
        return np.arange(-X, X + self.dx / 2, self.dx)

    def profile(self, x, cores):
        layers = [(xc - w / 2, xc + w / 2, self.n_s) for xc, w in cores]
        return slab.averaged_profile(x, layers, self.n_cl, pol="TM")

    def modes(self, x, cores, guided_only=True):
        n = self.profile(x, cores)
        ne, v = slab.fd_modes(x, n, self.lam, pol="TM", nmin=self.n_cl if guided_only else None)
        return ne, v


def _overlap(a, b, dx):
    """Power overlap-based amplitude coefficient <a|b> for real-normalised fields (scalar approximation)."""
    return np.sum(np.conj(a) * b) * dx


def mmi_transmission(lat: Lateral, W, L, w_in, x_in, w_out, x_outs, dW=0.0):
    """Power transmission from one input guide (width w_in at x_in) to output guides at x_outs.

    L may be a scalar or an array (the modes are solved once). dW is a uniform width error applied to all
    features (lithography bias). Returns T with shape (len(x_outs),) or (len(L), len(x_outs)).
    """
    x = lat.grid(W / 2 + dW / 2 + 1.0)
    dx = x[1] - x[0]
    _, a_in = lat.modes(x, [(x_in, w_in + dW)])
    a0 = a_in[:, 0]
    ne, V = lat.modes(x, [(0.0, W + dW)])
    c = V.T.conj() @ a0 * dx  # expansion coefficients onto guided MMI modes
    beta = 2 * np.pi / lat.lam * ne
    outs = []
    for xo in x_outs:
        _, a_out = lat.modes(x, [(xo, w_out + dW)])
        outs.append(a_out[:, 0])
    P = np.array([(V.T.conj() @ o) * dx for o in outs])  # (n_out, n_modes): <out|mode>
    Ls = np.atleast_1d(np.asarray(L, float))
    ph = np.exp(-1j * np.outer(Ls, beta - beta[0]))  # (nL, n_modes)
    amp = (ph * c[None, :]) @ P.T  # (nL, n_out)
    T = np.abs(amp) ** 2
    info = {"n_modes": len(ne), "beat_Lpi_um": float(np.pi / (beta[0] - beta[1])) if len(ne) > 1 else np.nan,
            "guided_fraction": float(np.sum(np.abs(c) ** 2))}
    return (T[0] if np.ndim(L) == 0 else T), info


def optimise_length(lat, W, w_in, x_in, w_out, x_outs, L_guess, span=0.25, n=501, dW=0.0):
    Ls = np.linspace(L_guess * (1 - span), L_guess * (1 + span), n)
    T, _ = mmi_transmission(lat, W, Ls, w_in, x_in, w_out, x_outs, dW=dW)
    i = int(np.argmax(T.sum(axis=1)))
    return float(Ls[i]), float(T[i].sum())
