"""Inverse-taper edge couplers: overlap of the tip mode with a Gaussian fibre mode.

The facet coupling efficiency is the scalar overlap of the dominant transverse field of the chip mode at the
taper tip with a Gaussian fibre mode of 1/e^2 intensity radius w_f = MFD/2 (flat phase, index-matched gap):
    eta = |int E G|^2 / (int |E|^2 int |G|^2).
Fresnel reflection at an unmatched facet and taper conversion loss are not included; adiabatic taper length
must be confirmed by EME or 3D FDTD. Lengths in micrometres.
"""
from __future__ import annotations

import numpy as np

from . import modes as M


def _solve_tip(w_tip, t, lam, half, n):
    xs = M.XSection(cores=[(0, w_tip)], t=t, x_left=-half, x_right=half, clad_bottom=half, t_top=half,
                    res_core=min(0.03, w_tip / 3), res_max=0.5)
    ms = M.solve(xs, lam, num_modes=2)
    te = M.fundamental(ms)
    x = np.linspace(-half + 0.5, half - 0.5, n)
    y = np.linspace(-half + 0.5 + t / 2, half - 0.5 + t / 2, n)
    E = M.field_on_grid(te.mode, x, y, 0)
    return x, y, E, te


def tip_mode_field(w_tip, t, lam, half=None, n=None, decay_lengths=7.0, half_max=35.0):
    """|Ex| of the TE0 mode at the tip on a square grid, plus grid vectors and mode data.

    half=None: the closed computational window grows until it extends `decay_lengths` evanescent decay lengths,
    lambda / (2 pi sqrt(n_eff^2 - n_clad^2)), beyond the core (at most half_max). Near cutoff the mode tail reaches
    several micrometres, and a fixed 7 um window clipped it (the HI1060 overlap of a 100 nm tip moved from -0.48 to
    -0.86 dB between 7 and 15 um). te.window_converged is False when half_max was not enough, or the mode is at or
    below the cladding index (a box mode, not a guided one).
    """
    from . import materials as mat
    n_clad = float(mat.n_sio2(lam))
    if half is not None:
        x, y, E, te = _solve_tip(w_tip, t, lam, half, n or 281)
        te.window_converged = None
        return x, y, E, te
    h = 7.0
    while True:
        x, y, E, te = _solve_tip(w_tip, t, lam, h, n or int(40 * h) + 1)
        dn2 = te.neff.real**2 - n_clad**2
        if dn2 <= 1e-7:
            te.window_converged = False
            break
        need = max(w_tip, t) / 2 + decay_lengths * lam / (2 * np.pi * np.sqrt(dn2))
        if need <= h:
            te.window_converged = True
            break
        if h >= half_max:
            te.window_converged = False
            break
        h = min(half_max, max(need, 1.4 * h))
    te.window_half_um = h
    return x, y, E, te


def gaussian(x, y, w, x0=0.0, y0=0.0):
    X, Y = np.meshgrid(x, y, indexing="xy")
    return np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / w**2)


def overlap(x, y, E, G):
    num = np.abs(np.sum(E * np.conj(G))) ** 2
    return float(num / (np.sum(np.abs(E) ** 2) * np.sum(np.abs(G) ** 2)))


def best_overlap(x, y, E, w_f):
    """Maximise over the fibre-centre height (the mode centroid sits near mid-film)."""
    I = np.abs(E) ** 2
    yc = float(np.sum(I.sum(axis=1) * y) / I.sum())
    ys = yc + np.linspace(-0.5, 0.5, 21)
    vals = [overlap(x, y, E, gaussian(x, y, w_f, 0.0, y0)) for y0 in ys]
    i = int(np.argmax(vals))
    return vals[i], float(ys[i])


def mode_size(x, y, E):
    """1/e^2 intensity full widths along x and y through the peak."""
    I = np.abs(E) ** 2
    j, i = np.unravel_index(np.argmax(I), I.shape)
    lvl = I[j, i] * np.exp(-2)
    wx = np.ptp(x[I[j, :] > lvl]) if np.any(I[j, :] > lvl) else np.nan
    wy = np.ptp(y[I[:, i] > lvl]) if np.any(I[:, i] > lvl) else np.nan
    return float(wx), float(wy)
