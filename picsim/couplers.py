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


def _solve_tip(w_tip, t, lam, half, n, res_max=0.5):
    xs = M.XSection(cores=[(0, w_tip)], t=t, x_left=-half, x_right=half, clad_bottom=half, t_top=half,
                    res_core=min(0.03, w_tip / 3), res_max=res_max)
    ms = M.solve(xs, lam, num_modes=2)
    te = M.fundamental(ms)
    x = np.linspace(-half + 0.5, half - 0.5, n)
    y = np.linspace(-half + 0.5 + t / 2, half - 0.5 + t / 2, n)
    E = M.field_on_grid(te.mode, x, y, 0)
    return x, y, E, te


def tip_mode_field(w_tip, t, lam, half=7.0, n=281, res_max=0.5):
    """|Ex| of the TE0 mode at the tip on a square grid (closed window of half-width `half`), plus grid vectors and
    mode data. Near cutoff the mode tail reaches several micrometres and the 7 um default window clips it; use
    `window_convergence` to check a case (the eigen-solve cost grows quickly with the window)."""
    x, y, E, te = _solve_tip(w_tip, t, lam, half, n, res_max)
    return x, y, E, te


def window_convergence(w_tip, t, lam, mfds, halves=(7.0, 11.0, 15.0, 20.0)):
    """Fibre overlaps (dB) of one tip versus the window half-width, with the evanescent decay length
    lambda / (2 pi sqrt(n_eff^2 - n_clad^2)). Larger windows use a 1 um cladding mesh (0.5 um below 12 um)."""
    from . import materials as mat
    n_clad = float(mat.n_sio2(lam))
    rows = []
    for h in halves:
        x, y, E, te = tip_mode_field(w_tip, t, lam, half=h, n=min(int(40 * h), 700) + 1,
                                     res_max=0.5 if h < 12 else 1.0)
        dn2 = te.neff.real**2 - n_clad**2
        r = {"half_um": h, "neff": float(te.neff.real), "n_clad": n_clad,
             "decay_um": float(lam / (2 * np.pi * np.sqrt(dn2))) if dn2 > 0 else None, "aeff_um2": float(te.aeff)}
        for mfd in mfds:
            eta, _ = best_overlap(x, y, E, mfd / 2)
            r[f"dB_mfd_{mfd}"] = float(10 * np.log10(eta))
        rows.append(r)
    return rows


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
