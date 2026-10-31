"""1D slab-waveguide modes: analytic symmetric slab and a finite-difference multilayer solver.

Used for (i) the effective-index reduction of SiN channel guides to a 2D lateral problem,
(ii) eigenmode expansion (EME) of MMI couplers, and (iii) unit tests of the FEM solver.
Lengths in micrometres.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq


def symmetric_slab_te(n_core, n_clad, d, lam, m=0):
    """Analytic TE_m effective index of a symmetric slab of thickness d."""
    k0 = 2 * np.pi / lam
    V = k0 * d / 2 * np.sqrt(n_core**2 - n_clad**2)

    # u in (m pi/2, (m+1) pi/2)
    lo_u, hi_u = m * np.pi / 2 + 1e-12, min((m + 1) * np.pi / 2 - 1e-12, V - 1e-12)
    if lo_u >= V:
        return np.nan

    def g(u):
        w = np.sqrt(max(V**2 - u**2, 0.0))
        return u * np.tan(u - m * np.pi / 2) - w

    u = brentq(g, lo_u, hi_u)
    b = 1 - (u / V) ** 2
    return np.sqrt(n_clad**2 + b * (n_core**2 - n_clad**2))


def fd_modes(x, n, lam, pol="TE", num=None, nmin=None):
    """Finite-difference modes of a 1D index profile n(x) on a uniform grid with E=0 at the ends.

    pol='TE': E parallel to layers (scalar Helmholtz).
    pol='TM': H parallel to layers; uses the 1/eps-weighted operator (symmetrised).
    Returns (neff [descending], fields [len(x), M]) with fields normalised to sum |E|^2 dx = 1.
    """
    x = np.asarray(x, float)
    dx = x[1] - x[0]
    k0 = 2 * np.pi / lam
    eps = np.asarray(n, float) ** 2
    if pol == "TE":
        d = -2.0 / dx**2 + k0**2 * eps
        e = np.ones(len(x) - 1) / dx**2
        w, v = eigh_tridiagonal(d, e)
    else:
        # d/dx (1/eps dH/dx) + k0^2 H = beta^2 / eps H  ->  symmetrise with S = diag(sqrt(eps))
        inv_mid = 1.0 / (0.5 * (eps[1:] + eps[:-1]))
        inv_l = np.concatenate([[1 / eps[0]], inv_mid])
        inv_r = np.concatenate([inv_mid, [1 / eps[-1]]])
        # generalised form: (D) H = beta^2 diag(1/eps) H ; transform G = H/sqrt(eps)
        s = np.sqrt(eps)
        diag = (-(inv_l + inv_r) / dx**2 + k0**2) * eps
        off = inv_mid / dx**2 * s[1:] * s[:-1]
        w, v = eigh_tridiagonal(diag, off)
        v = v * s[:, None]  # back to H
    order = np.argsort(-w)
    w, v = w[order], v[:, order]
    neff = np.sqrt(np.maximum(w, 0)) / k0
    if nmin is not None:
        keep = neff > nmin
        neff, v = neff[keep], v[:, keep]
    if num is not None:
        neff, v = neff[:num], v[:, :num]
    v = v / np.sqrt(np.sum(np.abs(v) ** 2, axis=0) * dx)
    return neff, v


def averaged_profile(x, layers, background, sub=16, pol="TE"):
    """Index profile on grid x with cell-averaged permittivity.

    layers: list of (x_lo, x_hi, n); later entries override earlier ones. background: index outside.
    Returns sqrt(<eps>) per cell for TE (E tangential to interfaces) and sqrt(1/<1/eps>) for TM
    (normal E), which removes the O(dx) staircase error.
    """
    dx = x[1] - x[0]
    xs = x[:, None] + (np.arange(sub)[None, :] + 0.5) / sub * dx - dx / 2
    eps = np.full(xs.shape, complex(background) ** 2)
    for lo, hi, n in layers:
        eps = np.where((xs >= lo) & (xs < hi), complex(n) ** 2, eps)
    e = eps.mean(axis=1) if pol == "TE" else 1.0 / (1.0 / eps).mean(axis=1)
    return np.sqrt(e.real) if np.allclose(e.imag, 0) else np.sqrt(e)


def slab_neff(t, lam, n_core, n_clad_lo, n_clad_hi=None, pol="TE", npts=2001, span=None):
    """Fundamental effective index of an (a)symmetric three-layer slab of thickness t."""
    n_clad_hi = n_clad_lo if n_clad_hi is None else n_clad_hi
    span = span or (t + 6 * lam)
    x = np.linspace(-span / 2, span / 2, npts)
    n = averaged_profile(x, [(-np.inf, -t / 2, n_clad_lo), (-t / 2, t / 2, n_core), (t / 2, np.inf, n_clad_hi)],
                         n_clad_lo, pol=pol)
    ne, _ = fd_modes(x, n, lam, pol=pol, num=1)
    return float(ne[0])
