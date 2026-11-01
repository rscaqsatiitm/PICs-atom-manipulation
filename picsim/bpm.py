"""2D wide-angle finite-difference beam propagation (FD-BPM) for planar SiN circuits.

The vertical dimension is removed by the effective-index method: the 400 nm SiN film is replaced by its TE slab
index n_s and the etched regions by SiO2. Because the film thickness is uniform and the etch is full, the
vertical mode profile is the same in every guide and in the slab, so the lateral problem carries all of the
physics of a star coupler (as in standard AWG design). The chip TE mode has its electric field in-plane and
normal to the lateral interfaces; hence the lateral operator is the TM-type one acting on the continuous H
field,
    P H = n^2 d/dx (1/n^2 dH/dx) + k0^2 (n^2 - nr^2) H.

Propagation (field = H exp(-i k_r z), k_r = k0 nr) uses the Pade(1,1) wide-angle equation
    (1 + P/(4 k_r^2)) dH/dz = -i P/(2 k_r) H,
discretised by Crank-Nicolson:
    [1 + P/(4 k_r^2) + i dz P/(4 k_r)] H^{n+1} = [1 + P/(4 k_r^2) - i dz P/(4 k_r)] H^n.
Several independent launches are propagated at once (multiple right-hand sides). Lateral edges carry graded
absorbers. Lengths in micrometres.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded


class BPM2D:
    def __init__(self, x, lam, n_ref, absorber=20.0, k_abs=0.02):
        self.x = np.asarray(x, float)
        self.dx = self.x[1] - self.x[0]
        self.lam = lam
        self.k0 = 2 * np.pi / lam
        self.nr = n_ref
        self.kr = self.k0 * n_ref
        X = self.x
        L = X[-1] - X[0]
        s = np.clip(np.maximum((X[0] + absorber) - X, X - (X[-1] - absorber)) / absorber, 0, 1)
        self.n_imag = k_abs * s**2  # n = n' - i n''  -> decay for exp(-i k0 n z)
        del L

    def _bands(self, n2):
        """Tridiagonal coefficients of P for the TM-type operator with complex n^2 (absorbers included)."""
        dx2 = self.dx**2
        inv = 1.0 / n2
        inv_half = 0.5 * (inv[1:] + inv[:-1])
        lower = np.zeros_like(n2)
        upper = np.zeros_like(n2)
        upper[:-1] = n2[:-1] * inv_half / dx2  # coefficient of H_{i+1} in row i
        lower[1:] = n2[1:] * inv_half / dx2  # coefficient of H_{i-1} in row i
        diag = -(upper + lower) + self.k0**2 * (n2 - self.nr**2)
        # boundary rows: Dirichlet (field vanishes beyond the absorbers)
        diag[0] = -2 * n2[0] * inv[0] / dx2 + self.k0**2 * (n2[0] - self.nr**2)
        diag[-1] = -2 * n2[-1] * inv[-1] / dx2 + self.k0**2 * (n2[-1] - self.nr**2)
        return lower, diag, upper

    def step(self, H, n_real, dz):
        n = n_real - 1j * self.n_imag
        n2 = n**2
        lo, d, up = self._bands(n2)
        a = 1 / (4 * self.kr**2)
        b = 1j * dz / (4 * self.kr)
        # left operator L = 1 + (a + b) P, right operator R = 1 + (a - b) P
        ab = np.zeros((3, len(d)), dtype=complex)
        ab[0, 1:] = (a + b) * up[:-1]
        ab[1, :] = 1 + (a + b) * d
        ab[2, :-1] = (a + b) * lo[1:]
        H2 = H if H.ndim == 2 else H[:, None]
        rhs = H2 + (a - b) * (d[:, None] * H2)
        rhs[1:] += (a - b) * lo[1:, None] * H2[:-1]
        rhs[:-1] += (a - b) * up[:-1, None] * H2[1:]
        out = solve_banded((1, 1), ab, rhs, overwrite_ab=True, overwrite_b=True, check_finite=False)
        return out if H.ndim == 2 else out[:, 0]


def tm_power(H, n):
    """Power of a TM-type field, proportional to int |H|^2 / n^2 dx (per column)."""
    w = 1.0 / np.asarray(n) ** 2
    return np.sum(np.abs(H) ** 2 * (w[:, None] if np.ndim(H) == 2 else w), axis=0)


def tm_overlap(H, mode, n):
    """Amplitude coupled into a normalised mode (TM-type inner product with weight 1/n^2)."""
    w = 1.0 / np.asarray(n) ** 2
    m = mode / np.sqrt(np.sum(np.abs(mode) ** 2 * w))
    return np.conj(m * w) @ H


def coverage(x, a, b):
    """Fraction of each grid cell (centres x, uniform spacing) covered by non-overlapping intervals [a, b).

    Vectorised over intervals with a difference array for fully covered cells and explicit fractions for the
    two boundary cells of each interval.
    """
    dx = x[1] - x[0]
    n = len(x)
    lo = (np.asarray(a, float) - (x[0] - dx / 2)) / dx
    hi = (np.asarray(b, float) - (x[0] - dx / 2)) / dx
    keep = (hi > 0) & (lo < n) & (hi > lo)
    lo, hi = np.clip(lo[keep], 0, n), np.clip(hi[keep], 0, n)
    out = np.zeros(n + 1)
    il, ih = np.floor(lo).astype(int), np.floor(hi).astype(int)
    same = il == ih
    # intervals inside a single cell
    np.add.at(out, il[same], hi[same] - lo[same])
    d = ~same
    il, ih, lo_, hi_ = il[d], ih[d], lo[d], hi[d]
    np.add.at(out, il, (il + 1) - lo_)  # partial first cell
    np.add.at(out, ih, hi_ - ih)  # partial last cell
    full = np.zeros(n + 2)
    np.add.at(full, il + 1, 1.0)  # full cells il+1 ... ih-1
    np.add.at(full, ih, -1.0)
    run = np.cumsum(full)[: n + 1]
    out += np.where(run > 0, run, 0)
    return np.clip(out[:n], 0, 1)


def lateral_profile(x, cores, n_core, n_clad):
    """TM-type (harmonic-mean) index profile of rectangular cores [(centre, width), ...] on grid x."""
    f = coverage(x, [c - w / 2 for c, w in cores], [c + w / 2 for c, w in cores])
    return np.sqrt(1.0 / (f / n_core**2 + (1 - f) / n_clad**2))


def eim_neff(width, lam, n_core, n_clad, half=6.0, dx=0.01):
    from .slab import fd_modes

    x = np.arange(-half - width / 2, half + width / 2, dx)
    ne, _ = fd_modes(x, lateral_profile(x, [(0.0, width)], n_core, n_clad), lam, pol="TM", num=1)
    return float(ne[0])


def eim_coupling_length(width, gap, lam, n_core, n_clad, dx=0.01):
    from .slab import fd_modes

    xc = (width + gap) / 2
    x = np.arange(-xc - width / 2 - 6, xc + width / 2 + 6, dx)
    ne, _ = fd_modes(x, lateral_profile(x, [(-xc, width), (xc, width)], n_core, n_clad), lam, pol="TM",
                     num=2, nmin=n_clad)
    return float(lam / (2 * (ne[0] - ne[1])))


def calibrate_eim(lam, targets_neff, targets_lc, n0):
    """Fit (n_core, n_clad) of the lateral model to vector-FEM effective indices and coupling lengths.

    targets_neff: {width: neff}; targets_lc: {(width, gap): Lc}. Returns (n_core, n_clad, residuals).
    """
    from scipy.optimize import minimize

    def cost(p):
        nc, ncl = p
        e = [(eim_neff(w, lam, nc, ncl) - v) / 1e-3 for w, v in targets_neff.items()]
        e += [np.log(eim_coupling_length(w, g, lam, nc, ncl) / v) / 0.02 for (w, g), v in targets_lc.items()]
        return float(np.sum(np.square(e)))

    res = minimize(cost, n0, method="Nelder-Mead", options={"xatol": 1e-5, "fatol": 1e-4, "maxiter": 400})
    nc, ncl = res.x
    resid = {"neff": {w: eim_neff(w, lam, nc, ncl) - v for w, v in targets_neff.items()},
             "Lc_rel": {f"{w}/{g}": eim_coupling_length(w, g, lam, nc, ncl) / v - 1 for (w, g), v in targets_lc.items()}}
    return float(nc), float(ncl), resid
