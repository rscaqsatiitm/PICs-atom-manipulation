"""Error-budget links from optical specifications to atom-array metrics.

Formulae reproduced from the proposal (Sec. 4, footnote 4) are marked USER; others are DERIVED here.
"""
from __future__ import annotations

import numpy as np
from scipy.special import j1

from . import transport as tr


def idle_error_from_extinction(ER_db, both_beams=False):
    """Rotation error of an idle qubit during a pi pulse on the target (proposal footnote 4).

    One Raman beam addressed (the other global): leaked Rabi frequency ~ sqrt(eps), error (pi^2/4) eps.
    Both beams addressed through the same switch: leaked Rabi frequency ~ eps, error (pi^2/4) eps^2.
    """
    eps = 10 ** (-np.asarray(ER_db, float) / 10)
    return (np.pi**2 / 4) * (eps**2 if both_beams else eps)


def psf_intensity(r_um, lam_um, NA, trunc=None, n=4000):
    """Focal-plane intensity I(r)/I(0) for a circular pupil with Gaussian illumination (scalar, Debye).

    trunc: ratio of the Gaussian 1/e^2 radius to the pupil radius (None = uniform pupil, Airy pattern).
    Uses the sine-condition pupil coordinate rho = sin(theta)/NA in [0, 1].
    """
    r = np.atleast_1d(np.asarray(r_um, float))
    k = 2 * np.pi / lam_um
    if trunc is None:
        v = k * NA * r
        out = np.where(v == 0, 1.0, (2 * j1(v) / np.where(v == 0, 1, v)) ** 2)
        return out
    rho = (np.arange(n) + 0.5) / n
    P = np.exp(-(rho / trunc) ** 2)  # amplitude for 1/e^2 intensity radius trunc
    from scipy.special import j0

    E = np.array([np.sum(P * j0(k * NA * rr * rho) * rho) / n for rr in r])
    E0 = np.sum(P * rho) / n
    return (E / E0) ** 2


def neighbour_gate_error(crosstalk_intensity, both_beams=False):
    return idle_error_from_extinction(-10 * np.log10(crosstalk_intensity), both_beams)


def depth_mK(P_W, w0_um, alpha_au):
    return tr.trap_depth_K(P_W, w0_um, alpha_au) * 1e3


def position_noise_psd_limit(quanta_per_s, f_trap):
    """One-sided position-noise PSD (m^2/Hz) at f_trap giving the stated heating rate."""
    return quanta_per_s / tr.heating_rate_position_noise(1.0, f_trap)


def rin_limit(gamma_per_s, f_trap):
    """One-sided fractional intensity-noise PSD (1/Hz) at 2 f_trap giving parametric rate gamma."""
    return gamma_per_s / (np.pi**2 * f_trap**2)
