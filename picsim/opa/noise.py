"""Effect of channel phase and amplitude errors on tweezer position and intensity.

Analytic reference (DERIVED): for independent phase errors of rms sigma on N equally spaced channels, the
least-squares tilt of the phase front has rms sigma * sqrt(12 / (N (N^2 - 1))) per channel. The main lobe
shifts by d s = (lambda / (2 pi p)) * tilt, hence at the atoms
    delta x_rms = sigma * sqrt(12 / (N (N^2 - 1))) * N lambda / (4 pi NA F_fill)
               ~ sigma * sqrt(12 / N) * lambda / (4 pi NA F_fill)        (N >> 1).
"""
from __future__ import annotations

import numpy as np

from .array import OPA, spot_metrics


def analytic_position_rms(sigma, N, lam, NA, fill=1.0):
    return sigma * np.sqrt(12.0 / (N * (N**2 - 1))) * N * lam / (4 * np.pi * NA * fill)


def monte_carlo(opa: OPA, s0, sigma_phase, sigma_amp=0.0, trials=200, seed=1, half_window=None):
    """Monte Carlo of the steered spot under random static channel errors.

    Returns arrays of peak position shift (um) and relative peak-intensity change.
    """
    rng = np.random.default_rng(seed)
    base = opa.steer_phases(s0)
    x0 = s0 * opa.x_scale
    hw = half_window or 1.5 * opa.spot_pitch
    xa = np.linspace(x0 - hw, x0 + hw, 241)
    ref = spot_metrics(xa, opa.atom_intensity(base, xa))
    apod0 = opa.apod.copy()
    dx, dI = [], []
    try:
        for _ in range(trials):
            opa.apod = apod0 * (1 + sigma_amp * rng.standard_normal(opa.N))
            ph = base + sigma_phase * rng.standard_normal(opa.N)
            m = spot_metrics(xa, opa.atom_intensity(ph, xa))
            dx.append(m["x_peak"] - ref["x_peak"])
            dI.append(m["peak"] / ref["peak"] - 1)
    finally:
        opa.apod = apod0
    return np.array(dx), np.array(dI)


def quantisation_sigma(bits, span=2 * np.pi):
    """rms phase error of uniform quantisation over 'span' with 'bits' bits."""
    return span / 2**bits / np.sqrt(12)
