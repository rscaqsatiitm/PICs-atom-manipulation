"""Refractive indices (wavelength in micrometres). Sources in picsim.params."""
from __future__ import annotations

import numpy as np

from . import params as prm


def _sellmeier(lam_um, terms, a0=0.0):
    l2 = np.asarray(lam_um, dtype=float) ** 2
    n2 = 1.0 + a0
    for b, c in terms:
        n2 = n2 + b * l2 / (l2 - c**2)
    return np.sqrt(n2)


def n_sin(lam_um, dT=0.0):
    """LPCVD Si3N4, Luke et al. 2015 (0.31-5.5 um), plus thermo-optic shift dT [K]."""
    return _sellmeier(lam_um, [(3.0249, 0.1353406), (40314.0, 1239.842)]) + prm.DNDT_SIN.value * dT


def n_sio2(lam_um, dT=0.0):
    """Fused silica, Malitson 1965 (0.21-6.7 um), plus thermo-optic shift dT [K]."""
    return _sellmeier(lam_um, [(0.6961663, 0.0684043), (0.4079426, 0.1162414), (0.8974794, 9.896161)]) \
        + prm.DNDT_SIO2.value * dT


def n_aln(lam_um):
    """AlN ordinary index, Pastrnak & Roskovcova 1966 (0.22-5 um). Formula-1 with a0=2.1399."""
    return _sellmeier(lam_um, [(1.3786, 0.1715), (3.861, 15.03)], a0=2.1399)


def n_si(lam_um):
    """Crystalline silicon near 1 um (RECALLED constant, used only as substrate for leakage)."""
    return prm.N_SI_1064.value + 0.0 * np.asarray(lam_um)
