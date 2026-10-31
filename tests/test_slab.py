import numpy as np
from scipy.optimize import brentq

from picsim import materials as mat
from picsim import slab


def _tm_analytic(nc, ncl, d, lam):
    k0 = 2 * np.pi / lam
    V = k0 * d / 2 * np.sqrt(nc**2 - ncl**2)
    g = lambda u: u * np.tan(u) * (ncl / nc) ** 2 - np.sqrt(V**2 - u**2)
    u = brentq(g, 1e-9, min(np.pi / 2 - 1e-9, V - 1e-12))
    beta2 = (k0 * nc) ** 2 - (2 * u / d) ** 2
    return np.sqrt(beta2) / k0


def test_fd_te_matches_analytic():
    lam = 1.0644
    nc, ncl = mat.n_sin(lam), mat.n_sio2(lam)
    for t in (0.2, 0.4):
        a = slab.symmetric_slab_te(nc, ncl, t, lam)
        b = slab.slab_neff(t, lam, nc, ncl, npts=4001)
        assert abs(a - b) < 2e-4


def test_fd_tm_matches_analytic():
    lam = 0.8946
    nc, ncl = mat.n_sin(lam), mat.n_sio2(lam)
    for t in (0.2, 0.4):
        a = _tm_analytic(nc, ncl, t, lam)
        b = slab.slab_neff(t, lam, nc, ncl, pol="TM", npts=4001)
        assert abs(a - b) < 5e-4


def test_materials_reference_values():
    # Spot checks against the Sellmeier formulas themselves (Luke 2015, Malitson 1965)
    assert abs(mat.n_sio2(1.0644) - 1.4496) < 2e-4
    assert abs(mat.n_sin(1.0644) - 2.0112) < 2e-4
