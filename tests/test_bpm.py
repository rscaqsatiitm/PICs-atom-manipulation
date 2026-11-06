import numpy as np

from picsim import slab
from picsim.bpm import BPM2D, coverage, lateral_profile, tm_overlap, tm_power

LAM, NC, NCL = 1.0644, 1.8423, 1.3551


def test_coverage_fractions():
    x = np.arange(0, 2, 0.1) + 0.05
    f = coverage(x, [0.23, 1.0], [0.71, 1.05])
    assert np.allclose(f[2:8], [0.7, 1, 1, 1, 1, 0.1]) and abs(f[10] - 0.5) < 1e-12


def test_straight_guide_power_and_phase():
    x = np.arange(-10, 10, 0.025)
    n = lateral_profile(x, [(0.0, 0.6)], NC, NCL)
    ne, V = slab.fd_modes(x, n, LAM, pol="TM", num=1)
    H = V[:, 0].astype(complex)
    b = BPM2D(x, LAM, NC, absorber=3.0)
    p0 = tm_power(H, n)
    L, dz = 300.0, 0.125
    for _ in range(int(L / dz)):
        H = b.step(H, n, dz)
    assert abs(tm_power(H, n) / p0 - 1) < 2e-3
    ph = np.angle(tm_overlap(H, V[:, 0], n))
    expected = -(2 * np.pi / LAM) * (ne[0] - NC) * L
    # effective-index error of the propagator below 3e-4 for a guide 0.12 below the reference index
    dn = np.angle(np.exp(1j * (ph - expected))) / (2 * np.pi / LAM * L)
    assert abs(dn) < 3e-4


def test_coupled_guides_match_supermodes():
    x = np.arange(-8, 8, 0.025)
    xc = 1.0
    n2 = lateral_profile(x, [(-xc, 1.7), (xc, 1.7)], NC, NCL)
    ne, _ = slab.fd_modes(x, n2, LAM, pol="TM", num=2, nmin=NCL)
    Lc = LAM / (2 * (ne[0] - ne[1]))
    _, V1 = slab.fd_modes(x, lateral_profile(x, [(-xc, 1.7)], NC, NCL), LAM, pol="TM", num=1)
    _, VR = slab.fd_modes(x, lateral_profile(x, [(xc, 1.7)], NC, NCL), LAM, pol="TM", num=1)
    H = V1[:, 0].astype(complex)
    p0 = tm_power(H, n2)
    b = BPM2D(x, LAM, NC, absorber=3.0)
    dz, L = 0.125, Lc / 2
    for _ in range(int(L / dz)):
        H = b.step(H, n2, dz)
    # half a coupling length transfers half the power
    assert abs(abs(tm_overlap(H, VR[:, 0], n2)) ** 2 / p0 - 0.5) < 0.03
