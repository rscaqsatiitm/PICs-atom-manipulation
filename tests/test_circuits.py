import numpy as np

from picsim import circuits as c
from picsim import transport as tr


def test_mzi_closed_form_matches_matrix():
    for k1, k2, a1, phi in ((0.5, 0.5, 1.0, 0.3), (0.47, 0.52, 0.98, 2.0)):
        assert np.isclose(abs(c.mzi(phi, k1, k2, a1)[1, 0]) ** 2, c.mzi_cross_T(phi, k1, k2, a1))


def test_extinction_vs_imbalance():
    # power imbalance of 0.05 dB in one coupler -> ~50.8 dB (DERIVED estimate)
    a = 10 ** (-0.05 / 20)
    on, off = c.mzi_on_off(1 / (1 + a**2), 0.5)
    assert abs(c.extinction_db(on, off) - 50.8) < 0.1


def test_balanced_phase_error():
    on, off = c.mzi_on_off(0.5, 0.5, dphi_off=0.01)
    assert abs(c.extinction_db(on, off) + 20 * np.log10(0.005)) < 0.01


def test_staircase_matches_time_integration():
    f = 77e3
    t, x = tr.trajectory(100e-6, 4e-6, f_update=100e3, tau=0.4e-6, dt=2e-9)
    assert abs(tr.n_bar(t, x, f) / tr.n_bar_staircase(100e-6, 4e-6, 100e3, f, 0.4e-6) - 1) < 0.01


def test_trap_depth():
    # 6.6 mW, w0 = 1.04 um, 1160 a.u. -> ~1.01 mK (DERIVED)
    assert abs(tr.trap_depth_K(6.6e-3, 1.04, 1160) * 1e3 - 1.014) < 0.005
