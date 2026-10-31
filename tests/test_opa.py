import numpy as np

from picsim.opa import noise
from picsim.opa.array import OPA, channel_response, gaussian_element, spot_metrics


def _opa(N=32, p=2.0):
    return OPA(N=N, pitch=p, lam=1.0644, element=gaussian_element(0.6), NA=0.5)


def test_energy_conservation():
    o = _opa()
    ph = o.steer_phases(0.1)
    xa = np.linspace(-1.3 * o.x_range, 1.3 * o.x_range, 8001)
    assert abs(np.trapezoid(o.atom_intensity(ph, xa), xa) - o.transmission(ph)) < 5e-3


def test_steering_position_and_range():
    o = _opa(N=64)
    assert abs(o.x_range - 64 * 1.0644 / (4 * 0.5)) < 1e-9
    for s0 in (0.05, 0.15):
        x0 = s0 * o.x_scale
        xa = np.linspace(x0 - 2, x0 + 2, 801)
        assert abs(spot_metrics(xa, o.atom_intensity(o.steer_phases(s0), xa))["x_peak"] - x0) < 0.01


def test_array_factor_first_null():
    o = _opa(N=32)
    s_null = o.lam / (o.N * o.pitch)
    af = np.abs(o.array_factor(o.steer_phases(0.0), np.array([0.0, s_null])))
    assert af[1] / af[0] < 1e-6


def test_channel_response_matches_intensity():
    o = _opa(N=16)
    ph = np.random.default_rng(0).uniform(0, 2 * np.pi, 16)
    xa = np.linspace(-5, 5, 11)
    G = channel_response(o, xa)
    p_em = o._powers(ph)[0]
    I = np.abs((o.apod * np.exp(1j * ph)) @ G) ** 2 / (o.lam * o.x_scale * p_em)
    assert np.allclose(I, o.atom_intensity(ph, xa))


def test_pointing_noise_matches_analytic():
    o = _opa(N=64)
    dx, _ = noise.monte_carlo(o, 0.1, 0.1, trials=150)
    ana = noise.analytic_position_rms(0.1, 64, 1.0644, 0.5)
    assert abs(dx.std() / ana - 1) < 0.2
