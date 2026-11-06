import numpy as np

from picsim.opa.star import StarCoupler, design_radius


def _sc(w_arm=1.7):
    R = design_radius(1.0644, 0.4, 32, 2.0, 0.6, w_arm, edge=np.exp(-2))
    return StarCoupler(1.0644, 0.4, 32, 2.0, R, 0.6, w_arm)


def test_power_conserved_on_arc():
    sc = _sc()
    p = np.sum(np.abs(sc.field(0.0)) ** 2) * sc.dy
    assert 0.93 < p <= 1.0 + 1e-6  # edge at 1/e^2: ~95% inside the arm window


def test_capture_bounded_and_falls_with_tilt():
    sc = _sc()
    c0, c1 = sc.capture(0.0), sc.capture(10.0)
    assert 0 < c1 < c0 <= 1.0


def test_linear_phase_matches_geometry():
    sc = _sc()
    d = 3.0
    ph = np.unwrap(np.angle(sc.arm_amplitudes(d)))
    slope = np.polyfit(np.arange(sc.N), ph, 1)[0]
    assert abs(slope + sc.k_s * d * sc.p_a / sc.R) < 1e-3


def test_design_radius_sets_edge():
    R = design_radius(1.0644, 0.4, 32, 2.0, 0.6, 1.7, edge=np.exp(-2))
    assert abs(StarCoupler(1.0644, 0.4, 32, 2.0, R, 0.6, 1.7).envelope_edge() - np.exp(-2)) < 0.03


def test_dense_ports_map_to_resolution_cells():
    """With R = N p_a p_port / lambda_s, adjacent ports land one resolution cell lambda/(2 NA) apart."""
    from picsim.opa.array import OPA, channel_response, gaussian_element, spot_metrics

    lam, n, p_a, p_port, p_out = 1.0644, 16, 2.0, 2.0, 1.5
    sc0 = StarCoupler(lam, 0.4, 2, p_a, 100.0, 1.7, 1.7)
    R = n * p_a * p_port / (lam / sc0.n_s)
    sc = StarCoupler(lam, 0.4, n, p_a, R, 1.7, 1.7)
    opa = OPA(N=n, pitch=p_out, lam=lam, element=gaussian_element(0.37), NA=0.5)
    x = np.arange(-opa.x_range, opa.x_range, 0.01)
    G = channel_response(opa, x)
    xs = [spot_metrics(x, np.abs(sc.arm_amplitudes(k * p_port) @ G) ** 2)["x_peak"] for k in (-1, 0, 1)]
    assert np.allclose(np.abs(np.diff(xs)), lam / (2 * 0.5), atol=0.01)
