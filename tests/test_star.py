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
