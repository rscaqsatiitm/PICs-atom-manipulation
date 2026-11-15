"""Helpers of modules 29 (concurrency), 30 (planar tile) and 31 (addressing)."""
import sys

import numpy as np

from picsim.provenance import RESULTS

sys.path.insert(0, str(RESULTS.parent / "scripts"))
import m29_concurrency as m29  # noqa: E402
import m30_planar_tile as m30  # noqa: E402
import m31_addressing as m31  # noqa: E402


def test_zone_holds_at_least_its_sites():
    for N in (3.9e3, 1.5625e4, 1e6):
        R, C = m29.zone_dims(N)
        comp = sum(1 for r in range(R) if r % m29.B < m29.BC) * C
        assert N <= comp < 1.02 * N + C


def test_rectangle_bounds_on_known_patterns():
    assert m31.rect_bounds(np.ones((6, 6), int)) == (1, 1)  # all on: one rectangle
    lo, up = m31.rect_bounds(np.eye(5, dtype=int))
    assert lo == up == 5  # identity: rank n
    M = np.zeros((6, 6), int)
    M[:2, :3] = 1
    M[3:, 3:] = 1
    assert m31.rect_bounds(M) == (2, 2)
    assert m31.rect_bounds(np.zeros((4, 4), int)) == (0, 0)


def test_rectangle_bounds_bracket_random_patterns():
    rng = np.random.default_rng(0)
    for _ in range(30):
        M = (rng.random((7, 7)) < 0.4).astype(int)
        lo, up = m31.rect_bounds(M)
        assert lo <= up <= M.any(axis=1).sum()


def test_third_order_products():
    assert m31.imd3_hits([3, 5]) == {1, 7}
    assert m31.imd3_hits([0, 1, 2]) == {-2, -1, 3, 4}  # 2*1-2 = 0 and 2*1-0 = 2 are lit sites


def test_y_imaging_ideal_and_tilt():
    u = np.linspace(-3000.0, 3000.0, 4001)
    A = np.exp(-(u / 1200.0) ** 2)
    y, I0, tr0 = m30._line(A, u)
    assert abs(y[np.argmax(I0)]) < 0.02 and 0 < tr0 <= 1.0
    assert abs(np.sum(I0) * (y[1] - y[0]) - tr0) < 1e-9  # normalisation: integral = transmitted fraction
    peaks = [m30._line(A, u, dtheta=t)[1].max() for t in (0.0, 3e-5, 6e-5)]
    assert peaks[0] > peaks[1] > peaks[2]
