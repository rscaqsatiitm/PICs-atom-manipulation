import numpy as np

from picsim import raytrace as rt
from picsim.opa import collimator as C

N = 1.4496


def test_fresnel_normal_incidence_and_flat_refraction():
    r, _ = rt.trace(0.0, 0.0, np.array([0.0, 0.3]), [rt.Flat(10.0)], [1.0, N], 20.0)
    assert abs(r.T[0] - (1 - ((N - 1) / (N + 1)) ** 2)) < 1e-12
    assert abs(r.dy[1] - np.sin(0.3) / N) < 1e-12  # Snell


def test_immersed_ellipse_collimates_exactly():
    d = C.ellipse_immersed(N, 175.0)
    L = d.lensmap(1.0644, 0.9)
    assert np.ptp(L.W) < 1e-6
    assert abs(L.g[1] / L.s[1] - 1) < 1e-3  # paraxial focal length


def test_weierstrass_aplanat_obeys_sine_condition():
    d = C.weierstrass_aplanat(N, 175.0, 80.0, 10.0, 250.0, n2=1.8)
    L = d.lensmap(1.0644, 0.9)
    assert L.s_max >= 0.9 - 1e-9
    assert np.max(np.abs(L.g[1:] / L.s[1:] - 1)) < 1e-4
    assert np.ptp(L.W) < 1e-4
