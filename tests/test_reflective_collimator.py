"""All-reflective collimator (picsim/opa/mirror_collimator.py): exact collimation at every x angle, mapping,
footprint, polarisation conventions, and that the m26 helpers leave the refractive results unchanged."""
import json
import sys

import numpy as np
import pytest

from picsim.opa import collimator as C
from picsim.opa import mirror_collimator as MC
from picsim.provenance import RESULTS

N = 1.4496261908688477  # fused silica at 1064.4 nm (materials.n_sio2)
LAM = 1.0644


@pytest.fixture(scope="module")
def design():
    return MC.reflective_aplanat(N, 175.0, 50.0, -390.0, -70.0, 1)


def test_collimates_and_obeys_sine_condition(design):
    s = np.linspace(-0.899, 0.899, 301)
    r = design.trace(0.0, s)
    W = (r["opl"] - r["opl"][150]) / LAM
    assert np.ptp(W) < 1e-6
    g = design.sign * (r["y"] - design.y_c) / design.f
    assert np.max(np.abs(g - s)) < 1e-6
    assert np.nanmax(np.abs(r["d"][1])) < 1e-8


@pytest.mark.parametrize("sx", [0.15, 0.3, 0.355])
def test_exact_for_skew_rays_with_scaled_focal_length(design, sx):
    """3D vector trace: the k_x = k0 sx component stays collimated; focal length scales by n / n_eff."""
    ne = np.sqrt(N**2 - sx**2)
    s = np.linspace(-1, 1, 201) * ne * np.sin(design.theta_m) * 0.999
    r = design.trace(sx, s)
    assert np.all(np.isfinite(r["y"]))
    W = (r["opl"] - sx * r["x"]) / LAM
    assert np.ptp(W) < 1e-6
    g = design.sign * (r["y"] - design.y_c) / design.f
    big = np.abs(s) > 0.05
    assert np.max(np.abs(g[big] / (s[big] * N / ne) - 1)) < 1e-6
    # the 2D effective-index lens map used by the wave model agrees with the 3D trace
    L = MC.skewed(design, sx).lensmap(LAM, 0.9)
    m = L.s <= s.max()
    assert np.max(np.abs(np.interp(L.s[m], s, g) - L.g[m])) < 1e-6
    assert np.max(np.abs(L.W)) < 1e-6
    # the rim of M1 fixes the aperture in the glass: s_max = n_eff sin(theta_m)
    assert L.s_max == pytest.approx(ne * np.sin(design.theta_m), rel=1e-6)


def test_footprint_fits_pitch_without_obscuration(design):
    ok, g = design.ok(500.0, 10.0)
    assert ok
    assert g["max_cross_section_um"] < 480.0
    assert g["light_behind_facet_um2"] < 1.0


def test_jones_convention_perfect_conductor():
    rng = np.random.default_rng(3)
    d = rng.normal(size=(3, 7))
    d /= np.linalg.norm(d, axis=0)
    Nv = rng.normal(size=(2, 7))
    Nv /= np.linalg.norm(Nv, axis=0)
    E = rng.normal(size=(3, 7)) + 1j * rng.normal(size=(3, 7))
    E -= np.sum(E * d, axis=0) * d
    a = MC._reflect(d, Nv[0], Nv[1], E, None)[1]
    b = MC._reflect(d, Nv[0], Nv[1], E, lambda c: MC.fresnel_metal(1.45, 1e8j, c))[1]
    assert np.max(np.abs(a - b)) < 1e-6
    # energy is conserved by total internal reflection
    c = MC._reflect(d, Nv[0], Nv[1], E, lambda c: MC.fresnel_metal(1.45, 1.0 + 0j, c))[1]
    ci = np.abs(np.sum(d * np.stack([0 * Nv[0], Nv[0], Nv[1]]), axis=0))
    tir = ci < np.cos(np.arcsin(1 / 1.45))
    assert np.allclose(np.linalg.norm(c[:, tir], axis=0), np.linalg.norm(E[:, tir], axis=0))


def test_no_cross_polarisation_without_x_angle(design):
    s = np.linspace(-0.89, 0.89, 101)
    tir = lambda c: MC.fresnel_metal(N, 1.0 + 0j, c)  # noqa: E731
    for coat in (None, [lambda c: MC.fresnel_metal(N, 0.04 + 7.6j, c), tir]):
        r = design.trace(0.0, s, pol=True, coat=coat)
        assert np.max(np.abs(r["E"][1])) < 1e-9


def test_m26_helpers_reproduce_refractive_geometry():
    sys.path.insert(0, str(RESULTS.parent / "scripts"))
    import m26_stack_combine as m26

    d = m26.fac_adopted()
    assert m26._bfp(d) == C.back_focal_plane(d)
    geo = json.loads((RESULTS / "m26_lens.json").read_text())["outputs"]["geometry"]
    sig, zb, shift = m26._sigma_from_facet_um(d)
    assert sig == pytest.approx(geo["sigma_from_facet_um"], rel=1e-9)
    assert m26._bfp(d) == pytest.approx(geo["S_from_facet_um"], rel=1e-9)
