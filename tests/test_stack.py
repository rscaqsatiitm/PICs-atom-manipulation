import numpy as np

from picsim.opa.focal_plane import FocalPlane, LinearSynth
from picsim.opa.stack import FacetY, Focus2D, RelayY, fwhm


def gaussian_relay(beta, NA_FAC=0.99):
    return RelayY(FacetY(1.064435, gaussian_w=1.0), beta, NA=0.5, NA_FAC=NA_FAC)


def test_y_line_conserves_power_and_matches_gaussian_imaging():
    # Gaussian facet (w = 1 um, far-field NA 0.34), pupil much larger than its image: no loss
    r = gaussian_relay(beta=3.0)
    y = np.arange(-15, 15, 0.005)
    I = r.atom_line(y)
    assert abs(np.trapezoid(I, y) - r.efficiency()) < 2e-3
    assert r.efficiency() > 0.999
    # facet far field -> S (sine condition) -> atoms: the y spot is the scaled far field, w_a = beta NA_f / NA
    w_a = 3.0 * (1.064435 / (np.pi * 1.0)) / 0.5
    assert abs(fwhm(y, I) / (w_a * np.sqrt(2 * np.log(2))) - 1) < 0.03


def test_decentre_shifts_the_facet_image_in_the_pupil():
    r = gaussian_relay(beta=2.0)
    P0 = np.abs(r.pupil())
    P1 = np.abs(r.pupil(delta=0.4))
    eta = r.eta
    c0 = np.sum(eta * P0**2) / np.sum(P0**2)
    c1 = np.sum(eta * P1**2) / np.sum(P1**2)
    assert abs((c1 - c0) - 0.4 / 2.0) < 0.01


def test_2d_field_is_separable_without_the_circle():
    fp = FocalPlane(N_ports=16, N_arms=24, model="paraxial", dx_atoms=0.02)
    ls = LinearSynth(fp, x_half=7.0, step=0.04)
    w = ls.weights([1.3])
    r = gaussian_relay(beta=1.0)
    F = Focus2D(fp, r, stride=1)
    x = np.arange(0.0, 2.6, 0.02)
    y = np.array([0.0, 0.5])
    I2 = F.intensity(w, x, y, circle=False)
    Ix = np.interp(x, fp.x, fp.intensity(w))
    Iy = r.atom_line(y)
    assert np.max(np.abs(I2 / Iy[:, None] - Ix[None, :])) / Ix.max() < 2e-3
    # the circular pupil can only remove light
    assert np.all(F.intensity(w, x, y)[0].max() <= I2[0].max() * (1 + 1e-9))


def test_vector_focus_matches_scalar_and_separates_polarisations():
    fp = FocalPlane(N_ports=16, N_arms=24, model="paraxial", dx_atoms=0.02)
    ls = LinearSynth(fp, x_half=7.0, step=0.04)
    w = ls.weights([1.3])
    r = gaussian_relay(beta=1.0)
    F = Focus2D(fp, r, stride=1)
    x = np.arange(0.0, 2.6, 0.02)
    y = np.array([0.0])
    Es = F.field(w, x, y)
    Ev = F.field_vector(w, x, y, pol="x")
    Iv = np.sum(np.abs(Ev) ** 2, axis=0)
    # same peak within the aplanatic apodisation factor (a few %) and negligible cross-polarised field on axis
    assert abs(Iv.max() / (np.abs(Es) ** 2).max() - 1) < 0.08
    assert np.abs(Ev[1]).max() / np.abs(Ev[0]).max() < 0.02
