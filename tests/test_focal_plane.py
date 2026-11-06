import numpy as np

from picsim.opa.focal_plane import FocalPlane, LinearSynth


def test_least_squares_tweezers_land_on_target_and_switch_off():
    fp = FocalPlane(N_ports=16, N_arms=24, model="paraxial", dx_atoms=0.02)
    ls = LinearSynth(fp, x_half=7.0, step=0.04)
    xt = np.array([-3.1, 2.7])
    I = fp.intensity(ls.weights(xt, phases=[0.0, 1.0]))
    for x0 in xt:
        assert abs(fp.peak_near(I, x0)["x_peak"] - x0) < 0.03
    # switch the second tweezer off: its site goes dark, the first is unchanged
    I1 = fp.intensity(ls.weights(xt[:1], phases=[0.0]))
    p0 = fp.peak_near(I, xt[0])["peak"]
    assert np.interp(xt[1], fp.x, I1) / p0 < 1e-3
    assert abs(fp.peak_near(I1, xt[0])["peak"] / p0 - 1) < 0.01
