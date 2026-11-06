import numpy as np

from picsim import materials as mat
from picsim import modes as M
from picsim import slab


def test_wide_guide_approaches_slab():
    """A 6 um wide, 400 nm thick guide: TE0 lies just below the TE slab index (lateral correction small)."""
    lam = 1.0644
    n_s = slab.symmetric_slab_te(mat.n_sin(lam), mat.n_sio2(lam), 0.4, lam)
    xs = M.XSection(cores=[(0, 6.0)], t=0.4, x_left=-5, x_right=5, res_core=0.05, res_max=0.5)
    te = M.fundamental(M.solve(xs, lam, num_modes=2))
    assert 0 < n_s - te.neff.real < 0.01


def test_te_tm_degenerate_for_square_core():
    lam = 1.0644
    ms = M.solve(M.XSection(cores=[(0, 0.4)], t=0.4), lam, num_modes=2)
    assert abs(ms[0].neff.real - ms[1].neff.real) < 2e-3
