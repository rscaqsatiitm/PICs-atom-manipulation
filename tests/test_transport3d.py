import numpy as np
from scipy.constants import k as kB

from picsim.transport3d import GridTrap, Sequence, StaticArray, run


def gaussian_grid():
    x = np.linspace(-3e-6, 3e-6, 121)
    y = np.linspace(-3e-6, 3e-6, 121)
    z = np.linspace(-8e-6, 8e-6, 65)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    w, zr = 1.2e-6, 4e-6
    g = np.exp(-2 * (X**2 + Y**2) / w**2 - (Z / zr) ** 2)
    return GridTrap(x, y, z, g), w


def test_grid_gradient_matches_analytic():
    tr, w = gaussian_grid()
    r = np.array([[0.4e-6, -0.3e-6, 1.0e-6]])
    g, dg = tr.eval(r)
    ga = np.exp(-2 * (0.4e-6**2 + 0.3e-6**2) / w**2 - (1e-6 / 4e-6) ** 2)
    assert abs(g[0] / ga - 1) < 2e-3
    assert abs(dg[0, 0] / (ga * -4 * 0.4e-6 / w**2) - 1) < 2e-2


def test_static_hold_keeps_atoms_and_energy():
    """Static array alone (moving tweezer off): atoms stay in their site and the integrator conserves energy."""
    st = StaticArray(np.array([[0.0, 0.0, 0.0], [4e-6, 0.0, 0.0]]), 0.85e-6, 1.075e-6, 0.2e-3 * kB)
    tr, _ = gaussian_grid()
    seq = Sequence(src=np.zeros(3), dst=np.zeros(3), depth_move=0.0, t_ramp=50e-6, t_move=50e-6, t_hold=50e-6)
    r = run(st, tr, seq, T_atom=10e-6, n=300, gravity=False)
    assert r["delivered"] == 1.0
    assert abs(r["dE_uK"]) < 0.3
