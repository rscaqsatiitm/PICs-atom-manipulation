import numpy as np

from picsim.abcd import aug_free, aug_lens, chain, trace
from picsim.megasort import Block, Layout, Params, Sim


def test_no_losses_fill_and_conservation():
    L = Layout(60, 45, 3, 1)
    P = Params(tau=1e12, eps_img=0.0, p_move=0.0, f_load=0.5, t_reload=1e9)
    s = Sim(L, Block(12, 1, 45), P, seed=2)
    n0 = s.occ.sum()
    s.run(2.0)
    assert s.occ.sum() == n0  # no losses, no reloads: atoms conserved
    assert s.defect_fraction() == 0.0  # 1/3 computation rows, reservoir 2:1 at 50% loading fills every hole


def test_engine_block_positions_cover_every_segment():
    for Mx in (1, 4):
        b = Block(16, Mx, 45)
        pos = b.positions(64, 900)
        segs = {(k, g) for k, gs in pos for g in gs}
        assert segs == {(k, g) for k in range(0, 64, 16) for g in range(20)}


def test_lenslet_array_and_common_lens_overlap_in_pupil():
    f1, fc = 3.0, 180.0
    for X in (0.0, 4.4, 8.8):
        M = chain(aug_free(f1), aug_lens(f1, X), aug_free(f1), aug_free(fc), aug_lens(fc), aug_free(fc))
        y, _ = trace(M, np.full(5, X), np.linspace(-0.3, 0.3, 5))
        assert np.ptp(y) < 1e-9 and abs(y[0]) < 1e-9  # every engine's emitter centre lands on the pupil centre


def test_heating_is_one_move_per_refill():
    L = Layout(60, 45, 3, 1)
    P = Params(tau=1e12, eps_img=0.0, p_move=0.0, f_load=0.5, t_reload=1e9, heat_slow=0.3)
    s = Sim(L, Block(12, 1, 45), P, seed=5)
    s.run(2.0)
    st = s.heat_stats()
    assert st["max_uK"] <= 0.3 + 1e-6 and 0 < st["moved_fraction"] < 1


def test_last_block_is_full_height():
    b = Block(32, 1, 45)
    ks = sorted({k for k, _ in b.positions(1155, 90)})
    assert ks[-1] == 1155 - 32 and all(k + 32 <= 1155 for k in ks)


def test_time_average_between_pre_and_post_image():
    from picsim.megasort import Layout, Params, Sim
    P = Params(tau=300.0, eps_img=1e-3, T_img=0.1, f_load=0.85)
    sim = Sim(Layout(80, 92, 4, 3), Block(16, 1, 46), P, seed=1)
    sim.run(3.0)
    d = np.array(sim.log)
    late = d[d[:, 0] >= 1.5]
    m = sim.mean_defect(1.5)
    assert late[:, 2].mean() <= late[:, 1].mean()  # the image itself removes atoms
    assert 0.5 * late[:, 2].mean() <= m <= 1.5 * late[:, 1].mean()


def test_time_average_vacuum_only():
    """No repair (no reservoir), no imaging loss: the mean over one long idle period is T/(2 tau)."""
    from picsim.megasort import Layout, Params, Sim
    P = Params(tau=10.0, eps_img=0.0, p_move=0.0, T_img=1.0, t_img=0.0, f_load=1.0, t_reload=1e9)
    sim = Sim(Layout(200, 200, 1, 1), Block(10, 1, 50), P, seed=3)
    sim.run(1.0 + 1e-9)
    m = sim.mean_defect(0.0)
    assert abs(m - (1 - np.exp(-0.5 / 10.0))) < 0.01
