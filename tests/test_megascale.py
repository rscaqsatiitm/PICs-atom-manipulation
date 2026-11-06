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
