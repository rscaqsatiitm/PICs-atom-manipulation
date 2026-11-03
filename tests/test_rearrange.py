import numpy as np

from picsim.rearrange import schedule


def test_scheduler_conserves_atoms_and_fills_block():
    rng = np.random.default_rng(3)
    occ = rng.random((24, 20)) < 0.55
    res, out, moves = schedule(occ, 12, 12, 8)
    assert out.sum() == occ.sum()
    assert res.target_filled
    assert max(res.moves_per_atom) <= 3
    assert res.time_s > 0


def test_order_preserving_match():
    from picsim.rearrange import _monotone_match
    pairs = _monotone_match(np.array([1, 5, 9]), np.array([2, 3, 4, 6, 10]))
    src = [a for a, _ in pairs]
    dst = [b for _, b in pairs]
    assert src == sorted(src) and dst == sorted(dst) and len(set(dst)) == 3


def test_repair_fills_holes_one_move_each():
    from picsim.rearrange import repair
    rng = np.random.default_rng(4)
    occ = rng.random((24, 20)) < 0.55
    _, out, _ = schedule(occ, 12, 12, 8)
    holes = [(7, 6), (10, 9), (14, 13)]  # inside the 12 x 12 block (rows 6-17, columns 4-15)
    for h in holes:
        out[h] = False
    new, rep, npass, t, pairs = repair(out, 12, 12, 8)
    assert new.sum() == out.sum()
    assert len(pairs) == 3 and rep.sum() == 3
    assert all(new[h] for h in holes)
