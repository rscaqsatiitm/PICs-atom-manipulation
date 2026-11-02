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
