"""The doubled-resolution comparator of module 32 passes identical values and fails perturbed ones."""
import sys

import numpy as np

from picsim.provenance import REPO

sys.path.insert(0, str(REPO / "scripts"))
import m32_meep_convergence as m32  # noqa: E402


def test_stored_cases_resolve():
    s = m32.stored()
    assert len(s["m3"]) == len(m32.M3_CASES) and len(s["m14"]) == 2 and len(s["m10"]) == 1
    assert abs(s["m3"][m32.M3_CASES[0][0]]["directionality_c"] - 0.715) < 5e-4  # the regression case


def test_selftest_detects_a_planted_shift():
    out = m32.selftest()
    for module, r in out.items():
        assert r["identical_all_pass"], module
        assert r["perturbed_any_fail"], module


def test_tolerance_edges():
    low = {"R": 0.10, "T": 0.80, "T_within_s0.9": 0.70, "fraction_within_s0.9": 0.9, "NA_1e2": 1.0}
    inside = dict(low, R=0.1099)
    outside = dict(low, R=0.1101)
    assert all(r["verdict"] == "PASS" for r in m32.compare("m14", low, inside))
    assert any(r["verdict"] == "FAIL" for r in m32.compare("m14", low, outside))
    a = {"T_re": [1.0, 0.5], "T_im": [0.0, 0.2], "capture": 0.4}
    assert m32.fidelity(np.array([1, 0.5 + 0.2j]), np.array([1, 0.5 + 0.2j])) > 1 - 1e-12
    assert all(r["verdict"] == "PASS" for r in m32.compare("m10", a, a))


def test_cheap_m10_pair_resolves():
    lo, hi = m32.M10_CHEAP["port 0"]
    assert (lo, hi) == (15, 30)
    assert (REPO / "results" / f"m10_meep_port0_res{hi}.json").exists()
    assert 2000 < m32._cheap_predicted_seconds() < 4000  # 833 s at 10/um scaled by 1.5^3


def test_pdg_rounding_rule():
    """The rounding rule fixed before the adopted uncertainties were computed."""
    cases = {0.0134: (0.013, 3), 0.0163: (0.016, 3), 0.0108: (0.011, 3), 0.382: (0.4, 1), 0.0083: (0.008, 3),
             0.354: (0.35, 2), 0.355: (0.4, 1), 0.96: (1.0, 1), 1.08: (1.1, 1), 153.1: (150.0, -1)}
    for u, (val, dec) in cases.items():
        got = m32.pdg_round(u)
        assert got[1] == dec and abs(got[0] - val) < 1e-12, (u, got)
