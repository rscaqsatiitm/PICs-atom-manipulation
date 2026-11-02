import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location("m15", Path(__file__).resolve().parent.parent / "scripts" / "m15_coherence.py")
m15 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m15)


def test_msequence_and_zadoff_chu_autocorrelation():
    s = m15.mseq(7)
    assert len(s) == 127
    ph = np.repeat(s * np.pi, 4).astype(float)
    assert abs(m15.coherence_residual(ph, 4 * 5) - 1 / 127) < 1e-12
    zc = np.repeat(m15.unwrap_levels(m15.zadoff_chu(127)), 4)
    assert m15.coherence_residual(zc, 4 * 9) < 1e-10


def test_gaussian_phase_noise_decorrelates():
    rng = np.random.default_rng(0)
    ph = m15.gauss_noise_phase(rng, 2**14, 8, 3.0)
    assert m15.coherence_residual(ph, 8 * 6) < 0.03
    assert m15.coherence_residual(ph, 1) > 0.5  # still coherent over a small fraction of the correlation time
