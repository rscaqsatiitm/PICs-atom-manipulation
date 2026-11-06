import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("m26", Path(__file__).resolve().parent.parent / "scripts" / "m26_stack_combine.py")
m26 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m26)

CY = m26.CombY()
BASE = CY.metrics(CY.pupil(CY.field_S()))


def test_fft_route_reproduces_relay():
    ref = CY.metrics(CY.Py_ref)
    assert abs(BASE["transmission"] - ref["transmission"]) < 1e-4
    assert abs(BASE["xt_4um"] / ref["xt_4um"] - 1) < 0.05
    assert CY.scale_check < 1e-2


def test_window_far_outside_changes_nothing():
    m = CY.metrics(CY.pupil(CY.field_S(20e3, 2000.0)))
    assert abs(m["transmission"] / BASE["transmission"] - 1) < 1e-4
    assert abs(m["xt_4um"] / BASE["xt_4um"] - 1) < 0.02


def test_loss_grows_with_distance():
    loss = [1 - CY.metrics(CY.pupil(CY.field_S(z, 250.0)))["transmission"] / BASE["transmission"]
            for z in (20e3, 30e3, 60e3)]
    assert loss[0] < loss[1] < loss[2] and loss[2] > 0.01
