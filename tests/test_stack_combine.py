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


def test_effective_index_matches_vector_snell():
    """One cylindrical surface (normal in the y-z plane): 3D vector refraction of a ray with x direction sine sx
    equals the 2D trace with indices sqrt(n^2 - sx^2), in projected direction and in optical path."""
    import numpy as np

    from picsim import raytrace as rt
    n1, n2, sx, R = 1.45, 1.0, 0.3, 80.0
    surf = rt.Conic(z_v=60.0, R=-R, k=0.0, y_max=0.99 * R)
    for sy in (0.2, 0.5, 0.8):
        # 3D: unit direction d in medium n1 with n1 d_x = sx, n1 d_y = sy
        dx, dy = sx / n1, sy / n1
        dz = np.sqrt(1 - dx**2 - dy**2)
        # projected 2D trace with effective indices
        e1, e2 = np.sqrt(n1**2 - sx**2), np.sqrt(n2**2 - sx**2)
        th = np.arcsin(sy / e1)
        r2, _ = rt.trace(0.0, 0.0, np.array([th]), [surf], [e1, e2], 200.0)
        # 3D refraction at the same intersection point (the projected path is the same line)
        r1, _ = rt.launch(0.0, 0.0, np.array([th]), e1)
        r1, _ = rt.refract(r1, surf, e1, e2)
        m = surf.dsag(r1.y[0])
        N = np.array([0.0, -m, 1.0]) / np.hypot(m, 1.0)
        d = np.array([dx, dy, dz])
        ci = d @ N
        mu = n1 / n2
        ct = np.sqrt(1 - mu**2 * (1 - ci**2))
        d2 = mu * d + (ct - mu * ci) * N
        assert abs(n2 * d2[0] - sx) < 1e-12  # k_x conserved
        p = n2 * d2[1:] / np.hypot(*(n2 * d2[1:]))
        assert abs(p[0] - r2.dy[0]) < 1e-9 and abs(p[1] - r2.dz[0]) < 1e-9
        # optical path to the surface: the 3D phase k0 n1 L minus the k_x x part equals the projected k0 e1 L_yz
        L_yz = np.hypot(r1.y[0], r1.z[0])  # projected path from the source at the origin
        L3 = L_yz / np.sqrt(1 - dx**2)
        assert abs((n1 * L3 - sx * L3 * dx) - r1.opl[0]) < 1e-9


def test_skew_zero_reproduces_lensmap():
    import numpy as np
    d = m26.fac_adopted()
    L0 = d.lensmap(m26.LAM, m26.NA_FAC, coated=True)
    L1 = m26.skewed(d, 0.0).lensmap(m26.LAM, m26.NA_FAC, coated=True)
    assert np.allclose(L0.g, L1.g) and np.allclose(L0.W, L1.W)


def test_lens_baseline_matches_relay():
    cy = m26.CombY(design=m26.fac_adopted())
    b = cy.metrics(cy.pupil(cy.field_S()))
    r = cy.metrics(cy.Py_ref)
    assert abs(b["xt_4um"] / r["xt_4um"] - 1) < 0.06
    assert abs(b["transmission"] / r["transmission"] - 1) < 1e-3
