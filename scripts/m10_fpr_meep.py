"""Module 10b: 2D FDTD (MEEP) reference for the star-coupler BPM, on a reduced 8-port, 8-arm FPR.

Same geometry (FPRLayout) and the same calibrated lateral indices as the BPM, so the comparison tests the
BPM's wide-angle propagation and discretisation, not the effective-index reduction (validated separately
against vector FEM). Run with the MEEP interpreter:
    PYTHONPATH=. <meep-env>/bin/python scripts/m10_fpr_meep.py <port> <resolution>
Writes results/m10_meep_port<k>_res<r>.json with complex arm amplitudes in the BPM phase convention.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import meep as mp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from picsim.opa.fpr_bpm import narrow_mode, tilted  # noqa: E402
from picsim.opa.fpr_layout import FPRLayout  # noqa: E402

LAM = 1.064435
N_CORE, N_CLAD = 1.84228, 1.35515  # calibrated lateral indices (m10_fpr.json)


def small_layout():
    N, p = 8, 2.0
    R = N * p * p / (LAM / 1.8423)
    return FPRLayout(N, N, p, p, R, L_taper=15, L_port=25, L_arm=30, z0=30)


def run(port, res):
    lay = small_layout()
    dpml = 2.0
    X = lay.x_slab + 6.0
    zl, zm = lay.z_launch, lay.z_measure
    zmin, zmax = zl - 3.0, zm + 3.0
    sx, sy = 2 * X + 2 * dpml, (zmax - zmin) + 2 * dpml
    yc = 0.5 * (zmin + zmax)
    core = mp.Medium(index=N_CORE)
    geom = []
    for poly in lay.polygons():
        verts = [mp.Vector3(px, pz - yc) for px, pz in poly]
        geom.append(mp.Prism(verts, height=mp.inf, material=core))
    xc, al = lay.port_axis(port, zl)
    fcen = 1 / LAM
    src = mp.EigenModeSource(mp.GaussianSource(fcen, fwidth=0.05 * fcen), center=mp.Vector3(xc, zl - yc),
                             size=mp.Vector3(3.0, 0), direction=mp.NO_DIRECTION,
                             eig_kpoint=mp.Vector3(np.sin(al), np.cos(al)), eig_band=1,
                             eig_parity=mp.EVEN_Z, eig_match_freq=True)
    sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), geometry=geom, sources=[src], resolution=res,
                        boundary_layers=[mp.PML(dpml)], default_material=mp.Medium(index=N_CLAD))
    z_in = zl + 1.5
    d_in = sim.add_dft_fields([mp.Hz], fcen, 0, 1, center=mp.Vector3(0, z_in - yc), size=mp.Vector3(2 * X, 0))
    d_out = sim.add_dft_fields([mp.Hz], fcen, 0, 1, center=mp.Vector3(0, zm - yc), size=mp.Vector3(2 * X, 0))
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Hz, mp.Vector3(0, zm - yc), 1e-6))
    h_in = np.conj(sim.get_dft_array(d_in, mp.Hz, 0))  # conjugate: MEEP exp(+ikz) -> BPM exp(-ikz)
    h_out = np.conj(sim.get_dft_array(d_out, mp.Hz, 0))
    x_in = np.linspace(-X, X, len(h_in))
    x_out = np.linspace(-X, X, len(h_out))
    xi, m, ne = narrow_mode(lay.w_narrow, LAM, N_CORE, N_CLAD)
    beta = 2 * np.pi / LAM * ne

    def ov(h, x, xc_, al_):
        md = tilted(x, xc_, al_, xi, m, beta)
        nprof = np.where(np.abs(x - xc_) < lay.w_narrow / 2 / np.cos(al_), N_CORE, N_CLAD)
        w = 1 / nprof**2
        md = md / np.sqrt(np.sum(np.abs(md) ** 2 * w))
        return np.sum(np.conj(md) * w * h)

    xc_in, al_in = lay.port_axis(port, z_in)
    a_in = ov(h_in, x_in, xc_in, al_in)
    T = []
    for n in range(lay.N_arms):
        xa, aa = lay.arm_axis(n, zm)
        T.append(ov(h_out, x_out, xa, aa) / a_in)
    T = np.array(T)
    out = {"port": port, "resolution": res, "T_re": T.real.tolist(), "T_im": T.imag.tolist(),
           "capture": float(np.sum(np.abs(T) ** 2)), "N_CORE": N_CORE, "N_CLAD": N_CLAD,
           "z_in_measure": z_in, "z_out_measure": zm}
    (REPO / "results" / f"m10_meep_port{port}_res{res}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({"port": port, "res": res, "capture": out["capture"], "absT": np.abs(T).round(4).tolist()}))


if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 20)
