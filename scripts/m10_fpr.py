"""Module 10: beam-propagation design of the dense-port star coupler of route E (docs/independent_steering.md).

Parts:
  calib    lateral effective-index model calibrated to vector-FEM (femwell) indices and coupling lengths
  small    BPM of a reduced 8 x 8 FPR (same geometry rules) with grid convergence, for the FDTD comparison
  compare  BPM vs 2D FDTD (MEEP, scripts/m10_fpr_meep.py) on the reduced FPR
  full     BPM of the 128-port, 128-arm FPR (R = 886 um) for calibrations A and B
  route    route-E performance with the BPM transfer matrix (window study, independence, power)

Outputs: results/m10_*.json, results/m10_T_full_*.npz, figures.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from picsim import modes as M
from picsim import params as prm
from picsim.bpm import calibrate_eim, eim_coupling_length, eim_neff
from picsim.mmi import Lateral
from picsim.opa.array import spot_metrics
from picsim.opa.fpr_bpm import transfer_matrix
from picsim.opa.fpr_layout import FPRLayout
from picsim.plotting import fig
from picsim.provenance import RESULTS, P, save_result

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m9_independent as m9  # noqa: E402

SCRIPT = "scripts/m10_fpr.py"
LAM = prm.LAM_TRANSPORT.value
T_FILM = prm.T_SIN.value


# ------------------------------------------------------------------ calibration
def femwell_targets():
    neff = {}
    for w in (0.6, 1.0, 1.7):
        neff[w] = M.fundamental(M.solve(M.XSection(cores=[(0, w)], t=T_FILM, x_left=-4 - w / 2, x_right=4 + w / 2),
                                        LAM, num_modes=2)).neff.real
    lc = {}
    for w, g in ((1.7, 0.3), (1.4, 0.49), (0.6, 1.0)):
        xc = (w + g) / 2
        ms = M.solve(M.XSection(cores=[(-xc, w), (xc, w)], t=T_FILM, x_left=-xc - w / 2 - 4, x_right=xc + w / 2 + 4),
                     LAM, num_modes=8)
        te = [m for m in ms if m.te_fraction > 0.5][:2]
        lc[(w, g)] = LAM / (2 * (te[0].neff.real - te[1].neff.real))
    return neff, lc


def calib():
    neff, lc = femwell_targets()
    lat = Lateral(LAM, T_FILM)
    out = {"femwell_neff": {str(k): v for k, v in neff.items()},
           "femwell_Lc_um": {f"{w}/{g}": v for (w, g), v in lc.items()},
           "uncalibrated": {"n_core": lat.n_s, "n_clad": lat.n_cl,
                            "neff": {str(w): eim_neff(w, LAM, lat.n_s, lat.n_cl) for w in neff},
                            "Lc_um": {f"{w}/{g}": eim_coupling_length(w, g, LAM, lat.n_s, lat.n_cl) for (w, g) in lc}}}
    for name, tl in (("A", {(1.7, 0.3): lc[(1.7, 0.3)]}),
                     ("B", {(1.7, 0.3): lc[(1.7, 0.3)], (0.6, 1.0): lc[(0.6, 1.0)]})):
        nc, ncl, res = calibrate_eim(LAM, neff, tl, [lat.n_s, lat.n_cl])
        out[name] = {"n_core": nc, "n_clad": ncl, "fit_targets": [f"{w}/{g}" for (w, g) in tl],
                     "neff_resid": {str(k): v for k, v in res["neff"].items()},
                     "Lc_um": {f"{w}/{g}": eim_coupling_length(w, g, LAM, nc, ncl) for (w, g) in lc}}
    save_result("m10_eim_calibration", script=SCRIPT,
                inputs={"femwell": "vector FEM, order 2 (picsim.modes), 400 nm SiN in SiO2, 1064.4 nm, TE",
                        "lateral_model": "TM-type operator, harmonic-mean averaging (picsim.bpm)"},
                outputs=out,
                notes="A: fitted to indices of 0.6/1.0/1.7 um guides and the 1.7 um / 0.3 um-gap coupling length "
                      "(dense ports). B: additionally the 0.6 um / 1.0 um-gap coupling length (compromise). The "
                      "spread between A and B is carried as model uncertainty.")
    return out


def cal(name="A"):
    d = json.loads((RESULTS / "m10_eim_calibration.json").read_text())["outputs"][name]
    return d["n_core"], d["n_clad"]


# ------------------------------------------------------------------ small FPR and FDTD comparison
def small_layout():
    N, p = 8, 2.0
    R = N * p * p / (LAM / 1.8423)
    return FPRLayout(N, N, p, p, R, L_taper=15, L_port=25, L_arm=30, z0=30)


def small():
    nc, ncl = cal("A")
    lay = small_layout()
    runs = {}
    for dx, dz in ((0.04, 0.25), (0.02, 0.25), (0.04, 0.125), (0.02, 0.125)):
        T, _ = transfer_matrix(lay, LAM, nc, ncl, dx=dx, dz=dz)
        runs[f"dx{dx}_dz{dz}"] = {"T_re": T.real.tolist(), "T_im": T.imag.tolist(),
                                  "capture": np.sum(np.abs(T) ** 2, axis=1).tolist()}
    save_result("m10_bpm_small", script=SCRIPT, inputs={"N": 8, "R_um": lay.R, "calibration": "A"}, outputs=runs)
    return runs


def compare():
    b = json.loads((RESULTS / "m10_bpm_small.json").read_text())["outputs"]["dx0.02_dz0.125"]
    Tb = np.array(b["T_re"]) + 1j * np.array(b["T_im"])
    rows = []
    for f in sorted(RESULTS.glob("m10_meep_port*_res*.json")):
        d = json.loads(f.read_text())
        k = d["port"]
        Tm = np.array(d["T_re"]) + 1j * np.array(d["T_im"])
        tb = Tb[k]
        # remove a common phase (reference planes differ by a constant), compare amplitudes and phases
        c = np.vdot(tm := Tm, tb)
        ph0 = np.angle(c)
        dphi = np.angle(tb * np.exp(-1j * ph0) * np.conj(tm))
        w = np.abs(tb) ** 2 / np.sum(np.abs(tb) ** 2)
        rows.append({"port": k, "res": d["resolution"], "capture_meep": d["capture"],
                     "capture_bpm": float(np.sum(np.abs(tb) ** 2)),
                     "amp_rel_err_rms": float(np.sqrt(np.sum(w * (np.abs(tb) / np.abs(tm) - 1) ** 2))),
                     "phase_err_rms_rad": float(np.sqrt(np.sum(w * dphi**2))),
                     "fidelity": float(np.abs(np.vdot(tm, tb)) ** 2 / (np.vdot(tm, tm).real * np.vdot(tb, tb).real))})
    save_result("m10_bpm_vs_fdtd", script=SCRIPT,
                inputs={"bpm": "dx 0.02 um, dz 0.125 um", "fdtd": "MEEP 1.35, 2D, Hz, resolution as listed"},
                outputs={"rows": rows},
                notes="fidelity = |<T_fdtd|T_bpm>|^2 / (|T_fdtd|^2 |T_bpm|^2): 1 means identical arm amplitude and "
                      "phase distributions up to a constant.")
    return rows


# ------------------------------------------------------------------ full 128 x 128 FPR
def full_layout():
    N, p = 128, 2.0
    R = N * p * p / (LAM / 1.8423)  # 886 um, one port per resolvable direction (as m9)
    return FPRLayout(N, N, p, p, R, L_taper=30, L_port=55, L_arm=70, z0=60)


def full(calname="A"):
    nc, ncl = cal(calname)
    lay = full_layout()
    t0 = time.time()
    T, info = transfer_matrix(lay, LAM, nc, ncl, dx=0.04, dz=0.125, progress=True)
    np.savez_compressed(RESULTS / f"m10_T_full_{calname}.npz", T=T)
    cap = np.sum(np.abs(T) ** 2, axis=1)
    # arm-phase calibration with the two central ports (flat phase for the centre of the field)
    c = np.exp(-1j * np.angle(T[63] + T[64]))
    Tc = T * c[None, :]
    resid, slope = [], []
    n = np.arange(lay.N_arms)
    for k in range(lay.N_ports):
        ph = np.unwrap(np.angle(Tc[k]))
        w = np.abs(Tc[k]) ** 2
        pf = np.polyfit(n, ph, 1, w=np.sqrt(w))
        r = ph - np.polyval(pf, n)
        slope.append(pf[0])
        resid.append(float(np.sqrt(np.sum(w * r**2) / np.sum(w))))
    out = {"calibration": calname, "n_core": nc, "n_clad": ncl, "R_um": lay.R, "runtime_s": time.time() - t0,
           "capture": cap.tolist(), "capture_centre": float(cap[63:65].mean()), "capture_edge": float(cap[[0, -1]].mean()),
           "phase_slope_per_port": float(np.polyfit(np.arange(lay.N_ports), slope, 1)[0]),
           "phase_residual_rms_rad": resid, "phase_residual_max_rad": float(max(resid)),
           "strehl_like_min": float(np.exp(-max(resid) ** 2))}
    save_result(f"m10_fpr_full_{calname}", script=SCRIPT,
                inputs={"layout": {"N": 128, "p_port": 2.0, "p_a": 2.0, "R": lay.R, "w_wide": 1.7, "w_narrow": 0.6,
                                   "L_taper": 30.0}, "calibration": calname, "dx_um": 0.04, "dz_um": 0.125},
                outputs=out,
                notes="T saved in m10_T_full_<cal>.npz. Arm phases calibrated with the two central ports; the "
                      "residual is the rms deviation of each port's arm phases from a straight line (aberration).")
    return out


# ------------------------------------------------------------------ route E with the BPM matrix
class BPMFocalPlane(m9.FocalPlaneOPA):
    def __init__(self, T):
        super().__init__()
        c = np.exp(-1j * np.angle(T[63] + T[64]))
        A = T * c[None, :]
        # orient: port index increasing -> +x at the atoms (same convention as the paraxial model)
        xs = [spot_metrics(self.x, np.abs(A[k] @ self.G) ** 2)["x_peak"] for k in (10, 117)]
        self.A = A if xs[1] > xs[0] else A[::-1]


def route(calname="A"):
    T = np.load(RESULTS / f"m10_T_full_{calname}.npz")["T"]
    fp = BPMFocalPlane(T)
    ports_x = [spot_metrics(fp.x, np.abs(fp.A[k] @ fp.G) ** 2)["x_peak"] for k in range(0, 128, 8)]
    slope = float(np.polyfit(np.arange(0, 128, 8), ports_x, 1)[0])
    rows, field = m9.window_study(fp)
    calt = m9.calibrated(fp, 0.9)
    ind = m9.independence(fp, cal=calt)
    us = np.arange(-60, 61, 1.0)
    ind["power_factor_mean_over_field"] = float(np.mean([calt[1](u) ** 2 for u in us]))
    ind["power_factor_edge"] = float(calt[1](60.0) ** 2)
    noise = m9.amplitude_noise(fp)
    s0 = spot_metrics(fp.x, fp.intensity(fp.weights(0.0, 0.9)))
    par = json.loads((RESULTS / "m9_E_power.json").read_text())["outputs"]
    # scale the paraxial power result by the ratio of centre peak intensity per unit port power
    fp_par = m9.FocalPlaneOPA()
    s_par = spot_metrics(fp_par.x, fp_par.intensity(fp_par.weights(0.0, 0.9)))
    ratio = s_par["peak"] / s0["peak"]
    power = {"P_centre_mW": par["P_fibre_per_tweezer_centre_mW"] * ratio,
             "P_field_mean_mW": par["P_fibre_per_tweezer_centre_mW"] * ratio * ind["power_factor_mean_over_field"],
             "P_edge_mW": par["P_fibre_per_tweezer_centre_mW"] * ratio * ind["power_factor_edge"],
             "P_20_W": 20 * par["P_fibre_per_tweezer_centre_mW"] * ratio * ind["power_factor_mean_over_field"] / 1e3,
             "peak_ratio_paraxial_over_bpm": ratio}
    out = {"calibration": calname, "port_spacing_um": slope, "window": rows, "field": field, "independence": ind,
           "noise": noise, "power": power, "centre_fwhm_um": s0["fwhm"]}
    save_result(f"m10_route_E_{calname}", script=SCRIPT,
                inputs={"T": f"results/m10_T_full_{calname}.npz", "sigma_ports": 0.9,
                        "power_reference": "m9_E_power (paraxial), scaled by the centre peak ratio"},
                outputs=out)
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["calib", "small", "full", "route"]
    if "calib" in parts:
        print(json.dumps(calib(), indent=1, default=float), flush=True)
    if "small" in parts:
        r = small()
        for k, v in r.items():
            print(k, np.round(v["capture"], 4).tolist(), flush=True)
    if "compare" in parts:
        for r in compare():
            print(r, flush=True)
    for cn in ("A", "B"):
        if f"full{cn}" in parts or ("full" in parts and cn == "A"):
            o = full(cn)
            print({k: (v if not isinstance(v, list) else f"[{len(v)}]") for k, v in o.items()}, flush=True)
        if f"route{cn}" in parts or ("route" in parts and cn == "A"):
            o = route(cn)
            print(json.dumps({k: v for k, v in o.items() if k not in ("window", "field")}, indent=1, default=float),
                  flush=True)
