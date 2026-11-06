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
def full_layout(n_ports=128, n_arms=128, rp_ratio=1.0, w_port=1.7, w_arm=1.7):
    p = 2.0
    R = n_arms * p * p / (LAM / 1.8423)  # one port per resolvable direction
    return FPRLayout(n_ports, n_arms, p, p, R, w_wide=w_port, w_wide_arm=w_arm, L_taper=30, L_port=55, L_arm=70,
                     z0=60, R_p=rp_ratio * R)


def full(calname="A", tag=None, **layout):
    nc, ncl = cal(calname)
    lay = full_layout(**layout)
    tag = tag or f"full_{calname}"
    t0 = time.time()
    T, info = transfer_matrix(lay, LAM, nc, ncl, dx=0.04, dz=0.125, progress=True)
    np.savez_compressed(RESULTS / f"m10_T_{tag}.npz", T=T)
    cap = np.sum(np.abs(T) ** 2, axis=1)
    # arm-phase calibration with the two central ports (flat phase for the centre of the field)
    c0 = lay.N_ports // 2
    c = np.exp(-1j * np.angle(T[c0 - 1] + T[c0]))
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
           "capture": cap.tolist(), "capture_centre": float(cap[c0 - 1:c0 + 1].mean()),
           "capture_edge": float(cap[[0, -1]].mean()),
           "phase_slope_per_port": float(np.polyfit(np.arange(lay.N_ports), slope, 1)[0]),
           "phase_residual_rms_rad": resid, "phase_residual_max_rad": float(max(resid)),
           "strehl_like_min": float(np.exp(-max(resid) ** 2))}
    save_result(f"m10_fpr_{tag}", script=SCRIPT,
                inputs={"layout": {"N_ports": lay.N_ports, "N_arms": lay.N_arms, "p_port": 2.0, "p_a": 2.0, "R": lay.R,
                                   "R_p": lay.R_p, "w_port_wide": lay.w_wide, "w_arm_wide": lay.w_wide_arm, "w_narrow": 0.6, "L_taper": 30.0},
                        "calibration": calname, "dx_um": 0.04, "dz_um": 0.125},
                outputs=out,
                notes=f"T saved in m10_T_{tag}.npz. Arm phases calibrated with the two central ports; the "
                      "residual is the rms deviation of each port's arm phases from a straight line (aberration).")
    return out


# ------------------------------------------------------------------ route E with the BPM matrix
def route(calname="A", tag=None, fp_kwargs=None):
    """Route-E metrics with a BPM transfer matrix (static arm and port trims applied, picsim.opa.focal_plane)."""
    from picsim.opa.focal_plane import FocalPlane

    tag = tag or f"full_{calname}"
    T = np.load(RESULTS / f"m10_T_{tag}.npz")["T"]
    fp = FocalPlane(A=T, **(fp_kwargs or {}))
    ports_x = [spot_metrics(fp.x, np.abs(fp.A[k] @ fp.G) ** 2)["x_peak"] for k in range(0, fp.N_ports, 8)]
    slope = float(np.polyfit(np.arange(0, fp.N_ports, 8), ports_x, 1)[0])
    rows, field = m9.window_study(fp)
    calt = m9.calibrated(fp, 0.9)
    ind = m9.independence(fp, cal=calt)
    us = np.arange(-60, 61, 1.0)
    ind["power_factor_mean_over_field"] = float(np.mean([calt[1](u) ** 2 for u in us]))
    ind["power_factor_edge"] = float(calt[1](60.0) ** 2)
    noise = m9.amplitude_noise(fp)
    scan = fp.field_scan(0.9, x_half=60.0, n=13)
    pk = np.array([q["peak"] for q in scan])
    out = {"tag": tag, "port_spacing_um": slope, "window": rows, "field": field, "independence": ind, "noise": noise,
           "peak_centre": float(pk[len(pk) // 2]), "P_eff": float(1 / np.mean(1 / pk)),
           "edge_over_centre": float(pk[0] / pk[len(pk) // 2]), "centre_fwhm_um": float(scan[len(scan) // 2]["fwhm"])}
    save_result(f"m10_route_E_{tag}", script=SCRIPT,
                inputs={"T": f"results/m10_T_{tag}.npz", "sigma_ports": 0.9, "focal_plane": fp_kwargs or {}},
                outputs=out,
                notes="P_eff: harmonic mean over |x| <= 60 um of the peak intensity per unit port power; compare "
                      "with m11 rows of the same geometry.")
    return out


def route_ls(tag="full_A", fp_kwargs=None, seed=0):
    """Route E with least-squares port synthesis on the calibrated response (picsim.opa.focal_plane.LinearSynth)."""
    from picsim.opa.focal_plane import FocalPlane, LinearSynth

    T = np.load(RESULTS / f"m10_T_{tag}.npz")["T"]
    fp = FocalPlane(A=T, **(fp_kwargs or {}))
    ls = LinearSynth(fp)
    rng = np.random.default_rng(seed)
    while True:
        xt = np.sort(rng.uniform(-60, 60, 20))
        if np.min(np.diff(xt)) >= 4.0:
            break
    th = rng.uniform(0, 2 * np.pi, 20)
    w0 = ls.weights(xt, phases=th)
    I0 = fp.intensity(w0)
    base = [fp.peak_near(I0, x) for x in xt]
    pk0 = np.array([b["peak"] for b in base])
    out = {"tag": tag, "static_pos_err_max_nm": float(max(abs(b["x_peak"] - x) for b, x in zip(base, xt)) * 1e3),
           "static_peak_spread": float(np.ptp(pk0) / pk0.mean()),
           "peak_per_unit_port_power_mean": float(pk0.mean() / (np.sum(np.abs(w0) ** 2) / 20)), "moves": []}
    j = int(np.argmax(np.diff(xt)))
    for d in (0.1, 0.5, 1.0, 2.0, 4.0):
        xs = xt.copy()
        xs[j] += d
        I = fp.intensity(ls.weights(xs, phases=th))
        mv = fp.peak_near(I, xs[j])
        oth = [fp.peak_near(I, x) for i, x in enumerate(xs) if i != j]
        ob = [b for i, b in enumerate(base) if i != j]
        out["moves"].append({"move_um": d, "moved_pos_err_nm": float((mv["x_peak"] - xs[j]) * 1e3),
                             "others_pos_change_max_nm": float(max(abs(o["x_peak"] - b["x_peak"]) for o, b in zip(oth, ob)) * 1e3),
                             "others_peak_change_max": float(max(abs(o["peak"] / b["peak"] - 1) for o, b in zip(oth, ob)))})
    xo, tho = np.delete(xt, j), np.delete(th, j)
    I = fp.intensity(ls.weights(xo, phases=tho))
    ob = [b for i, b in enumerate(base) if i != j]
    out["switch_off"] = {"residual_at_site": float(np.interp(xt[j], fp.x, I) / base[j]["peak"]),
                         "others_peak_change_max": float(max(abs(fp.peak_near(I, x)["peak"] / b["peak"] - 1)
                                                             for x, b in zip(xo, ob)))}
    # single tweezer across the field: peak per unit port power (sets laser power with depth equalised)
    xs1 = np.linspace(-60, 60, 13)
    pk1 = []
    for x0 in xs1:
        w = ls.weights([x0])
        I = fp.intensity(w / np.sqrt(np.sum(np.abs(w) ** 2)))
        pk1.append(fp.peak_near(I, x0)["peak"])
    pk1 = np.array(pk1)
    out.update({"single_peak_centre": float(pk1[6]), "single_P_eff": float(1 / np.mean(1 / pk1)),
                "single_edge_over_centre": float(pk1[0] / pk1[6]),
                "single_fwhm_centre_um": float(fp.peak_near(fp.intensity(ls.weights([0.0])), 0.0)["fwhm"])})
    # actuator errors on the synthesised port amplitudes
    noise = []
    w1 = ls.weights([12.3])
    s1 = fp.peak_near(fp.intensity(w1), 12.3)
    for sa in (0.003, 0.01, 0.03):
        dx = []
        for _ in range(150):
            w = w1 * (1 + sa * rng.standard_normal(len(w1))) * np.exp(1j * sa * rng.standard_normal(len(w1)))
            dx.append(fp.peak_near(fp.intensity(w), 12.3)["x_peak"] - s1["x_peak"])
        noise.append({"sigma_amp_and_phase": sa, "pos_rms_nm": float(np.std(dx) * 1e3)})
    out["noise"] = noise
    save_result(f"m10_route_E_LS_{tag}", script=SCRIPT,
                inputs={"T": f"results/m10_T_{tag}.npz", "synthesis": "LinearSynth: Gaussian targets, FWHM 1.68 um, "
                        "Tikhonov 1e-4 of the mean eigenvalue", "focal_plane": fp_kwargs or {}},
                outputs=out)
    return out


if __name__ == "__main__":
    if "routels" in sys.argv:
        i = sys.argv.index("routels")
        tag = sys.argv[i + 1]
        kw = json.loads(sys.argv[i + 2]) if len(sys.argv) > i + 2 else {}
        print(json.dumps(route_ls(tag, kw), indent=1, default=float), flush=True)
        sys.exit(0)
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
    if "fullopt" in parts:
        cfg = json.loads(sys.argv[sys.argv.index("fullopt") + 1])
        o = full("A", tag="opt_A", **cfg)
        print({k: (v if not isinstance(v, list) else f"[{len(v)}]") for k, v in o.items()}, flush=True)
    if "routeopt" in parts:
        cfg = json.loads(sys.argv[sys.argv.index("routeopt") + 1])
        o = route("A", tag="opt_A", fp_kwargs=cfg)
        print(json.dumps({k: v for k, v in o.items() if k not in ("window", "field")}, indent=1, default=float),
              flush=True)
    for cn in ("A", "B"):
        if f"full{cn}" in parts or ("full" in parts and cn == "A"):
            o = full(cn)
            print({k: (v if not isinstance(v, list) else f"[{len(v)}]") for k, v in o.items()}, flush=True)
        if f"route{cn}" in parts or ("route" in parts and cn == "A"):
            o = route(cn)
            print(json.dumps({k: v for k, v in o.items() if k not in ("window", "field")}, indent=1, default=float),
                  flush=True)
