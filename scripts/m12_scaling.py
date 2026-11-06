"""Module 12: how far the focal-plane OPA scales.

Parts:
  packing     tweezers per row: fill a +/-60 um row at spacing s with least-squares synthesis on the BPM matrix
              of design D1 (alternating 0 / pi/2 phases); uniformity, efficiency, barrier between neighbours,
              position error, independence of a moved / switched-off tweezer
  aberration  star-coupler aberration vs size (exact path lengths): N_ports = 128 ... 1024, 1.5x oversampled
              arms; Rowland-type circle and a circle + quartic port curve
  table       derived scaling of field, tweezers per row, slab size, actuators, data rate, power

Outputs: results/m12_packing.json, m12_aberration.json, m12_table.json
"""
from __future__ import annotations

import json
import sys

import numpy as np
from scipy.optimize import minimize

from picsim import params as prm
from picsim.opa.focal_plane import FocalPlane, LinearSynth
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m12_scaling.py"
LAM = prm.LAM_TRANSPORT.value
D1 = {"N_ports": 128, "N_arms": 192, "w_emit": 1.2, "p_out": 1.5, "w_port": 1.7, "w_arm": 1.7}


# ------------------------------------------------------------------ packing
def packing():
    T = np.load(RESULTS / "m10_T_opt_A.npz")["T"]
    fp = FocalPlane(A=T, **D1)
    ls = LinearSynth(fp)
    rows = []
    for s in (6.0, 4.0, 3.0, 2.5, 2.0, 1.75):
        xt = np.arange(-60.0, 60.0 + 1e-9, s)
        K = len(xt)
        th = np.where(np.arange(K) % 2 == 0, 0.0, np.pi / 2)
        w = ls.weights(xt, phases=th)
        I = fp.intensity(w)
        win = min(0.6, s / 3)
        pk = [fp.peak_near(I, x, win) for x in xt]
        pkv = np.array([p["peak"] for p in pk])
        mids = 0.5 * (xt[1:] + xt[:-1])
        barrier = np.interp(mids, fp.x, I) / pkv.mean()
        port_pow = np.sum(np.abs(w) ** 2)
        # independence: move the middle tweezer by +0.25 um; switch it off
        j = K // 2
        xs = xt.copy()
        xs[j] += 0.25
        I2 = fp.intensity(ls.weights(xs, phases=th))
        oth = [fp.peak_near(I2, x, win) for i, x in enumerate(xs) if i != j]
        ob = [p for i, p in enumerate(pk) if i != j]
        I3 = fp.intensity(ls.weights(np.delete(xt, j), phases=np.delete(th, j)))
        oth3 = [fp.peak_near(I3, x, win) for x in np.delete(xt, j)]
        rows.append({
            "spacing_um": s, "K": K, "peak_spread": float(np.ptp(pkv) / pkv.mean()),
            "peak_per_tweezer_per_unit_power": float(pkv.mean() / (port_pow / K)),
            "barrier_max_rel": float(barrier.max()), "pos_err_max_nm": float(max(abs(p["x_peak"] - x) for p, x in zip(pk, xt)) * 1e3),
            "fwhm_mean_um": float(np.mean([p["fwhm"] for p in pk])),
            "move_others_pos_nm": float(max(abs(o["x_peak"] - b["x_peak"]) for o, b in zip(oth, ob)) * 1e3),
            "move_others_peak": float(max(abs(o["peak"] / b["peak"] - 1) for o, b in zip(oth, ob))),
            "off_residual": float(np.interp(xt[j], fp.x, I3) / pkv.mean()),
            "off_others_peak": float(max(abs(o["peak"] / b["peak"] - 1) for o, b in zip(oth3, ob))),
        })
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    f, ax = fig(1, 2, h=3.0)
    S = [r["spacing_um"] for r in rows]
    ax[0, 0].plot(S, [r["peak_per_tweezer_per_unit_power"] for r in rows], "-o", ms=3)
    ax[0, 0].set_xlabel("tweezer spacing (um)")
    ax[0, 0].set_ylabel("peak per tweezer per unit port power")
    ax[0, 1].semilogy(S, [r["barrier_max_rel"] for r in rows], "-o", ms=3, label="intensity between neighbours")
    ax[0, 1].semilogy(S, [max(r["move_others_peak"], 1e-6) for r in rows], "-s", ms=3, label="others' change, 0.25 um move")
    ax[0, 1].set_xlabel("tweezer spacing (um)")
    ax[0, 1].legend(fontsize=7)
    save_result("m12_packing", script=SCRIPT, fig=f,
                inputs={"design": D1, "T": "results/m10_T_opt_A.npz (BPM)", "row_um": "+/-60",
                        "phases": "alternating 0, pi/2", "synthesis": "LinearSynth, Gaussian targets FWHM 1.68 um"},
                outputs={"rows": rows},
                notes="barrier: intensity midway between neighbours relative to the mean tweezer peak.")
    return rows


# ------------------------------------------------------------------ aberration vs size
def residuals(N_ports, N_arms, p=2.0, rp_ratio=0.555, c4=0.0, weight_edge=np.exp(-6)):
    n_s = 1.8423
    lam_s = LAM / n_s
    k_s = 2 * np.pi / lam_s
    R = N_arms * p * p / lam_s
    phi = (np.arange(N_arms) - (N_arms - 1) / 2) * p / R
    Ax, Az = R * np.sin(phi), R * np.cos(phi)
    n = np.arange(N_arms)
    # arm-envelope weight (Gaussian, edge level weight_edge) as produced by the window/least-squares synthesis
    y = (n - (N_arms - 1) / 2) / ((N_arms - 1) / 2)
    wgt = np.exp(np.log(weight_edge) * y**2)
    Rp = rp_ratio * R
    out = []
    used = np.arange(int(0.06 * N_ports), N_ports - int(0.06 * N_ports))  # |x| <= ~88% of the port field
    for k in used:
        s = (k - (N_ports - 1) / 2) * p
        a = s / Rp
        Px = Rp * np.sin(a)
        Pz = Rp - Rp * np.cos(a) + c4 * Px**4 / R**3
        ph = k_s * np.hypot(Ax - Px, Az - Pz)
        cf = np.polyfit(n, ph, 1, w=np.sqrt(wgt))
        r = ph - np.polyval(cf, n)
        r -= np.sum(wgt * r) / np.sum(wgt)
        out.append(np.sqrt(np.sum(wgt * r**2) / np.sum(wgt)))
    return np.array(out), R


def aberration():
    rows = []
    for Np in (128, 192, 256, 384, 512, 768, 1024):
        Na = int(round(1.5 * Np))
        conf, R = residuals(Np, Na, rp_ratio=1.0)
        best = minimize(lambda q: residuals(Np, Na, rp_ratio=q[0])[0].max(), [0.555], method="Nelder-Mead",
                        options={"xatol": 1e-4, "fatol": 1e-5})
        rw, _ = residuals(Np, Na, rp_ratio=best.x[0])
        best2 = minimize(lambda q: residuals(Np, Na, rp_ratio=q[0], c4=q[1])[0].max(), [best.x[0], 0.0],
                         method="Nelder-Mead", options={"xatol": 1e-5, "fatol": 1e-6, "maxiter": 400})
        rq, _ = residuals(Np, Na, rp_ratio=best2.x[0], c4=best2.x[1])
        rows.append({"N_ports": Np, "N_arms": Na, "R_mm": R / 1e3, "confocal_rms_max_rad": float(conf.max()),
                     "rowland_rp_ratio": float(best.x[0]), "rowland_rms_max_rad": float(rw.max()),
                     "quartic_rp_ratio": float(best2.x[0]), "quartic_c4": float(best2.x[1]),
                     "quartic_rms_max_rad": float(rq.max()),
                     "strehl_rowland": float(np.exp(-rw.max() ** 2)), "strehl_quartic": float(np.exp(-rq.max() ** 2))})
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    f, ax = fig(1, 1, w=4.5, h=3.0)
    N = [r["N_ports"] for r in rows]
    for key, lab in (("confocal_rms_max_rad", "confocal"), ("rowland_rms_max_rad", "Rowland-type circle"),
                     ("quartic_rms_max_rad", "circle + quartic")):
        ax[0, 0].loglog(N, [r[key] for r in rows], "-o", ms=3, label=lab)
    ax[0, 0].axhline(0.2, color="k", lw=0.6)
    ax[0, 0].set_xlabel("ports (arms = 1.5 x ports)")
    ax[0, 0].set_ylabel("worst-port rms arm-phase error (rad)")
    ax[0, 0].legend(fontsize=7)
    save_result("m12_aberration", script=SCRIPT, fig=f,
                inputs={"geometry": "arms on a circle of radius R; ports on a circle of radius R_p through the slab "
                        "vertex, optionally + c4 x^4 / R^3; exact path lengths; arm weights Gaussian with edge "
                        "exp(-6)", "ports_used": "central 88%"},
                outputs={"rows": rows},
                notes="Check against BPM for D1 (N = 128 ports, 192 arms): BPM worst-port residual 0.079 rad rms.")
    return rows


# ------------------------------------------------------------------ derived scaling table
def table(pack, aber):
    p_ref = json.loads((RESULTS / "m10_route_E_LS_opt_A.json").read_text())["outputs"]
    P_tw = 47.95 * 0.2281 / p_ref["single_P_eff"]  # mW per 1.2 mK tweezer, field mean (m9_E_power scaled)
    res = LAM / (2 * prm.NA_OBJ.value)
    # practical packing: barrier between neighbours <= 25% of the peak (>= 75% of the depth separates them),
    # moving one tweezer changes the others by < 0.3% and < 5 nm
    s_min = min(r["spacing_um"] for r in pack if r["barrier_max_rel"] <= 0.25 and r["move_others_peak"] < 3e-3
                and r["move_others_pos_nm"] < 5.0)
    rows = []
    for Np in (128, 192, 256, 384, 512, 1024):
        Na = int(round(1.5 * Np))
        ab = [a for a in aber if a["N_ports"] == Np]
        field = Np * res
        R = Na * 2.0 * 2.0 / (LAM / 1.8423)
        rows.append({
            "N_ports": Np, "field_um": field, "tweezers_per_row_at_4um": int(field * 0.88 // 4.0) + 1,
            "tweezers_per_row_at_min_spacing": int(field * 0.88 // s_min) + 1, "min_spacing_um": s_min,
            "slab_R_mm": R / 1e3, "slab_area_mm2": float(R * (Na * 2.0) / 1e6),
            "emitter_aperture_mm": Na * 1.5 / 1e3,
            "fast_actuators": 2 * Np - 1, "static_trims": Na,
            "data_rate_Gbit_s_at_200kHz_16bit": (2 * Np - 1) * 200e3 * 16 / 1e9,
            "worst_aberration_rms_rad_rowland": ab[0]["rowland_rms_max_rad"] if ab else None,
            "worst_aberration_rms_rad_quartic": ab[0]["quartic_rms_max_rad"] if ab else None,
            "power_W_full_row_at_4um_1.2mK": (int(field * 0.88 // 4.0) + 1) * P_tw / 1e3,
            "power_W_full_row_at_4um_0.5mK": (int(field * 0.88 // 4.0) + 1) * P_tw * 0.5 / 1.2 / 1e3,
        })
    out = {"P_tweezer_mW_1.2mK": P_tw, "resolution_cell_um": res, "rows": rows}
    save_result("m12_table", script=SCRIPT,
                inputs={"P_tweezer": "m10_route_E_LS_opt_A (D1, BPM) with the m9_E_power loss chain",
                        "field_used": "88% of the port field", "oversampling": "1.5", "pitches": "2.0 um slab, 1.5 um emitters"},
                outputs=out,
                notes="Power assumes the D1 efficiency holds at larger N (oversampling and facets unchanged); "
                      "aberration from m12_aberration.")
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["packing", "aberration", "table"]
    pack = packing() if "packing" in parts else None
    aber = aberration() if "aberration" in parts else None
    if "table" in parts:
        pack = pack or json.loads((RESULTS / "m12_packing.json").read_text())["outputs"]["rows"]
        aber = aber or json.loads((RESULTS / "m12_aberration.json").read_text())["outputs"]["rows"]
        print(json.dumps(table(pack, aber), indent=1, default=float))
