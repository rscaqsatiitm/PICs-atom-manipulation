"""Module 11: optimisation of emitter and slab geometry for the sliding-window focal-plane OPA (route E).

Figure of merit: the field-averaged laser power per tweezer for equal depth, i.e. the harmonic mean over the
field of use (|x| <= 60 um) of the peak intensity per unit port power, P_eff = 1 / <1/peak(x)>. Constraints:
neighbour intensity at 4 um <= 1e-4; gaps >= 0.3 um; the fast-actuator count (128 ports) is fixed.

Stage 1: oversampling N_arms (128 ... 256 for 128 ports) x emitter facet width.
Stage 2: arm and port widths at the FPR at the best stage-1 point; window width sigma for the crosstalk limit.
Coupling limits for the wide sections use the calibrated lateral model (m10_eim_calibration, A).
"""
from __future__ import annotations

import json
import sys

import numpy as np

from picsim.bpm import eim_coupling_length
from picsim.opa.focal_plane import FocalPlane
from picsim.parallel import pmap
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m11_optimise.py"
LAM = 1.064435
X_HALF = 60.0


def cal():
    d = json.loads((RESULTS / "m10_eim_calibration.json").read_text())["outputs"]["A"]
    return d["n_core"], d["n_clad"]


def evaluate(cfg):
    sigma = cfg.pop("sigma", 0.9)
    fp = FocalPlane(lam=LAM, **cfg)
    scan = fp.field_scan(sigma, x_half=X_HALF, n=13)
    pk = np.array([s["peak"] for s in scan])
    c = scan[len(scan) // 2]
    # decomposition at centre and edge: FPR capture of the window vs relay transmission
    dec = {}
    for name, u in (("centre", 0.0), ("edge", X_HALF / fp.dx_port)):
        w = fp.weights(u, sigma)
        capt = float(np.sum(np.abs(w @ fp.A) ** 2))
        dec[name] = {"fpr_capture": capt, "delivered": fp.spot(u, sigma)["delivered"]}
        dec[name]["relay_transmission"] = dec[name]["delivered"] / capt
    return {**cfg, "sigma": sigma, "R_um": fp.R, "aperture_um": fp.N_arms * cfg.get("p_out", 1.5),
            "peak_centre": float(c["peak"]), "P_eff": float(1 / np.mean(1 / pk)), "edge_over_centre": float(pk[0] / c["peak"]),
            "fwhm_um": float(c["fwhm"]), "nb4_max": float(max(s["nb4"] for s in scan)), "decomposition": dec,
            "peaks": pk.tolist()}


def coupling_limits(widths, pitch):
    nc, ncl = cal()
    out = {}
    for w in widths:
        g = pitch - w
        if g < 0.3:
            out[w] = None
            continue
        Lc = eim_coupling_length(w, g, LAM, nc, ncl)
        out[w] = {"gap": g, "Lc_um": Lc, "max_parallel_um_for_-20dB": 0.064 * Lc}
    return out


RP_ROWLAND = 0.555  # port-arc radius / R minimising aberration over the field (exact path lengths, m11 notes)


def stage1():
    cfgs = [dict(N_arms=na, w_emit=we, model="rs", rp_ratio=RP_ROWLAND) for na in (128, 160, 192, 256)
            for we in (0.6, 0.9, 1.2)]
    cfgs += [dict(N_arms=128, w_emit=0.6, model="rs", rp_ratio=1.0), dict(N_arms=128, w_emit=0.6, model="paraxial")]
    rows = pmap(evaluate, cfgs)
    base = [r for r in rows if r["N_arms"] == 128 and r["w_emit"] == 0.6 and r.get("model") == "rs"
            and r.get("rp_ratio") == 1.0][0]
    for r in rows:
        r["P_eff_rel_baseline"] = r["P_eff"] / base["P_eff"]
    lim = {"emitters_p1.5": coupling_limits((0.6, 0.9, 1.2), 1.5)}
    f, ax = fig(1, 2, h=3.0)
    for we, c in zip((0.6, 0.9, 1.2), ("C0", "C1", "C2")):
        r = sorted([x for x in rows if x["w_emit"] == we and x.get("rp_ratio") == RP_ROWLAND], key=lambda x: x["N_arms"])
        ax[0, 0].plot([x["N_arms"] for x in r], [x["P_eff_rel_baseline"] for x in r], "-o", color=c, ms=3,
                      label=f"facet {we} um")
        ax[0, 1].plot([x["N_arms"] for x in r], [x["edge_over_centre"] for x in r], "-o", color=c, ms=3)
    ax[0, 0].set_xlabel("arms / emitters (128 ports)")
    ax[0, 0].set_ylabel("field-averaged peak per W (rel. baseline)")
    ax[0, 1].set_xlabel("arms / emitters (128 ports)")
    ax[0, 1].set_ylabel("edge (60 um) / centre peak")
    ax[0, 0].legend(fontsize=7)
    save_result("m11_stage1", script=SCRIPT, fig=f,
                inputs={"N_ports": 128, "p_port": 2.0, "p_a": 2.0, "w_port": 1.7, "w_arm": 1.7, "p_out": 1.5,
                        "sigma": 0.9, "field_of_use_um": X_HALF,
                        "model": "picsim.opa.focal_plane with the exact-path star coupler (star_rs), Rowland-type "
                                 "port arc R_p = 0.555 R; baseline = confocal (R_p = R); paraxial row for reference"},
                outputs={"rows": rows, "coupling_limits": lim})
    return rows, lim


def stage2(n_arms, w_emit):
    base = dict(N_arms=n_arms, w_emit=w_emit, model="rs", rp_ratio=RP_ROWLAND)
    cfgs = [dict(base, w_arm=wa, w_port=wp) for wa in (1.2, 1.5, 1.7) for wp in (1.4, 1.7)]
    rows = pmap(evaluate, cfgs)
    best = max(rows, key=lambda r: r["P_eff"])
    sig = pmap(evaluate, [dict(base, w_arm=best["w_arm"], w_port=best["w_port"], sigma=s)
                          for s in (0.8, 0.85, 0.9, 1.0)])
    lim = {"arms_p2.0": coupling_limits((1.2, 1.5, 1.7), 2.0), "ports_p2.0": coupling_limits((1.4, 1.7), 2.0)}
    save_result("m11_stage2", script=SCRIPT,
                inputs={"N_arms": n_arms, "w_emit": w_emit, "p_out": 1.5, "N_ports": 128},
                outputs={"rows": rows, "sigma_scan": sig, "coupling_limits": lim})
    return rows, sig, lim


if __name__ == "__main__":
    parts = sys.argv[1:] or ["stage1"]
    if "stage1" in parts:
        rows, lim = stage1()
        for r in sorted(rows, key=lambda r: (r["N_arms"], r["w_emit"])):
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                   if k not in ("peaks", "decomposition")}, flush=True)
            print("    ", {k: {kk: round(vv, 3) for kk, vv in v.items()} for k, v in r["decomposition"].items()}, flush=True)
        print(json.dumps(lim, indent=1, default=float), flush=True)
    if "stage2" in parts:
        na, we = int(sys.argv[sys.argv.index("stage2") + 1]), float(sys.argv[sys.argv.index("stage2") + 2])
        rows, sig, lim = stage2(na, we)
        for r in rows + sig:
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                   if k not in ("peaks", "decomposition")}, flush=True)
        print(json.dumps(lim, indent=1, default=float), flush=True)
