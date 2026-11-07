"""Module 2: edge couplers and MMI splitters.

Outputs (results/):
  m2_couplers   inverse-taper tip sweep: mode size, overlap with lensed and standard fibres, misalignment
  m2_mmi        1x2 and 2x2 MMI design at 1064.4 and 894.6 nm: length, excess loss, bandwidth,
                fabrication sensitivity, and split-ratio statistics used by the MZI Monte Carlo (Module 4)
"""
from __future__ import annotations

import sys

import numpy as np

from picsim import couplers as cp
from picsim import params as prm
from picsim.mmi import Lateral, mmi_transmission, optimise_length
from picsim.parallel import pmap
from picsim.plotting import fig
from picsim.provenance import save_result

SCRIPT = "scripts/m2_couplers_mmi.py"
LAMS = {"1064": prm.LAM_TRANSPORT.value, "895": prm.LAM_ADDRESS.value}
TIPS = [0.10, 0.125, 0.15, 0.20, 0.25, 0.30, 0.40]
MFDS = {"1064": {"lensed_2.5": prm.MFD_LENSED.value, "UHNA_4.0": 4.0, "HI1060_6.2": prm.MFD_HI1060.value},
        "895": {"lensed_2.5": prm.MFD_LENSED.value, "UHNA_4.0": 4.0, "PM780_5.3": prm.MFD_PM780.value}}


def _tip_job(args):
    lk, t, w = args
    lam = LAMS[lk]
    x, y, E, te = cp.tip_mode_field(w, t, lam)
    wx, wy = cp.mode_size(x, y, E)
    r = {"lam": lk, "t": t, "w_tip": w, "neff": te.neff.real, "aeff_um2": te.aeff, "mfd_x": wx, "mfd_y": wy,
         "window_half_um": te.window_half_um, "window_converged": bool(te.window_converged)}
    for name, mfd in MFDS[lk].items():
        eta, y0 = cp.best_overlap(x, y, E, mfd / 2)
        r[f"eta_{name}"] = eta
        r[f"y0_{name}"] = y0
        # 1 dB lateral misalignment tolerance
        offs = np.linspace(0, 2.0, 41)
        etas = [cp.overlap(x, y, E, cp.gaussian(x, y, mfd / 2, dx, y0)) for dx in offs]
        lim = offs[np.array(etas) >= eta * 10 ** (-0.1)]
        r[f"tol_1dB_x_{name}"] = float(lim.max())
    return r


def couplers():
    jobs = [(lk, t, w) for lk in LAMS for t in (0.2, 0.4) for w in TIPS]
    rows = pmap(_tip_job, jobs)
    best = {}
    f, ax = fig(1, 2, h=3.0)
    for j, lk in enumerate(LAMS):
        for t, ls in ((0.2, "--"), (0.4, "-")):
            r = [x for x in rows if x["lam"] == lk and x["t"] == t]
            for (name, _), c in zip(MFDS[lk].items(), ("C0", "C1", "C2")):
                eta = np.array([x[f"eta_{name}"] if x["window_converged"] else np.nan for x in r])
                ax[0, j].plot([x["w_tip"] for x in r], 10 * np.log10(eta), ls + "o", color=c, ms=3,
                              label=f"{name}, t={t*1e3:.0f} nm")
                k = int(np.nanargmax(eta))
                best[f"{lk}_{t}_{name}"] = {"w_tip": r[k]["w_tip"], "dB": float(10 * np.log10(eta[k])),
                                            "tol_1dB_x_um": r[k][f"tol_1dB_x_{name}"]}
        ax[0, j].set_title(f"lambda = {LAMS[lk]*1e3:.1f} nm")
        ax[0, j].set_xlabel("taper tip width (um)")
        ax[0, j].set_ylim(-12, 0)
        ax[0, j].axhline(-2, color="k", lw=0.8)
    ax[0, 0].set_ylabel("mode overlap (dB)")
    ax[0, 0].legend(fontsize=5)
    out = {"rows": rows, "best": best,
           "best_lensed_1064_dB": max(v["dB"] for k, v in best.items() if k.startswith("1064") and "lensed" in k)}
    # facet intensity at 1.4 W into the best 1064 tip (both thicknesses)
    for t in (0.2, 0.4):
        b = best[f"1064_{t}_lensed_2.5"]
        rr = [x for x in rows if x["lam"] == "1064" and x["t"] == t and x["w_tip"] == b["w_tip"]][0]
        out[f"facet_intensity_MW_cm2_at_1.4W_t{t}"] = prm.P_CHIP.value / (rr["aeff_um2"] * 1e-8) / 1e6
    save_result("m2_couplers", script=SCRIPT, fig=f,
                inputs={"MFDs": {"lensed": prm.MFD_LENSED, "HI1060": prm.MFD_HI1060, "PM780": prm.MFD_PM780,
                                 "UHNA": "CHOICE 4.0 um"}, "target": prm.FACET_LOSS,
                        "method": "femwell TE0 at the tip, scalar overlap with a Gaussian, fibre height optimised"},
                outputs=out,
                notes="Excludes Fresnel loss (index-matched or AR-coated facet assumed) and taper conversion loss.")
    return out


# ------------------------------------------------------------------ MMI
def mmi_design(lk, W, wa):
    lam = LAMS[lk]
    lat = Lateral(lam, prm.T_SIN.value)
    x = lat.grid(W / 2 + 1)
    ne, _ = lat.modes(x, [(0, W)])
    Lpi = lam / (2 * (ne[0] - ne[1]))
    # 1x2 symmetric interference, outputs at +/- W/4
    L12, T12 = optimise_length(lat, W, wa, 0.0, wa, [-W / 4, W / 4], 3 * Lpi / 8)
    # 2x2 paired interference (3 dB), ports at +/- W/6; maximise total, then report the split ratio
    L22g = Lpi / 2
    Ls = np.linspace(0.8 * L22g, 1.2 * L22g, 801)
    T, _ = mmi_transmission(lat, W, Ls, wa, -W / 6, wa, [-W / 6, W / 6])
    score = T.sum(axis=1) - np.abs(T[:, 0] - T[:, 1])
    i = int(np.argmax(score))
    L22, T22 = float(Ls[i]), T[i]
    return lat, Lpi, L12, T12, L22, T22


def _mmi_case(args):
    lk, W, wa = args
    lat, Lpi, L12, T12, L22, T22 = mmi_design(lk, W, wa)
    lam = LAMS[lk]
    r = {"lam": lk, "W": W, "w_access": wa, "Lpi": Lpi, "L_1x2": L12, "excess_1x2_dB": -10 * np.log10(T12),
         "L_2x2": L22, "excess_2x2_dB": -10 * np.log10(T22.sum()),
         "split_2x2_cross": float(T22[1] / T22.sum())}
    # bandwidth of the 1x2 (excess loss at +/- 10 nm)
    for dl in (-0.01, 0.01):
        lat2 = Lateral(lam + dl, prm.T_SIN.value)
        T, _ = mmi_transmission(lat2, W, L12, wa, 0.0, wa, [-W / 4, W / 4])
        r[f"excess_1x2_dB_dlam{dl*1e3:+.0f}nm"] = -10 * np.log10(T.sum())
    # fabrication: width bias and thickness error at fixed design length
    for dW in (-0.02, 0.02):
        T, _ = mmi_transmission(lat, W, L12, wa, 0.0, wa, [-W / 4, W / 4], dW=dW)
        r[f"excess_1x2_dB_dW{dW*1e3:+.0f}nm"] = -10 * np.log10(T.sum())
        T, _ = mmi_transmission(lat, W, L22, wa, -W / 6, wa, [-W / 6, W / 6], dW=dW)
        r[f"split_2x2_cross_dW{dW*1e3:+.0f}nm"] = float(T[1] / T.sum())
    for dt in (-0.01, 0.01):
        latt = Lateral(lam, prm.T_SIN.value, dt=dt)
        T, _ = mmi_transmission(latt, W, L12, wa, 0.0, wa, [-W / 4, W / 4])
        r[f"excess_1x2_dB_dt{dt*1e3:+.0f}nm"] = -10 * np.log10(T.sum())
        T, _ = mmi_transmission(latt, W, L22, wa, -W / 6, wa, [-W / 6, W / 6])
        r[f"split_2x2_cross_dt{dt*1e3:+.0f}nm"] = float(T[1] / T.sum())
    # 1x2 imbalance from input offset (asymmetric perturbation)
    for d in (0.005, 0.01, 0.02, 0.05):
        T, _ = mmi_transmission(lat, W, L12, wa, d, wa, [-W / 4, W / 4])
        r[f"imbalance_1x2_dB_offset{d*1e3:.0f}nm"] = float(10 * np.log10(T[0] / T[1]))
    return r


def _mmi_mc(args):
    lk, W, wa, L12, L22, seed, n, sW, st, sd = args
    rng = np.random.default_rng(seed)
    lam = LAMS[lk]
    out = []
    for _ in range(n):
        dW, dt, d = sW * rng.standard_normal(), st * rng.standard_normal(), sd * rng.standard_normal()
        lat = Lateral(lam, prm.T_SIN.value, dt=dt)
        T1, _ = mmi_transmission(lat, W, L12, wa, d, wa, [-W / 4, W / 4], dW=dW)
        T2, _ = mmi_transmission(lat, W, L22, wa, -W / 6 + d, wa, [-W / 6, W / 6], dW=dW)
        out.append((T1[0] / T1.sum(), T1.sum(), T2[1] / T2.sum(), T2.sum()))
    return out


def mmi():
    cases = [(lk, W, wa) for lk in LAMS for W, wa in ((3.0, 1.0), (4.0, 1.2), (6.0, 1.5))]
    rows = pmap(_mmi_case, cases)
    # Monte Carlo of split ratios for the chosen 4 um designs (feeds Module 4)
    sW, st, sd = 0.010, 0.005, 0.005  # CHOICE: 1-sigma width bias, thickness error, lateral offset (um)
    mc = {}
    for lk in LAMS:
        d = [r for r in rows if r["lam"] == lk and r["W"] == 4.0][0]
        res = pmap(_mmi_mc, [(lk, 4.0, 1.2, d["L_1x2"], d["L_2x2"], s, 50, sW, st, sd) for s in range(4)])
        arr = np.array([x for chunk in res for x in chunk])
        mc[lk] = {"kappa_1x2_mean": float(arr[:, 0].mean()), "kappa_1x2_std": float(arr[:, 0].std()),
                  "T_1x2_mean": float(arr[:, 1].mean()),
                  "kappa_2x2_mean": float(arr[:, 2].mean()), "kappa_2x2_std": float(arr[:, 2].std()),
                  "T_2x2_mean": float(arr[:, 3].mean()), "n": len(arr)}
    f, ax = fig(1, 1, w=4.5, h=3.0)
    for lk, mk in (("1064", "o"), ("895", "s")):
        r = [x for x in rows if x["lam"] == lk]
        ax[0, 0].plot([x["W"] for x in r], [x["excess_1x2_dB"] for x in r], "-" + mk, label=f"1x2, {lk} nm")
        ax[0, 0].plot([x["W"] for x in r], [x["excess_2x2_dB"] for x in r], "--" + mk, label=f"2x2, {lk} nm")
    ax[0, 0].set_xlabel("MMI width (um)")
    ax[0, 0].set_ylabel("excess loss (dB)")
    ax[0, 0].legend(fontsize=7)
    d1064 = [r for r in rows if r["lam"] == "1064" and r["W"] == 4.0][0]
    d895 = [r for r in rows if r["lam"] == "895" and r["W"] == 4.0][0]
    out = {"rows": rows, "mc": mc, "mc_sigmas_um": {"width": sW, "thickness": st, "offset": sd},
           "excess_loss_1064_dB": d1064["excess_1x2_dB"], "excess_loss_895_dB": d895["excess_1x2_dB"],
           "tree_6_stage_1064_dB": 6 * d1064["excess_1x2_dB"], "tree_4_stage_895_dB": 4 * d895["excess_1x2_dB"]}
    save_result("m2_mmi", script=SCRIPT, fig=f,
                inputs={"t": prm.T_SIN, "method": "effective-index + guided-mode EME (picsim.mmi)",
                        "fab_sigmas": "CHOICE: width bias 10 nm, thickness 5 nm, lateral offset 5 nm (1 sigma)",
                        "tree_target": prm.TREE_LOSS},
                outputs=out,
                notes="2D effective-index model; excess loss from unguided content only (no sidewall scattering, "
                      "no back-reflection). Absolute loss should be confirmed by 3D EME/FDTD.")
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["couplers", "mmi"]
    if "couplers" in parts:
        o = couplers()
        for k, v in o["best"].items():
            print(k, v, flush=True)
        print({k: v for k, v in o.items() if k.startswith("facet") or k.startswith("best_l")}, flush=True)
    if "mmi" in parts:
        o = mmi()
        for r in o["rows"]:
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
        print(o["mc"], flush=True)
