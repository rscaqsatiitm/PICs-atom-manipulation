"""Module 25: imaging cadence and loss in steady-state refilling.

The steady-state defect fraction of `picsim/megasort.py` for a 1e5-site array (4 x 4 um, 85% loading, 32 x 4 block,
intermediate steps) over the imaging period T_img, the loss per image eps_img, the imaging time t_img and the
vacuum lifetime tau; a few 1e6 points. Fit rho = a eps_img + b (T_img + t_img)/tau and locate where the repair
capacity, not the loss rate, sets the floor.

Output: results/m25_imaging.json (+ .png)
"""
from __future__ import annotations

import itertools
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from picsim.megasort import Block, Layout, Params, Sim
from picsim.plotting import fig
from picsim.provenance import save_result
from picsim.rearrange import Timing

SCRIPT = "scripts/m25_imaging.py"
MID = dict(t_ramp=100e-6, t_lane=200e-6, t_x40=300e-6, t_min=150e-6)
FAST = dict(t_ramp=50e-6, t_lane=100e-6, t_x40=200e-6, t_min=100e-6)


def one(cfg):
    N = cfg["N"]
    B, Bc = 4, 3
    C = int(np.ceil(np.sqrt(N)))
    R = int(np.ceil(N / C * B / Bc))
    tm = Timing(**(FAST if cfg["timing"] == "fast" else MID))
    P = Params(tau=cfg["tau"], eps_img=cfg["eps"], t_img=cfg["t_img"], T_img=cfg["T_img"], f_load=0.85, timing=tm)
    sim = Sim(Layout(R, C, B, Bc), Block(32, 4, 46), P, seed=6)
    t_end = cfg["t_end"]
    sim.run(t_end)
    d = np.array(sim.log)
    late = d[d[:, 0] >= t_end / 2]
    return {**cfg, "steady_defect_fraction": float(late[:, 1].mean()), "steady_std": float(late[:, 1].std()),
            "pre_image_defect_fraction": float(late[:, 2].mean()), "mean_defect_fraction": sim.mean_defect(t_end / 2),
            "fills_per_s": float(sim.n_fills / t_end), "comp_sites": int(sim.comp.sum())}


def main():
    cfgs = []
    for T_img, eps, tau in itertools.product((0.03, 0.1, 0.3, 1.0), (1e-5, 1e-4, 1e-3), (300.0, 3000.0, 30000.0)):
        cfgs.append({"N": 1e5, "timing": "intermediate", "T_img": T_img, "eps": eps, "tau": tau, "t_img": 0.01,
                     "t_end": max(8.0, 8 * T_img)})
    for t_img in (0.003, 0.02):
        for T_img in (0.03, 0.1):
            cfgs.append({"N": 1e5, "timing": "intermediate", "T_img": T_img, "eps": 1e-4, "tau": 3000.0,
                         "t_img": t_img, "t_end": 8.0})
    for T_img in (0.03, 0.1, 0.3):
        for tm in ("intermediate", "fast"):
            cfgs.append({"N": 1e6, "timing": tm, "T_img": T_img, "eps": 1e-4, "tau": 3000.0, "t_img": 0.01,
                         "t_end": 12.0})
    cfgs.sort(key=lambda c: -c["N"])
    with ProcessPoolExecutor(4) as ex:
        rows = list(ex.map(one, cfgs))
    # fit on the 1e5, t_img = 10 ms points: rho = a eps + b (T_img + t_img)/tau
    sel = [r for r in rows if r["N"] == 1e5 and r["t_img"] == 0.01]
    X = np.array([[r["eps"], (r["T_img"] + r["t_img"]) / r["tau"]] for r in sel])
    y = np.array([r["steady_defect_fraction"] for r in sel])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    rel = np.abs(X @ coef - y) / np.maximum(y, 1e-9)
    fit = {"a": float(coef[0]), "b": float(coef[1]), "median_rel_error": float(np.median(rel)),
           "max_rel_error": float(rel.max())}
    for r in sorted(rows, key=lambda r: (r["N"], r["timing"], r["t_img"], r["eps"], r["tau"], r["T_img"])):
        print(int(r["N"]), r["timing"][:4], "t_img", r["t_img"], "eps %.0e tau %5.0f T %.2f" % (r["eps"], r["tau"], r["T_img"]),
              "| rho %.2e pre %.2e mean %.2e" % (r["steady_defect_fraction"], r["pre_image_defect_fraction"],
                                                 r["mean_defect_fraction"]), "| fills/s %.0f" % r["fills_per_s"], flush=True)
    print("fit", fit, flush=True)
    f, ax = fig(1, 1, h=3.0)
    for eps in (1e-5, 1e-4, 1e-3):
        for tau in (300.0, 30000.0):
            s = sorted([r for r in sel if r["eps"] == eps and r["tau"] == tau], key=lambda r: r["T_img"])
            ax[0, 0].loglog([r["T_img"] for r in s], [max(r["steady_defect_fraction"], 1e-7) for r in s], "-o", ms=3,
                            label=f"eps {eps:g}, tau {tau:g} s")
    ax[0, 0].set_xlabel("imaging period (s)")
    ax[0, 0].set_ylabel("steady-state defect fraction")
    ax[0, 0].legend(fontsize=6)
    save_result("m25_imaging", script=SCRIPT, fig=f,
                inputs={"layout": "4 x 4 um, 3:1 computation:reservoir rows", "block": [32, 4], "f_load": 0.85,
                        "timings": {"intermediate": MID, "fast": FAST}, "p_move": 1e-4, "t_reload_s": 1.0},
                outputs={"rows": rows, "fit_1e5": fit})
    return rows


def latency_fit():
    """Fit the time-averaged defect fraction of results/m25_imaging.json to the hole-lifetime model

        rho_mean = eps (T + d)/T + (T/2 + d)/tau,   T = T_img + t_img (imaging period),

    in which an imaging loss is unseen for one period and then repaired after a mean delay d, and a vacuum loss
    waits T/2 on average for the next image and then d. One d per (N, timing) group, least squares in log rho."""
    import json

    from scipy.optimize import minimize_scalar

    from picsim.provenance import RESULTS
    rows = json.load(open(RESULTS / "m25_imaging.json"))["outputs"]["rows"]

    def model(r, d):
        T = r["T_img"] + r["t_img"]
        return r["eps"] * (T + d) / T + (T / 2 + d) / r["tau"]

    out = {}
    for key in sorted({(r["N"], r["timing"]) for r in rows}):
        g = [r for r in rows if (r["N"], r["timing"]) == key]
        y = np.log([r["mean_defect_fraction"] for r in g])
        res = minimize_scalar(lambda d: np.sum((np.log([model(r, d) for r in g]) - y) ** 2), bounds=(0, 2), method="bounded")
        rel = [model(r, res.x) / r["mean_defect_fraction"] - 1 for r in g]
        d0 = [model(r, 0.0) / r["mean_defect_fraction"] - 1 for r in g]
        out[f"{int(key[0])} {key[1]}"] = {"d_s": float(res.x), "n_points": len(g),
                                           "median_abs_rel_error": float(np.median(np.abs(rel))),
                                           "max_abs_rel_error": float(np.max(np.abs(rel))),
                                           "median_rel_error_with_d0": float(np.median(d0))}
        print(key, out[f"{int(key[0])} {key[1]}"], flush=True)
    save_result("m25_imaging_fit", script=SCRIPT,
                inputs={"source": "results/m25_imaging.json", "model": latency_fit.__doc__.strip()},
                outputs={"fits": out})
    return out


if __name__ == "__main__":
    {"all": main, "fit": latency_fit}[sys.argv[1] if len(sys.argv) > 1 else "all"]()
