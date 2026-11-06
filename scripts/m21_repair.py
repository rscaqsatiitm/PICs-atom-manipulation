"""Module 21: repair rounds after the rearrangement cycle, and scheduler variants.

Parts:
  modes   pass-selection modes of the balancing stage (greedy, rate, lookahead) on 40 loadings
  repair  stochastic losses (per move, vacuum during every interval, imaging-induced), then rounds of
          (image, direct repair from reservoir atoms); defect-free probability and time per round

Outputs: results/m21_modes.json, results/m21_repair.json
"""
from __future__ import annotations

import itertools
import sys

import numpy as np

from picsim.provenance import save_result
from picsim.rearrange import Timing, cycle_with_repair, schedule

SCRIPT = "scripts/m21_repair.py"
R, C = 50, 45
R_T, C_T, S = 32, 32, 16
T_IMG = 10e-3  # s per image (CHOICE)


def modes():
    out = {}
    for mode in ("greedy", "rate", "lookahead"):
        P, T = [], []
        for sd in range(40):
            occ = np.random.default_rng(sd).random((R, C)) < 0.5
            res, _, _ = schedule(occ, R_T, C_T, S, mode=mode)
            P.append(res.passes_vertical + res.passes_horizontal)
            T.append(res.time_s * 1e3)
        out[mode] = {"passes_mean": float(np.mean(P)), "time_ms_mean": float(np.mean(T)), "time_ms_max": float(np.max(T))}
        print(mode, out[mode], flush=True)
    save_result("m21_modes", script=SCRIPT, inputs={"array": [R, C], "target": [R_T, C_T], "S": S, "loadings": 40,
                                                    "p_load": 0.5},
                outputs=out, notes="Lower bound on balancing passes about 4 (Sec. 2 of docs/repair.md); compaction 2.")
    return out


def repair_study():
    loads = [np.random.default_rng(100 + sd).random((R, C)) < 0.5 for sd in range(40)]
    scheds = [schedule(o, R_T, C_T, S) for o in loads]
    rows = []
    for p_move, tau, eps in itertools.product((1e-3, 1e-4), (30.0, 300.0), (1e-3, 1e-4)):
        rng = np.random.default_rng(7)
        D, NP, TR, TC = [], [], [], []
        for o, sc in zip(loads, scheds):
            for _ in range(25):
                r = cycle_with_repair(o, R_T, C_T, S, rng, p_move, tau, eps, t_img=T_IMG, rounds=2, sched=sc)
                D.append(r["defects"])
                NP.append(r["repair_passes"])
                TR.append(r["t_round_ms"])
                TC.append(r["t_cycle_ms"])
        D, NP, TR = np.array(D), np.array(NP), np.array(TR)
        row = {"p_move": p_move, "tau_s": tau, "eps_img": eps,
               "P_defect_free": [float((D[:, i] == 0).mean()) for i in range(3)],
               "defects_mean": [float(D[:, i].mean()) for i in range(3)],
               "repair_passes_mean": [float(NP[:, i].mean()) for i in range(2)],
               "t_round_ms_mean": [float(TR[:, i].mean()) for i in range(2)], "t_cycle_ms_mean": float(np.mean(TC)),
               "t_total_ms_mean": float(np.mean(TC) + TR.sum(axis=1).mean())}
        rows.append(row)
        print({k: (np.round(v, 3).tolist() if isinstance(v, list) else v) for k, v in row.items()}, flush=True)
    save_result("m21_repair", script=SCRIPT,
                inputs={"array": [R, C], "target": [R_T, C_T], "S": S, "loadings": 40, "loss_draws_per_loading": 25,
                        "t_img_s": T_IMG, "timing": Timing().__dict__,
                        "model": "detection precedes imaging-induced loss; repairs move one reservoir atom each"},
                outputs={"rows": rows},
                notes="Index 0: after the rearrangement cycle; 1 and 2: after each image-and-repair round. No final "
                      "verification image.")
    return rows


if __name__ == "__main__":
    for p in sys.argv[1:] or ["modes", "repair"]:
        {"modes": modes, "repair": repair_study}[p]()
