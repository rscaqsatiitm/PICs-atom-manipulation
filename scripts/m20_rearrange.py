"""Module 20: rearrangement cycles for a stacked focal-plane OPA serving a 50 x 45 static array.

Monte Carlo over random loadings: the scheduler of picsim/rearrange.py builds a centred defect-free block using
vertical-balancing passes (common row shift per pass, independent x per atom) and horizontal compaction passes.
Per configuration: passes, cycle time, moves per atom, accumulated heating, and the probability that the block is
complete after per-move loss and vacuum loss. Heating and loss per move come from modules 18 and 19.

Output: results/m20_rearrange.json (+ .png)
"""
from __future__ import annotations

import json
import sys

import numpy as np

from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result
from picsim.rearrange import Timing, schedule

SCRIPT = "scripts/m20_rearrange.py"
R, C = 50, 45  # rows x sites per row: +/-100 um field rows, 192-port chips (docs/scaling_beyond.md)


def run(S, R_t, C_t, p_load, n_mc=40, timing=Timing(), seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_mc):
        occ = rng.random((R, C)) < p_load
        res, _, _ = schedule(occ, R_t, C_t, S, timing)
        m = np.array(res.moves_per_atom) if res.moves_per_atom else np.zeros(1)
        out.append({"filled": res.target_filled, "deficit": res.deficit, "pv": res.passes_vertical,
                    "ph": res.passes_horizontal, "t_ms": res.time_s * 1e3, "moves_mean": float(m.mean()),
                    "moves_max": int(m.max()), "moves_total": int(m.sum()),
                    "frac_0": float((m == 0).mean()), "frac_ge2": float((m >= 2).mean())})
    return out


def summarise(rows, R_t, C_t):
    f = np.array([r["filled"] for r in rows])
    ok = [r for r in rows if r["filled"]] or rows
    t = np.array([r["t_ms"] for r in ok])
    mt = np.array([r["moves_total"] for r in ok])
    mm = np.array([r["moves_mean"] for r in ok])
    mx = np.array([r["moves_max"] for r in ok])
    s = {"filled_fraction_before_loss": float(f.mean()), "deficit_mean_when_unfilled": float(np.mean([r["deficit"] for r in rows if not r["filled"]])) if (~f).any() else 0.0,
         "passes_vertical_mean": float(np.mean([r["pv"] for r in ok])), "passes_horizontal_mean": float(np.mean([r["ph"] for r in ok])),
         "time_ms_mean": float(t.mean()), "time_ms_max": float(t.max()), "moves_per_atom_mean": float(mm.mean()),
         "moves_per_atom_max": int(mx.max()), "moves_total_mean": float(mt.mean()),
         "frac_atoms_moved_twice_or_more": float(np.mean([r["frac_ge2"] for r in ok]))}
    # defect-free probability from transport alone (every block atom survives its moves), and with vacuum loss of
    # every block atom during the cycle for two lifetimes
    for p in (1e-3, 3e-4, 1e-4):
        s[f"P_defect_free_moves_p{p:g}"] = float(np.mean(np.exp(-p * mt)) * f.mean())
        s[f"expected_defects_moves_p{p:g}"] = float(np.mean(p * mt))
    for tau in (30.0, 300.0):
        s[f"expected_defects_vacuum_tau{tau:g}s"] = float(np.mean(R_t * C_t * t * 1e-3 / tau))
        s[f"P_defect_free_p1e-4_tau{tau:g}s"] = float(np.mean(np.exp(-1e-4 * mt - R_t * C_t * t * 1e-3 / tau)) * f.mean())
    return s


def main():
    cases = []
    for S, label in ((16, "one stack (16 rows)"), (32, "two stacks (32 rows)")):
        for (R_t, C_t), p in (((30, 30), 0.5), ((32, 32), 0.5), ((32, 32), 0.6), ((36, 36), 0.6)):
            rows = run(S, R_t, C_t, p)
            s = summarise(rows, R_t, C_t)
            cases.append({"stack_rows": S, "label": label, "target": [R_t, C_t], "p_load": p, **s})
            print(label, R_t, C_t, p, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()}, flush=True)
    # timing sensitivity for the base case (one stack, 32 x 32, p = 0.5)
    sens = []
    for name, tm in (("base (ramp 100, lane 400, x 400 us/40 um)", Timing()),
                     ("fast lane 100 us", Timing(t_lane=100e-6)),
                     ("fast all: ramp 50, lane 100, x 200 us/40 um", Timing(t_ramp=50e-6, t_lane=100e-6, t_x40=200e-6,
                                                                        t_min=100e-6)),
                     ("slow: ramp 200, lane 400", Timing(t_ramp=200e-6)),
                     ("scanner settle 1 ms", Timing(t_reposition=1e-3))):
        s = summarise(run(16, 32, 32, 0.5, n_mc=20, timing=tm), 32, 32)
        sens.append({"timing": name, "time_ms_mean": s["time_ms_mean"],
                     "expected_defects_vacuum_tau30s": s["expected_defects_vacuum_tau30s"]})
        print(sens[-1], flush=True)
    # heating: per-move energy gain (m19 quantum, slow lane; m18 classical, fast lane)
    heat = {"slow lane (ramp 200, lane 400 us), 0.1/0.3 mK, quantum": 0.12 + 0.17,
            "fast lane (100/100 us), 0.1/0.5 mK, quantum": 0.84 + 0.37,
            "fast lane (100/100 us), 0.1/0.3 mK, classical 3D": 10.0}
    base = [c for c in cases if c["stack_rows"] == 16 and c["target"] == [32, 32] and c["p_load"] == 0.5][0]
    heating = {k: {"mean_uK": v * base["moves_per_atom_mean"], "worst_atom_uK": v * base["moves_per_atom_max"]}
               for k, v in heat.items()}
    f, ax = fig(1, 1, h=3.0)
    for S in (16, 32):
        sel = [c for c in cases if c["stack_rows"] == S]
        ax[0, 0].plot(range(len(sel)), [c["time_ms_mean"] for c in sel], "-o", ms=3, label=f"{S} rows addressed")
    ax[0, 0].set_xticks(range(4))
    ax[0, 0].set_xticklabels(["30x30, p .5", "32x32, p .5", "32x32, p .6", "36x36, p .6"], fontsize=7)
    ax[0, 0].set_ylabel("cycle time (ms)")
    ax[0, 0].legend(fontsize=7)
    save_result("m20_rearrange", script=SCRIPT, fig=f,
                inputs={"array": [R, C], "monte_carlo": 40, "timing": Timing().__dict__, "tau_vac_s": [30, 300],
                        "loss_per_move": [1e-3, 3e-4, 1e-4], "heating_per_move_uK": heat},
                outputs={"cases": cases, "timing_sensitivity": sens, "heating_base_case": heating},
                notes="One rearrangement cycle, no repair round. P_defect_free includes per-move loss and vacuum loss "
                      "during the cycle, not imaging or scattering.")
    return cases


if __name__ == "__main__":
    main()
