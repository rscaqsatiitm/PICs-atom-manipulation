"""Module 23: assembly and steady-state refilling of 1e4, 1e5 and 1e6-site arrays (picsim/megasort.py).

For each configuration (array size, pitch, loading and reservoir layout, engine block, vacuum lifetime, timing):
time to reach 1e-2 and 1e-3 defect fraction from a random loading, steady-state defect fraction (mean over the
last 3 s), refills per second, reservoir consumption, and the latency estimate it implies. A zoned layout is the
same simulation per zone (independent zones), so its time is that of the smaller array.

Output: results/m23_steady.json (+ .png)
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from picsim.megasort import Block, Layout, Params, Sim
from picsim.plotting import fig
from picsim.provenance import save_result
from picsim.rearrange import Timing

SCRIPT = "scripts/m23_steady.py"
LAYOUTS = {0.5: (3, 1), 0.85: (4, 3)}  # loading -> band period, computation rows per band
TIMINGS = {"slow (ramps 200, lanes 400 us)": Timing(t_ramp=200e-6, t_lane=400e-6),
           "fast (ramps 50, lanes 100, x 200 us/40 um)": Timing(t_ramp=50e-6, t_lane=100e-6, t_x40=200e-6, t_min=100e-6)}


def one(cfg):
    N, pitch, f, (S, Mx), tau, tname = cfg["N"], cfg["pitch"], cfg["f"], cfg["block"], cfg["tau"], cfg["timing"]
    B, Bc = LAYOUTS[f]
    C = int(np.ceil(np.sqrt(N)))
    R = int(np.ceil(N / C * B / Bc))
    seg = int(180.0 // pitch) + 1
    tm = TIMINGS[tname]
    tm = Timing(**{**tm.__dict__, "pitch": pitch * 1e-6})
    sim = Sim(Layout(R, C, B, Bc), Block(S, Mx, seg), Params(tau=tau, f_load=f, timing=tm), seed=cfg.get("seed", 1))
    t_end = cfg.get("t_end", 12.0)
    t3 = sim.run(t_end, target=1e-3)
    d = np.array(sim.log)
    t2 = next((t for t, x in d if x <= 1e-2), None)
    late = d[d[:, 0] >= t_end - 3.0, 1]
    # refill rate over the last 3 s from move counts is not logged per time; use fills / (time) after assembly as
    # an upper estimate, and the analytic generation rate N (1/tau + eps/T_img)
    return {**{k: v for k, v in cfg.items() if k != "seed"}, "rows": R, "cols": C, "comp_sites": int(sim.comp.sum()),
            "seg_sites": seg, "positions": len(sim.positions), "t_1e-2_s": t2, "t_1e-3_s": t3,
            "steady_defect_fraction": float(late.mean()) if len(late) else None,
            "steady_defect_fraction_std": float(late.std()) if len(late) else None,
            "fills_total": sim.n_fills, "moves_total": sim.n_moves, "consumed_total": sim.consumed,
            "visits": int(getattr(sim, "n_visits", 0)), "t_end_s": t_end,
            "generation_per_s": float(sim.comp.sum() * (1 / tau + sim.P.eps_img / sim.P.T_img)),
            "curve": d[:: max(1, len(d) // 120)].tolist()}


def main():
    cfgs = []
    for N in (1e4, 1e5, 1e6):
        for f in (0.5, 0.85):
            for blk in ((16, 1), (32, 4)):
                for tau in (300.0, 3000.0):
                    cfgs.append({"N": N, "pitch": 4.0, "f": f, "block": blk, "tau": tau,
                                 "timing": "slow (ramps 200, lanes 400 us)", "t_end": 12.0 if N < 1e6 else 20.0})
    for N in (1e5, 1e6):  # 3 um pitch and fast timing, best block
        cfgs.append({"N": N, "pitch": 3.0, "f": 0.85, "block": (32, 4), "tau": 3000.0,
                     "timing": "slow (ramps 200, lanes 400 us)", "t_end": 12.0 if N < 1e6 else 20.0})
        cfgs.append({"N": N, "pitch": 4.0, "f": 0.85, "block": (32, 4), "tau": 3000.0,
                     "timing": "fast (ramps 50, lanes 100, x 200 us/40 um)", "t_end": 12.0 if N < 1e6 else 20.0})
    cfgs.sort(key=lambda c: -c["N"])  # long runs first
    with ProcessPoolExecutor(4) as ex:
        rows = list(ex.map(one, cfgs))
    rows.sort(key=lambda r: (r["N"], r["pitch"], r["f"], r["block"], r["tau"], r["timing"]))
    for r in rows:
        print(int(r["N"]), r["pitch"], r["f"], r["block"], int(r["tau"]), r["timing"][:4], "| t1e-2", r["t_1e-2_s"],
              "t1e-3", r["t_1e-3_s"], "| steady %.2e" % r["steady_defect_fraction"], "| visits", r["visits"], flush=True)
    f, ax = fig(1, 1, h=3.2)
    for r in rows:
        if r["pitch"] == 4.0 and r["f"] == 0.85 and r["tau"] == 3000.0 and r["timing"].startswith("slow"):
            c = np.array(r["curve"])
            ax[0, 0].semilogy(c[:, 0], np.maximum(c[:, 1], 1e-6), label=f"N={int(r['N']):.0e}, block {r['block']}")
    ax[0, 0].set_xlabel("time (s)")
    ax[0, 0].set_ylabel("defect fraction")
    ax[0, 0].legend(fontsize=6)
    save_result("m23_steady", script=SCRIPT, fig=f,
                inputs={"layouts": {str(k): {"band_period": v[0], "comp_rows_per_band": v[1]} for k, v in LAYOUTS.items()},
                        "params": {"eps_img": 1e-4, "p_move": 1e-4, "t_img_s": 0.01, "T_img_s": 0.1,
                                   "t_reload_s": 1.0, "max_passes_per_visit": 4},
                        "timings": {k: v.__dict__ for k, v in TIMINGS.items()}, "segment_um": 180.0},
                outputs={"rows": rows},
                notes="Single large-field objective assumed (all sites reachable by the block). Zones behave as "
                      "independent smaller arrays.")
    return rows


if __name__ == "__main__":
    main()
