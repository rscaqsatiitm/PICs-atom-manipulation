"""Module 28: the PIC rearrangement engine against crossed AODs, a single moving tweezer and an all-parallel ceiling.

The same steady-state simulator (picsim/megasort.py), layout, loading, timing physics, losses and imaging are used
for every method; only the move constraint and the field change (Block.mode):
  PIC 16 x 1, PIC 32 x 4   stacked PIC engines (one stack; 128 engines), 180 um segments (megascale.md)
  AOD row-wise             crossed AODs over a square field: one row tone and many column tones per pass (lower
                           bound on AOD throughput: multi-row passes with a shared column pattern are not used)
  single tweezer           one moving spot in the same field, one atom per pass
  all-parallel             every known hole matched in one pass with independent 2D moves (SLM-like ceiling)

Parts:
  main      assembly time, steady-state time-averaged defect fraction, refills, tweezers per pass, heating; 1e3, 1e4,
            1e5 sites; tau 3000 / 300 s; cold steps; fast steps; AOD field sensitivity; several seeds
  hardware  laser power (peak simultaneous moving tweezers x power per tweezer), fast control channels, field
            positions (DERIVED from main and the stated per-tweezer powers)
  figure    defect fraction vs time at 1e4 sites; steady-state defects and assembly time per method

Outputs: results/m28_main.json, results/m28_hardware.json, results/m28_figure.png
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from picsim.megasort import Block, Layout, Params, Sim
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result
from picsim.rearrange import Timing

SCRIPT = "scripts/m28_aod_baseline.py"
B, BC, F_LOAD = 4, 3, 0.85  # band period, computation rows per band, loading (m23 layout for 85%)
PITCH = 4e-6
TIMINGS = {"slow": Timing(t_ramp=200e-6, t_lane=400e-6, pitch=PITCH),  # m23 cold steps (0.3 uK per move)
           "fast": Timing(t_ramp=50e-6, t_lane=100e-6, t_x40=200e-6, t_min=100e-6, pitch=PITCH)}
P_MOVE = {"slow": 1e-4, "fast": 4e-3}  # m23 (fast: m18_scan SIM; slow: CHOICE)
HEAT = {"slow": 0.3, "fast": 10.0}  # uK per move (m23)
# AOD field: 1064 nm TeO2 deflectors list 90 (3 mm) and 150 (5 mm) resolvable spots (45035-x-6.5DEG-1.06, SOURCE:
# sintec.sg/ao/497.html, read 2026-11-14); one resolution cell at NA 0.5 is 1.06 um, so 150 spots span about
# 160 um = 40 sites at 4 um. 40 is the base field; 24 (90 spots) and 64 / 128 are the sensitivity cases.
AOD_FIELD = 40
SIZES = {1e3: {"t_end": 200.0, "seeds": 4}, 1e4: {"t_end": 60.0, "seeds": 3}, 1e5: {"t_end": 40.0, "seeds": 2}}


def geometry(N):
    C = int(np.ceil(np.sqrt(N)))
    R = int(np.ceil(N / C * B / BC))
    return R, C


def block_for(method, R, C):
    if method == "PIC 16 x 1":
        return Block(16, 1, 46)
    if method == "PIC 32 x 4":
        return Block(32, 4, 46)
    if method.startswith("AOD row-wise"):
        n = int(method.split("field ")[1])
        return Block(n, 1, n, mode="aod_row")
    if method.startswith("single tweezer"):
        n = int(method.split("field ")[1])
        return Block(n, 1, n, mode="single")
    if method == "all-parallel":
        return Block(R, 1, C, mode="parallel")
    raise ValueError(method)


def one(cfg):
    N, method, tau, tname, seed = cfg["N"], cfg["method"], cfg["tau"], cfg["timing"], cfg["seed"]
    R, C = geometry(N)
    t_end = cfg["t_end"]
    P = Params(tau=tau, f_load=F_LOAD, timing=TIMINGS[tname], p_move=P_MOVE[tname], heat_slow=HEAT[tname],
               skip_reposition_same=True)
    sim = Sim(Layout(R, C, B, BC), block_for(method, R, C), P, seed=seed)
    t0 = t_end / 2
    t3 = sim.run(t0, target=1e-3)
    n_pass0, fills0, t_mid = len(sim.pass_tweezers), sim.n_fills, sim.t
    sim.run(t_end)
    d = np.array(sim.log)
    t2 = next((r[0] for r in d if r[1] <= 1e-2), None)
    if t3 is None:
        t3 = next((r[0] for r in d if r[1] <= 1e-3), None)
    pt = np.array(sim.pass_tweezers) if sim.pass_tweezers else np.zeros(1)
    pt_late = pt[n_pass0:] if len(pt) > n_pass0 else np.zeros(1)
    late = d[d[:, 0] >= t0]
    hs = sim.heat_stats()
    return {**cfg, "rows": R, "cols": C, "comp_sites": int(sim.comp.sum()), "positions": len(sim.positions),
            "t_1e-2_s": t2, "t_1e-3_s": t3, "assembled_before_window": bool(t3 is not None and t3 <= t0),
            "steady_mean_defect": sim.mean_defect(t0), "steady_at_image_defect": float(late[:, 1].mean()),
            "steady_pre_image_defect": float(late[:, 2].mean()), "final_defect": float(d[-1, 1]),
            "refills_per_s_window": (sim.n_fills - fills0) / max(sim.t - t_mid, 1e-9),
            "passes": int(len(sim.pass_tweezers)), "tweezers_per_pass_max": int(pt.max()),
            "tweezers_per_pass_p99": float(np.percentile(pt, 99)), "tweezers_per_pass_mean": float(pt.mean()),
            "tweezers_per_pass_mean_window": float(pt_late.mean()), "moves_total": sim.n_moves,
            "fills_total": sim.n_fills, "consumed_total": sim.consumed, "heat_mean_uK": hs["mean_uK"],
            "heat_p99_uK": hs["p99_uK"], "moved_fraction": hs["moved_fraction"],
            "generation_per_s": float(sim.comp.sum() * (1 / tau + P.eps_img / P.T_img)),
            "curve": d[:: max(1, len(d) // 200)].tolist()}


METHODS = ["PIC 16 x 1", "PIC 32 x 4", f"AOD row-wise, field {AOD_FIELD}", "AOD row-wise, field 128",
           f"single tweezer, field {AOD_FIELD}", "all-parallel"]


def configs():
    cfgs = []
    for N, s in SIZES.items():
        for method in METHODS:
            for tau in (3000.0, 300.0):
                for seed in range(1, s["seeds"] + 1):
                    cfgs.append({"set": "cold", "N": N, "method": method, "tau": tau, "timing": "slow", "seed": seed,
                                 "t_end": s["t_end"]})
    for N in (1e4, 1e5):
        for method in METHODS:
            for seed in (1, 2):
                cfgs.append({"set": "fast", "N": N, "method": method, "tau": 3000.0, "timing": "fast", "seed": seed,
                             "t_end": SIZES[N]["t_end"]})
        for n in (24, 64):
            for seed in (1, 2):
                cfgs.append({"set": "aod field", "N": N, "method": f"AOD row-wise, field {n}", "tau": 3000.0,
                             "timing": "slow", "seed": seed, "t_end": SIZES[N]["t_end"]})
    return cfgs


def _cost(c):  # longest first
    w = {"single": 50.0, "all-p": 3.0}.get(c["method"][:6], 1.0)
    return -c["N"] * c["t_end"] * w


def summarise(rows):
    groups = {}
    for r in rows:
        groups.setdefault((r["set"], r["N"], r["method"], r["tau"], r["timing"]), []).append(r)
    out = []
    for (st, N, m, tau, tn), rs in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][3], kv[0][2])):
        def ms(key):
            v = [x[key] for x in rs if x[key] is not None]
            return (float(np.mean(v)) if v else None, float(np.std(v)) if len(v) > 1 else None, len(v))
        out.append({"set": st, "N": N, "method": m, "tau": tau, "timing": tn, "seeds": len(rs),
                    "comp_sites": rs[0]["comp_sites"], "positions": rs[0]["positions"],
                    "t_1e-2_s": ms("t_1e-2_s"), "t_1e-3_s": ms("t_1e-3_s"),
                    "assembled_runs": sum(1 for x in rs if x["t_1e-3_s"] is not None),
                    "steady_mean_defect": ms("steady_mean_defect"), "steady_at_image_defect": ms("steady_at_image_defect"),
                    "final_defect": ms("final_defect"), "refills_per_s_window": ms("refills_per_s_window"),
                    "tweezers_per_pass_max": max(x["tweezers_per_pass_max"] for x in rs),
                    "tweezers_per_pass_mean": ms("tweezers_per_pass_mean"),
                    "tweezers_per_pass_mean_window": ms("tweezers_per_pass_mean_window"),
                    "heat_mean_uK": ms("heat_mean_uK"), "generation_per_s": rs[0]["generation_per_s"]})
    return out


def main():
    cfgs = sorted(configs(), key=_cost)
    with ProcessPoolExecutor(4) as ex:
        rows = list(ex.map(one, cfgs))
    summ = summarise(rows)
    for s in summ:
        f = lambda v: "-" if v[0] is None else ("%.3g" % v[0])  # noqa: E731
        print(s["set"], "%.0e" % s["N"], "tau", int(s["tau"]), s["method"][:30].ljust(30), "| t1e-3", f(s["t_1e-3_s"]),
              "(%d/%d)" % (s["assembled_runs"], s["seeds"]), "| steady", f(s["steady_mean_defect"]),
              "+/-", "-" if s["steady_mean_defect"][1] is None else "%.1g" % s["steady_mean_defect"][1],
              "| tw max", s["tweezers_per_pass_max"], "| refills/s", f(s["refills_per_s_window"]), flush=True)
    save_result("m28_main", script=SCRIPT,
                inputs={"layout": {"band_period": B, "comp_rows_per_band": BC, "loading": F_LOAD, "pitch_um": 4.0},
                        "timings": {k: v.__dict__ for k, v in TIMINGS.items()}, "p_move": P_MOVE, "heat_uK_per_move": HEAT,
                        "params": {"eps_img": 1e-4, "t_img_s": 0.01, "T_img_s": 0.1, "t_reload_s": 1.0,
                                   "max_passes_per_visit": 4, "max_passes_single": 400,
                                   "skip_reposition_same": True, "t_reposition_s": 300e-6},
                        "aod_field_sites": {"base": AOD_FIELD, "sensitivity": [24, 64, 128],
                                            "source": "1064 nm TeO2 AODs with 90 / 150 resolvable spots (sintec.sg/ao/"
                                                      "497.html, read 2026-11-14); 1.06 um resolution cell at NA 0.5"},
                        "sizes": {str(k): v for k, v in SIZES.items()},
                        "steady_window": "second half of each run (t_end / 2 to t_end)"},
                outputs={"summary": summ, "runs": [{k: v for k, v in r.items() if k != "curve"} for r in rows],
                         "curves_1e4_tau3000_seed1": {r["method"]: r["curve"] for r in rows
                                                      if r["N"] == 1e4 and r["tau"] == 3000.0 and r["seed"] == 1
                                                      and r["set"] == "cold"}},
                notes="One objective field covering the array is assumed for every method. AOD row-wise is a lower "
                      "bound on crossed-AOD throughput; all-parallel is an idealised ceiling, not a hardware model.")
    return summ


# ---------------------------------------------------------------------------------------------- hardware
P_TW_PIC_W = 57.7e-3  # fibre input per 1.2 mK tweezer, reflective collimator, M1 silver (results/m27_wave.json, SIM)
P_TW_ATOMS_W = 6.6e-3 * 1.2  # proposal's 6.6 mW per tweezer gives 1.0 mK for w0 = 1.04 um (docs/findings.md), so
# 7.9 mW at the atoms for 1.2 mK (DERIVED from USER)
ETA_AOD = 0.65  # per deflector, mid-band, single tone, 1064 nm (45035-3-6.5DEG-1.06-XY: >65%, SOURCE as above)
ACT_ENGINE = 383  # fast actuators per 192-port engine (megascale.md Sec. 4); D1 (128 ports): 255


def hardware():
    summ = json.loads((RESULTS / "m28_main.json").read_text())["outputs"]["summary"]
    rows = []
    for s in summ:
        if s["set"] != "cold" or s["tau"] != 3000.0:
            continue
        m, n = s["method"], s["tweezers_per_pass_max"]
        if m.startswith("PIC"):
            engines = 16 if "16" in m else 128
            p_tw, chans = P_TW_PIC_W, engines * ACT_ENGINE
            note = f"{engines} engines x {ACT_ENGINE} fast actuators"
        elif m.startswith("AOD") or m.startswith("single"):
            p_tw, chans = P_TW_ATOMS_W / ETA_AOD**2, 2
            note = "2 RF channels (multi-tone synthesis); AOD pair efficiency 0.65^2 = 0.42 (single tone, mid-band)"
        else:
            p_tw, chans = None, None
            note = "SLM-like ceiling: not a hardware model"
        rows.append({"N": s["N"], "method": m, "peak_tweezers_per_pass": n,
                     "P_per_tweezer_mW": None if p_tw is None else p_tw * 1e3,
                     "P_transport_W_at_1.2mK": None if p_tw is None else n * p_tw,
                     "P_transport_W_at_0.3mK": None if p_tw is None else n * p_tw * 0.25,
                     "fast_control_channels": chans, "field_positions": s["positions"], "note": note})
        print(rows[-1], flush=True)
    sweep = {f"{e:.1f}": P_TW_ATOMS_W / e**2 * 1e3 for e in (0.55, 0.65, 0.75, 0.9)}
    save_result("m28_hardware", script=SCRIPT,
                inputs={"P_tweezer_PIC_mW": P_TW_PIC_W * 1e3, "P_tweezer_atoms_mW": P_TW_ATOMS_W * 1e3,
                        "eta_AOD_per_deflector": ETA_AOD, "actuators_per_engine": ACT_ENGINE,
                        "depth": "1.2 mK reference; the transport model moves at 0.3 mK (power scales linearly)",
                        "objective_and_relay": "transmission not included for either method"},
                outputs={"rows": rows, "AOD_P_per_tweezer_mW_vs_eta_per_deflector": sweep},
                notes="Laser power is provisioned for the peak number of simultaneously moving tweezers seen in the "
                      "simulation (assembly included). AOD multi-tone efficiency and intermodulation are not included "
                      "(per-tone efficiency falls with tone count; [VALUE?]).")
    return rows


# ---------------------------------------------------------------------------------------------- figure
def figure():
    o = json.loads((RESULTS / "m28_main.json").read_text())["outputs"]
    col = {m: f"C{i}" for i, m in enumerate(METHODS)}
    f, ax = fig(1, 2, h=3.2)
    for m in METHODS:
        c = np.array(o["curves_1e4_tau3000_seed1"][m])
        ax[0, 0].semilogy(c[:, 0], np.maximum(c[:, 1], 1e-5), label=m, lw=0.9, color=col[m])
    ax[0, 0].set_xlim(0, 4)
    ax[0, 0].set_xlabel("time (s)")
    ax[0, 0].set_ylabel("defect fraction at image (1e4 sites)")
    ax[0, 0].legend(fontsize=5.5)
    summ = [s for s in o["summary"] if s["set"] == "cold" and s["tau"] == 3000.0]
    Ns = sorted({s["N"] for s in summ})
    w = 0.13
    for i, m in enumerate(METHODS):
        y = [next((s["t_1e-3_s"][0] for s in summ if s["N"] == N and s["method"] == m), np.nan) for N in Ns]
        ax[0, 1].bar(np.arange(len(Ns)) + (i - 2.5) * w, [np.nan if v is None else v for v in y], w, color=col[m])
    ax[0, 1].axhline(0.12, color="0.5", lw=0.7)
    ax[0, 1].set_yscale("log")
    ax[0, 1].set_xticks(range(len(Ns)))
    ax[0, 1].set_xticklabels([f"{int(N):.0e}" for N in Ns])
    ax[0, 1].set_xlabel("sites")
    ax[0, 1].set_ylabel("assembly time to 1e-3 (s)")
    save_result("m28_figure", script=SCRIPT, fig=f, inputs={"source": "results/m28_main.json"}, outputs={},
                notes="Left: defect fraction just after each image, 1e4 sites, tau 3000 s, cold steps, seed 1. Right: "
                      "time to the first image with defect fraction <= 1e-3 from 85% random loading, mean over seeds, "
                      "tau 3000 s, cold steps; the grey line (0.12 s) is the first image after the start, the "
                      "resolution of this metric.")
    return f


if __name__ == "__main__":
    for part in sys.argv[1:] or ["main", "hardware", "figure"]:
        globals()[part]()
