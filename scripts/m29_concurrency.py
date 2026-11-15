"""Module 29: concurrency at 1e6 sites. K independent units, each serving its own zone of 1e6 / K sites.

Zones are independent (each holds its own reservoir rows and its own device), so K units on 1e6 sites are one zone
of N_z = 1e6 / K sites simulated with the m28 units and settings (picsim/megasort.py). The zone ladder runs from
one unit scanning the whole array (K = 1) to the fixed-field limit, where the zone equals what the unit reaches
without moving its field.

Parts:
  capacity  moves per second a unit sustains when demand far exceeds it (loss per image 0.1, reservoir refreshed
            every 0.1 s so supply does not bind), on 1e4 and 1e5-site zones; peak simultaneous tweezers when busy
  ladder    steady state from an assembled start (computation sites filled at t = 0): time-averaged defect fraction
            vs zone size, tau 3000 / 300 s, cold steps; assembly time from a random loading for zones <= 6.25e4
  hardware  vs K at 1e6 sites: defect fraction, capacity x K against demand, laser power, fast channels, units
            (separate beam paths into the objective); smallest K for 1.5e-4 and 2e-4
  figure    defect fraction vs K per unit at 1e6 sites

Outputs: results/m29_capacity.json, m29_ladder.json, m29_hardware.json, m29_figure.png
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from picsim.megasort import Block, Layout, Params, Sim
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

sys.path.insert(0, str(RESULTS.parent / "scripts"))
from m28_aod_baseline import (ACT_ENGINE, B, BC, ETA_AOD, F_LOAD, HEAT, P_MOVE, P_TW_ATOMS_W, P_TW_PIC_W,  # noqa: E402
                              TIMINGS)

SCRIPT = "scripts/m29_concurrency.py"
N_TOTAL = 1e6
UNITS = ["PIC 16 x 1", "PIC 32 x 4", "AOD 40", "AOD 128", "single 40"]
# fixed-field zones (rows x columns of sites, reservoir rows included), i.e. what the unit reaches without scanning
FIXED = {"PIC 16 x 1": (16, 46), "PIC 32 x 4": (32, 8 * 46), "AOD 40": (40, 40), "AOD 128": (128, 128),
         "single 40": (40, 40)}
LADDER = [1e6, 2.5e5, 6.25e4, 1.5625e4, 3.9e3]


def zone_dims(Nz):
    """Rows and columns of a square-ish zone holding Nz computation sites (band layout B, BC)."""
    C = int(np.ceil(np.sqrt(Nz)))
    R = int(np.ceil(Nz / C * B / BC))
    return R, C


def block(unit, R, C):
    if unit == "PIC 16 x 1":
        return Block(16, 1, 46)
    if unit == "PIC 32 x 4":
        return Block(32, 4, 46)
    if unit.startswith("AOD"):
        n = int(unit.split()[1])
        return Block(n, 1, n, mode="aod_row")
    if unit.startswith("single"):
        n = int(unit.split()[1])
        return Block(n, 1, n, mode="single")
    raise ValueError(unit)


def run(cfg):
    R, C = cfg["R"], cfg["C"]
    tn = cfg.get("timing", "slow")
    P = Params(tau=cfg["tau"], f_load=F_LOAD, timing=TIMINGS[tn], p_move=P_MOVE[tn], heat_slow=HEAT[tn],
               skip_reposition_same=True, eps_img=cfg.get("eps_img", 1e-4), t_reload=cfg.get("t_reload", 1.0))
    sim = Sim(Layout(R, C, B, BC), block(cfg["unit"], R, C), P, seed=cfg["seed"])
    if cfg.get("assembled", False):
        sim.occ[sim.comp] = True
    t_end = cfg["t_end"]
    t0 = cfg.get("t0", t_end / 3)
    t3 = sim.run(t0, target=1e-3)
    n_pass0, fills0, moves0, t_mid = len(sim.pass_tweezers), sim.n_fills, sim.n_moves, sim.t
    sim.run(t_end)
    d = np.array(sim.log)
    if t3 is None:
        t3 = next((r[0] for r in d if r[1] <= 1e-3), None)
    pt = np.array(sim.pass_tweezers[n_pass0:] or [0])
    return {**cfg, "comp_sites": int(sim.comp.sum()), "positions": len(sim.positions), "t_1e-3_s": t3,
            "steady_mean_defect": sim.mean_defect(t0), "steady_at_image_defect": float(d[d[:, 0] >= t0, 1].mean()),
            "refills_per_s": (sim.n_fills - fills0) / max(sim.t - t_mid, 1e-9),
            "moves_per_s": (sim.n_moves - moves0) / max(sim.t - t_mid, 1e-9),
            "demand_per_s": float(sim.comp.sum() * (1 / cfg["tau"] + P.eps_img / P.T_img)),
            "tweezers_per_pass_max_window": int(pt.max()), "tweezers_per_pass_mean_window": float(pt.mean())}


def _pool(cfgs):
    cfgs = sorted(cfgs, key=lambda c: -c["R"] * c["C"] * c["t_end"] * (30.0 if c["unit"].startswith("single") else 1.0))
    with ProcessPoolExecutor(4) as ex:
        return list(ex.map(run, cfgs))


# ---------------------------------------------------------------------------------------------- capacity
def capacity():
    cfgs = []
    for unit in UNITS:
        for Nz in (1e4, 1e5):
            R, C = zone_dims(Nz)
            for tn in ("slow", "fast"):
                cfgs.append({"unit": unit, "Nz": Nz, "R": R, "C": C, "tau": 3000.0, "timing": tn, "seed": 1,
                             "eps_img": 0.1, "t_reload": 0.1, "assembled": True, "t_end": 6.0, "t0": 2.0})
    rows = _pool(cfgs)
    for r in sorted(rows, key=lambda r: (r["unit"], r["Nz"], r["timing"])):
        print(r["unit"].ljust(11), "%.0e" % r["Nz"], r["timing"], "| moves %.0f /s, refills %.0f /s (demand %.0f) | tw max %d mean %.1f" % (
            r["moves_per_s"], r["refills_per_s"], r["demand_per_s"], r["tweezers_per_pass_max_window"],
            r["tweezers_per_pass_mean_window"]), flush=True)
    save_result("m29_capacity", script=SCRIPT,
                inputs={"eps_img": 0.1, "t_reload_s": 0.1, "start": "computation sites filled", "window": "2 to 6 s",
                        "units": UNITS, "timings": {k: v.__dict__ for k, v in TIMINGS.items()}},
                outputs={"rows": rows},
                notes="Demand (1 refill/s per site plus vacuum) exceeds every unit's capacity, so moves_per_s is the "
                      "throughput (capacity) of one unit, travel and field scanning included; refills_per_s is lower "
                      "because some reservoir atoms seen in the image are lost before they are moved.")
    return rows


# ---------------------------------------------------------------------------------------------- ladder
def _t_end(n_comp):
    return float(np.clip(600.0 / (n_comp * (1 / 3000.0 + 1e-3)), 15.0, 400.0))  # about 600 losses per run


def ladder():
    cfgs = []
    for unit in UNITS:
        sizes = [("ladder", Nz) + zone_dims(Nz) for Nz in LADDER] + [("fixed", None) + FIXED[unit]]
        for kind, Nz, R, C in sizes:
            if unit.startswith("single") and kind == "ladder" and Nz > 6.25e4:
                continue  # one moving spot cannot hold 1e5 or more sites (m28)
            n_comp = R * C * BC / B
            seeds = 2 if n_comp > 2e5 else (3 if n_comp > 2e4 else 6)
            for tau in (3000.0, 300.0):
                for seed in range(1, seeds + 1):
                    cfgs.append({"unit": unit, "kind": kind, "Nz": Nz, "R": R, "C": C, "tau": tau, "seed": seed,
                                 "assembled": True, "t_end": _t_end(n_comp) if n_comp < 2e5 else 15.0})
            if n_comp <= 7e4:  # assembly time from a random loading
                cfgs.append({"unit": unit, "kind": kind + " assembly", "Nz": Nz, "R": R, "C": C, "tau": 3000.0,
                             "seed": 1, "assembled": False, "t_end": 40.0 if unit.startswith("single") else 10.0})
    rows = _pool(cfgs)
    groups = {}
    for r in rows:
        groups.setdefault((r["unit"], r["kind"], r["R"], r["C"], r["tau"]), []).append(r)
    summ = []
    for (unit, kind, R, C, tau), rs in sorted(groups.items()):
        v = [x["steady_mean_defect"] for x in rs]
        K = N_TOTAL / rs[0]["comp_sites"]
        summ.append({"unit": unit, "kind": kind, "R": R, "C": C, "tau": tau, "comp_sites": rs[0]["comp_sites"],
                     "K_at_1e6": K, "positions": rs[0]["positions"], "seeds": len(rs),
                     "steady_mean_defect": float(np.mean(v)), "steady_std": float(np.std(v)) if len(v) > 1 else None,
                     "t_1e-3_s": rs[0]["t_1e-3_s"] if "assembly" in kind else None,
                     "refills_per_s": float(np.mean([x["refills_per_s"] for x in rs])),
                     "demand_per_s": rs[0]["demand_per_s"]})
        s = summ[-1]
        print(unit.ljust(11), kind.ljust(16), "tau %4d" % tau, "zone %4dx%-4d K %7.1f pos %4d | defect %.3e +/- %s | t1e-3 %s"
              % (R, C, K, s["positions"], s["steady_mean_defect"], "-" if s["steady_std"] is None else "%.1e" % s["steady_std"],
                 s["t_1e-3_s"]), flush=True)
    save_result("m29_ladder", script=SCRIPT,
                inputs={"units": UNITS, "ladder_comp_sites": LADDER, "fixed_zones": FIXED, "eps_img": 1e-4,
                        "T_img_s": 0.1, "timing": "slow (cold) steps", "start": "assembled (steady runs); random 85% "
                        "loading (assembly runs)", "window": "last two thirds of each run"},
                outputs={"summary": summ, "runs": rows},
                notes="K_at_1e6 = 1e6 / computation sites per zone. Zones are independent, so the time-averaged defect "
                      "fraction of one zone is that of 1e6 sites served by K such units.")
    return summ


# ---------------------------------------------------------------------------------------------- hardware
def hardware():
    lad = json.loads((RESULTS / "m29_ladder.json").read_text())["outputs"]["summary"]
    cap = json.loads((RESULTS / "m29_capacity.json").read_text())["outputs"]["rows"]
    p_tw = {"PIC": P_TW_PIC_W, "AOD": P_TW_ATOMS_W / ETA_AOD**2, "single": P_TW_ATOMS_W / ETA_AOD**2}
    chans = {"PIC 16 x 1": 16 * ACT_ENGINE, "PIC 32 x 4": 128 * ACT_ENGINE, "AOD 40": 2, "AOD 128": 2, "single 40": 2}
    rows, best = [], {}
    for s in lad:
        if "assembly" in s["kind"]:
            continue
        u = s["unit"]
        c = [r for r in cap if r["unit"] == u and r["timing"] == "slow"]
        cap_u = max(r["moves_per_s"] for r in c)  # moves/s per unit, largest of the 1e4 / 1e5 zones
        tw = max(r["tweezers_per_pass_max_window"] for r in c)  # peak when busy
        K = s["K_at_1e6"]
        P = K * tw * p_tw[u.split()[0]]
        row = {**{k: s[k] for k in ("unit", "kind", "tau", "R", "C", "comp_sites", "positions", "steady_mean_defect",
                                    "steady_std")},
               "K": K, "capacity_per_unit": cap_u, "capacity_total": K * cap_u,
               "demand_1e6": N_TOTAL * (1 / s["tau"] + 1e-4 / 0.1), "peak_tweezers_per_unit": tw,
               "P_transport_W_1.2mK": P, "P_transport_W_0.3mK": P / 4, "fast_channels": K * chans[u]}
        rows.append(row)
    for u in UNITS:
        for tau in (3000.0, 300.0):
            for tgt in (1.5e-4, 2e-4, 3e-4, 5e-4):
                ok = [r for r in rows if r["unit"] == u and r["tau"] == tau and r["steady_mean_defect"] <= tgt]
                if ok:
                    r = min(ok, key=lambda r: r["K"])
                    best[f"{u}|{tau:g}|{tgt:g}"] = {"K": r["K"], "zone": [r["R"], r["C"]], "defect": r["steady_mean_defect"],
                                                    "P_W_1.2mK": r["P_transport_W_1.2mK"], "fast_channels": r["fast_channels"]}
                else:
                    best[f"{u}|{tau:g}|{tgt:g}"] = None
    for k, v in best.items():
        print(k, v, flush=True)
    save_result("m29_hardware", script=SCRIPT,
                inputs={"P_tweezer_PIC_mW": P_TW_PIC_W * 1e3, "P_tweezer_AOD_mW": P_TW_ATOMS_W / ETA_AOD**2 * 1e3,
                        "fast_channels_per_unit": chans, "capacity": "results/m29_capacity.json (cold steps)",
                        "defects": "results/m29_ladder.json", "units": "one unit = one separate beam path into the "
                        "objective (stack, block, AOD pair, single spot)"},
                outputs={"rows": rows, "smallest_K": best},
                notes="Power is provisioned for K units each at its peak simultaneous tweezers when busy. The ladder "
                      "is coarse, so smallest_K is the smallest simulated K meeting the target, not an interpolation.")
    return rows, best


# ---------------------------------------------------------------------------------------------- figure
def figure():
    rows = json.loads((RESULTS / "m29_hardware.json").read_text())["outputs"]["rows"]
    f, ax = fig(1, 2, h=3.2)
    for j, tau in enumerate((3000.0, 300.0)):
        a = ax[0, j]
        for i, u in enumerate(UNITS):
            r = sorted([x for x in rows if x["unit"] == u and x["tau"] == tau], key=lambda x: x["K"])
            a.loglog([x["K"] for x in r], [x["steady_mean_defect"] for x in r], "o-", ms=3, color=f"C{i}", label=u)
        a.axhline(1e-4 + 0.05 / tau, color="0.5", lw=0.7)
        a.set_xlabel("units K serving 1e6 sites")
        a.set_ylabel("time-averaged defect fraction")
        a.set_title(f"tau = {tau:g} s", fontsize=8)
    ax[0, 0].legend(fontsize=6)
    save_result("m29_figure", script=SCRIPT, fig=f, inputs={"source": "results/m29_hardware.json"}, outputs={},
                notes="Steady-state defect fraction at 1e6 sites vs number of independent units (zones of 1e6/K "
                      "sites), cold steps, assembled start. Grey line: eps + T/(2 tau), the floor for zero latency.")
    return f


if __name__ == "__main__":
    for part in sys.argv[1:] or ["capacity", "ladder", "hardware", "figure"]:
        globals()[part]()
