"""Module 5: thermo-optic phase shifters (P_pi, speed, thermal crosstalk) and actuator summary.

Outputs (results/):
  m5_heaters    P_pi and 10-90% rise time for heater-to-core spacing and isolation options at 894.6 and
                1064.4 nm; crosstalk vs pitch; mesh-convergence check
  m5_actuators  lumped comparison of thermo-optic, AlN and PZT shifters for the OPA and router
"""
from __future__ import annotations

import json

import numpy as np

from picsim import heaters as h
from picsim import params as prm
from picsim.parallel import pmap
from picsim.plotting import fig
from picsim.provenance import RESULTS, P, save_result

SCRIPT = "scripts/m5_heaters.py"
PITCHES = [5, 10, 15, 20, 30, 50, 75, 100, 150]
CONFIGS = {
    "top1um": dict(t_top=1.0), "top2um": dict(t_top=2.0), "top3um": dict(t_top=3.0),
    "top2um_trench": dict(t_top=2.0, trench=6.0, w_trench=4.0),
    "top2um_trench_undercut": dict(t_top=2.0, trench=6.0, w_trench=4.0, undercut=8.0),
}


def dneff_dT():
    rows = json.loads((RESULTS / "m1_design_modes.json").read_text())["outputs"]["rows"]
    return {r["key"]: r for r in rows}


def _job(args):
    name, kw, w_core, t = args
    xs = h.HeaterXS(w_core=w_core, t=t, **kw)
    v, _ = h.steady(xs, x_probe=PITCHES)
    ts, y = h.transient(xs, t_end=5e-3, n=200)
    return {"config": name, "dT_per_W_per_m": float(v[0]), "xtalk": (v[1:] / v[0]).tolist(),
            "rise_s": h.rise_time(ts, y)}


def heaters():
    dm = dneff_dT()
    out = []
    for lk, lam in (("895", prm.LAM_ADDRESS.value), ("1064", prm.LAM_TRANSPORT.value)):
        d = dm[f"{lk}_0.4"]
        res = pmap(_job, [(n, kw, d["w"], d["t"]) for n, kw in CONFIGS.items()])
        k0 = 2 * np.pi / (lam * 1e-6)
        for r in res:
            r.update({"lam": lk, "dneff_dT": d["dneff_dT"],
                      "P_pi_mW": float(np.pi / (k0 * d["dneff_dT"] * r["dT_per_W_per_m"]) * 1e3)})
            out.append(r)
    # mesh convergence on the baseline at 1064 nm
    d = dm["1064_0.4"]
    base = h.HeaterXS(w_core=d["w"], t=d["t"], t_top=2.0)
    v1, _ = h.steady(base)
    fine = h.HeaterXS(w_core=d["w"], t=d["t"], t_top=2.0, d_si=200.0, half_width=300.0)
    v2, _ = h.steady(fine)
    conv = {"baseline_dT": float(v1[0]), "domain_x2_dT": float(v2[0]), "rel_change": float(v2[0] / v1[0] - 1)}
    f, ax = fig(1, 1, w=4.5, h=3.0)
    for r in out:
        if r["lam"] == "1064":
            ax[0, 0].semilogy(PITCHES, r["xtalk"], "-o", ms=3,
                              label=f"{r['config']}: P_pi {r['P_pi_mW']:.0f} mW, {r['rise_s']*1e6:.0f} us")
    ax[0, 0].set_xlabel("heater-to-neighbour pitch (um)")
    ax[0, 0].set_ylabel("neighbour dT / own dT")
    ax[0, 0].legend(fontsize=6)
    save_result("m5_heaters", script=SCRIPT, fig=f,
                inputs={"k_SiO2": prm.K_SIO2, "k_Si": prm.K_SI, "k_SiN": prm.K_SIN, "k_metal": prm.K_METAL,
                        "rhocp": {"SiO2": prm.RHO_CP_SIO2, "Si": prm.RHO_CP_SI}, "t_box": prm.T_BOX,
                        "heater": "CHOICE 2 um x 100 nm strip centred over the guide", "configs": CONFIGS,
                        "dneff_dT": "SIM m1_design_modes (uniform heating)"},
                outputs={"rows": out, "convergence": conv},
                notes="2D cross-section (long heater). Si substrate 100 um with the bottom at ambient; doubling the "
                      "domain changes dT by rel_change. Metal absorption of the optical mode is not computed.")
    return out, conv


def actuators(heat_rows):
    base = [r for r in heat_rows if r["lam"] == "1064" and r["config"] == "top2um"][0]
    tr = [r for r in heat_rows if r["lam"] == "1064" and r["config"] == "top2um_trench_undercut"][0]
    vpi_pzt_1064 = prm.PZT_VPI_532.value * prm.LAM_TRANSPORT.value / 0.532
    rows = {
        "thermo-optic (plain)": {"P_pi_mW": base["P_pi_mW"], "bandwidth_Hz": 0.35 / base["rise_s"],
                                 "64ch_static_power_W": 64 * base["P_pi_mW"] / 1e3, "hysteresis": "none"},
        "thermo-optic (trench + undercut)": {"P_pi_mW": tr["P_pi_mW"], "bandwidth_Hz": 0.35 / tr["rise_s"],
                                             "64ch_static_power_W": 64 * tr["P_pi_mW"] / 1e3, "hysteresis": "none"},
        "PZT stress-optic": {"V_pi_V_at_1064_lambda_scaled": vpi_pzt_1064, "bandwidth_Hz": prm.PZT_BW.value,
                             "holding_power": "tens of nW per actuator (SOURCE arXiv:2601.15695)",
                             "hysteresis": "yes (proposal Sec. 4); magnitude [VALUE?]"},
        "AlN piezo-optomechanical": {"V_pi_V_795": prm.ALN_VPI_795.value, "V_pi_V_1013": prm.ALN_VPI_1013.value,
                                     "rise_s": prm.ALN_RISE.value, "hysteresis": "negligible (AlN not ferroelectric)"},
    }
    save_result("m5_actuators", script=SCRIPT,
                inputs={"pzt": [prm.PZT_VPI_532, prm.PZT_BW], "aln": [prm.ALN_VPI_795, prm.ALN_VPI_1013, prm.ALN_RISE],
                        "PZT_Vpi_1064": P(vpi_pzt_1064, "V", "DERIVED",
                                          "2.8 V x 1064.4/532 assuming the same coil length and stress-optic "
                                          "coefficient; real value depends on mode overlap")},
                outputs=rows,
                notes="Thermo-optic bandwidth ~0.35/rise. 64-channel power assumes on average P_pi per channel "
                      "(phases uniformly spread over 0..2 pi).")
    return rows


if __name__ == "__main__":
    rows, conv = heaters()
    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k != "xtalk"},
              [round(x, 4) for x in r["xtalk"]], flush=True)
    print(conv)
    print(json.dumps(actuators(rows), indent=1))
