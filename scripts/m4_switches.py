"""Module 4: MZI switches and the 894.6 nm routing tree.

Outputs (results/):
  m4_extinction    Monte Carlo extinction for four switch architectures vs phase-setting accuracy, using the
                   2x2 MMI split-ratio spread from Module 2; probability of >= 50 dB
  m4_polarisation  on-chip TM suppression needed vs input polarisation extinction
  m4_dynamics      optical rise and time-to-50 dB for an AlN actuator modelled as a second-order resonance,
                   with and without input shaping (ZV, ZVD) and with resonance-frequency error
  m4_router        16-output tunable-splitter tree: leakage to unaddressed outputs, combined idle extinction,
                   insertion loss and input power
"""
from __future__ import annotations

import json
import sys

import numpy as np

from picsim import circuits as c
from picsim import params as prm
from picsim.plotting import fig
from picsim.provenance import RESULTS, P, save_result

SCRIPT = "scripts/m4_switches.py"
TARGET = prm.ER_TARGET.value


def kappa_sigma():
    f = RESULTS / "m2_mmi.json"
    if f.exists():
        mc = json.loads(f.read_text())["outputs"]["mc"]["895"]
        return float(mc["kappa_2x2_std"]), float(abs(mc["kappa_2x2_mean"] - 0.5)), "SIM m2_mmi (895 nm, 4 um 2x2)"
    return 0.01, 0.0, "CHOICE fallback 0.01"


def extinction():
    sk, bias, src = kappa_sigma()
    sk_eff = float(np.hypot(sk, bias))  # treat the systematic offset as part of the spread (conservative)
    archs = {"single": dict(stages=1, tunable=False), "cascade2": dict(stages=2, tunable=False),
             "single_trimmed": dict(stages=1, tunable=True), "cascade2_trimmed": dict(stages=2, tunable=True)}
    sig_phis = [0.001, 0.003, 0.01, 0.03]
    floors = {"none": 0.0, "-70dB": 1e-7}
    rows = []
    for name, kw in archs.items():
        for sp in sig_phis:
            for fl_name, fl in floors.items():
                er = c.mc_extinction(4000, sk_eff, sp, sigma_a=0.005, floor=fl, seed=1, **kw)
                rows.append({"arch": name, "sigma_phi_rad": sp, "floor": fl_name, "median_dB": float(np.median(er)),
                             "p5_dB": float(np.percentile(er, 5)), "P_ge_target": float(np.mean(er >= TARGET))})
    f, ax = fig(1, 1, w=5.0, h=3.2)
    for name, col in zip(archs, ("C0", "C1", "C2", "C3")):
        r = [x for x in rows if x["arch"] == name and x["floor"] == "none"]
        ax[0, 0].semilogx([x["sigma_phi_rad"] for x in r], [x["p5_dB"] for x in r], "-o", color=col, ms=3,
                          label=f"{name} (5th percentile)")
        ax[0, 0].semilogx([x["sigma_phi_rad"] for x in r], [x["median_dB"] for x in r], ":", color=col)
    ax[0, 0].axhline(TARGET, color="k", lw=0.8)
    ax[0, 0].set_xlabel("phase-setting error rms (rad)")
    ax[0, 0].set_ylabel("extinction (dB)")
    ax[0, 0].legend(fontsize=6)
    # analytic: balanced MZI, phase error only: T_off = sin^2(dphi/2) -> ER = -20 log10(dphi/2)
    phi_req = 2 * 10 ** (-TARGET / 20)
    save_result("m4_extinction", script=SCRIPT, fig=f,
                inputs={"sigma_kappa": P(sk_eff, "-", "SIM", src + " (std and mean offset combined)"),
                        "sigma_arm_amplitude": "CHOICE 0.005", "ER_target": prm.ER_TARGET,
                        "floors": floors, "trials": 4000},
                outputs={"rows": rows, "single_balanced_phase_error_for_target_rad": phi_req},
                notes="Trimmed: the input splitter is an MZI trimmed to match the output coupler, limited by the "
                      "same phase error. Floor: incoherent stray light relative to input.")
    return rows, sk_eff


def polarisation():
    rows = []
    for per in (15, 20, 25, 30):
        for T_tm in (1.0, 0.25):
            need = TARGET - per + 10 * np.log10(T_tm)
            rows.append({"input_PER_dB": per, "TM_transmission_at_TE_null": T_tm,
                         "onchip_TM_suppression_needed_dB": float(max(need, 0))})
    save_result("m4_polarisation", script=SCRIPT, inputs={"ER_target": prm.ER_TARGET},
                outputs={"rows": rows},
                notes="TM light is not nulled when the TE path is at its null. T_TM = 1 is the worst case for one "
                      "MZI; 0.25 is the mean of two cascaded MZIs with random TM phase. Required suppression "
                      "S = ER_target - PER + 10 log10(T_TM).")
    return rows


def dynamics():
    rows = []
    for f0 in (10e6, 20e6, 40e6):
        for Q in (3, 10, 30):
            for shape in (False, True, "zvd"):
                for err in (0.0, 0.02, 0.05):
                    m = c.switching_metrics(f0, Q, shaped=shape, f0_true=f0 * (1 + err), t_max=max(20 * Q / f0, 2e-6),
                                            ER_target_db=TARGET)
                    rows.append({"f0_MHz": f0 / 1e6, "Q": Q, "shaper": {False: "none", True: "ZV", "zvd": "ZVD"}[shape],
                                 "f0_error": err, "rise_ns": m["rise_10_90_s"] * 1e9,
                                 "t_to_50dB_ns": m["t_to_ER_s"] * 1e9})
    # half-wave voltage at 894.6 nm: bracket from the two published points (different device designs)
    v795, v1013 = prm.ALN_VPI_795.value, prm.ALN_VPI_1013.value
    vpi_scaled = v795 * prm.LAM_ADDRESS.value / 0.795
    vpi_interp = float(np.interp(prm.LAM_ADDRESS.value, [0.795, 1.013], [v795, v1013]))
    drive = {"Vpi_895_lambda_scaled_V": vpi_scaled, "Vpi_895_linear_interp_V": vpi_interp,
             "V_noise_rms_for_target_V": {f"Vpi={v:.0f}": float(2 * 10 ** (-TARGET / 20) * v / np.pi)
                                          for v in (vpi_scaled, vpi_interp)}}
    f, ax = fig(1, 1, w=5.0, h=3.2)
    for shaper, col in (("none", "C0"), ("ZV", "C1"), ("ZVD", "C2")):
        for err, ls in ((0.0, "-o"), (0.02, "--s")):
            r = [x for x in rows if x["shaper"] == shaper and x["f0_error"] == err and x["f0_MHz"] == 20]
            ax[0, 0].semilogy([x["Q"] for x in r], [x["t_to_50dB_ns"] for x in r], ls, color=col, ms=3,
                              label=f"{shaper}, f0 error {err:.0%}")
    ax[0, 0].axhline(prm.T_SWITCH.value * 1e9, color="k", lw=0.8)
    ax[0, 0].set_xlabel("actuator Q (f0 = 20 MHz)")
    ax[0, 0].set_ylabel("time to 50 dB after switch-off (ns)")
    ax[0, 0].legend(fontsize=6)
    save_result("m4_dynamics", script=SCRIPT, fig=f,
                inputs={"rise_ref": prm.ALN_RISE, "Vpi_795": prm.ALN_VPI_795, "Vpi_1013": prm.ALN_VPI_1013,
                        "model": "single mechanical mode, f0 and Q are CHOICE sweep values (not measured)"},
                outputs={"rows": rows, "drive": drive},
                notes="Vpi at 894.6 nm is bracketed, not known: lambda scaling of the 795 nm device and linear "
                      "interpolation to the 1013 nm device (a different design) differ.")
    return rows, drive


def router(er_rows):
    stages = int(np.log2(prm.N_ROUTER.value))
    rows = []
    for sp in (0.003, 0.01, 0.03, 0.1):
        P_ = c.router_leakage(stages, sp, n=2000, seed=3)
        leak = P_[:, 1:]  # target output 0
        rows.append({"sigma_phi_rad": sp, "router_leak_max_mean": float(leak.max(axis=1).mean()),
                     "router_leak_max_p95": float(np.percentile(leak.max(axis=1), 95)),
                     "router_leak_sum_mean": float(leak.sum(axis=1).mean())})
    # combine with output switch OFF state (cascade2, sigma_phi 0.003, 5th percentile)
    sw = [r for r in er_rows if r["arch"] == "cascade2" and r["sigma_phi_rad"] == 0.003 and r["floor"] == "none"][0]
    for r in rows:
        r["idle_site_total_dB"] = float(-10 * np.log10(r["router_leak_max_p95"]) + sw["p5_dB"])
    # insertion loss and input power
    m2 = json.loads((RESULTS / "m2_mmi.json").read_text())["outputs"] if (RESULTS / "m2_mmi.json").exists() else None
    m2c = json.loads((RESULTS / "m2_couplers.json").read_text())["outputs"] if (RESULTS / "m2_couplers.json").exists() else None
    ex22 = [x for x in m2["rows"] if x["lam"] == "895" and x["W"] == 4.0][0]["excess_2x2_dB"] if m2 else 0.1
    facet = max(v["dB"] for k, v in m2c["best"].items() if k.startswith("895_0.4")) if m2c else -1.5
    il = {"facet (SIM m2)": facet, "tree: 4 MZI x 2 MMI (SIM m2)": -2 * stages * ex22,
          "output switch: 2 MZI x 2 MMI (SIM m2)": -4 * ex22, "propagation 1 cm @0.3 dB/cm (CHOICE)": -0.3,
          "output facet + free-space delivery + meta-objective (CHOICE 0.5)": 10 * np.log10(0.5)}
    il_tot = sum(il.values())
    p_site = max(prm.P_SITE_ADDR.value)
    p_in = prm.N_ROUTER.value * p_site / 10 ** (il_tot / 10)
    out = {"rows": rows, "switch_used": sw, "insertion_items_dB": il, "insertion_total_dB": il_tot,
           "P_in_all_16_sites_at_20uW_W": p_in}
    save_result("m4_router", script=SCRIPT,
                inputs={"N_outputs": prm.N_ROUTER, "P_site": prm.P_SITE_ADDR},
                outputs=out,
                notes="Router leakage is coherent-free (powers add); the output switch extinction multiplies.")
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["extinction", "polarisation", "dynamics", "router"]
    er_rows = None
    if "extinction" in parts:
        er_rows, sk = extinction()
        print("sigma_kappa", sk)
        for r in er_rows:
            if r["floor"] == "none":
                print(r, flush=True)
    if "polarisation" in parts:
        for r in polarisation():
            print(r)
    if "dynamics" in parts:
        rows, drive = dynamics()
        for r in rows:
            if r["f0_MHz"] == 20:
                print(r)
        print(drive)
    if "router" in parts:
        if er_rows is None:
            er_rows = json.loads((RESULTS / "m4_extinction.json").read_text())["outputs"]["rows"]
        print(json.dumps(router(er_rows), indent=1))
