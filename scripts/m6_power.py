"""Module 6: on-chip power handling at 1064.4 nm (OPA input up to 1.4 W).

Output results/m6_power.json:
  self-heating of the input guide vs propagation loss (absorbed fraction 30% and 100%), its thermo-optic
  phase per cm, the per-channel figures after a 64-way split, facet and in-guide intensities, Kerr phase,
  and the number of input fibres needed to keep the facet intensity below a chosen limit.
"""
from __future__ import annotations

import json

import numpy as np

from picsim import heaters as h
from picsim import params as prm
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m6_power.py"

if __name__ == "__main__":
    dm = {r["key"]: r for r in json.loads((RESULTS / "m1_design_modes.json").read_text())["outputs"]["rows"]}
    d = dm["1064_0.4"]
    lam = prm.LAM_TRANSPORT.value
    xs = h.HeaterXS(w_core=d["w"], t=d["t"], t_top=2.0, source="core")
    dT_per, _ = h.steady(xs)
    dT_per = float(dT_per[0])  # K per (W/m) absorbed in the core
    k0 = 2 * np.pi / (lam * 1e-6)
    rows = []
    for P in (prm.P_CHIP.value, prm.P_CHIP.value / 64):
        for loss in (0.05, 0.1, 0.3, 0.5):
            for f_abs in (0.3, 1.0):
                Pabs = P * (1 - 10 ** (-loss / 10)) * f_abs * 100  # W/m (loss per cm -> per m)
                dT = Pabs * dT_per
                rows.append({"P_W": P, "loss_dB_cm": loss, "absorbed_fraction": f_abs, "P_abs_W_per_m": Pabs,
                             "dT_core_K": dT, "thermo_phase_rad_per_cm": k0 * d["dneff_dT"] * dT * 1e-2})
    m2 = json.loads((RESULTS / "m2_couplers.json").read_text())["outputs"] if (RESULTS / "m2_couplers.json").exists() else {}
    out = {"dT_core_per_W_per_m": dT_per, "rows": rows,
           "I_guide_MW_cm2_at_1.4W": d["I_peak_MW_cm2_at_1.4W"],
           "kerr_rad_per_cm_at_1.4W": d["kerr_phase_rad_per_cm_at_1.4W"],
           "facet_MW_cm2_at_1.4W_t0.4": m2.get("facet_intensity_MW_cm2_at_1.4W_t0.4")}
    # fibres needed if the facet intensity is held below a limit (limit is a CHOICE; damage data [VALUE?])
    for lim in (5.0, 20.0):
        f = out["facet_MW_cm2_at_1.4W_t0.4"]
        out[f"n_fibres_for_facet_below_{lim:.0f}MW_cm2"] = int(np.ceil(f / lim)) if f else None
    save_result("m6_power", script=SCRIPT,
                inputs={"P_chip": prm.P_CHIP, "dneff_dT": "SIM m1_design_modes 1064_0.4",
                        "thermal": "picsim.heaters with the core as heat source (2D, long guide)",
                        "loss_sweep": "CHOICE 0.05-0.5 dB/cm; absorbed fraction CHOICE 0.3 or 1.0",
                        "facet_limits": "CHOICE 5 and 20 MW/cm^2; CW facet damage threshold of SiN/SiO2 at 1 um "
                                        "[VALUE?] not found in session"},
                outputs=out,
                notes="Temperature rise is the steady state of a uniformly absorbing straight guide; the phase is "
                      "common-mode before the split and differential only through channel-to-channel loss spread.")
    print(json.dumps(out, indent=1))
