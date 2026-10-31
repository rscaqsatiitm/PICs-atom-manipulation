"""Module 7: error-budget links and the parameter ledger.

Outputs:
  results/m7_budget.json   extinction -> idle error; PSF crosstalk -> neighbour error; trap depth cross-check;
                           position- and intensity-noise limits at the trap frequency
  docs/parameters.md       every provenance-tagged input found in results/*.json, plus the key SIM outputs
"""
from __future__ import annotations

import json

import numpy as np

from picsim import budget as b
from picsim import params as prm
from picsim import transport as tr
from picsim.plotting import fig
from picsim.provenance import REPO, RESULTS, P, save_result

SCRIPT = "scripts/m7_budget.py"


def links():
    out = {}
    # 1. switch extinction -> idle-qubit error (proposal footnote 4)
    out["idle_error_vs_ER"] = [{"ER_dB": e, "one_beam_addressed": float(b.idle_error_from_extinction(e)),
                                "both_beams_addressed": float(b.idle_error_from_extinction(e, True))}
                               for e in (30, 40, 50, 60, 70)]
    # 2. optical crosstalk of the 894.6 nm addressing beam at the nearest site (envelope max over +/-0.5 um)
    lam, NA, a = prm.LAM_ADDRESS.value, prm.NA_OBJ.value, prm.SITE_PITCH.value
    rr = np.linspace(a - 0.5, a + 0.5, 101)
    rows = []
    for trunc in (None, 1.5, 1.0, 0.85, 0.7, 0.5):
        c = float(b.psf_intensity(rr, lam, NA, trunc).max())
        r_f = np.linspace(0, 1.5, 1501)
        I = b.psf_intensity(r_f, lam, NA, trunc)
        fwhm = 2 * r_f[np.argmax(I < 0.5)]
        rows.append({"pupil_gauss_trunc": "uniform" if trunc is None else trunc, "crosstalk_at_neighbour": c,
                     "neighbour_error_one_beam": float(b.neighbour_gate_error(c)), "fwhm_um": float(fwhm)})
    out["psf_crosstalk_895"] = rows
    # 3. trap depth cross-check (proposal: 6.6 mW -> ~1.2 mK at w0 = 1.04 um)
    w0, al = prm.W0_TWEEZER.value, prm.ALPHA_CS_1064.value
    out["depth_check"] = {"depth_mK_at_6.6mW": float(b.depth_mK(prm.P_TWEEZER_ATOMS.value, w0, al)),
                          "power_mW_for_1.2mK": float(tr.power_for_depth(prm.DEPTH_TRANSPORT.value, w0, al) * 1e3),
                          "waist_um_for_1.2mK_at_6.6mW": float(w0 * np.sqrt(b.depth_mK(prm.P_TWEEZER_ATOMS.value, w0, al)
                                                                            / (prm.DEPTH_TRANSPORT.value * 1e3)))}
    m3s = RESULTS / "m3_spot.json"
    if m3s.exists():
        s = json.loads(m3s.read_text())["outputs"]["uniform"]
        # 1D atom-plane profile: peak per unit transmitted power along x (1/um); y Gaussian with w0
        peak_x = s["peak"] / s["transmission"]
        peak_y = np.sqrt(2 / np.pi) / w0
        I_per_W = peak_x * peak_y  # 1/um^2
        I_gauss = 2 / (np.pi * w0**2)
        out["depth_check"]["opa_spot_peak_vs_gaussian_w0"] = float(I_per_W / I_gauss)
    # 4. position noise and intensity noise limits
    nl = []
    for U in (0.3e-3, prm.DEPTH_TRANSPORT.value):
        f_tr = tr.radial_freq_Hz(U, w0)
        for rate in (1.0, 10.0, 100.0):
            Sx = b.position_noise_psd_limit(rate, f_tr)
            nl.append({"U_mK": U * 1e3, "f_trap_kHz": f_tr / 1e3, "heating_quanta_per_s": rate,
                       "S_x_limit_m2_per_Hz": float(Sx), "x_rms_in_10kHz_band_nm": float(np.sqrt(Sx * 1e4) * 1e9),
                       "RIN_limit_1_per_Hz_for_gamma_1_per_s": float(b.rin_limit(1.0, f_tr))})
        # the proposal's 10 nm rms, if white over 1 MHz
        Sx_w = (prm.POS_NOISE.value * 1e-6) ** 2 / 1e6
        nl.append({"U_mK": U * 1e3, "f_trap_kHz": f_tr / 1e3, "case": "10 nm rms white over 1 MHz",
                   "S_x_m2_per_Hz": Sx_w, "heating_quanta_per_s": float(tr.heating_rate_position_noise(Sx_w, f_tr))})
    out["noise_limits"] = nl
    f, ax = fig(1, 1, w=4.5, h=3.0)
    for r in rows:
        ax[0, 0].plot(r["fwhm_um"], r["crosstalk_at_neighbour"], "o", label=str(r["pupil_gauss_trunc"]))
    ax[0, 0].set_yscale("log")
    ax[0, 0].axhline(1e-4, color="k", lw=0.8)
    ax[0, 0].set_xlabel("addressing spot FWHM (um)")
    ax[0, 0].set_ylabel(f"intensity at {a:.0f} um neighbour")
    ax[0, 0].legend(fontsize=6, title="Gaussian pupil trunc.")
    save_result("m7_budget", script=SCRIPT, fig=f,
                inputs={"NA": prm.NA_OBJ, "lam_addr": prm.LAM_ADDRESS, "site_pitch": prm.SITE_PITCH, "w0": prm.W0_TWEEZER,
                        "alpha": prm.ALPHA_CS_1064, "P_tweezer": prm.P_TWEEZER_ATOMS, "depth": prm.DEPTH_TRANSPORT,
                        "pos_noise": prm.POS_NOISE,
                        "formulas": P(0, "-", "RECALLED", "heating rates: Savard, O'Hara, Thomas, PRA 56, R1095 (1997)")},
                outputs=out,
                notes="PSF: scalar Debye, aberration-free; a Strehl of 0.8 adds a halo not included here.")
    return out


def ledger():
    """Collect provenance-tagged inputs and summary outputs from all results into docs/parameters.md."""
    seen = {}

    def walk(x, where):
        if isinstance(x, dict):
            if set(x.keys()) >= {"value", "unit", "kind", "ref"}:
                key = (json.dumps(x["value"]), x["unit"], x["ref"])
                seen.setdefault(key, {**x, "used_in": set()})["used_in"].add(where)
                return
            for v in x.values():
                walk(v, where)
        elif isinstance(x, list):
            for v in x:
                walk(v, where)

    files = sorted(RESULTS.glob("m*.json"))
    for fp in files:
        try:
            d = json.loads(fp.read_text())
        except Exception:
            continue
        if isinstance(d, dict):
            walk(d.get("inputs", {}), fp.stem)
    lines = ["# Parameter ledger", "",
             "Generated by `scripts/m7_budget.py` from the input records of every `results/*.json`. "
             "Types: USER (proposal text), SOURCE (read in session, pointer given), DERIVED, SIM, CHOICE (design or "
             "sweep choice), RECALLED (not verified; check before use). Simulation outputs are SIM and must be "
             "reviewed before submission.", "",
             "| value | unit | type | pointer | used in |", "|---|---|---|---|---|"]
    for v in sorted(seen.values(), key=lambda r: (r["kind"], r["ref"])):
        val = json.loads(json.dumps(v["value"]))
        if val == 0 and v["unit"] == "-":
            val = "(model)"
        lines.append(f"| {val} | {v['unit']} | {v['kind']} | {v['ref'].replace('|', '/')} | "
                     f"{', '.join(sorted(v['used_in']))} |")
    lines += ["", "## Result files", ""]
    for fp in files:
        try:
            d = json.loads(fp.read_text())
            if isinstance(d, dict):
                lines.append(f"- `results/{fp.name}`: {d.get('notes', '')[:200]}")
        except Exception:
            pass
    (REPO / "docs" / "parameters.md").write_text("\n".join(lines) + "\n")
    return len(seen)


if __name__ == "__main__":
    o = links()
    print(json.dumps(o, indent=1, default=float))
    print("ledger entries:", ledger())
