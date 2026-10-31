"""Module 1: SiN platform at 894.6 nm and 1064.4 nm.

Outputs (results/):
  m1_width_sweep     TE0/TE1/TM0 effective indices vs width for t = 200, 300, 400 nm; single-mode width
  m1_design_modes    design-width properties: neff, ng, confinement, Aeff, dneff/dT, peak intensity, Kerr phase
  m1_bend_loss       TE0 and TM0 bend loss per 90 deg vs radius (conformal map, graded absorber)
  m1_box_leakage     substrate leakage vs buried-oxide thickness
  m1_coupling        coupling length vs gap of two identical guides (sets OPA pitch and arm spacing)
"""
from __future__ import annotations

import json
import sys

import numpy as np

from picsim import materials as mat
from picsim import modes as M
from picsim import params as prm
from picsim.parallel import pmap
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

LAMS = {"895": prm.LAM_ADDRESS.value, "1064": prm.LAM_TRANSPORT.value}
TS = prm.T_SIN_SWEEP.value
WIDTHS = np.round(np.arange(0.3, 1.81, 0.1), 3)
SCRIPT = "scripts/m1_platform.py"


# ------------------------------------------------------------------ width sweep
def _width_job(args):
    t, lk, w = args
    lam = LAMS[lk]
    ms = M.solve(M.XSection(cores=[(0, w)], t=t, x_left=-4 - w / 2, x_right=4 + w / 2), lam, num_modes=6)
    nc = mat.n_sio2(lam)
    te = [m for m in ms if m.te_fraction > 0.5 and m.neff.real > nc + 1e-3]
    tm = [m for m in ms if m.te_fraction <= 0.5 and m.neff.real > nc + 1e-3]
    g = lambda lst, i: lst[i].neff.real if len(lst) > i else np.nan
    return {"t": t, "lam": lk, "w": w, "TE0": g(te, 0), "TE1": g(te, 1), "TM0": g(tm, 0), "TM1": g(tm, 1)}


TE1_MARGIN = 0.015  # TE1 counted as stripped if neff(TE1) < n_clad + margin (radiates in bends)
BIREF_MIN = 0.015  # minimum TE0-TM0 index splitting to avoid polarisation hybridisation


def select_design(rows):
    """Design width per (lambda, t): widest grid width with TE1 cut off or weakly guided and TE0-TM0 split
    >= BIREF_MIN. Rationale: widest width minimises sidewall-scattering loss and width sensitivity."""
    design = {}
    for lk in LAMS:
        nc = float(mat.n_sio2(LAMS[lk]))
        for t in TS:
            r = sorted([x for x in rows if x["t"] == t and x["lam"] == lk], key=lambda x: x["w"])
            ok = [x for x in r if (not np.isfinite(x["TE1"]) or x["TE1"] < nc + TE1_MARGIN)
                  and abs(x["TE0"] - x["TM0"]) >= BIREF_MIN]
            te1 = [x["w"] for x in r if np.isfinite(x["TE1"])]
            w_d = float(max(x["w"] for x in ok)) if ok else float("nan")
            sel = [x for x in r if x["w"] == w_d][0] if ok else None
            design[f"{lk}_{t}"] = {"lam": LAMS[lk], "t": t, "w_TE1_guided_at": min(te1) if te1 else None,
                                   "w_design": w_d,
                                   "TE1_margin_at_design": (sel["TE1"] - nc) if sel and np.isfinite(sel["TE1"]) else None,
                                   "biref_at_design": (sel["TE0"] - sel["TM0"]) if sel else None}
    return design


def plot_width(rows, design):
    f, ax = fig(1, 2, h=3.0)
    for j, lk in enumerate(LAMS):
        nc = mat.n_sio2(LAMS[lk])
        for t, c in zip(TS, ["C0", "C1", "C2"]):
            r = sorted([x for x in rows if x["t"] == t and x["lam"] == lk], key=lambda x: x["w"])
            w = np.array([x["w"] for x in r])
            for key, ls in (("TE0", "-"), ("TE1", "--"), ("TM0", ":")):
                ax[0, j].plot(w, [x[key] for x in r], ls, color=c, label=f"{key}, t={t*1e3:.0f} nm" if j == 0 else None)
            ax[0, j].axvline(design[f"{lk}_{t}"]["w_design"], color=c, lw=0.6, alpha=0.6)
        ax[0, j].axhline(nc, color="k", lw=0.8)
        ax[0, j].set_title(f"lambda = {LAMS[lk]*1e3:.1f} nm (vertical lines: design widths)")
        ax[0, j].set_xlabel("width (um)")
    ax[0, 0].set_ylabel("n_eff")
    ax[0, 0].legend(fontsize=6, ncol=1)
    return f


def width_sweep(rows=None):
    if rows is None:
        jobs = [(t, lk, w) for t in TS for lk in LAMS for w in WIDTHS]
        rows = pmap(_width_job, jobs)
    design = select_design(rows)
    f = plot_width(rows, design)
    save_result("m1_width_sweep", script=SCRIPT, fig=f,
                inputs={"thicknesses": prm.T_SIN_SWEEP, "lam_895": prm.LAM_ADDRESS, "lam_1064": prm.LAM_TRANSPORT,
                        "n_SiN": prm.N_SIN_MODEL, "n_SiO2": prm.N_SIO2_MODEL,
                        "widths_um": WIDTHS, "solver": "femwell 0.1.12 vector FEM, order 2, 0.03 um core mesh",
                        "TE1_margin": TE1_MARGIN, "birefringence_min": BIREF_MIN},
                outputs={"rows": rows, "design": design},
                notes="Design width: widest 0.1 um grid width with TE1 cut off or within TE1_margin of the "
                      "cladding index, and TE0-TM0 splitting >= birefringence_min (a square 400 x 400 nm guide "
                      "at 895 nm is polarisation-degenerate).")
    return design


# ------------------------------------------------------------------ design-width properties
def _design_job(args):
    key, d = args
    lam, t, w = d["lam"], d["t"], d["w_design"]
    xs = M.XSection(cores=[(0, w)], t=t)
    te = M.fundamental(M.solve(xs, lam, num_modes=2))
    out = {"key": key, "lam": lam, "t": t, "w": w, "neff": te.neff.real, "gamma_core": te.gamma_core,
           "aeff_um2": te.aeff}
    # group index from +/- 2 nm
    nep = M.fundamental(M.solve(xs, lam + 0.002, num_modes=2)).neff.real
    nem = M.fundamental(M.solve(xs, lam - 0.002, num_modes=2)).neff.real
    out["ng"] = te.neff.real - lam * (nep - nem) / 0.004
    # thermo-optic: uniform +10 K in core and cladding
    xs_h = M.XSection(cores=[(0, w)], t=t, dT_core=10.0, dT_clad=10.0)
    out["dneff_dT"] = (M.fundamental(M.solve(xs_h, lam, num_modes=2)).neff.real - te.neff.real) / 10.0
    # intensity and Kerr phase
    for P in (prm.P_CHIP.value, 0.066):
        out[f"I_peak_MW_cm2_at_{P}W"] = P / (te.aeff * 1e-12) / 1e10
        out[f"kerr_phase_rad_per_cm_at_{P}W"] = (2 * np.pi / (lam * 1e-6) * prm.N2_SIN.value * P
                                               / (te.aeff * 1e-12) * te.gamma_core * 1e-2)
    return out


def design_modes(design):
    rows = pmap(_design_job, list(design.items()))
    save_result("m1_design_modes", script=SCRIPT,
                inputs={"design_widths": design, "n2_SiN": prm.N2_SIN, "dndT_SiN": prm.DNDT_SIN,
                        "dndT_SiO2": prm.DNDT_SIO2, "P_chip": prm.P_CHIP},
                outputs={"rows": rows},
                notes="Kerr phase uses gamma_core x n2 x P/Aeff (upper bound: Aeff from transverse |E|^2). "
                      "dneff/dT for uniform heating of core and cladding.")
    return {r["key"]: r for r in rows}


# ------------------------------------------------------------------ bend loss
R_SETS = {0.2: [20, 30, 40, 60, 80, 120], 0.3: [8, 12, 16, 20, 30, 40], 0.4: [4, 6, 8, 10, 12, 16]}


def _bend_job(args):
    key, lam, t, w, R, pol, n_guess = args
    xr = R * (n_guess / mat.n_sio2(lam) - 1)
    absw = 4.0
    xs = M.XSection(cores=[(0, w)], t=t, x_left=-3 - w / 2, x_right=max(xr, 2) + 3 + absw, absorber=absw,
                    res_max=0.3)
    ms = M.solve(xs, lam, num_modes=4, radius=R, n_guess=n_guess, absorber_k=0.3)
    cand = [m for m in ms if (m.te_fraction > 0.5) == (pol == "TE")]
    if not cand:
        return {"key": key, "R": R, "pol": pol, "db_per_90": np.nan, "gamma": np.nan}
    m = max(cand, key=lambda m: m.gamma_core)
    return {"key": key, "R": R, "pol": pol, "db_per_90": max(m.loss_db_per_cm, 0) * 1e-4 * np.pi / 2 * R,
            "gamma": m.gamma_core, "neff": m.neff.real}


def bend_loss(design, straight):
    jobs = []
    for key, d in design.items():
        lam, t, w = d["lam"], d["t"], d["w_design"]
        st = M.solve(M.XSection(cores=[(0, w)], t=t), lam, num_modes=4)
        for pol in ("TE", "TM"):
            m0 = M.fundamental(st, pol)
            if m0 is None or m0.neff.real < mat.n_sio2(lam) + 1e-3:
                continue
            for R in R_SETS[t]:
                jobs.append((key, lam, t, w, R, pol, m0.neff.real))
    rows = pmap(_bend_job, jobs)
    f, ax = fig(1, 2, h=3.0)
    for j, lk in enumerate(LAMS):
        for t, c in zip(TS, ["C0", "C1", "C2"]):
            key = f"{lk}_{t}"
            for pol, ls in (("TE", "-o"), ("TM", "--s")):
                r = [x for x in rows if x["key"] == key and x["pol"] == pol]
                if r:
                    ax[0, j].semilogy([x["R"] for x in r], np.maximum([x["db_per_90"] for x in r], 1e-8), ls,
                                      color=c, ms=3, label=f"{pol}0, t={t*1e3:.0f} nm, w={design[key]['w_design']} um")
        ax[0, j].axhline(0.01, color="k", lw=0.8)
        ax[0, j].set_xscale("log")
        ax[0, j].set_xlabel("bend radius (um)")
        ax[0, j].set_title(f"lambda = {LAMS[lk]*1e3:.1f} nm")
        ax[0, j].legend(fontsize=5)
    ax[0, 0].set_ylabel("radiation loss (dB / 90 deg)")
    # minimum radius for < 0.01 dB / 90 deg, by interpolation in log space
    rmin = {}
    for key in design:
        for pol in ("TE", "TM"):
            r = sorted([x for x in rows if x["key"] == key and x["pol"] == pol], key=lambda x: x["R"])
            if not r:
                continue
            Rs = np.array([x["R"] for x in r], float)
            L = np.maximum(np.array([x["db_per_90"] for x in r]), 1e-12)
            ok = np.where(L < 0.01)[0]
            rmin[f"{key}_{pol}"] = float(Rs[ok[0]]) if len(ok) else f">{Rs[-1]}"
    save_result("m1_bend_loss", script=SCRIPT, fig=f,
                inputs={"design": design, "method": "femwell conformal bend, quadratic graded absorber 4 um, k_max 0.3; "
                        "mode chosen by maximum core power fraction; checked insensitive to absorber 4 vs 8 um "
                        "and k 0.3 vs 1.0 (R = 15-60 um, t = 400 nm)"},
                outputs={"rows": rows, "R_min_for_0.01dB_per_90": rmin},
                notes="Values below ~1e-6 dB/90 are at the numerical floor. Mode-mismatch loss at straight-bend "
                      "junctions is not included.")
    return rmin


# ------------------------------------------------------------------ substrate leakage
def _box_job(args):
    key, lam, t, w, tbox, pol = args
    xs = M.XSection(cores=[(0, w)], t=t, substrate=True, t_box=tbox, x_left=-5, x_right=5)
    ms = M.solve(xs, lam, num_modes=4, n_guess=None)
    cand = [m for m in ms if (m.te_fraction > 0.5) == (pol == "TE")]
    m = max(cand, key=lambda m: m.gamma_core) if cand else None
    return {"key": key, "t_box": tbox, "pol": pol, "db_per_cm": max(m.loss_db_per_cm, 0) if m else np.nan}


def box_leakage(design):
    jobs = [(k, d["lam"], d["t"], d["w_design"], tb, pol) for k, d in design.items()
            for tb in (1.0, 1.5, 2.0, 3.0, 4.0) for pol in ("TE", "TM")]
    rows = pmap(_box_job, jobs)
    f, ax = fig(1, 1, w=4.0, h=3.0)
    for k in design:
        for pol, ls in (("TE", "-o"), ("TM", "--s")):
            r = [x for x in rows if x["key"] == k and x["pol"] == pol]
            ax[0, 0].semilogy([x["t_box"] for x in r], np.maximum([x["db_per_cm"] for x in r], 1e-8), ls, ms=3,
                              label=f"{pol}0 {k}")
    ax[0, 0].axhline(0.01, color="k", lw=0.8)
    ax[0, 0].set_xlabel("buried oxide thickness (um)")
    ax[0, 0].set_ylabel("substrate leakage (dB/cm)")
    ax[0, 0].legend(fontsize=5, ncol=2)
    save_result("m1_box_leakage", script=SCRIPT, fig=f,
                inputs={"design": design, "n_Si": prm.N_SI_1064,
                        "method": "Si substrate with 1 um graded absorbing layer below 1.5 um of Si"},
                outputs={"rows": rows}, notes="Values below ~1e-5 dB/cm are at the numerical floor.")


# ------------------------------------------------------------------ directional coupling vs gap
GAPS = [0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]


def _cpl_job(args):
    key, lam, t, w, gap = args
    xc = (gap + w) / 2
    xs = M.XSection(cores=[(-xc, w), (xc, w)], t=t, x_left=-xc - w / 2 - 4, x_right=xc + w / 2 + 4)
    ms = M.solve(xs, lam, num_modes=6)
    te = [m for m in ms if m.te_fraction > 0.5][:2]
    dn = te[0].neff.real - te[1].neff.real
    return {"key": key, "gap": gap, "pitch": gap + w, "dn": dn, "Lc_um": lam / (2 * dn)}


def coupling(design):
    keys = [k for k in design if design[k]["t"] in (0.2, 0.4)]
    rows = pmap(_cpl_job, [(k, design[k]["lam"], design[k]["t"], design[k]["w_design"], g) for k in keys for g in GAPS])
    f, ax = fig(1, 1, w=4.2, h=3.0)
    for k in keys:
        r = [x for x in rows if x["key"] == k]
        ax[0, 0].semilogy([x["pitch"] for x in r], [x["Lc_um"] for x in r], "-o", ms=3,
                          label=f"{k} (w={design[k]['w_design']} um)")
    ax[0, 0].set_xlabel("centre-to-centre pitch (um)")
    ax[0, 0].set_ylabel("coupling length L_c (um)")
    ax[0, 0].legend(fontsize=6)
    save_result("m1_coupling", script=SCRIPT, fig=f,
                inputs={"design": design, "method": "even/odd TE supermodes, L_c = lambda/(2 dn)"},
                outputs={"rows": rows},
                notes="Power transfer over a parallel length L is sin^2(pi L / (2 L_c)); crosstalk below -30 dB "
                      "needs L < 0.02 L_c.")


if __name__ == "__main__":
    parts = sys.argv[1:] or ["width", "design", "bend", "box", "coupling"]
    if parts == ["rest"]:
        parts = ["design", "bend", "box", "coupling"]
    dfile = RESULTS / "m1_width_sweep.json"
    if "width" in parts:
        design = width_sweep()
    else:  # re-select from the stored sweep (cheap)
        design = width_sweep(json.loads(dfile.read_text())["outputs"]["rows"])
    print("design widths:", {k: v["w_design"] for k, v in design.items()}, flush=True)
    if "design" in parts:
        dm = design_modes(design)
        for k, v in dm.items():
            print(k, {kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in v.items()}, flush=True)
    if "bend" in parts:
        print("R_min:", bend_loss(design, None), flush=True)
    if "box" in parts:
        box_leakage(design)
        print("box leakage done", flush=True)
    if "coupling" in parts:
        coupling(design)
        print("coupling done", flush=True)
