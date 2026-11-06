"""Module 17: scaling beyond one row of ~41 tweezers without large penalties.

Parts:
  moves  how shallow can a transport tweezer be? Classical Monte Carlo of thermal Cs atoms in a 1D Gaussian
         potential (the steering axis, waist from the D1 spot) whose centre follows a minimum-jerk move, smooth
         or sampled-and-held at the 200 kHz update rate. Loss and energy gain per move vs depth, move time and
         distance; minimum depth and power per tweezer for given criteria
  table  scaling configurations (DERIVED from modules 12-16 and the moves part): simultaneous tweezers, sites
         served, chips, fast actuators, laser power, passes per rearrangement cycle

Outputs: results/m17_moves.json (+ .png), results/m17_table.json
"""
from __future__ import annotations

import json
import sys

import numpy as np
from scipy.constants import hbar, k as kB

from picsim import params as prm
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result
from picsim.transport import min_jerk

SCRIPT = "scripts/m17_scale.py"
M = prm.CS_MASS.value
W_X = 1.68 / np.sqrt(2 * np.log(2)) * 1e-6  # 1/e^2 intensity radius along x from the D1 spot (FWHM 1.68 um)
T_ATOM = 10e-6  # K, atom temperature before the move (CHOICE)
F_UPDATE = 200e3
LOSS_MAX = 1e-3  # per move (CHOICE)
DN_MAX = 1.0  # classical energy gain per move in units of hbar omega (CHOICE)


def omega(U):
    return np.sqrt(4 * U * kB / (M * W_X**2))


def simulate(U_K, d, T, staircase, n=4000, seed=0):
    rng = np.random.default_rng(seed)
    U = U_K * kB
    w = omega(U_K)
    # thermal start in the harmonic approximation, bound atoms only
    x = rng.normal(0, np.sqrt(kB * T_ATOM / (M * w**2)), n)
    v = rng.normal(0, np.sqrt(kB * T_ATOM / M), n)

    def force(x, xc):
        dx = x - xc
        return -U * 4 * dx / W_X**2 * np.exp(-2 * dx**2 / W_X**2)

    def energy(x, v, xc):
        return 0.5 * M * v**2 - U * np.exp(-2 * (x - xc) ** 2 / W_X**2)

    E0 = energy(x, v, 0.0)
    bound0 = E0 < 0
    dt = min(2 * np.pi / w / 80, 1 / F_UPDATE / 8)
    t_end = T * 1.3
    t = 0.0

    def xc_at(t):
        tt = np.floor(t * F_UPDATE) / F_UPDATE if staircase else t
        return min_jerk(tt, T, d)

    a = force(x, xc_at(0.0)) / M
    while t < t_end:
        v += 0.5 * dt * a
        x += dt * v
        t += dt
        a = force(x, xc_at(t)) / M
        v += 0.5 * dt * a
    E1 = energy(x, v, d)
    lost = (E1 >= 0) & bound0
    ok = bound0 & ~lost
    dE = (E1[ok] - E0[ok]).mean() if ok.any() else np.nan
    return {"loss": float(lost.sum() / bound0.sum()), "dE_uK": float(dE / kB * 1e6), "dn": float(dE / (hbar * w)),
            "f_trap_kHz": float(w / 2 / np.pi / 1e3)}


def moves():
    rows = []
    for d_um in (12.0, 40.0):
        for T in (50e-6, 100e-6, 200e-6, 400e-6):
            for U in (0.05e-3, 0.1e-3, 0.2e-3, 0.3e-3, 0.5e-3, 0.8e-3, 1.2e-3):
                for st in (False, True):
                    r = simulate(U, d_um * 1e-6, T, st)
                    rows.append({"d_um": d_um, "T_us": T * 1e6, "U_mK": U * 1e3, "staircase_200kHz": st, **r})
                    print(rows[-1], flush=True)
    # minimum depth meeting the criteria, per (d, T), staircase case
    mins = []
    for d_um in (12.0, 40.0):
        for T in (50.0, 100.0, 200.0, 400.0):
            ok = [r for r in rows if r["d_um"] == d_um and r["T_us"] == T and r["staircase_200kHz"]
                  and r["loss"] <= LOSS_MAX and r["dn"] <= DN_MAX]
            Umin = min((r["U_mK"] for r in ok), default=None)
            mins.append({"d_um": d_um, "T_us": T, "U_min_mK": Umin})
    f, ax = fig(1, 2, h=3.0)
    for k, d_um in enumerate((12.0, 40.0)):
        for T in (50.0, 100.0, 200.0, 400.0):
            sel = [r for r in rows if r["d_um"] == d_um and r["T_us"] == T and r["staircase_200kHz"]]
            ax[0, k].semilogy([r["U_mK"] for r in sel], [max(r["loss"], 1e-5) for r in sel], "-o", ms=3,
                              label=f"T = {T:.0f} us")
        ax[0, k].axhline(LOSS_MAX, color="k", lw=0.8, ls="--")
        ax[0, k].set_xlabel("depth (mK)")
        ax[0, k].set_ylabel(f"loss per {d_um:.0f} um move")
        ax[0, k].legend(fontsize=7)
    save_result("m17_moves", script=SCRIPT, fig=f,
                inputs={"model": "classical, 1D Gaussian potential along the move, velocity Verlet", "w_x_um": W_X * 1e6,
                        "T_atom_uK": T_ATOM * 1e6, "atoms": 4000, "trajectory": "minimum jerk, settle 0.3 T",
                        "update": "smooth or sample-and-hold at 200 kHz", "criteria": {"loss": LOSS_MAX, "dn": DN_MAX}},
                outputs={"rows": rows, "minimum_depth": mins},
                notes="1D model: transverse and axial motion and their anharmonic coupling are not included.")
    return rows, mins


def table():
    mv = json.loads((RESULTS / "m17_moves.json").read_text())["outputs"]["minimum_depth"]
    plan = json.loads((RESULTS / "m14_plan.json").read_text())["outputs"]
    p12 = plan["P_tweezer_adopted_W"]  # W per tweezer at 1.2 mK, glued collimator
    u200 = [m for m in mv if m["d_um"] == 12.0 and m["T_us"] == 200.0][0]["U_min_mK"]
    u100 = [m for m in mv if m["d_um"] == 12.0 and m["T_us"] == 100.0][0]["U_min_mK"]
    site = 4.0
    rows = []

    def cfg(name, ports, chips, stacks, rows_field, passes_note=""):
        row_len = 0.88 * ports * 1.0644  # um, 88% of the port field (m12)
        per_row = int(row_len // site) + 1  # sites including both ends, as in m12
        sim_rows = chips * stacks
        simult = per_row * sim_rows
        sites_rows = max(sim_rows, rows_field)
        passes = int(np.ceil(sites_rows / sim_rows))
        act = (2 * ports - 1) * chips * stacks
        rows.append({"configuration": name, "ports_per_chip": ports, "chips": chips * stacks, "row_length_um": row_len,
                     "tweezers_per_row": per_row, "simultaneous_tweezers": simult,
                     "sites_served": per_row * sites_rows, "passes_per_cycle": passes, "fast_actuators": act,
                     "actuators_per_simultaneous_tweezer": act / simult,
                     "power_W_1.2mK": simult * p12, f"power_W_{u200}mK_(200us moves)": simult * p12 * u200 / 1.2,
                     f"power_W_{u100}mK_(100us moves)": simult * p12 * u100 / 1.2, "note": passes_note})

    cfg("one D1 chip", 128, 1, 1, 1)
    cfg("one chip, 192 ports (objective field +/-100 um)", 192, 1, 1, 1)
    cfg("stack of 16 chips", 192, 16, 1, 16)
    cfg("stack of 16 + row select (y tilt at a pupil-conjugate plane)", 192, 16, 1, 50, "50 rows = +/-100 um at 4 um")
    cfg("two stacks combined by polarisation + row select", 192, 16, 2, 50, "second stack orthogonally polarised")
    out = {"rows": rows, "site_pitch_um": site, "P_tweezer_1.2mK_W": p12, "U_min_mK_200us_12um": u200,
           "U_min_mK_100us_12um": u100}
    for r in rows:
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    save_result("m17_table", script=SCRIPT,
                inputs={"site_pitch_um": site, "row_length": "88% of the port field (m12)", "P_tweezer": "m14_plan",
                        "depths": "m17_moves minimum depth for 12 um moves at 200 kHz updates"},
                outputs=out, notes="DERIVED. Row select moves all rows of a stack rigidly in y; atoms wait in the static "
                                   "array between passes.")
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["moves", "table"]
    for p in parts:
        {"moves": moves, "table": table}[p]()
