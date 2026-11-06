"""Module 22: scaling budget for 1e4, 1e5 and 1e6-atom arrays at 3-4 um pitch, and engine-tiling optics.

Parts:
  budget  geometry and space-bandwidth product; single large-field objective vs zones of fixed corrected field;
          static-array power (hold and loading/imaging depths); transport engine blocks (engines, actuators,
          simultaneous tweezers, laser power at 0.3 mK, field-select positions); loss rates and the repair
          throughput and reservoir they require; analytic steady-state defect density
  tiling  tiling engines along x: grating-lobe geometry, a shared slit plane with per-engine lenslets and a common
          Fourier lens (ABCD check: every engine fills the pupil, offsets map to atom segments), coverage

Outputs: results/m22_budget.json, results/m22_tiling.json
"""
from __future__ import annotations

import itertools
import json
import sys

import numpy as np

from picsim import params as prm
from picsim.abcd import aug_free, aug_lens, chain, trace
from picsim.provenance import RESULTS, save_result
from picsim.transport import trap_depth_K

SCRIPT = "scripts/m22_megascale.py"
LAM = prm.LAM_TRANSPORT.value  # um
NA = prm.NA_OBJ.value
CELL = LAM / (2 * NA)  # um, resolution cell at the atoms
W0_STATIC = 0.85  # um (CHOICE, as in m18)
ALPHA_STATIC = prm.ALPHA_CS_1075.value
ETA_STATIC = 0.3  # laser-to-atoms efficiency of the static array optics (CHOICE)
U_HOLD, U_LOAD = 0.1e-3, 1.0e-3  # K: hold depth during rearrangement (m18) and loading/imaging depth (CHOICE)
F_OBJ = prm.WD.value * np.tan(np.arcsin(NA)) / NA  # mm
A_EMIT = 192 * 1.5  # um, emitter aperture of a 192-arm engine
P_EMIT = 1.5  # um, emitter pitch
SEG = 180.0  # um, row segment served by one engine (88% of its 204 um reach)


def sites_per_segment(pitch):
    return int(SEG // pitch) + 1


# ------------------------------------------------------------------ budget
def budget():
    p_tw_12 = json.loads((RESULTS / "m14_plan.json").read_text())["outputs"]["P_tweezer_adopted_W"]
    p_tw = p_tw_12 * 0.3 / 1.2  # W per tweezer at 0.3 mK (m18/m19 operating point)
    p_site = lambda U: U / trap_depth_K(1.0, W0_STATIC, ALPHA_STATIC)  # W at the atoms per static site
    rows = []
    for N, pitch, f_load in itertools.product((1e4, 1e5, 1e6), (4.0, 3.0), (0.5, 0.85)):
        n_side = int(np.ceil(np.sqrt(N)))
        side = (n_side - 1) * pitch
        cells = side / CELL
        field_r = side / np.sqrt(2) / 2 * 2 / 2  # half-diagonal (um) = side / sqrt(2)
        field_r = side / np.sqrt(2)
        zones = {}
        for zr in (100.0, 250.0):
            zs = int((2 * zr / np.sqrt(2)) // pitch) + 1  # sites per side of the square inscribed in the field
            zones[f"+/-{int(zr)} um"] = {"sites_per_zone": zs * zs, "zones": int(np.ceil(N / (zs * zs)))}
        seg_sites = sites_per_segment(pitch)
        block = {}
        for S, Mx in ((16, 1), (16, 4), (32, 4)):
            tw = S * Mx * seg_sites
            pos = int(np.ceil(n_side / S) * np.ceil(n_side / (seg_sites * Mx)))
            if Mx > 1:
                pos *= 2  # shared slit plane: every other segment per position (tiling part)
            block[f"{S} rows x {Mx} segments"] = {"engines": S * Mx, "fast_actuators": 383 * S * Mx,
                                                  "simultaneous_tweezers": tw, "laser_W_0.3mK": tw * p_tw,
                                                  "field_positions": pos}
        # reservoir needed so that, after loading at f_load, reservoir atoms fill the computation holes:
        # r f >= (1 - r)(1 - f)  ->  r >= (1 - f)
        r_res = 1 - f_load
        loss = {f"tau {int(t)} s": N / t for t in (300.0, 3000.0, 30000.0)}
        img = {f"eps {e:g}": N * e for e in (1e-3, 1e-4, 1e-5)}
        rows.append({
            "N": N, "pitch_um": pitch, "f_load": f_load, "sites_per_side": n_side, "side_mm": side / 1e3,
            "cells_per_side": cells, "space_bandwidth_cells": cells**2, "single_objective_field_radius_mm": field_r / 1e3,
            "zones": zones, "static_W_at_atoms_hold": N / (1 - r_res) * p_site(U_HOLD) if False else
            N * p_site(U_HOLD) / (1 - r_res), "static_W_at_atoms_load": N * p_site(U_LOAD) / (1 - r_res),
            "static_laser_W_load": N * p_site(U_LOAD) / (1 - r_res) / ETA_STATIC,
            "reservoir_fraction_of_sites": r_res, "total_sites": N / (1 - r_res),
            "block": block, "vacuum_losses_per_s": loss, "imaging_losses_per_image": img})
        print(N, pitch, f_load, round(side / 1e3, 2), "mm", {k: round(v, 1) for k, v in rows[-1].items()
                                                           if isinstance(v, float) and "W" in k}, flush=True)
    out = {"p_site_W_per_mK": p_site(1e-3), "p_tweezer_0.3mK_W": p_tw, "rows": rows}
    save_result("m22_budget", script=SCRIPT,
                inputs={"static": {"w0_um": W0_STATIC, "alpha_au": ALPHA_STATIC, "eta_laser_to_atoms": ETA_STATIC,
                                   "U_hold_mK": U_HOLD * 1e3, "U_load_mK": U_LOAD * 1e3},
                        "engine": {"ports": 192, "segment_um": SEG, "actuators": 383},
                        "reservoir": "fraction r with r f >= (1 - r)(1 - f): r = 1 - f (one-shot fill)",
                        "objective_cell_um": CELL},
                outputs=out, notes="Static powers include the reservoir sites (total sites = N / (1 - r)).")
    return out


# ------------------------------------------------------------------ static potential: tweezers vs lattice
def lattice():
    """Laser power for the static potential: focused tweezers (one beam per site) against a 2D optical lattice of
    the same spacing (two orthogonal retro-reflected beams, added incoherently: peak 8 I_b; coherent four-beam: 16).
    The lattice beams must keep the depth within -20% over the array (tolerance from m18_uniformity): a Gaussian
    beam of waist w = r_c sqrt(2 / ln(1/0.8)) for corner radius r_c, or a flat-top beam of area 1.2 x the array.
    Axial confinement of the lattice needs a light sheet (not included; small, see notes)."""
    from picsim.transport import AU_POL, EPS0, C as C0
    I_per_K = lambda U: U * 1.380649e-23 * 2 * EPS0 * C0 / (ALPHA_STATIC * AU_POL)  # W/m^2 of peak intensity
    rows = []
    for N, (px, py), U in itertools.product((1e4, 1e5, 1e6), ((4.0, 4.0), (3.0, 4.0)), (0.1e-3, 1.0e-3)):
        n_tot = N / 0.75  # sites incl. a 25% reservoir (85% loading, 3:1 layout of m23)
        area = n_tot * px * py * 1e-12  # m^2
        L = np.sqrt(area)
        I0 = I_per_K(U)
        tweezer = n_tot * (np.pi * (W0_STATIC * 1e-6) ** 2 / 2) * I0 / ETA_STATIC
        r_c = L / np.sqrt(2)
        w = r_c * np.sqrt(2 / np.log(1 / 0.8))
        out = {"N": N, "pitch_um": [px, py], "U_mK": U * 1e3, "array_side_mm": L * 1e3, "tweezers_W": tweezer}
        for G, gname in ((8, "incoherent"), (16, "coherent")):
            out[f"lattice_gauss_{gname}_W"] = np.pi * w**2 / 2 * 2 * I0 / G / ETA_STATIC  # two beams, eta for delivery
            out[f"lattice_flattop_{gname}_W"] = 1.2 * area * 2 * I0 / G / ETA_STATIC
        for B in (10, 100):
            out[f"lattice_gauss_incoherent_cavity{B}_W"] = out["lattice_gauss_incoherent_W"] / B
        out["gaussian_waist_mm"] = w * 1e3
        rows.append(out)
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in out.items()}, flush=True)
    save_result("m22_lattice", script=SCRIPT,
                inputs={"alpha_au": ALPHA_STATIC, "w0_tweezer_um": W0_STATIC, "eta_delivery": ETA_STATIC,
                        "reservoir": "25% of sites", "uniformity": "depth >= 0.8 of peak over the array (m18: +/-20% "
                        "costs nothing)", "interference_gain": {"incoherent (two 1D standing waves)": 8, "coherent four-beam": 16},
                        "cavity_buildup": [10, 100]},
                outputs={"rows": rows},
                notes="Lattice per site ~ 15.4 a^2 I / G (Gaussian, -20% corners) or 1.2 a^2 I / G (flat top), against "
                      "pi w0^2 I / 2 for a tweezer. A light sheet for axial confinement adds roughly I x (pi/2) w_x w_z, "
                      "small for w_z of a few um; not included.")
    return rows


# ------------------------------------------------------------------ tiling
def tiling():
    x_scale = A_EMIT / (2 * NA)  # um at the atoms per unit s
    s_cut = LAM / (2 * P_EMIT)
    reach = s_cut * x_scale  # half-range of one engine at the atoms
    lobe = (LAM / P_EMIT) * x_scale  # first grating lobe offset at the atoms
    out = {"x_scale_um": x_scale, "reach_half_um": reach, "first_lobe_offset_um": lobe, "segment_um": SEG,
           "lobe_zone_um": [lobe - reach, lobe + reach],
           "neighbour_zone_um": [SEG - reach, SEG + reach]}
    # shared slit plane Sigma: slit zone of each engine 2 f1 s_cut wide; lobes occupy f1 (lambda/p -+ s_cut).
    # Opaque gaps must contain the lobes: engine pitch Pe >= f1 (lambda/p + s_cut) + f1 s_cut = 2 f1 lambda/p
    f1 = 3.12  # mm (m13 x relay)
    Pe_min = f1 * (LAM / P_EMIT + 2 * s_cut)
    out["shared_slit"] = {"f1_mm": f1, "slit_zone_mm": 2 * f1 * s_cut, "lobe_zone_mm": [f1 * (LAM / P_EMIT - s_cut),
                                                                                     f1 * (LAM / P_EMIT + s_cut)],
                          "engine_pitch_min_mm": Pe_min, "coverage_per_position": 2 * f1 * s_cut / Pe_min}
    # ABCD of the combiner: emitter plane -> lenslet i (f1, centred X_i) -> Sigma -> common lens F (f_c) -> pupil
    M_pupil = 2 * F_OBJ * NA / (A_EMIT * 1e-3)  # pupil diameter / emitter aperture
    f_c = M_pupil * f1  # common lens: pupil image of the emitter array has magnification f_c / f1
    seg_pitch_atoms = (F_OBJ / f_c) * Pe_min * 1e3  # um between adjacent engines' segment centres
    checks = []
    for i in range(4):
        X = i * Pe_min
        M = chain(aug_free(f1), aug_lens(f1, X), aug_free(f1), aug_free(f_c), aug_lens(f_c, 0.0), aug_free(f_c))
        # rays from the engine's emitter centre at angles theta: positions at the pupil (should all be 0)
        y, th = trace(M, np.full(5, X), np.linspace(-0.3, 0.3, 5))
        # rays from emitter offset xi (mm) parallel to the axis: pupil position (magnification) and angle
        y2, th2 = trace(M, np.array([X + 0.1]), np.array([0.0]))
        # atom position of the beam steered to s: chief ray from the emitter centre with angle s -> pupil angle
        _, th3 = trace(M, np.array([X]), np.array([0.2]))
        checks.append({"engine": i, "X_mm": X, "pupil_y_spread_mm": float(np.ptp(y)),
                       "pupil_magnification": float(y2[0] / 0.1), "atom_x_um_of_s0": float(F_OBJ * (-th[2]) * 1e3),
                       "atom_x_um_of_s0.2": float(F_OBJ * (-th3[0]) * 1e3)})
    out["combiner"] = {"f_c_mm": f_c, "pupil_magnification_target": M_pupil, "segment_pitch_at_atoms_um": seg_pitch_atoms,
                       "checks": checks}
    print(json.dumps(out, indent=1, default=float), flush=True)
    save_result("m22_tiling", script=SCRIPT,
                inputs={"engine": "192 emitters at 1.5 um (A = 288 um)", "f1_mm": f1, "NA": NA, "f_obj_mm": F_OBJ},
                outputs=out,
                notes="Shared slit plane: each engine's slit zone must be followed by an opaque gap that holds its "
                      "grating lobes, so one block position addresses every other segment.")
    return out


if __name__ == "__main__":
    for p in sys.argv[1:] or ["budget", "lattice", "tiling"]:
        {"budget": budget, "lattice": lattice, "tiling": tiling}[p]()
