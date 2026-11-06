"""Module 18: 3D transport with handoff between the static 1075 nm array and a moving 1064 nm tweezer.

Parts:
  field  3D intensity of one moving tweezer of the adopted design (D1 chip, least-squares synthesis; Weierstrass +
         aplanat collimator; circular NA 0.5 pupil; scalar focusing with defocus), tabulated around its focus
  scan   pick-up, move by 12 um (three sites), drop-off, for static hold depth, moving depth, ramp time and move
         time; on the row (through two empty static sites) or along a lane 2 um between rows
  tol    misalignment of the two foci (z), lateral offsets, atom temperature, update staircase, longer moves
Static traps stay on at a global hold depth throughout (a fixed array cannot switch single sites).

Outputs: results/m18_trap.npz, results/m18_field.json, results/m18_scan.json, results/m18_tol.json
"""
from __future__ import annotations

import itertools
import json
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.constants import k as kB

from picsim import materials as mat
from picsim import params as prm
from picsim.opa import collimator as C
from picsim.opa.focal_plane import FocalPlane, LinearSynth
from picsim.opa.stack import FacetY, Focus2D, RelayY, fwhm
from picsim.provenance import RESULTS, save_result
from picsim.transport3d import M_CS, GridTrap, Sequence, StaticArray, run

SCRIPT = "scripts/m18_transport3d.py"
LAM = prm.LAM_TRANSPORT.value
LAM_S = prm.LAM_STATIC.value
W0_S = 0.85e-6  # static-trap waist (CHOICE: NA 0.5 diffraction limit 0.68 um at 1075 nm, mildly apodised)
PITCH = 4e-6
D1 = {"N_ports": 128, "N_arms": 192, "w_emit": 1.2, "p_out": 1.5, "w_port": 1.7, "w_arm": 1.7}


# ------------------------------------------------------------------ moving-trap field
def field():
    n = float(mat.n_sio2(LAM))
    T = np.load(RESULTS / "m10_T_opt_A.npz")["T"]
    fp = FocalPlane(A=T, **D1)
    ls = LinearSynth(fp)
    r = RelayY(FacetY(LAM, 1.2, 0.4, n_m=n), 0.7, NA=0.5, NA_FAC=0.9, lens=C.fac_design(n).lensmap(LAM, 0.9, coated=True))
    F = Focus2D(fp, r)
    w = ls.weights([0.0])
    w = w / np.linalg.norm(w)
    x = np.arange(-3.5, 3.5001, 0.05)
    y = np.arange(-3.0, 3.0001, 0.05)
    z = np.arange(-10.0, 10.0001, 0.25)
    I = np.empty((len(x), len(y), len(z)), dtype=np.float32)
    for k, zz in enumerate(z):
        I[:, :, k] = F.intensity(w, x, y, dz_a=zz).T
    i, j, kk = np.unravel_index(np.argmax(I), I.shape)
    g = I / I.max()
    np.savez_compressed(RESULTS / "m18_trap.npz", x=(x - x[i]) * 1e-6, y=(y - y[j]) * 1e-6, z=(z - z[kk]) * 1e-6, g=g)
    ax = {"fwhm_x_um": fwhm(x, g[:, j, kk]), "fwhm_y_um": fwhm(y, g[i, :, kk]), "fwhm_z_um": fwhm(z, g[i, j, :]),
          "peak_offset_um": [float(x[i]), float(y[j]), float(z[kk])]}
    # harmonic frequencies at 1 mK from the curvature at the peak
    tr = GridTrap((x - x[i]) * 1e-6, (y - y[j]) * 1e-6, (z - z[kk]) * 1e-6, g.astype(float))
    h = 0.05e-6
    U = 1e-3 * kB
    freqs = []
    for e in np.eye(3):
        gp, _ = tr.eval(np.array([e * h]))
        gm, _ = tr.eval(np.array([-e * h]))
        g0, _ = tr.eval(np.zeros((1, 3)))
        curv = (gp[0] + gm[0] - 2 * g0[0]) / h**2
        freqs.append(float(np.sqrt(-U * curv / M_CS) / (2 * np.pi)))
    ax["freq_at_1mK_kHz_xyz"] = [f / 1e3 for f in freqs]
    print(ax, flush=True)
    save_result("m18_field", script=SCRIPT,
                inputs={"design": "D1 (BPM matrix), LinearSynth tweezer at x = 0; Weierstrass + aplanat; beta 0.7 um",
                        "grid_um": {"x": [-3.5, 3.5, 0.05], "y": [-3, 3, 0.05], "z": [-10, 10, 0.25]}},
                outputs=ax, notes="Grid stored in results/m18_trap.npz (normalised intensity, peak at the origin).")
    return ax


# ------------------------------------------------------------------ simulation helpers
_G = {}


def trap():
    if "trap" not in _G:
        d = np.load(RESULTS / "m18_trap.npz")
        _G["trap"] = GridTrap(d["x"], d["y"], d["z"], d["g"].astype(float))
    return _G["trap"]


def static_array(depth_mK, dz=0.0, src_factor=1.0, dst_factor=1.0, sigma=0.0, seed=0, d_um=12.0, pitch=PITCH,
                 rows=None):
    """Static sites at x pitch `pitch` (15 per row) on rows at y = `rows` (m); default rows at 0 and +/- pitch."""
    xs = np.arange(-2, 13) * pitch
    rows = (-pitch, 0.0, pitch) if rows is None else rows
    sites = np.array([[x, y, dz] for y in rows for x in xs])
    fac = np.ones(len(sites))
    if sigma:
        fac = np.clip(1 + sigma * np.random.default_rng(seed).standard_normal(len(sites)), 0.2, None)
    fac[np.argmin(np.linalg.norm(sites - np.array([0, 0, dz]), axis=1))] *= src_factor
    fac[np.argmin(np.linalg.norm(sites - np.array([d_um * 1e-6, 0, dz]), axis=1))] *= dst_factor
    return StaticArray(sites, W0_S, LAM_S * 1e-6, depth_mK * 1e-3 * kB, factors=fac)


def one(p):
    st = static_array(p["Us_mK"], p.get("static_dz_um", 0.0) * 1e-6, p.get("src_factor", 1.0), p.get("dst_factor", 1.0),
                      p.get("sigma", 0.0), p.get("seed", 0), p.get("d_um", 12.0), p.get("pitch_um", 4.0) * 1e-6,
                      None if p.get("rows_um") is None else tuple(np.array(p["rows_um"]) * 1e-6))
    seq = Sequence(src=np.zeros(3), dst=np.array([p.get("d_um", 12.0) * 1e-6, 0, 0]), depth_move=p["Um_mK"] * 1e-3 * kB,
                   t_ramp=p["t_ramp_us"] * 1e-6, t_move=p["t_move_us"] * 1e-6, lane_dy=p.get("lane_um", 0.0) * 1e-6,
                   f_update=p.get("f_update", 200e3), offset=np.array(p.get("offset_um", [0, 0, 0])) * 1e-6,
                   t_lane=p.get("t_lane_us", 100.0) * 1e-6)
    site = p.get("bystander_um")
    res = run(st, trap(), seq, T_atom=p.get("T_uK", 10.0) * 1e-6, n=p.get("n", 1500), seed=p.get("seed", 0),
              gravity=p.get("gravity", True), atom_site=None if site is None else np.array(site) * 1e-6)
    return {**p, **res}


def pool_map(params):
    trap()
    with ProcessPoolExecutor(4) as ex:
        return list(ex.map(one, params))


# ------------------------------------------------------------------ scan
def scan():
    P = []
    for Us, Um, tr, tm, lane in itertools.product((0.1, 0.2), (0.2, 0.3, 0.5), (50.0, 100.0), (200.0, 400.0), (0.0, 2.0)):
        if Um <= Us:
            continue
        P.append({"Us_mK": Us, "Um_mK": Um, "t_ramp_us": tr, "t_move_us": tm, "lane_um": lane})
    rows = pool_map(P)
    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    save_result("m18_scan", script=SCRIPT,
                inputs={"static": f"1075 nm Gaussian sites, w0 {W0_S*1e6} um, 4 um pitch, three rows, global hold depth",
                        "moving": "results/m18_trap.npz", "move": "x = 0 -> 12 um (two empty sites crossed on the row)",
                        "lane": "2 um in y via the row-select tilt, 100 us each way", "update": "200 kHz sample-and-hold "
                        "of position and depth", "atoms": "1500, 10 uK thermal in the source trap", "gravity": "along y"},
                outputs={"rows": rows},
                notes="delivered: bound in the target static trap after the moving trap is off; dE: mean energy above "
                      "the static trap bottom, final minus initial.")
    return rows


# ------------------------------------------------------------------ tolerances
def tol():
    P = []
    for Um in (0.3, 0.5):
        base = {"Us_mK": 0.1, "Um_mK": Um, "t_ramp_us": 100.0, "t_move_us": 400.0, "lane_um": 2.0}
        tag = f"Um {Um}: "
        P.append(dict(base, label=tag + "base (lane)"))
        for dz in (0.5, 1.0, 2.0, 4.0):
            P.append(dict(base, static_dz_um=dz, label=tag + f"focus offset {dz} um"))
        for dy in (0.1, 0.2, 0.4):
            P.append(dict(base, offset_um=[0, dy, 0], label=tag + f"y offset {dy} um"))
        for dx in (0.1, 0.2, 0.4):
            P.append(dict(base, offset_um=[dx, 0, 0], label=tag + f"x offset {dx} um"))
        for T in (5.0, 20.0, 40.0):
            P.append(dict(base, T_uK=T, label=tag + f"atom {T} uK"))
        P.append(dict(base, f_update=None, label=tag + "smooth (no staircase)"))
        P.append(dict(base, gravity=False, label=tag + "no gravity"))
        for d, tm in ((40.0, 400.0), (40.0, 800.0)):
            P.append(dict(base, d_um=d, t_move_us=tm, label=tag + f"{d} um in {tm} us"))
        for tr in (25.0, 200.0):
            P.append(dict(base, t_ramp_us=tr, label=tag + f"ramp {tr} us"))
    rows = pool_map(P)
    for r in rows:
        print(r["label"], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                           if k in ("delivered", "left_behind", "misplaced", "lost", "dE_uK")}, flush=True)
    save_result("m18_tol", script=SCRIPT, inputs={"base": "Us 0.1 mK, ramps 100 us, 12 um move in 400 us along a 2 um "
                                                          "lane; Um 0.3 and 0.5 mK", "static_w0_um": W0_S * 1e6},
                outputs={"rows": rows})
    return rows


# ------------------------------------------------------------------ atoms waiting next to the lane
def bystander():
    P = []
    for Um in (0.3, 0.5):
        for site, name in (([8.0, 0.0, 0.0], "same row, passed at 2 um"), ([8.0, 4.0, 0.0], "next row, passed at 2 um"),
                           ([12.0, 4.0, 0.0], "next row, beside the drop site")):
            P.append({"Us_mK": 0.1, "Um_mK": Um, "t_ramp_us": 100.0, "t_move_us": 400.0, "lane_um": 2.0,
                      "bystander_um": site, "label": f"Um {Um}: {name}"})
        P.append({"Us_mK": 0.1, "Um_mK": Um, "t_ramp_us": 100.0, "t_move_us": 400.0, "lane_um": 0.0,
                  "bystander_um": [8.0, 4.0, 0.0], "label": f"Um {Um}: next row, on-row transport (4 um away)"})
    rows = pool_map(P)
    for r in rows:
        print(r["label"], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                           if k in ("delivered", "lost", "misplaced", "dE_uK")}, flush=True)
    save_result("m18_bystander", script=SCRIPT,
                inputs={"move": "loaded tweezer from (0,0) to (12,0) um along a lane at y = 2 um (or on the row)",
                        "bystander": "atom in a static site that is not addressed", "hold_mK": 0.1},
                outputs={"rows": rows}, notes="'delivered' here means the bystander stayed in its own site.")
    return rows


# ------------------------------------------------------------------ hold-depth margin for bystanders
def margin():
    d = np.load(RESULTS / "m18_trap.npz")
    g, x, y, z = d["g"], d["x"] * 1e6, d["y"] * 1e6, d["z"] * 1e6
    i, k = np.argmin(abs(x)), np.argmin(abs(z))
    side = {f"{yy}": {"focal_plane": float(g[i, np.argmin(abs(y - yy)), k]),
                      "max_over_z": float(g[i, np.argmin(abs(y - yy)), :].max())} for yy in (1.5, 2.0, 2.5)}
    P = []
    for T, Us, Um in ((10.0, 0.1, 0.3), (10.0, 0.2, 0.6), (5.0, 0.1, 0.3), (5.0, 0.1, 0.5), (10.0, 0.2, 1.0)):
        for site in ([8.0, 0.0, 0.0], [12.0, 4.0, 0.0]):
            P.append({"Us_mK": Us, "Um_mK": Um, "t_ramp_us": 100.0, "t_move_us": 400.0, "lane_um": 2.0,
                      "bystander_um": site, "T_uK": T, "n": 3000, "label": f"{T} uK, {Us}/{Um} mK, bystander {site}"})
        P.append({"Us_mK": Us, "Um_mK": Um, "t_ramp_us": 100.0, "t_move_us": 400.0, "lane_um": 2.0, "T_uK": T,
                  "n": 3000, "label": f"{T} uK, {Us}/{Um} mK, transported"})
    rows = pool_map(P)
    for r in rows:
        print(r["label"], round(r["delivered"], 4), round(r["lost"], 4), r["dE_uK"], flush=True)
    save_result("m18_margin", script=SCRIPT,
                inputs={"atoms": 3000, "move": "12 um in 400 us along a 2 um lane", "bystanders": "(8,0) and (12,4) um"},
                outputs={"moving_trap_intensity_at_lane_distance": side, "rows": rows},
                notes="Bystander rows: 'delivered' = stayed in its own site.")
    return rows


# ------------------------------------------------------------------ slower lane steps (3D check of m19)
def slowlane():
    P = []
    for T, Us, Um in ((5.0, 0.1, 0.3), (5.0, 0.1, 0.5), (10.0, 0.2, 0.6)):
        for tr, tl in ((200.0, 400.0),):
            P.append({"Us_mK": Us, "Um_mK": Um, "t_ramp_us": tr, "t_lane_us": tl, "t_move_us": 400.0, "lane_um": 2.0,
                      "T_uK": T, "n": 3000, "label": f"{T} uK, {Us}/{Um} mK, ramp {tr}, lane {tl} us"})
            P.append({"Us_mK": Us, "Um_mK": Um, "t_ramp_us": tr, "t_lane_us": tl, "t_move_us": 400.0, "lane_um": 2.0,
                      "T_uK": T, "n": 3000, "bystander_um": [12.0, 4.0, 0.0],
                      "label": f"{T} uK, {Us}/{Um} mK, ramp {tr}, lane {tl} us, bystander beside drop"})
    rows = pool_map(P)
    for r in rows:
        print(r["label"], round(r["delivered"], 4), round(r["lost"], 4), r["dE_uK"], flush=True)
    save_result("m18_slowlane", script=SCRIPT, inputs={"atoms": 3000, "move": "12 um in 400 us"}, outputs={"rows": rows})
    return rows


# ------------------------------------------------------------------ anisotropic 3 x 4 um lattice and 3 um corridors
def aniso():
    layouts = {
        "3 x 4 um (rows 4 um apart, lane 2 um)": {"pitch_um": 3.0, "rows_um": [-4.0, 0.0, 4.0], "lane_um": 2.0,
                                                  "by": [[3.0, 0.0, 0.0], [-3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [9.0, 4.0, 0.0]]},
        "3 um with corridors (rows 6 um apart, lane 3 um)": {"pitch_um": 3.0, "rows_um": [-6.0, 0.0, 6.0], "lane_um": 3.0,
                                                             "by": [[3.0, 0.0, 0.0], [6.0, 0.0, 0.0], [9.0, 6.0, 0.0]]},
        "4 x 4 um reference (in-row neighbour)": {"pitch_um": 4.0, "rows_um": [-4.0, 0.0, 4.0], "lane_um": 2.0,
                                                  "by": [[4.0, 0.0, 0.0]], "d": 12.0},
    }
    P = []
    for name, L in layouts.items():
        d = L.get("d", 9.0)
        for T, Us, Um in ((5.0, 0.1, 0.3), (5.0, 0.1, 0.5), (10.0, 0.2, 0.6)):
            base = {"Us_mK": Us, "Um_mK": Um, "t_ramp_us": 200.0, "t_lane_us": 400.0, "t_move_us": 400.0, "T_uK": T,
                    "n": 3000, "pitch_um": L["pitch_um"], "rows_um": L["rows_um"], "lane_um": L["lane_um"], "d_um": d}
            tag = f"{name} | {T} uK, {Us}/{Um} mK: "
            P.append(dict(base, label=tag + f"transported {d} um"))
            for b in L["by"]:
                P.append(dict(base, bystander_um=b, label=tag + f"bystander at {b[:2]}"))
    rows = pool_map(P)
    for r in rows:
        print(r["label"], round(r["delivered"], 4), round(r["lost"], 4), r["dE_uK"], flush=True)
    save_result("m18_aniso", script=SCRIPT,
                inputs={"layouts": {k: {kk: vv for kk, vv in v.items()} for k, v in layouts.items()}, "ramps_us": 200,
                        "lane_steps_us": 400, "move_us": 400, "atoms": 3000},
                outputs={"rows": rows}, notes="Bystander rows: 'delivered' = stayed in its own site.")
    return rows


# ------------------------------------------------------------------ 3 um pitch (lane 1.5 um)
def pitch3():
    P = []
    for T, Us, Um in ((5.0, 0.1, 0.3), (5.0, 0.1, 0.5), (10.0, 0.2, 0.6), (10.0, 0.2, 1.0)):
        base = {"Us_mK": Us, "Um_mK": Um, "t_ramp_us": 200.0, "t_lane_us": 400.0, "t_move_us": 400.0,
                "lane_um": 1.5, "pitch_um": 3.0, "d_um": 9.0, "T_uK": T, "n": 3000}
        tag = f"{T} uK, {Us}/{Um} mK: "
        P.append(dict(base, label=tag + "transported 9 um"))
        P.append(dict(base, bystander_um=[6.0, 0.0, 0.0], label=tag + "bystander same row, passed at 1.5 um"))
        P.append(dict(base, bystander_um=[9.0, 3.0, 0.0], label=tag + "bystander next row, beside drop"))
    rows = pool_map(P)
    for r in rows:
        print(r["label"], round(r["delivered"], 4), round(r["lost"], 4), r["dE_uK"], flush=True)
    save_result("m18_pitch3", script=SCRIPT,
                inputs={"pitch_um": 3.0, "lane_um": 1.5, "move": "9 um (3 sites) in 400 us", "ramps_us": 200,
                        "lane_steps_us": 400, "atoms": 3000},
                outputs={"rows": rows}, notes="Bystander rows: 'delivered' = stayed in its own site.")
    return rows


# ------------------------------------------------------------------ static-array depth uniformity
def uniformity():
    P = []
    for T, Us, Um in ((5.0, 0.1, 0.3), (5.0, 0.1, 0.5), (10.0, 0.2, 0.6)):
        base = {"Us_mK": Us, "Um_mK": Um, "t_ramp_us": 100.0, "t_move_us": 400.0, "lane_um": 2.0, "T_uK": T, "n": 2000}
        tag = f"{T} uK, {Us}/{Um} mK: "
        for f in (1.2, 1.4, 1.6):
            P.append(dict(base, src_factor=f, label=tag + f"source site x{f}"))
        for f in (0.8, 0.6, 0.4):
            P.append(dict(base, dst_factor=f, label=tag + f"target site x{f}"))
        for sg in (0.05, 0.1, 0.2):
            for sd in (1, 2):
                P.append(dict(base, sigma=sg, seed=sd, label=tag + f"random {int(sg*100)}% rms, seed {sd}"))
    rows = pool_map(P)
    for r in rows:
        print(r["label"], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                           if k in ("delivered", "left_behind", "lost", "dE_uK")}, flush=True)
    save_result("m18_uniformity", script=SCRIPT,
                inputs={"move": "12 um in 400 us, lane 2 um, ramps 100 us", "atoms": 2000,
                        "random": "independent Gaussian depth factors on all 45 sites (clipped at 0.2)"},
                outputs={"rows": rows})
    return rows


# ------------------------------------------------------------------ depth vs atom temperature
def temp():
    P = []
    for T, Us, Ums in ((10.0, 0.1, (0.3, 0.5)), (20.0, 0.2, (0.6, 1.0)), (40.0, 0.4, (1.2,))):
        for Um in Ums:
            P.append({"Us_mK": Us, "Um_mK": Um, "t_ramp_us": 100.0, "t_move_us": 400.0, "lane_um": 2.0, "T_uK": T,
                      "n": 3000, "label": f"{T} uK, hold {Us} mK, moving {Um} mK"})
    rows = pool_map(P)
    for r in rows:
        print(r["label"], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                           if k in ("delivered", "lost", "dE_uK")}, flush=True)
    save_result("m18_temp", script=SCRIPT, inputs={"rule": "hold depth ~10 kT, moving depth ~3x hold", "atoms": 3000,
                                                   "move": "12 um in 400 us along a 2 um lane, ramps 100 us"},
                outputs={"rows": rows})
    return rows


if __name__ == "__main__":
    parts = sys.argv[1:] or ["field", "scan", "tol", "temp", "bystander", "margin", "uniformity", "slowlane", "pitch3", "aniso"]
    for p in parts:
        {"field": field, "scan": scan, "tol": tol, "temp": temp, "bystander": bystander, "margin": margin, "uniformity": uniformity, "slowlane": slowlane, "pitch3": pitch3, "aniso": aniso}[p]()
