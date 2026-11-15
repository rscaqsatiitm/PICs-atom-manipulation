"""Module 30: the planar surface-emitting tile as an alternative to the chip stack and its glued collimator.

Each engine lies flat on a carrier; its 192 emitters continue as weak gratings about 2 mm long and radiate upwards.
The 2 mm grating aperture is imaged onto the atoms' y axis with demagnification D = tile pitch / row pitch
(6 mm / 4 um = 1500, docs/scaling_beyond.md Sec. 5), so no fast-axis collimator is needed.

Parts:
  coupling        evanescent coupling of neighbouring 2 mm grating emitters (1.2 x 0.4 um SiN, 1064.4 nm): a
                  tweezer is close to one Bloch mode of the coupled array, whose propagation constant
                  beta0 + 2 kappa cos(q) depends on the x position of the tweezer (q = k0 s_x p), so the y emission
                  angle varies along the row by up to 4 kappa / k0 in sine. Supermodes vs gap (uniform array) and
                  for alternating widths (phase-mismatched array)
  directionality  2D FDTD (MEEP, run with the picsim-meep environment) of the weak grating with the Si substrate
                  as reflector (BOX sweep), a Si/SiO2 Bragg mirror under the BOX, and a metal mirror
  imaging         1D Fourier optics of the y relay: grating near field (uniform weak grating or Gaussian
                  apodization) -> NA 0.5 / D object-side pupil -> atoms; peak and crosstalk vs aperture, emission-
                  angle error and object defocus; sets the angle tolerance used in `coupling`
  power           power per 1.2 mK tweezer of the tile chain against the edge chain (m14 / m27), x pitch cost

Outputs: results/m30_coupling.json, m30_directionality.json, m30_imaging.json, m30_power.json (+ .png)
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from picsim import materials as mat
from picsim import params as prm
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m30_planar_tile.py"
LAM = prm.LAM_TRANSPORT.value
K0 = 2 * np.pi / LAM
NA = prm.NA_OBJ.value
ROW = 4.0  # um, row pitch at the atoms
TILE = 6000.0  # um, tile pitch (CHOICE, docs/scaling_beyond.md Sec. 5)
D = TILE / ROW
L_GRAT = 2000.0  # um, grating length (CHOICE)
W_EMIT, T_SIN = 1.2, 0.4  # um, D1 emitter
P_OUT = 1.5  # um, D1 emitter pitch
S_CUT_D1 = 0.355  # x direction sine at the slit edge (D1)


# ---------------------------------------------------------------------------------------------- coupling
def _pair(args):
    from picsim import modes as M
    w1, w2, gap = args
    xc1, xc2 = -(gap + w1) / 2, (gap + w2) / 2
    xs = M.XSection(cores=[(xc1, w1), (xc2, w2)], t=T_SIN, x_left=xc1 - w1 / 2 - 4, x_right=xc2 + w2 / 2 + 4)
    ms = M.solve(xs, LAM, num_modes=8)
    te = sorted([m.neff.real for m in ms if m.te_fraction > 0.5], reverse=True)
    return {"w1": w1, "w2": w2, "gap": gap, "pitch": gap + (w1 + w2) / 2, "n_top2": te[:2]}


def _single(w):
    from picsim import modes as M
    xs = M.XSection(cores=[(0.0, w)], t=T_SIN, x_left=-w / 2 - 4, x_right=w / 2 + 4)
    ms = M.solve(xs, LAM, num_modes=4)
    return max(m.neff.real for m in ms if m.te_fraction > 0.5)


def coupling():
    from picsim.parallel import pmap
    tol = json.loads((RESULTS / "m30_imaging.json").read_text())["outputs"]["angle_tolerance"]
    gaps = [0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3]
    rows = pmap(_pair, [(W_EMIT, W_EMIT, g) for g in gaps], processes=4)
    out_u = []
    for r in rows:
        dn = r["n_top2"][0] - r["n_top2"][1]
        kappa = K0 * dn / 2  # 1/um
        swing = 2 * dn  # Bloch band width in effective index (q = 0 to pi): 4 kappa / k0
        s_cut = LAM / (2 * r["pitch"])  # grating-lobe-free x half-range for this pitch
        out_u.append({**r, "dn_split": dn, "kappa_per_um": kappa, "Lc_mm": LAM / (2 * dn) * 1e-3,
                      "bloch_swing_sin": swing, "s_cut": s_cut, "x_field_fraction_of_D1": s_cut / S_CUT_D1})
        print("uniform gap %.2f pitch %.2f | dn %.2e Lc %.2f mm | swing %.2e (tol %.1e) | s_cut %.3f (%.2f of D1)" % (
            r["gap"], r["pitch"], dn, out_u[-1]["Lc_mm"], swing, tol["delta_sin_for_peak_0.97"], s_cut,
            out_u[-1]["x_field_fraction_of_D1"]), flush=True)
    # phase-mismatched (alternating widths) at the D1 pitch: band swing 4 kappa^2 / dbeta (dbeta >> kappa)
    out_a = []
    pairs = [(1.0, 1.4), (0.9, 1.5), (0.8, 1.6)]
    jobs = [(w1, w2, P_OUT - (w1 + w2) / 2) for w1, w2 in pairs] + [(w1, w2, 3.0) for w1, w2 in pairs]
    res = pmap(_pair, jobs, processes=3)
    for (w1, w2), r, far in zip(pairs, res[:3], res[3:]):
        gap = P_OUT - (w1 + w2) / 2
        split = r["n_top2"][0] - r["n_top2"][1]
        dn0 = far["n_top2"][0] - far["n_top2"][1]  # isolated index difference, same solver set-up at a 3 um gap
        n_iso = {w1: far["n_top2"][1], w2: far["n_top2"][0]}
        c = np.sqrt(max(split**2 - dn0**2, 0.0))  # = 2 kappa / k0 (two-guide coupled modes)
        kap = K0 * c / 2
        # alternating array, bands n(q) = nbar +/- sqrt((dn0/2)^2 + (c cos(q/2))^2): swing of the upper band
        sw = np.sqrt((dn0 / 2) ** 2 + c**2) - dn0 / 2
        out_a.append({"w1": w1, "w2": w2, "pitch": P_OUT, "gap": gap, "n_iso_at_3um_gap": [n_iso[w1], n_iso[w2]], "dn0": dn0,
                      "split": split, "kappa_per_um": kap, "bloch_swing_sin": float(sw),
                      "note": "two emitter types radiate at different angles unless their grating periods differ"})
        print("alternating %.1f/%.1f at %.2f | dn0 %.3e kappa %.2e /um | swing %.2e" % (w1, w2, P_OUT, dn0, kap, sw),
              flush=True)
    need = tol["delta_sin_for_peak_0.97"]
    ok = [r for r in out_u if r["bloch_swing_sin"] <= need]
    p_min = min(r["pitch"] for r in ok) if ok else None
    save_result("m30_coupling", script=SCRIPT,
                inputs={"emitter_um": [W_EMIT, T_SIN], "lam_um": LAM, "solver": "femwell (picsim.modes), two guides",
                        "tolerance": "angle error for peak >= 0.97 from m30_imaging"},
                outputs={"uniform": out_u, "alternating": out_a, "tolerance_delta_sin": need,
                         "min_pitch_uniform_um": p_min,
                         "x_field_fraction_at_min_pitch": None if p_min is None else LAM / (2 * p_min) / S_CUT_D1},
                notes="Bloch picture: a tweezer at x_a is formed by emitters with phase step q = k0 s_x p, close to one "
                      "Bloch mode of the coupled array; its effective index, and so the y emission sine, is "
                      "n0 + (dn_split/2)... the band spans 2 dn_split (uniform array) between q = 0 and pi. Mode-solver "
                      "noise is about 1e-7 in n_eff.")
    return out_u, out_a


# ---------------------------------------------------------------------------------------------- directionality (MEEP)
DBR = lambda pairs: [(LAM / (4 * 3.55), 3.55), (LAM / (4 * 1.4496), 1.4496)] * pairs  # noqa: E731  Si / SiO2


def _fdtd(args):
    sys.path.insert(0, str(RESULTS.parent / "scripts"))
    import m3_grating_fdtd as G
    etch, t_box, kind, pairs = args
    stack = DBR(pairs) if kind == "dbr" else None
    r = G.run_case(etch, t_box, stack=stack, substrate="metal" if kind == "metal" else "si")
    keep = ("etch_um", "t_box_um", "period_um", "L_um", "theta_air_deg", "alpha_per_um", "L_1e_power_um", "up_c",
            "down_c", "trans_c", "refl_c", "directionality_c")
    out = {"kind": kind, "dbr_pairs": pairs, **{k: r[k] for k in keep}}
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in out.items()}), flush=True)
    return out


def directionality():
    from multiprocessing import Pool
    cases = [(0.05, round(tb, 3), "si", 0) for tb in np.arange(2.0, 2.81, 0.05)]
    # mirrors: 0.1 um steps span one interference period (lambda / (2 n cos theta) = 0.37 um)
    cases += [(0.05, round(tb, 3), "dbr", p) for p in (2, 3) for tb in np.arange(1.8, 2.21, 0.1)]
    cases += [(0.05, round(tb, 3), "metal", 0) for tb in np.arange(1.8, 2.21, 0.1)]
    done = []
    for f in os.environ.get("FDTD_CACHE", "").split(","):  # completed cases from earlier logs (JSON lines)
        if f and os.path.exists(f):
            done += [json.loads(line) for line in open(f) if line.startswith('{"kind"')]
    key = lambda c: (c[0], round(c[1], 3), c[2], c[3])  # noqa: E731
    have = {key((r["etch_um"], r["t_box_um"], r["kind"], r["dbr_pairs"])): r for r in done}
    todo = [c for c in cases if key(c) not in have]
    with Pool(int(os.environ.get("NPROC", 4))) as pool:
        new = pool.map(_fdtd, todo, chunksize=1)
    rows = [have[key(c)] if key(c) in have else new[todo.index(c)] for c in cases]
    best = {}
    for kind in ("si", "dbr", "metal"):
        for p in sorted({r["dbr_pairs"] for r in rows if r["kind"] == kind}):
            rs = [r for r in rows if r["kind"] == kind and r["dbr_pairs"] == p]
            b = max(rs, key=lambda r: r["up_c"] / max(r["up_c"] + r["down_c"], 1e-12))
            best[f"{kind}{p if kind == 'dbr' else ''}"] = b
    with Pool(3) as pool:  # deeper etch at the best BOX of each kind (etch dependence of directionality)
        chk = pool.map(_fdtd, [(0.1, b["t_box_um"], b["kind"], b["dbr_pairs"]) for b in best.values()], chunksize=1)
    save_result("m30_directionality", script=SCRIPT,
                inputs={"tool": "MEEP 1.35 2D FDTD (scripts/m3_grating_fdtd.run_case), 40 px/um, 30 um uniform grating, "
                                "emission ~10 deg in air", "etch_um": 0.05, "dbr": "Si (n 3.55) / SiO2 quarter-wave pairs "
                                "under the BOX", "metal": "perfect electric conductor under the BOX"},
                outputs={"rows": rows, "best": best, "etch_0.1_at_best": chk},
                notes="directionality = up / (up + down) of the radiated power; for a long grating the transmitted "
                      "remainder vanishes and the up-fraction of guided power approaches it (minus back-reflection).")
    return rows


# ---------------------------------------------------------------------------------------------- imaging
def _line(A, u, dtheta=0.0, dz=0.0, n_eta=1201, y=np.arange(-12.0, 12.0001, 0.01)):
    """Atom-plane y line of a near field A(u) (um) imaged with demagnification D through the NA pupil. Normalised
    so that its integral is the fraction of the emitted power that passes the pupil."""
    na_o = NA / D
    eta = np.linspace(-1, 1, n_eta)
    s = eta * na_o
    du = u[1] - u[0]
    ph = np.exp(1j * K0 * dtheta * u)
    F = np.exp(-1j * K0 * np.outer(s, u)) @ (A * ph) * du  # angular spectrum
    F = F * np.exp(-1j * K0 * dz * s**2 / 2)  # object-side defocus (paraxial)
    P_tot = np.sum(np.abs(A) ** 2) * du
    P_in = np.sum(np.abs(F) ** 2) * (s[1] - s[0]) * K0 / (2 * np.pi)
    E = np.exp(1j * K0 * NA * np.outer(y, eta)) @ F * (eta[1] - eta[0])
    I = np.abs(E) ** 2
    I = I / (np.sum(I) * (y[1] - y[0])) * (P_in / P_tot)
    return y, I, P_in / P_tot


def _eps(y, I, d, half=0.5):
    m = np.abs(np.abs(y) - d) <= half
    return float(I[m].max() / I.max())


def imaging():
    u = np.linspace(-3500.0, 3500.0, 8001)
    profiles, length = {}, {}
    for L in (2000.0, 3000.0, 4000.0, 5000.0, 6000.0):
        for w in (600.0, 800.0, 1000.0, 1200.0, 1450.0, 1700.0):  # Gaussian apodization (1/e^2 radius), truncated at +/-L/2
            if w > 0.75 * L:
                continue
            k = f"L {L:g} um, gauss w {w:g} um"
            profiles[k] = np.where(np.abs(u) <= L / 2, np.exp(-(u / w) ** 2), 0.0)
            length[k] = L
    for Le in (500.0, 1000.0):  # uniform weak grating, 1/e power length Le, 2 mm long
        k = f"L 2000 um, uniform Le {Le:g} um"
        profiles[k] = np.where(np.abs(u) <= 1000.0, np.exp(-(u + 1000.0) / (2 * Le)), 0.0)
        length[k] = 2000.0
    base = {}
    for name, A in profiles.items():
        y, I, tr = _line(A, u)
        base[name] = {"grating_length_um": length[name], "peak": float(I.max()), "transmission": tr,
                      "eps3": _eps(y, I, 3.0), "eps4": _eps(y, I, 4.0), "fwhm_um": float(np.ptp(y[I > I.max() / 2]))}
        print(name, {k: round(v, 5) for k, v in base[name].items()}, flush=True)
    ok = {k: v for k, v in base.items() if v["eps4"] <= 1.2e-4}
    # design: shortest grating meeting the crosstalk of the edge emitter, then the highest peak
    design = (min(ok, key=lambda k: (ok[k]["grating_length_um"], -ok[k]["peak"])) if ok
              else min(base, key=lambda k: base[k]["eps4"]))
    A = profiles[design]
    p0 = base[design]["peak"]
    dths = np.array([0, 5, 10, 20, 30, 40, 60, 80, 100, 150]) * 1e-6
    tilt = []
    for dth in dths:
        y, I, tr = _line(A, u, dtheta=dth)
        cen = float(np.sum(y[I > I.max() / 2] * I[I > I.max() / 2]) / np.sum(I[I > I.max() / 2]))
        tilt.append({"dtheta_urad": dth * 1e6, "peak_rel": float(I.max() / p0), "eps4": _eps(y, I, 4.0),
                     "centroid_nm": cen * 1e3, "transmission": tr})
    foc = []
    for dz in (0.0, 100.0, 200.0, 400.0, 800.0, 1600.0):
        y, I, tr = _line(A, u, dz=dz)
        foc.append({"dz_obj_um": dz, "dz_atoms_um": dz / D**2, "peak_rel": float(I.max() / p0), "eps4": _eps(y, I, 4.0)})
    t97 = float(np.interp(0.97, [r["peak_rel"] for r in tilt][::-1], dths[::-1]))
    t99 = float(np.interp(0.99, [r["peak_rel"] for r in tilt][::-1], dths[::-1]))
    e = np.array([r["eps4"] for r in tilt])
    t_xt = float(np.interp(2e-4, e, dths)) if e.max() > 2e-4 else None  # eps4 grows monotonically with tilt
    tol = {"delta_sin_for_peak_0.97": t97, "delta_sin_for_peak_0.99": t99, "delta_sin_for_eps4_2e-4": t_xt,
           "estimate_scaling_beyond_urad": 30.0}
    for r in tilt:
        print("tilt %5.0f urad | peak %.4f | eps4 %.2e | centroid %.1f nm" % (r["dtheta_urad"], r["peak_rel"], r["eps4"],
                                                                           r["centroid_nm"]), flush=True)
    print("design", design, "tolerance", tol, flush=True)
    save_result("m30_imaging", script=SCRIPT,
                inputs={"D": D, "NA": NA, "object_NA": NA / D, "grating_lengths_um": [2000, 3000, 4000, 5000, 6000], "lam_um": LAM,
                        "model": "1D scalar Fourier optics, object-side pupil |s| <= NA / D, paraxial defocus"},
                outputs={"profiles": base, "design": design, "tilt": tilt, "focus": foc, "angle_tolerance": tol},
                notes="peak: line intensity (1/um) per unit emitted power, same normalisation as RelayY.atom_line. "
                      "Design: shortest grating with crosstalk at 4 um <= 1.2e-4 (the edge emitter gives 8e-5 to "
                      "1.1e-4), then the highest peak. The pupil-fill optimum needs a near field of 1/e^2 radius about "
                      "D x 0.97 um.")
    return tol


# ---------------------------------------------------------------------------------------------- power
def power():
    im = json.loads((RESULTS / "m30_imaging.json").read_text())["outputs"]
    di = json.loads((RESULTS / "m30_directionality.json").read_text())["outputs"]
    cp = json.loads((RESULTS / "m30_coupling.json").read_text())["outputs"]
    wv = json.loads((RESULTS / "m14_wave.json").read_text())["outputs"]
    fgl = json.loads((RESULTS / "m14_facet_meep_glass.json").read_text())["outputs"]
    m27 = json.loads((RESULTS / "m27_wave.json").read_text())["outputs"]
    corr = fgl["T_within_s0.9"] / fgl["kirchhoff"]["T_within_s0.9"]
    edge = wv["weierstrass_aplanat|AR|beta0.7"]["peak"] * wv["T_facet"]["SiO2"] * corr * \
        m27["mirror_pair_reflectance"]["adopted: M1 silver, M2 TIR"]  # line peak per unit guided power, edge chain
    p_edge = 57.7  # mW per 1.2 mK tweezer (m27, reflective collimator)
    tile_peak = im["profiles"][im["design"]]["peak"]
    rows = {}
    for k, b in di["best"].items():
        up = b["up_c"] / max(b["up_c"] + b["down_c"], 1e-12)  # long grating: all guided power radiates
        eff = tile_peak * up * (1 - b["refl_c"])
        rows[k] = {"t_box_um": b["t_box_um"], "directionality": up, "back_reflection": b["refl_c"],
                   "line_peak_per_guided": eff, "P_tweezer_mW": p_edge * edge / eff}
        print(k, rows[k], flush=True)
    xf = cp["x_field_fraction_at_min_pitch"]
    out = {"edge_line_peak_per_guided": edge, "tile_line_peak_per_emitted": tile_peak, "by_reflector": rows,
           "x_field_fraction_at_min_pitch_uniform": xf, "min_pitch_um": cp["min_pitch_uniform_um"],
           "y_relay_object_field_mm_16_tiles": 16 * TILE * 1e-3, "object_NA": NA / D,
           "lagrange_check": {"object_half_field_x_NA": 8 * TILE * NA / D, "atom_half_rows_x_NA": 8 * ROW * NA}}
    print(json.dumps({k: v for k, v in out.items() if k != "by_reflector"}), flush=True)
    save_result("m30_power", script=SCRIPT,
                inputs={"edge": "m14 wave (Weierstrass + aplanat, beta 0.7) x T_facet(SiO2) x FDTD/Kirchhoff x m27 "
                                "mirror pair (silver M1, TIR M2): 57.7 mW per 1.2 mK tweezer",
                        "tile": "m30 imaging design x FDTD directionality x (1 - back-reflection)",
                        "x_axis": "unchanged (same engine and x relay); the emitter pitch cost is in m30_coupling"},
                outputs=out,
                notes="P_tweezer scales the edge-chain power by the ratio of line peaks per unit guided power (y axis "
                      "only). If the emitter pitch must grow to suppress coupling, the x field per engine shrinks "
                      "by x_field_fraction (or the port count must grow by its inverse).")
    return out


if __name__ == "__main__":
    for part in sys.argv[1:] or ["imaging", "coupling", "power"]:
        globals()[part]()
