"""Module 27: an all-reflective fast-axis collimator that is insensitive to the x angle of the light.

The refractive collimator of module 14 detunes for tweezers away from the row centre (docs/megascale.md Sec. 12):
light of a tweezer at x_a has x direction sine s_x = x_a / 288 um, a cylindrical surface conserves k_x, and every
index acts as sqrt(n^2 - s_x^2), so the facet defocuses by -99 s_x^2 um and one engine serves |x_a| <= 45 um only.
Two cylindrical mirrors in one glass block (picsim/opa/mirror_collimator.py) collimate exactly at every s_x; with the
sine condition built in, only the focal length changes, by n / sqrt(n^2 - s_x^2).

Parts:
  scan       grid of constructions (r0, y_c, phi0, sign): no obscuration, light inside z >= 0, blocks of neighbouring
             chips (0.5 mm pitch, 10 um margin each) disjoint, radius of curvature; ranked by form-error sensitivity
  trace      independent 3D vector trace of the adopted design at s_x = 0 ... 0.355: wavefront, mapping, aperture
  wave       wave model at s_x = 0 against the adopted refractive collimator (m14); power per tweezer
  skew       m26 skew analysis with the reflective design: peak, crosstalk, y centroid vs x_a, 2D, comb
  tolerance  facet decentre and focus (first order and exact trace), mirror form errors (Monte Carlo)
  pol        Jones trace of the TE field through both mirrors: cross-polarised power and circular polarisation vs x_a
  xaxis      x slit-plane distance, comb start, glass-slab phase in x applied to D1 before and after recalibration
  figure     layout with neighbouring blocks; peak vs x_a, refractive and reflective

Outputs: results/m27_*.json (+ .png)
"""
from __future__ import annotations

import itertools
import json
import sys

import numpy as np

from picsim import materials as mat
from picsim import params as prm
from picsim.opa import mirror_collimator as MC
from picsim.parallel import pmap
from picsim.plotting import fig
from picsim.opa import collimator as C
from picsim.opa.stack import FacetY, Focus2D, RelayY
from picsim.provenance import RESULTS, save_result

sys.path.insert(0, str(RESULTS.parent / "scripts"))
import m14_collimator as m14  # noqa: E402
import m26_stack_combine as m26  # noqa: E402

SCRIPT = "scripts/m27_reflective_collimator.py"
LAM = prm.LAM_TRANSPORT.value
N_SIO2 = float(mat.n_sio2(LAM))
F_FAC = 175.0
NA_FAC = 0.9
PITCH = 500.0  # um, chip pitch (m13 CHOICE)
MARGIN = 10.0  # um, clearance of each block around its light (CHOICE): 20 um between neighbouring blocks
R_MIN = 20.0  # um, smallest radius of curvature accepted (CHOICE, fabrication)
BETA = 0.7
NA = prm.NA_OBJ.value
# metal optical constants at 1064.4 nm: Johnson & Christy, Phys. Rev. B 6, 4370 (1972), doi:10.1103/PhysRevB.6.4370,
# tables retrieved from refractiveindex.info (database/data/main/{Ag,Au}/nk/Johnson.yml) on 2026-11-09; linear
# interpolation between 0.9840 and 1.0880 um (Ag 0.04 + 6.992i, 0.04 + 7.795i; Au 0.22 + 6.350i, 0.27 + 7.150i)
N_AG = complex(0.04, 6.992 + (1.0644 - 0.984) / (1.088 - 0.984) * (7.795 - 6.992))
N_AU = complex(0.22 + (1.0644 - 0.984) / (1.088 - 0.984) * 0.05, 6.350 + (1.0644 - 0.984) / (1.088 - 0.984) * 0.80)
XA = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)


def _one(p):
    r0, y_c, phi0, sign = p
    try:
        d = MC.reflective_aplanat(N_SIO2, F_FAC, r0, y_c, phi0, sign)
        ok, g = d.ok(PITCH, MARGIN)
    except Exception as e:  # noqa: BLE001  construction or trace failure
        return {"params": p, "valid": False, "error": str(e)[:80]}
    g = {k: v for k, v in g.items()}
    valid = bool(ok and min(g["min_radius_of_curvature_um"]) >= R_MIN and g["exit_tilt_max_rad"] < 1e-6)
    return {"params": p, "valid": valid, **g}


SELECTION = ("adopted: block length <= 500 um (CHOICE, comparable to the 0.40 mm refractive collimator); among those, "
             "form sensitivity within 1% of the lowest; then the smallest y extent of the light; then the largest "
             "smallest radius of curvature. TIR variant: both mirrors used beyond the critical angle over the whole "
             "aperture (no coating); same tie-breaks.")


def _pick(cands):
    if not cands:
        return None
    smin = min(r["form_sensitivity_rms"] for r in cands)
    near = [r for r in cands if r["form_sensitivity_rms"] <= 1.01 * smin]
    near.sort(key=lambda r: (round(np.ptp(r["region_y_extent_um"])), -min(r["min_radius_of_curvature_um"])))
    return near[0]


def select(valid):
    return {"adopted": _pick([r for r in valid if r["block_length_um"] <= 500.0]),
            "tir": _pick([r for r in valid if r["TIR_fraction_M1"] == 1.0 and r["TIR_fraction_M2"] == 1.0])}


def _design(which="adopted", R=1.0):
    p = json.loads((RESULTS / "m27_scan.json").read_text())["outputs"]["adopted" if which == "adopted" else "tir_variant"]
    r0, y_c, phi0, sign = p["params"]
    return MC.reflective_aplanat(N_SIO2, F_FAC, r0, y_c, phi0, int(sign), R=R)


def scan():
    grid = list(itertools.product((40.0, 50.0, 60.0, 75.0, 90.0, 105.0, 120.0, 140.0, 160.0),
                                  (-160.0, -180.0, -200.0, -220.0, -240.0, -260.0, -280.0, -300.0, -330.0, -360.0,
                                   -390.0, -420.0, -450.0, -480.0, -510.0),
                                  (-30.0, -40.0, -50.0, -60.0, -65.0, -70.0, -75.0, -80.0, -90.0, -100.0, -110.0,
                                   -120.0, -135.0, -150.0),
                                  (-1, 1)))
    rows = pmap(_one, grid, processes=4)
    valid = [r for r in rows if r["valid"]]
    valid.sort(key=lambda r: (round(r["form_sensitivity_rms"], 2), r["block_length_um"]))
    print(len(rows), "designs,", len(valid), "valid", flush=True)
    for r in valid[:25]:
        print(r["params"], "sens %.3f len %.0f xs %.0f Rc %s inc1 %s inc2 %s" % (
            r["form_sensitivity_rms"], r["block_length_um"], r["max_cross_section_um"],
            [round(x) for x in r["min_radius_of_curvature_um"]], [round(x) for x in r["incidence_M1_deg"]],
            [round(x) for x in r["incidence_M2_deg"]]), flush=True)
    sel = select(valid)
    print("adopted", sel["adopted"]["params"], "TIR variant", sel["tir"]["params"] if sel["tir"] else None, flush=True)
    reasons = {"construction": sum(1 for r in rows if "error" in r),
               "obscured": sum(1 for r in rows if "error" not in r and r["obscured_length_um"] >= 0.5),
               "pitch": sum(1 for r in rows if "error" not in r and r["pitch_overlap_area_um2"] >= 1.0),
               "behind_facet": sum(1 for r in rows if "error" not in r and r["light_behind_facet_um2"] >= 1.0),
               "curvature": sum(1 for r in rows if "error" not in r and min(r["min_radius_of_curvature_um"]) < R_MIN)}
    save_result("m27_scan", script=SCRIPT,
                inputs={"n": N_SIO2, "f_um": F_FAC, "NA_FAC": NA_FAC, "pitch_um": PITCH, "margin_um": MARGIN,
                        "R_min_um": R_MIN, "grid": "r0 x y_c x phi0 x sign (see script)"},
                outputs={"n_designs": len(rows), "n_valid": len(valid), "rejected": reasons,
                         "adopted": sel["adopted"], "tir_variant": sel["tir"], "selection": SELECTION, "valid": valid},
                notes="form_sensitivity_rms: rms over the pupil of sqrt(sum_j (2 n cos i_j)^2), optical path per unit "
                      "normal sag error of the two mirrors. Ranking: sensitivity, then block length.")
    return valid


# ---------------------------------------------------------------------------------------------- trace
def _sxs():
    fp, _, _, _ = m26.setup()
    return fp.opa.x_scale, fp.opa.s_cut


def trace():
    """Independent 3D vector trace (k_x conserved by construction of the rays, not assumed) at several s_x."""
    x_scale, s_cut = _sxs()
    out = {"x_scale_um": x_scale, "s_cut": s_cut}
    ref = m26.fac_adopted()
    for which in ("adopted", "tir"):
        d = _design(which)
        rows = []
        for sx in (0.0, 0.1, 0.2, 0.312, s_cut):
            ne = float(np.sqrt(d.n**2 - sx**2))
            s_rim = ne * np.sin(d.theta_m)
            sy = np.linspace(-s_rim, s_rim, 801) * 0.9995
            r = d.trace(sx, sy)
            ok = np.isfinite(r["y"])
            W = r["opl"] - sx * r["x"]  # phase of the plane-wave component k_x = k0 sx, relative to the axial ray
            W = (W - W[len(W) // 2]) / LAM
            g = d.sign * (r["y"] - d.y_c) / d.f
            target = sy * d.n / ne
            big = np.abs(sy) > 0.05
            Le = MC.skewed(d, sx).lensmap(LAM, NA_FAC)
            cm = Le.s <= sy.max()  # common range of the two grids
            W3 = np.interp(Le.s[cm], sy, W)
            g3 = np.interp(Le.s[cm], sy, g)
            Lr = m26.skewed(ref, sx).lensmap(LAM, NA_FAC, coated=True)
            rows.append({"sx": sx, "x_a_um": sx * x_scale, "n_eff": ne, "s_rim": float(s_rim), "rays_ok": int(ok.sum()),
                         "W_pv_waves": float(np.nanmax(W) - np.nanmin(W)),
                         "g_over_target_max_dev": float(np.nanmax(np.abs(g[big] / target[big] - 1))),
                         "f_scale": d.n / ne, "antisymmetry_um": float(np.nanmax(np.abs(g + g[::-1])) * d.f),
                         "exit_tilt_max_rad": float(np.nanmax(np.abs(r["d"][1] / r["d"][2]))),
                         "stripe_edge_um": float(d.f * np.nanmax(np.abs(g))),
                         "effective_index_vs_3d": {"W_waves": float(np.max(np.abs(W3 - Le.W[cm]))),
                                                   "g": float(np.max(np.abs(g3 - Le.g[cm])))},
                         "refractive_W_pv_waves": float(np.ptp(Lr.W))})
            print(which, {k: (float("%.3g" % v) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
        out[which] = {"params": d.construction, "rows": rows}
    save_result("m27_trace", script=SCRIPT,
                inputs={"method": "3D rays from the facet in fused silica, direction (sx, sy, .)/n, vector law of "
                        "reflection on the spline-fitted mirror curves (cylinders along x), optical path n x 3D length to "
                        "the exit plane; W = OPL - sx x_exit (phase of the k_x = k0 sx component)",
                        "refractive": "m14 adopted collimator in the effective-index trace (m26.skewed)"},
                outputs=out,
                notes="g_over_target: exit height / (f sy n / n_eff); f_scale = n / n_eff. effective_index_vs_3d: largest "
                      "difference between the 2D effective-index lens map (used by the wave model) and the 3D trace.")
    return out


# ---------------------------------------------------------------------------------------------- wave
def _facet_corr():
    fgl = json.loads((RESULTS / "m14_facet_meep_glass.json").read_text())["outputs"]
    fair = json.loads((RESULTS / "m14_facet_meep_air.json").read_text())["outputs"]
    return {"SiO2": fgl["T_within_s0.9"] / fgl["kirchhoff"]["T_within_s0.9"],
            "air": fair["T_within_s0.9"] / fair["kirchhoff"]["T_within_s0.9"]}


def _metal(nm, n=N_SIO2):
    return lambda c: MC.fresnel_metal(n, nm, c)


def _tir(n=N_SIO2):
    return lambda c: MC.fresnel_metal(n, 1.0 + 0j, c)


def _mirror_R(d, coats, sx=0.0):
    """Power reflectance of the pair, weighted by the facet far field (TE launch)."""
    ne = float(np.sqrt(d.n**2 - sx**2))
    fa = FacetY(LAM, 1.2, 0.4, n_m=ne)
    sy = np.linspace(-1, 1, 401) * ne * np.sin(d.theta_m) * 0.999
    r = d.trace(sx, sy, pol=True, coat=coats)
    a2 = np.interp(sy, fa.s, np.abs(fa.A) ** 2 * fa.cos)
    P = np.sum(np.abs(r["E"]) ** 2, axis=0)
    return float(np.sum(a2 * P) / np.sum(a2))


def wave():
    y = np.arange(-12, 12, 0.005)
    fa_air, fa_im = FacetY(LAM, 1.2, 0.4), FacetY(LAM, 1.2, 0.4, n_m=N_SIO2)
    ref = m14.line_metrics(RelayY(fa_air, BETA, NA=NA, NA_FAC=NA_FAC), y)["peak"] * fa_air.T_facet
    designs = {"reflective (adopted)": _design(), "reflective, TIR variant": _design("tir"),
               "refractive (m14 adopted)": C.fac_design(N_SIO2, n2=1.8), "ideal lens, immersed": None}
    out = {"T_facet": {"air": fa_air.T_facet, "SiO2": fa_im.T_facet}}
    for name, d in designs.items():
        L = None if d is None else d.lensmap(LAM, NA_FAC, coated=True)
        for beta in (0.6, 0.7, 0.8):
            m = m14.line_metrics(RelayY(fa_im, beta, NA=NA, NA_FAC=NA_FAC, lens=L), y)
            m["peak_rel_ideal"] = m["peak"] * fa_im.T_facet / ref
            out[f"{name}|beta{beta}"] = m
            print(name, beta, {k: round(v, 5) for k, v in m.items()}, flush=True)
    # power per 1.2 mK tweezer (method of m14 plan), with the mirror reflectance of real coatings
    m13 = json.loads((RESULTS / "m13_power.json").read_text())["outputs"]["P_new_mW"]
    corr = _facet_corr()
    T_air = fa_air.T_facet
    ad, tv = designs["reflective (adopted)"], designs["reflective, TIR variant"]
    Rs = {"adopted: M1 silver, M2 TIR": _mirror_R(ad, [_metal(N_AG), _tir()]),
          "adopted: M1 gold, M2 TIR": _mirror_R(ad, [_metal(N_AU), _tir()]),
          "adopted: both silver": _mirror_R(ad, [_metal(N_AG), _metal(N_AG)]),
          "TIR variant: both TIR": _mirror_R(tv, [_tir(), _tir()])}
    power = {}
    for name in designs:
        rel = out[f"{name}|beta0.7"]["peak_rel_ideal"] * T_air * corr["SiO2"]
        power[name] = m13 / rel
    power_coated = {k: power["reflective (adopted)" if k.startswith("adopted") else "reflective, TIR variant"] / R
                    for k, R in Rs.items()}
    out.update({"mirror_pair_reflectance": Rs, "P_tweezer_mW_ideal_coatings": power, "P_tweezer_mW_real_mirrors": power_coated,
                "n_Ag": [N_AG.real, N_AG.imag], "n_Au": [N_AU.real, N_AU.imag]})
    print(json.dumps({k: out[k] for k in ("mirror_pair_reflectance", "P_tweezer_mW_ideal_coatings",
                                          "P_tweezer_mW_real_mirrors")}, indent=1), flush=True)
    save_result("m27_wave", script=SCRIPT,
                inputs={"beta_um": "0.6-0.8", "NA_FAC": NA_FAC, "facet": "D1 emitter radiating into SiO2 (glued block)",
                        "coatings": "AR exit face (t = 1); mirrors ideal in the line metrics, real in the power",
                        "power": "m13 power (47 mW) scaled as in m14 plan: peak_rel_ideal x T_facet(air) x FDTD/Kirchhoff "
                                 "(SiO2), divided by the mirror-pair reflectance",
                        "metals": "Johnson & Christy 1972 via refractiveindex.info, interpolated to 1.0644 um (SOURCE)"},
                outputs=out, notes="1D y line, s_x = 0. peak_rel_ideal: per unit guided power, relative to the ideal "
                                   "lens in air at beta 0.7 (as m14_wave).")
    return out


# ---------------------------------------------------------------------------------------------- skew
def _contiguous(rows, key, ok):
    last = None
    for r in sorted(rows, key=lambda r: r["x_a_um"]):
        if not ok(r[key]):
            break
        last = r["x_a_um"]
    return last


D1_REACH = 67.0  # um: the 2D model's engine (D1, 128 ports) cannot form a tweezer beyond about 67 um (m3/m9 range)


def skew():
    x_scale, s_cut = _sxs()
    m26.skew(design=_design(), name="m27_skew", xas=XA + (s_cut * x_scale,))
    return skew_summary()


def skew_summary():
    """Usable ranges from results/m27_skew.json and m26_skew.json; 2D values only within the D1 reach."""
    x_scale, s_cut = _sxs()
    d = _design()
    rows = json.loads((RESULTS / "m27_skew.json").read_text())["outputs"]["rows"]
    ref = json.loads((RESULTS / "m26_skew.json").read_text())["outputs"]["rows"]
    ya = np.arange(-12.0, 12.0001, 0.01)
    comp = {}
    for name, dd in (("reflective", d), ("refractive", m26.fac_adopted())):
        pk, cen = [], []
        for xa in XA + (s_cut * x_scale,):
            cy = m26.CombY(design=dd, sx=xa / x_scale)
            I = cy.r.atom_line(ya)
            pk.append(float(I.max()))
            w = I > 0.5 * I.max()
            cen.append(float(np.sum(ya[w] * I[w]) / np.sum(I[w]) * 1e3))
        comp[name] = {"x_a_um": list(XA + (s_cut * x_scale,)), "peak_rel_own_centre": [p / pk[0] for p in pk],
                      "centroid_shift_nm": [c - cen[0] for c in cen]}
        print(name, comp[name], flush=True)
    summary = {"contiguous_x_a_um": {
        "reflective peak >= 0.97": _contiguous(rows, "peak_rel_ideal", lambda v: v >= 0.97),
        "reflective 2D peak >= 0.97 (within D1 reach)": _contiguous([r for r in rows if r["x_a_um"] <= D1_REACH],
                                                                     "twod_peak_rel_ideal", lambda v: v >= 0.97),
        "reflective xt4 <= 2e-4": _contiguous(rows, "xt_4um", lambda v: v <= 2e-4),
        "refractive peak >= 0.97": _contiguous(ref, "peak_rel_ideal", lambda v: v >= 0.97),
        "refractive xt4 <= 2e-4": _contiguous(ref, "xt_4um", lambda v: v <= 2e-4)},
        "own_centre": comp, "twod_valid_up_to_x_a_um": D1_REACH,
        "note_2d": "2D rows with x_a > 67 um are outside the steering range of the D1 engine of the 2D model (its "
                   "synthesis cannot place a tweezer there, for any collimator); they are not collimator results"}
    print(summary["contiguous_x_a_um"], flush=True)
    save_result("m27_skew_summary", script=SCRIPT,
                inputs={"rows": "results/m27_skew.json (reflective) and results/m26_skew.json (refractive)",
                        "own_centre": "1D line peak and centroid relative to the same collimator at x_a = 0 (facet "
                                      "spectrum in the effective source index)"},
                outputs=summary)
    return summary


# ---------------------------------------------------------------------------------------------- tolerance
def exact_offaxis(d, dy0=0.0, dz0=0.0, n=601):
    s = np.linspace(-NA_FAC, NA_FAC, n) * 0.999
    r0 = d.trace(0.0, s)
    r1 = d.trace(0.0, s, y0=dy0, z0=dz0)
    ok = np.isfinite(r0["y"]) & np.isfinite(r1["y"])
    o = np.argsort(r0["y"][ok])
    y0s, o0s, ss = r0["y"][ok][o], r0["opl"][ok][o], s[ok][o]
    m = ok & np.isfinite(r1["y"])
    o0 = np.interp(r1["y"][m], y0s, o0s)
    s_at = np.interp(r1["y"][m], y0s, ss)
    exact = r1["opl"][m] - o0
    model = -s_at * dy0 - d.n * np.sqrt(1 - (s_at / d.n) ** 2) * dz0
    res = (exact - model) - np.mean(exact - model)
    return float(np.sqrt(np.mean(res**2)) / LAM), float(np.ptp(res) / LAM)


def _arc(M, u):
    ug = np.linspace(M.u[0], M.u[-1], 4001)
    sp = np.hypot(M.cy(ug, 1), M.cz(ug, 1))
    a = np.concatenate([[0.0], np.cumsum(0.5 * (sp[1:] + sp[:-1]) * np.diff(ug))])
    return np.interp(u, ug, a), a[-1]


def tolerance():
    y = np.arange(-12, 12, 0.005)
    fa = FacetY(LAM, 1.2, 0.4, n_m=N_SIO2)
    out = {}
    rng = np.random.default_rng(11)
    for which in ("adopted", "tir"):
        d = _design(which)
        L = d.lensmap(LAM, NA_FAC, coated=True)
        r = RelayY(fa, BETA, NA=NA, NA_FAC=NA_FAC, lens=L)
        nom = m14.line_metrics(r, y)
        rows = {"nominal": nom}
        for dl in (0.05, 0.1, 0.2, 0.3):
            m = m14.line_metrics(r, y, delta=dl)
            rows[f"decentre {dl}"] = {"peak_rel": m["peak"] / nom["peak"], "eps4": m["eps4"]}
        for dz in (0.25, 0.5, 1.0):
            m = m14.line_metrics(r, y, dz=dz)
            rows[f"focus {dz}"] = {"peak_rel": m["peak"] / nom["peak"], "eps4": m["eps4"]}
        rows["exact_check_waves_rms_pv"] = {"decentre 0.2": exact_offaxis(d, dy0=0.2),
                                            "decentre 0.5": exact_offaxis(d, dy0=0.5),
                                            "focus 0.5": exact_offaxis(d, dz0=0.5)}
        # form errors: normal displacement delta_j of mirror j changes the path by 2 n cos(i_j) delta_j
        sg = np.linspace(-NA_FAC, NA_FAC, 1201) * 0.999
        tr = d.trace(0.0, sg)
        arcs = [_arc(M, u) for M, u in zip((d.M1, d.M2), tr["u"])]
        form = []
        for corr in (20.0, 60.0):
            for sigma in (0.005, 0.01, 0.02, 0.05):
                pk, e4, wr = [], [], []
                for _ in range(24):
                    W = np.zeros_like(sg)
                    for (a, a_len), ci in zip(arcs, tr["cos_i"]):
                        grid = np.linspace(0, a_len, 2001)
                        prof = m14.random_sag(rng, grid, sigma, corr)
                        W += 2 * d.n * ci * np.interp(a, grid, prof) / LAM
                    Ws = (lambda s, sg=sg, W=W: np.interp(s, sg, W))
                    m = m14.line_metrics(r, y, Ws=Ws)
                    pk.append(m["peak"] / nom["peak"])
                    e4.append(m["eps4"])
                    wr.append(float(np.std(W)))
                form.append({"corr_um": corr, "sigma_nm": sigma * 1e3, "peak_rel_mean": float(np.mean(pk)),
                             "peak_rel_p10": float(np.percentile(pk, 10)), "eps4_p90": float(np.percentile(e4, 90)),
                             "W_rms_waves_mean": float(np.mean(wr))})
                print(which, form[-1], flush=True)
        rows["form"] = form
        rows["mirror_arc_lengths_um"] = [a[1] for a in arcs]
        out[which] = rows
        print(which, json.dumps({k: v for k, v in rows.items() if k != "form"}, default=float), flush=True)
    ref = json.loads((RESULTS / "m14_form.json").read_text())["outputs"]["weierstrass_aplanat"]["rows"]
    out["refractive_form_m14"] = ref
    save_result("m27_tolerance", script=SCRIPT,
                inputs={"beta_um": BETA, "decentre": "facet in y, um", "focus": "facet along z in the glass, um",
                        "form": "Gaussian random normal displacement along each mirror (arc length), Gaussian "
                                "correlation, independent per mirror, 24 samples; path change 2 n cos(i) delta",
                        "refractive": "m14_form (Weierstrass + aplanat), same statistics"},
                outputs=out,
                notes="exact_check: residual between exactly traced path change of a displaced source and the "
                      "first-order model of the wave model (waves rms, peak-to-valley), |s| <= 0.9.")
    return out


# ---------------------------------------------------------------------------------------------- polarisation
_FA = {}


def _facet(ne):
    key = round(ne, 12)
    if key not in _FA:
        _FA[key] = FacetY(LAM, 1.2, 0.4, n_m=ne)
    return _FA[key]


def _pol_rows(d, coats, x_scale, xas, n_rays=601):
    rows = []
    for xa in xas:
        sx = xa / x_scale
        ne = float(np.sqrt(d.n**2 - sx**2))
        fa = _facet(ne)
        sy = np.linspace(-1, 1, n_rays) * ne * np.sin(d.theta_m) * 0.999
        r = d.trace(sx, sy, pol=True, coat=coats)
        a = np.interp(sy, fa.s, fa.A.real) + 1j * np.interp(sy, fa.s, fa.A.imag)
        a = a * np.sqrt(np.interp(sy, fa.s, fa.cos))
        Ex, Ey = r["E"][0], r["E"][1]
        w = np.abs(a) ** 2
        inc = float(np.sum(w * np.abs(Ey) ** 2) / np.sum(w * (np.abs(Ex) ** 2 + np.abs(Ey) ** 2)))
        cx, cyy = np.sum(a * Ex), np.sum(a * Ey)
        tot = abs(cx) ** 2 + abs(cyy) ** 2
        rows.append({"x_a_um": xa, "cross_pol_power_pupil": inc, "cross_pol_focus": float(abs(cyy) ** 2 / tot),
                     "V_focus": float(2 * np.imag(np.conj(cx) * cyy) / tot),
                     "V_pupil_rms": float(np.sqrt(np.sum(w * (2 * np.imag(np.conj(Ex) * Ey)) ** 2) /
                                                  np.sum(w * (np.abs(Ex) ** 2 + np.abs(Ey) ** 2) ** 2)))})
    return rows


def pol():
    x_scale, s_cut = _sxs()
    xas = XA + (s_cut * x_scale,)
    ad, tv = _design(), _design("tir")
    cases = {"adopted, perfect conductors": (ad, [None, None]),
             "adopted, M1 silver, M2 TIR": (ad, [_metal(N_AG), _tir()]),
             "adopted, M1 gold, M2 TIR": (ad, [_metal(N_AU), _tir()]),
             "adopted, both silver": (ad, [_metal(N_AG), _metal(N_AG)]),
             "TIR variant, both TIR": (tv, [_tir(), _tir()])}
    for dlt in (2.0, 5.0, 10.0, 20.0):  # lossless coating whose s-p phase departs from a perfect mirror by dlt
        ret = (lambda c, dlt=dlt: (-np.ones_like(c) + 0j, np.exp(1j * np.radians(dlt)) * np.ones_like(c)))
        cases[f"adopted, both mirrors excess retardance {dlt:g} deg"] = (ad, [ret, ret])
    out = {}
    for name, (d, coats) in cases.items():
        out[name] = _pol_rows(d, coats, x_scale, xas)
        print(name, [(r["x_a_um"], "%.1e %.1e %+.1e" % (r["cross_pol_focus"], r["cross_pol_power_pupil"], r["V_focus"]))
                     for r in out[name]], flush=True)
    save_result("m27_pol", script=SCRIPT,
                inputs={"launch": "TE: x unit vector projected transverse to each ray", "weights": "facet far field "
                        "amplitude (FacetY, effective source index) times sqrt(cos)", "metals": "Johnson & Christy 1972 "
                        "via refractiveindex.info, interpolated to 1.0644 um (SOURCE)", "TIR": "fused silica to air",
                        "excess retardance": "lossless coating r_s = -1, r_p = exp(i delta) (perfect mirror: delta = 0), "
                                             "same delta at every angle of incidence"},
                outputs=out,
                notes="cross_pol_focus: |E_y|^2 / |E|^2 of the coherent pupil sum (field at the centre of the tweezer, "
                      "before the objective); V_focus: its degree of circular polarisation; cross_pol_power_pupil: "
                      "incoherent pupil average. The objective's own vector effects are not included.")
    return out


def _poltrade_one(p):
    r0, y_c, phi0, sign = p
    d = MC.reflective_aplanat(N_SIO2, F_FAC, r0, y_c, phi0, int(sign))
    out = {}
    for name, coats in (("both silver", [_metal(N_AG), _metal(N_AG)]), ("both perfect", [None, None])):
        rr = _pol_rows(d, coats, 288.0, (90.0,), n_rays=241)[0]
        out[name] = {"V_focus": rr["V_focus"], "cross_pol_focus": rr["cross_pol_focus"]}
    return out


def poltrade():
    """Polarisation at x_a = 90 um with two silvered mirrors, for every valid layout of the scan, against its form
    sensitivity: lower incidence angles mean less metal retardance but a more sensitive (steeper) fold."""
    v = json.loads((RESULTS / "m27_scan.json").read_text())["outputs"]["valid"]
    res = [_poltrade_one(tuple(r["params"])) for r in v]  # serial: ~0.1 s per layout once the facet is cached
    rows = []
    for r, q in zip(v, res):
        rows.append({"params": r["params"], "form_sensitivity_rms": r["form_sensitivity_rms"],
                     "block_length_um": r["block_length_um"], "incidence_M1_deg": r["incidence_M1_deg"],
                     "incidence_M2_deg": r["incidence_M2_deg"], **{k: q[k] for k in q}})
    rows.sort(key=lambda r: abs(r["both silver"]["V_focus"]))
    for r in rows[:8]:
        print(r["params"], "sens %.2f len %.0f |V| %.3f cross %.1e" % (r["form_sensitivity_rms"], r["block_length_um"],
              abs(r["both silver"]["V_focus"]), r["both silver"]["cross_pol_focus"]), flush=True)
    # Pareto front: |V| vs sensitivity
    front = []
    for r in sorted(rows, key=lambda r: r["form_sensitivity_rms"]):
        if not front or abs(r["both silver"]["V_focus"]) < abs(front[-1]["both silver"]["V_focus"]):
            front.append(r)
    for r in front:
        print("front", r["params"], "sens %.2f |V| %.3f" % (r["form_sensitivity_rms"], abs(r["both silver"]["V_focus"])))
    save_result("m27_poltrade", script=SCRIPT,
                inputs={"x_a_um": 90.0, "coatings": "both mirrors silver (Johnson & Christy), and perfect conductors",
                        "designs": "all valid layouts of m27_scan"},
                outputs={"rows": rows, "pareto_V_vs_sensitivity": front})
    return rows


# ---------------------------------------------------------------------------------------------- x axis
def _glass_sections(d):
    if isinstance(d, MC.ReflectiveDesign):
        return [(d.chief_glass_length(), d.n)]
    zs = [0.0] + [sf.z_v for sf in d.surfaces]
    return [(zs[j + 1] - zs[j], n) for j, n in enumerate(d.indices[:-1]) if n > 1]


def xaxis():
    from picsim.opa.focal_plane import FocalPlane, LinearSynth
    T = np.load(RESULTS / "m10_T_opt_A.npz")["T"]
    fp = FocalPlane(A=T, **m14.D1)
    ls = LinearSynth(fp)
    s_cut = fp.opa.s_cut
    xt = np.array([-60.0, -30.0, 0.0, 30.0, 60.0])
    out = {}
    for name, d in (("reflective", _design()), ("refractive", m26.fac_adopted())):
        sig, zb, shift = m26._sigma_from_facet_um(d)
        zS = m26._bfp(d)
        z_c = (sig - zS) * 1e-3 + 0.5
        secs = _glass_sections(d)

        def extra(s, xn=None, secs=secs):
            s = np.asarray(s)
            return 2 * np.pi / LAM * sum(t * (np.sqrt(nn**2 - s**2) - nn + s**2 / (2 * nn)) for t, nn in secs)

        ss = np.linspace(-s_cut, s_cut, 401)
        eff = m14.x_effect(fp, ls, extra, xt)
        out[name] = {"glass_sections_um_n": secs, "S_from_facet_um": zS, "sigma_from_facet_um": sig,
                     "glass_x_shift_um": shift, "z_c_mm": z_c, "far_row_z_mm_32_rows": z_c + 16.0,
                     "slab_phase_pv_waves_after_defocus": float(np.ptp(extra(ss)) / (2 * np.pi)), **eff}
        print(name, json.dumps(out[name], default=float), flush=True)
    save_result("m27_xaxis", script=SCRIPT,
                inputs={"x lens": "m14 aplanat SiO2 back focal distance (results/m14_xlens.json)",
                        "slab": "glass path of the axial ray along the unfolded path; phase k0 t sqrt(n^2 - s^2) minus "
                                "its defocus (absorbed by placing the x lens), applied to every emitter of D1",
                        "comb": "starts 0.5 mm beyond the x slit plane (CHOICE, m26)", "tweezers": xt.tolist()},
                outputs=out)
    return out


# ---------------------------------------------------------------------------------------------- figure
def figure():
    from shapely.affinity import translate
    d = _design()
    reg = d.light_region()
    f, ax = fig(1, 2, h=3.4)
    a = ax[0, 0]
    for k in (-1, 0, 1):
        g = translate(reg.buffer(MARGIN), yoff=k * PITCH)
        for poly in getattr(g, "geoms", [g]):
            zz, yy = poly.exterior.xy
            a.fill(zz, yy, color="0.85" if k else "0.75", lw=0)
        a.add_patch(__import__("matplotlib").patches.Rectangle((-150, k * PITCH - PITCH + 5), 150, PITCH - 10,
                                                                color="0.5"))
    th = np.linspace(-d.theta_m, d.theta_m, 13)
    for t in th:
        p1 = np.array(d.M1.point(t))
        dy, dz = d._dir_after_m1(np.array([t]))
        u2, tt = d.M2.intersect(np.array([p1[0]]), np.array([p1[1]]), dy, dz)
        p2 = np.array(d.M2.point(u2[0]))
        a.plot([0, p1[1], p2[1], d.z_e + 100], [0, p1[0], p2[0], p2[0]], color="C3", lw=0.5)
    tt = np.linspace(-d.theta_m, d.theta_m, 200)
    for M in (d.M1, d.M2):
        yy, zz = M.point(tt)
        a.plot(zz, yy, color="k", lw=1.2)
    a.axvline(d.z_e, color="k", lw=0.8)
    a.set_aspect("equal")
    a.set_xlim(-150, d.z_e + 120)
    a.set_ylim(-800, 400)
    a.set_xlabel("z from facet (um)")
    a.set_ylabel("y (um)")
    rf = json.loads((RESULTS / "m26_skew.json").read_text())["outputs"]["rows"]
    rr = json.loads((RESULTS / "m27_skew.json").read_text())["outputs"]["rows"]
    b = ax[0, 1]
    for rows, lab, c in ((rf, "refractive (m14)", "C0"), (rr, "reflective (m27)", "C3")):
        b.plot([r["x_a_um"] for r in rows], [r["peak_rel_ideal"] for r in rows], "o-", ms=3, color=c, label=lab + ", 1D")
        r2 = [r for r in rows if r["x_a_um"] <= D1_REACH]
        b.plot([r["x_a_um"] for r in r2], [r["twod_peak_rel_ideal"] for r in r2], "s--", ms=3, color=c,
               label=lab + ", 2D (D1 reach)")
    b.axhline(0.97, color="0.5", lw=0.7)
    b.set_xlabel("tweezer position x_a (um)")
    b.set_ylabel("peak / ideal collimator")
    b.legend(fontsize=6)
    save_result("m27_figure", script=SCRIPT, fig=f, inputs={"design": d.name},
                outputs={"layout": d.geometry(PITCH, MARGIN)},
                notes="Left: y-z section of the adopted block (rays every 6.4 deg, mirrors black, exit face), with the "
                      "light regions of the blocks of the neighbouring chips (grey, 10 um margin) and the chips (dark "
                      "grey, z < 0). Right: peak vs x_a, refractive (m26_skew) and reflective (m27_skew); 2D only within the "
                      "67 um reach of the D1 engine of the 2D model.")
    return f


if __name__ == "__main__":
    parts = sys.argv[1:] or ["scan"]
    for part in parts:
        globals()[part]()
