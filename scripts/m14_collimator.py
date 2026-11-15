"""Module 14: real fast-axis collimators for the stacked-chip relay, the relay's first x lens, and the plan.

Parts:
  designs    ray-traced collimators at NA 0.9, f = 175 um: immersed ellipse (one printed/glued surface), plano-
             asphere in air (n = 1.8), Weierstrass block + aplanat (adopted); fused-silica and plano-asphere NA
             limits; geometry, incidence angles, sine-condition error
  wave       each design in the wave model of module 13 (facet in air or immersed): transmission, y spot, row
             crosstalk, peak relative to the ideal lens; beta scan
  tolerance  decentre and focus with the real mapping; exact off-axis traces check the first-order model
  form       Monte Carlo of surface form errors (Gaussian random sag, correlation 20 and 60 um)
  xlens      first x lens of the relay (f = 3.12 mm, NA 0.355, field +/-144 um): off-axis aberration of an aplanat
             and of a plano-asphere, their effect on D1 tweezers before and after recalibration
  plan       per-chip power, input fibres, frequency offsets for 16 chips (DERIVED)

Outputs: results/m14_*.json (+ .png)
"""
from __future__ import annotations

import copy
import json
import sys

import numpy as np

from picsim import materials as mat
from picsim import params as prm
from picsim import raytrace as rt
from picsim.opa import collimator as C
from picsim.opa.array import channel_response
from picsim.opa.focal_plane import FocalPlane, LinearSynth
from picsim.opa.stack import FacetY, Focus2D, RelayY, fwhm, peak2d
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m14_collimator.py"
LAM = prm.LAM_TRANSPORT.value
NA = prm.NA_OBJ.value
N_SIO2 = float(mat.n_sio2(LAM))
N_HI = 1.8  # CHOICE: high-index glass class (dense flint / lanthanum glass) at 1064 nm; select a real glass later
F_FAC = 175.0
BETA = 0.7
NA_FAC = 0.9
D1 = {"N_ports": 128, "N_arms": 192, "w_emit": 1.2, "p_out": 1.5, "w_port": 1.7, "w_arm": 1.7}


def designs_dict():
    return {"ellipse_immersed": C.ellipse_immersed(N_SIO2, F_FAC),
            "plano_asphere_n1.8": C.plano_asphere(N_HI, F_FAC, 30.0),
            "weierstrass_aplanat": C.fac_design(N_SIO2, n2=N_HI)}


def facet_for(d):
    return FacetY(LAM, 1.2, 0.4, n_m=d.n_m if d is not None else 1.0)


def eps_at(y, I, d, half=0.5):
    m = np.abs(np.abs(y) - d) <= half
    return float(I[m].max() / I.max())


def line_metrics(r, y, **err):
    Py = r.pupil(**err)
    I = r.atom_line(y, Py)
    return {"eta_y": r.efficiency(Py), "peak": float(I.max()), "fwhm_um": fwhm(y, I),
            "eps3": eps_at(y, I, 3.0), "eps4": eps_at(y, I, 4.0), "eps6": eps_at(y, I, 6.0), "eps8": eps_at(y, I, 8.0)}


# ------------------------------------------------------------------ designs
def designs():
    out = {"n_SiO2": N_SIO2, "n_hi": N_HI,
           "limits": {"plano-hyperbolic in air (curved side to source), NA max sqrt(1 - 1/n^2)":
                      {"SiO2": float(np.sqrt(1 - 1 / N_SIO2**2)), "n=1.8": float(np.sqrt(1 - 1 / N_HI**2))}}}
    rows = {}
    f, ax = fig(1, 2, h=3.2)
    for name, d in designs_dict().items():
        L = d.lensmap(LAM, NA_FAC)
        g = d.geometry(min(NA_FAC, L.s_max))
        rows[name] = {"s_max": L.s_max, "g_over_s": {str(x): float(np.interp(x, L.s, L.g) / x) for x in (0.3, 0.6, 0.9)},
                      "W_pv_waves": float(np.ptp(L.W)), **g, "construction": getattr(d, "construction", None)}
        ax[0, 0].plot(L.s, L.g / np.maximum(L.s, 1e-9), label=name)
        ax[0, 1].plot(L.s, L.t**2, label=name)
        # geometry sketch
        print(name, json.dumps(rows[name], default=float)[:400], flush=True)
    for nn, wd in ((N_SIO2, 30.0), (N_SIO2, 60.0)):
        d = C.plano_asphere(nn, F_FAC, wd)
        L = d.lensmap(LAM, NA_FAC)
        g = d.geometry(min(NA_FAC, L.s_max))
        rows[f"plano_asphere_SiO2_wd{int(wd)}"] = {"s_max": L.s_max, "max_slope_deg": g["max_slope_deg"],
                                                   "max_incidence_deg": g["max_incidence_deg"],
                                                   "g_over_s_at_smax": float(L.g[-1] / L.s[-1])}
    ax[0, 0].set_xlabel("s = sin(theta) at the facet (n_m sin in the medium)")
    ax[0, 0].set_ylabel("y_S / (f s)  (1: sine condition)")
    ax[0, 0].legend(fontsize=7)
    ax[0, 1].set_xlabel("s")
    ax[0, 1].set_ylabel("Fresnel transmission, uncoated")
    out["rows"] = rows
    save_result("m14_designs", script=SCRIPT, fig=f,
                inputs={"f_um": F_FAC, "NA_FAC": NA_FAC, "n_SiO2": N_SIO2, "n_hi": "CHOICE 1.8",
                        "weierstrass_aplanat": "R 80 um, gap 10 um, aplanat T 250 um (n 1.8), block SiO2",
                        "ellipse_immersed": "SiO2, e = 1/n, vertex n f from the facet",
                        "plano_asphere": "wd 30 um, flat side to facet, exit asphere for exact collimation"},
                outputs=out, notes="2D ray trace (picsim/raytrace.py), s-polarised Fresnel, no coatings.")
    return out


# ------------------------------------------------------------------ scan of the Weierstrass + aplanat design
def _scan_one(p):
    R, gap, T2, n2 = p
    try:
        d = C.weierstrass_aplanat(N_SIO2, F_FAC, R, gap, T2, n2=n2, n_steps=2500)
        if d.construction["s_max_system"] < NA_FAC:
            return None
        g = d.geometry(NA_FAC)
        return {"R_um": R, "gap_um": gap, "T2_um": T2, "n2": n2, "max_incidence_deg": g["max_incidence_deg"],
                "max_slope_deg": g["max_slope_deg"], "T_uncoated_min": g["T_uncoated_min"],
                "length_um": g["z_vertices_um"][-1], "heights_um": g["ray_height_at_surfaces_um"]}
    except Exception:
        return None


def scan():
    import itertools
    from concurrent.futures import ProcessPoolExecutor
    P = list(itertools.product((50.0, 65.0, 80.0, 100.0), (10.0, 30.0, 60.0), (120.0, 160.0, 200.0, 250.0),
                               (N_SIO2, N_HI)))
    with ProcessPoolExecutor(4) as ex:
        res = [r for r in ex.map(_scan_one, P) if r]
    res.sort(key=lambda r: max(r["max_incidence_deg"]))
    print(len(res), "of", len(P), res[0], flush=True)
    save_result("m14_scan", script=SCRIPT,
                inputs={"R_um": [50, 65, 80, 100], "gap_um": [10, 30, 60], "T2_um": [120, 160, 200, 250],
                        "n2": [N_SIO2, N_HI], "criterion": "reach NA 0.9; ranked by the largest incidence angle"},
                outputs={"feasible": len(res), "total": len(P), "rows": res})
    return res


# ------------------------------------------------------------------ wave model
def wave():
    """Each design in the wave model. peak_rel_ideal compares per unit guided power: the line peak (per unit power
    radiated into the source medium) times the facet transmission T_facet of that medium, relative to the ideal
    sine-condition lens in air at beta = 0.7."""
    y = np.arange(-12, 12, 0.005)
    out = {}
    fa_air, fa_im = FacetY(LAM, 1.2, 0.4), FacetY(LAM, 1.2, 0.4, n_m=N_SIO2)
    out["T_facet"] = {"air": fa_air.T_facet, "SiO2": fa_im.T_facet}
    ref = line_metrics(RelayY(fa_air, BETA, NA=NA, NA_FAC=NA_FAC), y)["peak"] * fa_air.T_facet
    for key, fa in (("ideal_air", fa_air), ("ideal_immersed", fa_im)):
        for beta in (0.6, 0.7, 0.8):
            m = line_metrics(RelayY(fa, beta, NA=NA, NA_FAC=NA_FAC), y)
            m["peak_rel_ideal"] = m["peak"] * fa.T_facet / ref
            out[f"{key}|beta{beta}"] = m
            print(key, beta, {k: round(v, 5) for k, v in m.items()}, flush=True)
    for name, d in designs_dict().items():
        fa = fa_im if d.n_m > 1 else fa_air
        for coated in (False, True):
            L = d.lensmap(LAM, NA_FAC, coated=coated)
            for beta in (0.6, 0.7, 0.8, 0.9):
                m = line_metrics(RelayY(fa, beta, NA=NA, NA_FAC=NA_FAC, lens=L), y)
                m["peak_rel_ideal"] = m["peak"] * fa.T_facet / ref
                out[f"{name}|{'AR' if coated else 'bare'}|beta{beta}"] = m
                print(name, coated, beta, {k: round(v, 5) for k, v in m.items()}, flush=True)
    save_result("m14_wave", script=SCRIPT,
                inputs={"beta_um": "0.6-0.9", "NA_FAC": NA_FAC, "facet": "D1 emitter, femwell; immersed designs radiate "
                        "into SiO2", "AR": "coated: Fresnel transmission of the lens surfaces set to 1",
                        "T_facet": "Kirchhoff + normal-incidence Fresnel (FacetY.T_facet); checked by FDTD in m14_facet_meep"},
                outputs=out, notes="1D y line. peak_rel_ideal: peak per unit guided power relative to the ideal "
                                   "sine-condition lens in air, beta 0.7 (power per tweezer scales as its inverse).")
    return out


# ------------------------------------------------------------------ tolerances
def exact_offaxis(d, dy0=0.0, dz0=0.0, s_max=0.9, n=601):
    """Exact trace from a displaced source: optical path change at the exit plane vs exit height, compared with
    the first-order model (-s dy0 for decentre, -n_m cos(theta) dz0 for focus)."""
    s = np.linspace(-s_max, s_max, n)
    th = np.arcsin(s / d.n_m)
    r0, _ = rt.trace(0.0, 0.0, th, d.surfaces, d.indices, d.z_out)
    r1, _ = rt.trace(dy0, dz0, th, d.surfaces, d.indices, d.z_out)
    ok = np.isfinite(r0.y) & np.isfinite(r1.y)
    o0 = np.interp(r1.y[ok], r0.y[ok], r0.opl[ok])
    s_at = np.interp(r1.y[ok], r0.y[ok], s[ok])
    exact = r1.opl[ok] - o0
    model = -s_at * dy0 - d.n_m * np.sqrt(1 - (s_at / d.n_m) ** 2) * dz0
    res = (exact - model) - np.mean(exact - model)
    return float(np.sqrt(np.mean(res**2)) / LAM), float(np.ptp(res) / LAM)


def tolerance():
    y = np.arange(-12, 12, 0.005)
    out = {}
    for name, d in [("ideal", None)] + list(designs_dict().items()):
        L = None if d is None else d.lensmap(LAM, NA_FAC, coated=True)
        r = RelayY(facet_for(d), BETA, NA=NA, NA_FAC=NA_FAC, lens=L)
        nom = line_metrics(r, y)
        rows = {"nominal": nom}
        for dl in (0.05, 0.1, 0.2, 0.3):
            m = line_metrics(r, y, delta=dl)
            rows[f"decentre {dl}"] = {"peak_rel": m["peak"] / nom["peak"], "eps4": m["eps4"], "eta_y": m["eta_y"]}
        for dz in (0.25, 0.5, 1.0):
            m = line_metrics(r, y, dz=dz)
            rows[f"focus {dz}"] = {"peak_rel": m["peak"] / nom["peak"], "eps4": m["eps4"], "eta_y": m["eta_y"]}
        if d is not None:
            rows["exact_check_waves_rms_pv"] = {"decentre 0.2": exact_offaxis(d, dy0=0.2),
                                                "decentre 0.5": exact_offaxis(d, dy0=0.5),
                                                "focus 0.5": exact_offaxis(d, dz0=0.5)}
        out[name] = rows
        print(name, json.dumps(rows, default=float), flush=True)
    save_result("m14_tolerance", script=SCRIPT,
                inputs={"beta_um": BETA, "NA_FAC": NA_FAC, "coatings": "AR (Fresnel transmission 1)",
                        "decentre": "facet relative to the collimator axis (y), um", "focus": "facet along z, um "
                        "(in the source medium)"},
                outputs=out,
                notes="exact_check: residual between exactly traced optical-path change and the first-order model used "
                      "in the wave model (waves rms, peak-to-valley), over |s| <= 0.9.")
    return out


# ------------------------------------------------------------------ form errors
def random_sag(rng, y, sigma, corr):
    """Gaussian random profile with Gaussian correlation of length corr (um) and rms sigma (um)."""
    yy = np.arange(y.min() - 3 * corr, y.max() + 3 * corr, corr / 8)
    w = rng.standard_normal(len(yy))
    k = np.exp(-((np.arange(-24, 25) / 8) ** 2) * 2)
    p = np.convolve(w, k, mode="same")
    p = (p - p.mean()) / p.std() * sigma
    return np.interp(y, yy, p)


def form():
    rng = np.random.default_rng(7)
    y = np.arange(-12, 12, 0.005)
    out = {}
    for name, d in designs_dict().items():
        L = d.lensmap(LAM, NA_FAC, coated=True)
        r = RelayY(facet_for(d), BETA, NA=NA, NA_FAC=NA_FAC, lens=L)
        nom = line_metrics(r, y)
        sg = np.linspace(0, NA_FAC, 901)
        F, Y = d.sensitivity(sg)
        rows = []
        for corr in (20.0, 60.0):
            for sigma in (0.01, 0.02, 0.05, 0.1):
                pk, e4, wr = [], [], []
                for _ in range(24):
                    W = np.zeros_like(sg)
                    for Fj, Yj in zip(F, Y):
                        yy = np.concatenate([-Yj[::-1], Yj])
                        prof = random_sag(rng, yy, sigma, corr)
                        # symmetric fan: use the +s branch; a random profile need not be symmetric, so draw both sides
                        Wp = Fj * np.interp(Yj, yy, prof) / LAM
                        W += Wp
                    sw = np.concatenate([-sg[::-1], sg])
                    Wn = np.concatenate([W[::-1], W])
                    # independent errors on the -s side: redraw
                    W2 = np.zeros_like(sg)
                    for Fj, Yj in zip(F, Y):
                        yy = np.concatenate([-Yj[::-1], Yj])
                        prof = random_sag(rng, yy, sigma, corr)
                        W2 += Fj * np.interp(-Yj, yy, prof) / LAM
                    Wn[: len(sg)] = W2[::-1]
                    Ws = (lambda s, sw=sw, Wn=Wn: np.interp(s, sw, Wn))
                    m = line_metrics(r, y, Ws=Ws)
                    pk.append(m["peak"] / nom["peak"])
                    e4.append(m["eps4"])
                    wr.append(float(np.sqrt(np.mean((Wn - Wn.mean()) ** 2))))
                rows.append({"corr_um": corr, "sigma_nm": sigma * 1e3, "peak_rel_mean": float(np.mean(pk)),
                             "peak_rel_p10": float(np.percentile(pk, 10)), "eps4_p90": float(np.percentile(e4, 90)),
                             "W_rms_waves_mean": float(np.mean(wr))})
                print(name, rows[-1], flush=True)
        out[name] = {"nominal": nom, "rows": rows}
    save_result("m14_form", script=SCRIPT,
                inputs={"samples": 24, "profile": "Gaussian random sag, Gaussian correlation, independent per surface "
                        "and per side", "beta_um": BETA, "coatings": "AR"},
                outputs=out,
                notes="Wavefront from sag errors to first order: dOPL = (n_in dz_in - n_out dz_out) dsag per surface.")
    return out


# ------------------------------------------------------------------ first x lens of the relay
def x_effect(fp, ls, phase_fn, xt, w0=None, I0=None, b0=None, p0=None):
    """D1 tweezers at xt with an extra x-pupil phase phase_fn(s, x_emitter) (rad): peak and position before and
    after least-squares recalibration, and the change of the others when the middle tweezer moves by 0.25 um."""
    if w0 is None:
        w0 = ls.weights(xt)
        I0 = fp.intensity(w0)
        b0 = [fp.peak_near(I0, x) for x in xt]
        p0 = np.sum(np.abs(w0) ** 2)
    fp2 = copy.copy(fp)
    fp2.G = channel_response(fp.opa, fp.x, phase_fn=phase_fn)
    Iu = fp2.intensity(w0)
    bu = [fp2.peak_near(Iu, x) for x in xt]
    ls2 = LinearSynth(fp2)
    w2 = ls2.weights(xt)
    Ic = fp2.intensity(w2)
    bc = [fp2.peak_near(Ic, x) for x in xt]
    p2 = np.sum(np.abs(w2) ** 2)
    # independence after recalibration: move the middle tweezer
    im = len(xt) // 2
    xm = np.array(xt, float).copy()
    xm[im] += 0.25
    Im = fp2.intensity(ls2.weights(xm))
    oth = [fp2.peak_near(Im, x) for i, x in enumerate(xt) if i != im]
    obc = [b for i, b in enumerate(bc) if i != im]
    return {"uncalibrated_peak_rel": [u["peak"] / b["peak"] for u, b in zip(bu, b0)],
            "uncalibrated_pos_nm": [(u["x_peak"] - x) * 1e3 for u, x in zip(bu, xt)],
            "recalibrated_peak_per_power_rel": [(c["peak"] / p2) / (b["peak"] / p0) for c, b in zip(bc, b0)],
            "recalibrated_pos_nm": [(c["x_peak"] - x) * 1e3 for c, x in zip(bc, xt)],
            "recalibrated_move_others": {"peak": float(max(abs(o["peak"] / b["peak"] - 1) for o, b in zip(oth, obc))),
                                         "pos_nm": float(max(abs(o["x_peak"] - b["x_peak"]) for o, b in zip(oth, obc)) * 1e3)}}


def xlens():
    T = np.load(RESULTS / "m10_T_opt_A.npz")["T"]
    fp = FocalPlane(A=T, **D1)
    ls = LinearSynth(fp)
    f1 = 3120.0
    s_cut = fp.opa.s_cut
    # r0 2.0 mm: the ray construction then covers |s| <= 0.42 on axis and up to 0.36 at the 144 um field edge, beyond
    # the slit edge s_cut = 0.355 (r0 1.5 mm stopped at 0.33, so rays near the slit edge were missing)
    lenses = {"aplanat SiO2": C.aplanat(N_SIO2, f1, 2000.0, 1500.0, s_max=0.5),
              "plano-asphere SiO2": C.plano_asphere(N_SIO2, f1, 1500.0, s_max=0.45)}
    out = {"f1_um": f1, "s_cut": s_cut}
    xt = np.array([-60.0, -30.0, 0.0, 30.0, 60.0])
    w0 = ls.weights(xt)
    I0 = fp.intensity(w0)
    b0 = [fp.peak_near(I0, x) for x in xt]
    p0 = np.sum(np.abs(w0) ** 2)
    for name, d in lenses.items():
        zb = C.back_focal_plane(d)
        d.z_out = zb
        sg = np.linspace(-0.4, 0.4, 801)
        th = np.arcsin(sg)
        r0, _ = rt.trace(0.0, 0.0, th, d.surfaces, d.indices, zb)
        xs_tab = np.array([0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 144.0])  # 144 um: edge of the 288 um emitter row
        Wtab, cover = [], []
        for x0 in xs_tab:
            r1, _ = rt.trace(x0, 0.0, th, d.surfaces, d.indices, zb)
            okc = np.isfinite(r1.y)
            cover.append([float(sg[okc].min()), float(sg[okc].max())])
            # plane-wave component s_eff = (exit angle of the bundle) ~ position at the back focal plane / f
            s_eff = r1.y / f1
            ok = np.isfinite(s_eff)
            ideal = -s_eff[ok] * x0
            o = r1.opl[ok] - np.interp(s_eff[ok], r0.y[np.isfinite(r0.y)] / f1, r0.opl[np.isfinite(r0.y)])
            W = o - ideal
            W -= np.interp(0.0, s_eff[ok], W)
            Wtab.append(np.interp(sg, s_eff[ok], W))
        Wtab = np.array(Wtab)  # (x0, s) in um
        m = np.abs(sg) <= s_cut
        rms_edge = float(np.sqrt(np.mean((Wtab[-1][m] - np.polyval(np.polyfit(sg[m], Wtab[-1][m], 1), sg[m])) ** 2)) / LAM)

        def phase_fn(s, xn, Wtab=Wtab):
            ax_ = abs(xn)
            Wx = np.array([np.interp(ax_, xs_tab, Wtab[:, j]) for j in range(len(sg))])
            Wx = np.interp(np.sign(xn) * s if xn != 0 else s, sg, Wx)
            return 2 * np.pi / LAM * Wx

        out[name] = {"back_focal_z_um": zb, "W_rms_waves_field_144um_after_tilt": rms_edge,
                     "launch_s_covered_by_field": dict(zip([f"{x:g}" for x in xs_tab], cover)),
                     "W_pv_um_by_field": [float(np.ptp(Wt[m])) for Wt in Wtab], "field_um": xs_tab.tolist(),
                     **x_effect(fp, ls, phase_fn, xt, w0=w0, I0=I0, b0=b0, p0=p0),
                     "geometry": d.geometry(0.36)}
        print(name, json.dumps(out[name], default=float)[:900], flush=True)
    save_result("m14_xlens", script=SCRIPT,
                inputs={"f1_um": f1, "aplanat": "SiO2, r0 2.0 mm, T 1.5 mm", "plano_asphere": "SiO2, wd 1.5 mm",
                        "tweezers": xt.tolist(), "design": D1},
                outputs=out,
                notes="Off-axis wavefront of the first x lens applied to each emitter's angular spectrum; the second "
                      "lens and the objective are ideal. Recalibrated: least-squares synthesis on the perturbed response.")
    return out


# ------------------------------------------------------------------ plan
def plan():
    """Per-chip power with the facet and collimator of this module, input fibres, frequency offsets (DERIVED)."""
    m13 = json.loads((RESULTS / "m13_power.json").read_text())["outputs"]["P_new_mW"] * 1e-3  # 47 mW, T_facet = 1
    wv = json.loads((RESULTS / "m14_wave.json").read_text())["outputs"]
    fair = json.loads((RESULTS / "m14_facet_meep_air.json").read_text())["outputs"]
    fgl = json.loads((RESULTS / "m14_facet_meep_glass.json").read_text())["outputs"]
    corr = {"air": fair["T_within_s0.9"] / fair["kirchhoff"]["T_within_s0.9"],
            "SiO2": fgl["T_within_s0.9"] / fgl["kirchhoff"]["T_within_s0.9"]}
    T_air = wv["T_facet"]["air"]
    cases = {"ideal lens in air": ("ideal_air|beta0.7", "air"),
             "Weierstrass + aplanat (immersed, AR)": ("weierstrass_aplanat|AR|beta0.7", "SiO2"),
             "ellipse printed on facet (AR, beta 0.9)": ("ellipse_immersed|AR|beta0.9", "SiO2"),
             "plano-asphere n 1.8 in air (AR)": ("plano_asphere_n1.8|AR|beta0.7", "air")}
    per = {}
    for name, (key, med) in cases.items():
        rel = wv[key]["peak_rel_ideal"] * T_air * corr[med]  # peak per guided power relative to the m13 model
        per[name] = {"peak_vs_m13": rel, "P_tweezer_mW": m13 / rel * 1e3, "eps4": wv[key]["eps4"]}
    p_tw = per["Weierstrass + aplanat (immersed, AR)"]["P_tweezer_mW"] * 1e-3
    facet = json.loads((RESULTS / "m6_power.json").read_text())["outputs"]["facet_MW_cm2_at_1.4W_t0.4"]
    p_fibre_max = 1.4 * 20.0 / facet  # W per input facet at 20 MW/cm^2 (CHOICE limit, m6), lensed-fibre tip
    m2 = json.loads((RESULTS / "m2_couplers.json").read_text())["outputs"]
    hi = dict([r for r in m2["rows"] if r["lam"] == "1064" and r["t"] == 0.4 and r["w_tip"] == 0.1][0])
    hi_db = m2["best"]["1064_0.4_HI1060_6.2"]["dB"]
    conv = RESULTS / "m2_tip_convergence.json"
    if conv.exists():  # the 7 um window of the m2 sweep clips the 100 nm tip mode; use the converged window
        last = json.loads(conv.read_text())["outputs"]["1064_0.4_0.1"][-1]
        hi["aeff_um2"], hi_db = last["aeff_um2"], last["dB_mfd_6.2"]
    p_fibre_hi = 20e6 * hi["aeff_um2"] * 1e-8  # HI1060-matched tip (w_tip 0.1 um), same 20 MW/cm^2
    rows = []
    for per_row, depth in ((45, 1.2e-3), (60, 1.2e-3), (45, 0.5e-3), (60, 0.5e-3)):
        P = per_row * p_tw * depth / 1.2e-3
        rows.append({"tweezers_per_row": per_row, "depth_mK": depth * 1e3, "P_chip_W": P,
                     "fibres_per_chip_lensed_tip": int(np.ceil(P / p_fibre_max)),
                     "fibres_per_chip_HI1060_tip": int(np.ceil(P / p_fibre_hi)),
                     "facet_MW_cm2_HI1060_tip": P / (hi["aeff_um2"] * 1e-8) / 1e6, "P_16_chips_W": 16 * P})
    offs = 80e6 + 1.13e6 * np.arange(16)  # e.g. one acousto-optic shifter per chip feed (CHOICE)
    diffs = np.abs(offs[:, None] - offs[None, :])[np.triu_indices(16, 1)]
    out = {"facet_fdtd_correction": corr, "per_design": per, "P_tweezer_adopted_W": p_tw, "P_fibre_max_W": p_fibre_max,
           "P_fibre_max_HI1060_tip_W": p_fibre_hi, "HI1060_tip_coupling_dB": hi_db, "HI1060_tip_aeff_um2": hi["aeff_um2"],
           "rows": rows, "offsets_MHz": (offs / 1e6).tolist(), "min_pair_difference_MHz": float(diffs.min() / 1e6),
           "max_pair_difference_MHz": float(diffs.max() / 1e6), "2nu_trap_kHz": 2 * 84.0}
    print(json.dumps(out, indent=1), flush=True)
    save_result("m14_plan", script=SCRIPT,
                inputs={"P_m13": "m13_power (47 mW at 1.2 mK, facet transmission 1)", "facet": "m14_facet_meep (2D FDTD) "
                        "scaling the Kirchhoff estimate", "facet_limit": "20 MW/cm^2 (CHOICE, m6)",
                        "trap": "84 kHz radial at 1.2 mK"},
                outputs=out)
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["designs", "scan", "wave", "tolerance", "form", "xlens", "plan"]
    for p in parts:
        {"designs": designs, "scan": scan, "wave": wave, "tolerance": tolerance, "form": form, "xlens": xlens, "plan": plan}[p]()
