"""Module 13: stacked-chip anamorphic relay for a 2D array of independently steered rows.

Each chip is design D1 (one row of independently steerable tweezers along x). Chips are stacked in y at pitch d_c,
each with its own fast-axis collimator (FAC); a common y-relay and the objective image the stack onto the atoms
with demagnification D = d_c / Delta. Model: picsim/opa/stack.py (scalar wave optics; x from the BPM matrix of
D1, y from the femwell facet mode; circular pupil).

Parts:
  facet      vertical facet mode of the emitters (400 / 300 / 200 nm film): far-field divergence, FAC capture
  fill       y pupil fill beta (pupil radius referred to the facet) x FAC aperture: efficiency, y spot,
             row crosstalk vs distance, power relative to the Gaussian (w0 = 1.04 um) assumption used so far
  design     first-order anamorphic layout (ABCD check of both axes) for 16 chips at 0.5 mm, Delta = 4 um
  crosstalk  2D: a full row (41 tweezers at 3 um) of chip m seen by the tweezers of chip m+1; coherent (same
             laser frequency) and incoherent (frequency-offset chips); dynamic: one row-m tweezer moves
  tolerance  per-chip errors: FAC decentre, focus, roll, smile, aberration; chip piston (x axis); row position
  power      2D peak per unit power -> laser power per tweezer, replacing the Gaussian cross-axis assumption

Outputs: results/m13_*.json (+ .png)
"""
from __future__ import annotations

import copy
import json
import sys

import numpy as np
from numpy.polynomial import legendre as L

from picsim import modes as Mo
from picsim import params as prm
from picsim.abcd import abcd_free, abcd_lens, chain
from picsim.opa.array import channel_response
from picsim.opa.focal_plane import FocalPlane, LinearSynth
from picsim.opa.stack import FacetY, Focus2D, RelayY, fwhm, peak2d, raised_cosine
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m13_stack.py"
LAM = prm.LAM_TRANSPORT.value
NA = prm.NA_OBJ.value
D1 = {"N_ports": 128, "N_arms": 192, "w_emit": 1.2, "p_out": 1.5, "w_port": 1.7, "w_arm": 1.7}
W0_OLD = prm.W0_TWEEZER.value  # cross-axis Gaussian assumed in m9/m10 power estimates
D_CHIP = 500.0  # um, chip pitch in the stack (CHOICE: 0.5 mm die incl. bond line)
F_OBJ = prm.WD.value * np.tan(np.arcsin(NA)) / NA  # mm: pupil radius WD tan(asin NA) = f_obj NA (DERIVED)
BETA = 0.7  # design y fill (set by the 'fill' part; see docs/stacked_relay.md)
NA_FAC = 0.9  # FAC numerical aperture (sine of the largest accepted angle in air)


def eps_at(y, I, d, half=0.5):
    """Largest intensity within +/-half of distance d (both sides), relative to the peak."""
    m = (np.abs(np.abs(y) - d) <= half)
    return float(I[m].max() / I.max())


# ------------------------------------------------------------------ facet
def facet():
    rows = []
    f, ax = fig(1, 2, h=3.0)
    for t in (0.4, 0.2, 0.15, 0.1):
        fa = FacetY(LAM, 1.2, t)
        # coupling length of two 1.2 um facets at 0.3 um gap (1.5 um emitter pitch): cost of a thin facet section
        ms = Mo.solve(Mo.XSection(cores=[(-0.75, 1.2), (0.75, 1.2)], t=t, x_left=-6, x_right=6), LAM, num_modes=4)
        te = sorted([m.neff.real for m in ms if m.te_fraction > 0.5], reverse=True)[:2]
        rows.append({"t_um": t, "neff": fa.neff, "NA_1e2": fa.NA_1e2, "w_equiv_um": fa.w_equiv,
                     "Lc_pair_0p3gap_um": LAM / (2 * (te[0] - te[1])),
                     "z_R_equiv_um": np.pi * fa.w_equiv**2 / LAM,
                     "capture": {str(n): fa.fraction_within(n) for n in (0.5, 0.6, 0.7, 0.8, 0.9)}})
        ax[0, 0].plot(fa.y, np.abs(fa.Y), label=f"t = {t} um")
        d = np.abs(fa.A) ** 2 * fa.cos
        ax[0, 1].plot(fa.s, d / d.max(), label=f"t = {t} um")
        print(rows[-1], flush=True)
    ax[0, 0].set_xlim(-1.5, 1.5)
    ax[0, 0].set_xlabel("y at the facet (um)")
    ax[0, 0].set_ylabel("|E| (x-integrated)")
    ax[0, 1].set_xlabel("sin(theta_y) in air")
    ax[0, 1].set_ylabel("far-field intensity")
    ax[0, 0].legend(fontsize=7)
    save_result("m13_facet", script=SCRIPT, fig=f,
                inputs={"emitter": "1.2 um wide SiN facet, oxide clad, TE0, 1064.435 nm", "medium": "air"},
                outputs={"rows": rows},
                notes="Far field from the x-integrated Ex facet profile (Kirchhoff), power density |A|^2 cos(theta).")
    return rows


# ------------------------------------------------------------------ fill
def fill():
    fa = FacetY(LAM, 1.2, 0.4)
    y = np.arange(-12, 12, 0.005)
    peak_old = np.sqrt(2 / np.pi) / W0_OLD
    rows = []
    f, ax = fig(1, 2, h=3.2)
    cases = [(0.7, None, "hard"), (0.8, None, "hard"), (0.9, None, "hard"), (0.95, None, "hard"),
             (0.9, raised_cosine(0.5, 0.9), "graded 0.5-0.9"), (0.95, raised_cosine(0.6, 0.95), "graded 0.6-0.95")]
    for naf, ap, edge in cases:
        for beta in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2):
            r = RelayY(fa, beta, NA=NA, NA_FAC=naf, apod=ap)
            Py = r.pupil()
            I = r.atom_line(y, Py)
            row = {"NA_FAC": naf, "edge": edge, "beta_um": beta, "efficiency": r.efficiency(Py), "fwhm_um": fwhm(y, I),
                   "line_peak": float(I.max()), "power_factor_vs_gauss": float(peak_old / I.max()),
                   "eps": {str(d): eps_at(y, I, d) for d in (3.0, 4.0, 5.0, 6.0, 8.0)}}
            rows.append(row)
            print(naf, edge, beta, round(row["efficiency"], 4), round(row["fwhm_um"], 3),
                  round(row["power_factor_vs_gauss"], 3), {k: f"{v:.1e}" for k, v in row["eps"].items()}, flush=True)
            if naf == NA_FAC and edge == "hard" and beta in (0.5, 0.7, 1.0):
                ax[0, 0].semilogy(y, I / I.max(), lw=0.8, label=f"beta = {beta} um")
    sel = [r for r in rows if r["NA_FAC"] == NA_FAC and r["edge"] == "hard"]
    ax[0, 0].set_xlim(0, 10)
    ax[0, 0].set_ylim(1e-8, 1.2)
    ax[0, 0].set_xlabel("y at the atoms (um)")
    ax[0, 0].set_ylabel("intensity / peak")
    ax[0, 0].legend(fontsize=7)
    ax[0, 1].semilogy([r["power_factor_vs_gauss"] for r in sel], [r["eps"]["4.0"] for r in sel], "-o", ms=3,
                      label="4 um")
    ax[0, 1].semilogy([r["power_factor_vs_gauss"] for r in sel], [r["eps"]["3.0"] for r in sel], "-s", ms=3,
                      label="3 um")
    ax[0, 1].set_xlabel("power per tweezer relative to Gaussian w0 = 1.04 um")
    ax[0, 1].set_ylabel("row crosstalk (intensity)")
    ax[0, 1].legend(fontsize=7)
    save_result("m13_fill", script=SCRIPT, fig=f,
                inputs={"facet": "D1 emitter, 1.2 x 0.4 um", "NA": NA, "beta": "pupil radius referred to the facet",
                        "eps": "max intensity within +/-0.5 um of the given distance / peak (1D, y only)"},
                outputs={"rows": rows})
    return rows


# ------------------------------------------------------------------ first-order design
def design(delta=4.0, n_chips=16, beta=BETA, m_t=25.0, f_a=100.0):
    """First-order layout (mm). y: facet -> FAC -> S -> telescope (f_a, f_b) -> Fourier lens f_2 -> pupil;
    x: emitter plane -> f_1x -> slit -> f_2x -> pupil (4f). Lengths matched so both pupils coincide."""
    D = D_CHIP / delta
    F_y = D * F_OBJ  # mm
    f_fac = beta * D / NA * 1e-3  # mm
    f_b = f_a / m_t
    f_2 = F_y / m_t
    Ly = 2 * f_fac + 2 * f_a + 2 * f_b + 2 * f_2
    A_x = D1["N_arms"] * D1["p_out"] * 1e-3  # mm, emitter aperture
    M_x = 2 * F_OBJ * NA / A_x
    f_1x = Ly / (2 * (1 + M_x))
    f_2x = M_x * f_1x
    # ABCD checks (mm, rad)
    y_facet_pupil = chain(abcd_free(f_fac), abcd_lens(f_fac), abcd_free(f_fac), abcd_free(f_a), abcd_lens(f_a),
                          abcd_free(f_a + f_b), abcd_lens(f_b), abcd_free(f_b + f_2), abcd_lens(f_2), abcd_free(f_2))
    y_S_pupil = chain(abcd_free(f_a), abcd_lens(f_a), abcd_free(f_a + f_b), abcd_lens(f_b), abcd_free(f_b + f_2),
                      abcd_lens(f_2), abcd_free(f_2))
    x_chip_pupil = chain(abcd_free(f_1x), abcd_lens(f_1x), abcd_free(f_1x + f_2x), abcd_lens(f_2x), abcd_free(f_2x))
    obj = chain(abcd_free(F_OBJ), abcd_lens(F_OBJ), abcd_free(F_OBJ))  # pupil -> atoms
    y_S_atoms = obj @ y_S_pupil
    stack_h = n_chips * D_CHIP * 1e-3
    out = {
        "delta_um": delta, "n_chips": n_chips, "D_demag": D, "f_obj_mm": F_OBJ, "F_y_mm": F_y, "beta_um": beta,
        "f_FAC_um": f_fac * 1e3, "collimated_beam_half_height_um": f_fac * 1e3 * NA_FAC,
        "collimated_beam_1e2_radius_um": f_fac * 1e3 * FacetY(LAM, 1.2, 0.4).NA_1e2,
        "y_telescope": {"m": m_t, "f_a_mm": f_a, "f_b_mm": f_b, "NA_after": 0.0},
        "f_2_mm": f_2, "length_facet_to_pupil_mm": Ly, "x_relay": {"M_x": M_x, "f_1x_mm": f_1x, "f_2x_mm": f_2x,
                                                                   "slit_from_facet_mm": 2 * f_1x},
        "stack_height_mm": stack_h, "field_at_f_a_half_mm": stack_h / 2,
        "rows_at_atoms_um": [(i - (n_chips - 1) / 2) * delta for i in (0, n_chips - 1)],
        "abcd": {"y_S_to_atoms": y_S_atoms.tolist(), "y_facet_to_pupil": y_facet_pupil.tolist(),
                 "x_chip_to_pupil": x_chip_pupil.tolist()},
        "checks": {"y S->atoms is imaging (B=0)": abs(y_S_atoms[0, 1]) < 1e-9 * F_y,
                   "y magnification": float(y_S_atoms[0, 0]), "1/D": 1 / D,
                   "y facet->pupil is imaging (B=0)": abs(y_facet_pupil[0, 1]) < 1e-6,
                   "y facet->pupil magnification": float(y_facet_pupil[0, 0]),
                   "x chip->pupil is imaging (B=0)": abs(x_chip_pupil[0, 1]) < 1e-9,
                   "x magnification": float(x_chip_pupil[0, 0]),
                   "y telecentric at S (C of S->pupil = -1/F_y, A = 0)": abs(y_S_pupil[0, 0]) < 1e-9},
    }
    # collimated beam NA from each chip and NA after the telescope
    na_c = (1.0 / np.pi) * LAM * 1e-3 / (out["collimated_beam_1e2_radius_um"] * 1e-3)
    out["NA_collimated"] = na_c
    out["y_telescope"]["NA_after"] = na_c * m_t
    out["lagrange_check"] = {"stack_half_height_x_NA_mm": stack_h / 2 * na_c,
                             "atom_half_field_x_NA_mm": (n_chips / 2) * delta * 1e-3 * NA * (na_c * D / NA)}
    return out


def design_all():
    rows = [design(d) for d in (3.0, 4.0, 5.0, 6.0)]
    for r in rows:
        print({k: r[k] for k in ("delta_um", "D_demag", "F_y_mm", "f_FAC_um", "length_facet_to_pupil_mm")},
              r["x_relay"], r["checks"], flush=True)
    save_result("m13_design", script=SCRIPT,
                inputs={"d_chip_um": D_CHIP, "f_obj_mm": F_OBJ, "NA": NA, "beta_um": BETA, "NA_FAC": NA_FAC,
                        "chips": 16, "y_telescope_m": 25, "f_a_mm": 100},
                outputs={"rows": rows},
                notes="Thin cylindrical lenses, first order. f_obj from the pupil radius WD tan(asin NA) = f_obj NA.")
    return rows


# ------------------------------------------------------------------ shared 2D setup
_FP = {}


def setup(beta=BETA, naf=NA_FAC):
    if "fp" not in _FP:
        T = np.load(RESULTS / "m10_T_opt_A.npz")["T"]
        _FP["fp"] = FocalPlane(A=T, **D1)
        _FP["ls"] = LinearSynth(_FP["fp"])
        _FP["fa"] = FacetY(LAM, 1.2, 0.4)
    r = RelayY(_FP["fa"], beta, NA=NA, NA_FAC=naf)
    return _FP["fp"], _FP["ls"], r, Focus2D(_FP["fp"], r)


def unit(w):
    return w / np.linalg.norm(w)


# ------------------------------------------------------------------ 2D crosstalk between rows
def crosstalk(beta=BETA):
    fp, ls, r, F = setup(beta)
    rng = np.random.default_rng(3)
    xt = np.arange(-60.0, 60.0 + 1e-9, 3.0)
    th = np.where(np.arange(len(xt)) % 2 == 0, 0.0, np.pi / 2)
    w = ls.weights(xt, phases=th)
    x = np.arange(-63, 63, 0.04)
    y = np.arange(-1.6, 17.0, 0.04)
    E = F.field(w, x, y)  # chip m, row at y = 0
    I = np.abs(E) ** 2
    pk = np.mean([peak2d(I[:, np.abs(x - x0) < 0.6], x[np.abs(x - x0) < 0.6], y)[0] for x0 in xt])
    prof = I.max(axis=1) / pk  # worst over x, vs y
    eps2d = {str(d): float(prof[np.abs(np.abs(y) - d) <= 0.5].max()) for d in (3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 16.0)}
    # single tweezer, for reference
    E1 = F.field(unit(ls.weights([0.0])), x, y)
    I1 = np.abs(E1) ** 2
    p1 = I1.max()
    eps_single = {str(d): float((I1.max(axis=1) / p1)[np.abs(np.abs(y) - d) <= 0.5].max()) for d in (3.0, 4.0, 6.0, 8.0)}
    print("row eps", eps2d, "single", eps_single, flush=True)

    # perturbation of chip m+1 tweezers (row at y = Delta) by chip m's full row
    res = {}
    for delta in (3.0, 4.0, 6.0, 8.0, 12.0, 16.0):
        yy = np.arange(delta - 1.2, delta + 1.2, 0.02)
        xs_n = xt[15:26] + 1.5  # 11 tweezers of chip m+1 near the centre, offset by half a spacing
        xx = np.arange(xs_n[0] - 1.5, xs_n[-1] + 1.5, 0.02)
        En = F.field(ls.weights(xs_n, phases=np.where(np.arange(len(xs_n)) % 2 == 0, 0.0, np.pi / 2)), xx, yy - delta)
        Em = F.field(w, xx, yy)
        xs_m2 = xt.copy()
        xs_m2[20] += 0.25  # one tweezer of chip m moves (it sits under xs_n[5] +/- 1.5 um)
        Em2 = F.field(ls.weights(xs_m2, phases=th), xx, yy)

        def tweezers(I):
            out = []
            for x0 in xs_n:
                m = np.abs(xx - x0) < 0.6
                out.append(peak2d(I[:, m], xx[m], yy))
            return np.array(out)

        base = tweezers(np.abs(En) ** 2)
        rows = {}
        for mode in ("incoherent", "coherent"):
            dpk, dpos, dmove_pk, dmove_pos = 0.0, 0.0, 0.0, 0.0
            phis = [0.0] if mode == "incoherent" else np.linspace(0, 2 * np.pi, 8, endpoint=False)
            for ph in phis:
                if mode == "incoherent":
                    Ia = np.abs(En) ** 2 + np.abs(Em) ** 2
                    Ib = np.abs(En) ** 2 + np.abs(Em2) ** 2
                else:
                    Ia = np.abs(En + np.exp(1j * ph) * Em) ** 2
                    Ib = np.abs(En + np.exp(1j * ph) * Em2) ** 2
                a, b = tweezers(Ia), tweezers(Ib)
                dpk = max(dpk, np.max(np.abs(a[:, 0] / base[:, 0] - 1)))
                dpos = max(dpos, np.max(np.hypot(a[:, 1] - base[:, 1], a[:, 2] - base[:, 2])))
                dmove_pk = max(dmove_pk, np.max(np.abs(b[:, 0] / a[:, 0] - 1)))
                dmove_pos = max(dmove_pos, np.max(np.hypot(b[:, 1] - a[:, 1], b[:, 2] - a[:, 2])))
            rows[mode] = {"static_peak_change": float(dpk), "static_pos_nm": float(dpos * 1e3),
                          "move_peak_change": float(dmove_pk), "move_pos_nm": float(dmove_pos * 1e3)}
        res[str(delta)] = rows
        print(delta, rows, flush=True)

    f, ax = fig(1, 2, h=3.2)
    m = (x > -12) & (x < 12)
    ax[0, 0].imshow(np.log10(I[:, m] / pk + 1e-9), extent=[x[m][0], x[m][-1], y[0], y[-1]], origin="lower",
                    vmin=-8, vmax=0, aspect="auto", cmap="magma")
    ax[0, 0].set_xlabel("x (um)")
    ax[0, 0].set_ylabel("y (um)")
    ax[0, 0].set_title("chip m, 41 tweezers at 3 um, log10 I/peak", fontsize=8)
    ax[0, 1].semilogy(y, prof, lw=0.8, label="full row, max over x")
    ax[0, 1].semilogy(y, I1.max(axis=1) / p1, lw=0.8, label="single tweezer")
    ax[0, 1].set_xlabel("y (um)")
    ax[0, 1].set_ylim(1e-9, 1.5)
    ax[0, 1].legend(fontsize=7)
    save_result("m13_crosstalk", script=SCRIPT, fig=f,
                inputs={"beta_um": beta, "NA_FAC": NA_FAC, "row_m": "41 tweezers, 3 um, +/-60 um, phases 0/pi/2",
                        "row_m+1": "11 tweezers at 3 um offset by 1.5 um in x", "move": "row-m tweezer at 0 um "
                        "moves +0.25 um", "coherent": "worst of 8 relative phases between the chips"},
                outputs={"eps_row_max_over_x": eps2d, "eps_single": eps_single, "perturbation": res,
                         "row_peak_mean": float(pk)})
    return res


# ------------------------------------------------------------------ tolerances
def metrics(F, w, x0=0.0, Py=None, shear=None, **kw):
    """2D peak, peak position, y FWHM and crosstalk at 3, 4 and 6 um for a single tweezer at x0."""
    x = np.arange(x0 - 1.2, x0 + 1.2, 0.02)
    y = np.arange(-7, 7, 0.02)
    I = F.intensity(w, x, y, Py=Py, shear=shear, **kw)
    p, xp, yp = peak2d(I, x, y)
    prof = I.max(axis=1) / p
    j = int(np.argmin(np.abs(x - xp)))
    return {"peak": p, "x_nm": (xp - x0) * 1e3, "y_nm": yp * 1e3, "fwhm_y": fwhm(y, I[:, j]),
            "eps3": float(prof[np.abs(np.abs(y) - 3) <= 0.5].max()), "eps4": float(prof[np.abs(np.abs(y) - 4) <= 0.5].max()),
            "eps6": float(prof[np.abs(np.abs(y) - 6) <= 0.5].max())}


def best_focus(F, w, Py=None, shear=None, zs=np.arange(-3.0, 3.01, 0.1)):
    x = np.arange(-0.8, 0.8, 0.02)
    y = np.arange(-1.5, 1.5, 0.02)
    pk = [F.intensity(w, x, y, Py=Py, shear=shear, dz_a=z).max() for z in zs]
    i = int(np.argmax(pk))
    return float(zs[i]), float(pk[i])


def tolerance(beta=BETA):
    fp, ls, r, F = setup(beta)
    w0 = unit(ls.weights([0.0]))
    w60 = unit(ls.weights([60.0]))
    nom = metrics(F, w0)
    nom60 = metrics(F, w60, 60.0)
    z0, p0 = best_focus(F, w0)
    out = {"nominal": nom, "nominal_60um": nom60, "nominal_best_focus": [z0, p0]}
    print("nominal", nom, z0, p0, flush=True)

    def rel(m, ref):
        return {"peak_rel": m["peak"] / ref["peak"], "x_nm": m["x_nm"], "y_nm": m["y_nm"], "fwhm_y": m["fwhm_y"],
                "eps3": m["eps3"], "eps4": m["eps4"], "eps6": m["eps6"]}

    rows = []
    for d in (0.05, 0.1, 0.2, 0.3, 0.5):
        Py = r.pupil(delta=d)
        m = rel(metrics(F, w0, Py=Py), nom)
        rows.append({"error": "FAC decentre", "value": d, "unit": "um", "eta_y": r.efficiency(Py), **m})
    for dz in (0.25, 0.5, 1.0, 2.0, -1.0):
        Py = r.pupil(dz=dz)
        m = rel(metrics(F, w0, Py=Py), nom)
        zb, pb = best_focus(F, w0, Py=Py)
        rows.append({"error": "FAC focus", "value": dz, "unit": "um", "eta_y": r.efficiency(Py), **m,
                     "best_focus_um": zb, "peak_best_rel": pb / p0})
    half = fp.opa.A / 2
    for a in (0.25e-3, 0.5e-3, 1e-3, 2e-3):
        sh = (lambda xi, a=a: a * half * xi / beta)
        m = rel(metrics(F, w0, shear=sh), nom)
        m60 = rel(metrics(F, w60, 60.0, shear=sh), nom60)
        rows.append({"error": "FAC roll", "value": a * 1e3, "unit": "mrad", **m, "edge_60um": m60})
    for b in (0.05, 0.1, 0.2, 0.4):
        sh = (lambda xi, b=b: b * (xi**2 - 1 / 3) / beta)
        m = rel(metrics(F, w0, shear=sh), nom)
        m60 = rel(metrics(F, w60, 60.0, shear=sh), nom60)
        rows.append({"error": "smile (bow, peak-to-valley)", "value": b, "unit": "um", **m, "edge_60um": m60})
    for n, name in ((2, "FAC aberration, defocus-like (P2)"), (3, "FAC aberration, coma-like (P3)"),
                    (4, "FAC aberration, spherical-like (P4)")):
        for a in (0.05, 0.1, 0.25):
            c = np.zeros(n + 1)
            c[n] = a * np.sqrt(2 * n + 1)
            Wf = (lambda u, c=c: L.legval(u, c))
            Py = r.pupil(W=Wf)
            m = rel(metrics(F, w0, Py=Py), nom)
            zb, pb = best_focus(F, w0, Py=Py)
            rows.append({"error": name, "value": a, "unit": "waves rms over the FAC aperture",
                         "eta_y": r.efficiency(Py), **m, "best_focus_um": zb, "peak_best_rel": pb / p0})
    for row in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if k != "edge_60um"}, flush=True)
    out["rows"] = rows

    # chip piston (emitter plane displaced along z from the x-relay conjugate): x axis only, D1 1D model
    chipz = []
    xs_t = [0.0, 30.0, 60.0]
    I0 = fp.intensity(ls.weights(xs_t))
    b0 = [fp.peak_near(I0, x0) for x0 in xs_t]
    for dz in (10.0, 30.0, 100.0, 300.0):
        fp2 = copy.copy(fp)
        fp2.G = channel_response(fp.opa, fp.x, dz=dz)
        Iu = fp2.intensity(ls.weights(xs_t))
        ls2 = LinearSynth(fp2)
        Ic = fp2.intensity(ls2.weights(xs_t))
        bu = [fp2.peak_near(Iu, x0) for x0 in xs_t]
        bc = [fp2.peak_near(Ic, x0) for x0 in xs_t]
        p_u = np.sum(np.abs(ls.weights(xs_t)) ** 2)
        p_c = np.sum(np.abs(ls2.weights(xs_t)) ** 2)
        chipz.append({"dz_um": dz,
                      "uncalibrated_peak_rel": [u["peak"] / b["peak"] for u, b in zip(bu, b0)],
                      "uncalibrated_pos_nm": [(u["x_peak"] - x0) * 1e3 for u, x0 in zip(bu, xs_t)],
                      "recalibrated_peak_per_power_rel": [(c["peak"] / p_c) / (b["peak"] / p_u) for c, b in zip(bc, b0)],
                      "recalibrated_pos_nm": [(c["x_peak"] - x0) * 1e3 for c, x0 in zip(bc, xs_t)]})
        print("chip z", chipz[-1], flush=True)
    out["chip_piston_x"] = chipz
    save_result("m13_tolerance", script=SCRIPT,
                inputs={"beta_um": beta, "NA_FAC": NA_FAC, "tweezer": "single, LS synthesis, x = 0 (and 60 um)",
                        "aberration": "Legendre P_n over the FAC aperture u = s/NA_FAC, rms a (uniform weighting)",
                        "roll": "facet decentre a * x along the 288 um emitter row", "smile": "b (xi^2 - 1/3)"},
                outputs=out,
                notes="peak_rel: 2D peak at the nominal atom plane relative to no error; y_nm: row (y) position shift; "
                      "eps: worst intensity within +/-0.5 um of 3, 4, 6 um from the row, relative to the peak.")
    return out


# ------------------------------------------------------------------ power
def power(beta=BETA):
    fp, ls, r, F = setup(beta)
    old = json.loads((RESULTS / "m10_route_E_LS_opt_A.json").read_text())["outputs"]
    rows = []
    for x0 in (0.0, 30.0, 60.0):
        w = unit(ls.weights([x0]))
        I1 = fp.peak_near(fp.intensity(w), x0)["peak"]
        gauss = I1 * np.sqrt(2 / np.pi) / W0_OLD
        x = np.arange(x0 - 1.0, x0 + 1.0, 0.02)
        y = np.arange(-1.5, 1.5, 0.02)
        Ic = F.intensity(w, x, y)
        Ip = F.intensity(w, x, y, circle=False)
        rows.append({"x_um": x0, "peak_gauss_assumption": gauss, "peak_2d": float(Ic.max()),
                     "peak_2d_square_pupil": float(Ip.max()), "factor_vs_assumption": float(gauss / Ic.max()),
                     "circle_loss": float(1 - Ic.max() / Ip.max())})
        print(rows[-1], flush=True)
    fac = float(np.mean([r_["factor_vs_assumption"] for r_ in rows]))
    p_old = 43.0  # mW per 1.2 mK tweezer, field mean (docs/fpr_design.md, DERIVED from m9_E_power, m10 LS)
    out = {"rows": rows, "mean_factor": fac, "P_old_mW": p_old, "P_new_mW": p_old * fac,
           "FAC_capture": r.f.fraction_within(NA_FAC), "eta_y": r.efficiency(),
           "single_P_eff_old": old["single_P_eff"]}
    print(out, flush=True)
    save_result("m13_power", script=SCRIPT,
                inputs={"beta_um": beta, "NA_FAC": NA_FAC, "assumption_replaced": "cross axis Gaussian w0 = 1.04 um, "
                        "no cross-axis loss (m9_E_power)"},
                outputs=out,
                notes="factor_vs_assumption multiplies the laser power per tweezer estimated earlier (43 mW at 1.2 mK).")
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["facet", "fill", "design", "power", "crosstalk", "tolerance"]
    for p in parts:
        {"facet": facet, "fill": fill, "design": design_all, "crosstalk": crosstalk, "tolerance": tolerance,
         "power": power}[p]()
