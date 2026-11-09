"""Module 26: combining several 16-chip stacks into one engine block of 32 to 64 rows.

The plane S (back focal planes of the per-chip collimators, FAC) is conjugate to the atoms with magnification 1/D
(D = 125 at 4 um rows). At S each chip's beam is a stripe of +/- f_FAC NA_FAC = +/-158 um on the 0.5 mm chip pitch,
so neighbouring beams leave 184 um gaps, and the beams stay nearly collimated (Rayleigh range ~66 mm). A slotted
45 deg mirror ("comb") in these gaps interleaves the beams of two stacks built at twice the chip pitch: stack A
passes through the slots, stack B is reflected by the teeth. The cost is edge diffraction: the FAC aperture is a
hard edge, its Fresnel tail spreads into the gap as sqrt(lambda z), the comb edge clips it, and since S is conjugate
to the atoms the clipped tail sits midway between rows.

Parts:
  edge    1D wave optics (y): field at the FAC exit (aperture plane, f_FAC before S) -> propagate to the comb at
          distance z from S -> window (slot |y| < a, or a single edge y < a) -> propagate back to S -> pupil ->
          atoms (picsim.opa.stack.RelayY). Transmission, peak, crosstalk at 3 / 4 / 8 um, row shift vs z and a.
  twod    2D atom-plane check through the circular pupil (m13 Focus2D): single tweezers at x = 0 / 60 um, a full
          row of 41 tweezers (crosstalk max over x), and the neighbour row (11 tweezers) static and with a moving
          tweezer, incoherent and coherent, with the comb-windowed y pupil of both chips; comb roll and edge
          straightness by a local window a(X) (Sigma space is telecentric in x)
  lens    the same with the adopted traced collimator (glued Weierstrass + aplanat, docs/collimator.md), on axis,
          with the real x slit-plane distance (m14 x lens)
  skew    x dependence through the cylindrical collimator: light of a tweezer at x_a has k_x = k0 x_a / x_scale,
          and the y-z problem is the 2D trace with effective indices sqrt(n^2 - sx^2)
  layout  first-order layouts for 32 / 48 / 64 rows (m13 design), comb z-range, path equality, tilt tolerances,
          chip pitch and I/O per chip.

Output: results/m26_edge.json (+ .png), results/m26_layout.json
"""
from __future__ import annotations

import copy
import sys

import numpy as np

from picsim import materials as mat
from picsim.opa import collimator as C
from picsim.opa.stack import FacetY, Focus2D, RelayY, raised_cosine
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

sys.path.insert(0, str(RESULTS.parent / "scripts"))
from m13_stack import BETA, D_CHIP, F_OBJ, LAM, NA, NA_FAC, design, eps_at, metrics, setup, unit  # noqa: E402
from picsim.opa.stack import peak2d  # noqa: E402

SCRIPT = "scripts/m26_stack_combine.py"
N_SIO2 = float(mat.n_sio2(LAM))
N_HI = 1.8  # aplanat glass of the adopted collimator (m14, CHOICE)
DELTA = 4.0  # um, row pitch at the atoms
D = D_CHIP / DELTA
F_FAC = BETA * D / NA  # um
DY = 0.5  # um, grid at S
NY = 2**16  # 32.8 mm span


def skewed(d, sx):
    """The cylindrical collimator d seen by light with k_x = k0 sx: the normals of its surfaces have no x component,
    so k_x is conserved and the y-z problem is the 2D trace with effective indices sqrt(n^2 - sx^2) (exact)."""
    dd = copy.deepcopy(d)
    dd.indices = [float(np.sqrt(n**2 - sx**2)) for n in d.indices]
    return dd


def _bfp(d):
    """S of a traced collimator: its own method (folded reflective designs, picsim.opa.mirror_collimator) or the
    chief-ray construction of picsim.opa.collimator."""
    own = getattr(d, "back_focal_plane", None)
    return own() if callable(own) else C.back_focal_plane(d)


def fac_adopted():
    """Adopted collimator of docs/collimator.md: Weierstrass block (fused silica, glued) + 10 um gap + 250 um
    aplanat of index 1.8."""
    return C.fac_design(N_SIO2, n2=N_HI)


class CombY:
    """1D field at S of one chip, with an optional comb window at distance z from S.

    design None: ideal thin sine-condition collimator in air (aperture f_FAC before S). design given: the traced
    collimator (facet in its first medium), exit field at the plane d.z_out from its lens map, S at the nominal back
    focal plane; sx > 0: the same lens seen by light of x direction sine sx (effective indices)."""

    def __init__(self, apod=None, clip_at_S=False, design=None, sx=0.0):
        self.sx = sx
        if design is None:
            self.fa = FacetY(LAM, 1.2, 0.4)
            self.r = RelayY(self.fa, BETA, NA=NA, NA_FAC=NA_FAC, apod=apod)
            self.z_ap = 0.0 if clip_at_S else -F_FAC  # exit-field plane relative to S
            n_air = 1.0
        else:
            dd = skewed(design, sx)
            self.L = dd.lensmap(LAM, NA_FAC, coated=True)
            self.fa = FacetY(LAM, 1.2, 0.4, n_m=dd.indices[0])
            self.r = RelayY(self.fa, BETA, NA=NA, NA_FAC=NA_FAC, apod=apod, lens=self.L)
            self.z_ap = design.z_out - _bfp(design)
            n_air = dd.indices[-1]
        self.y = (np.arange(NY) - NY // 2) * DY
        s, amp, kern = self.r.spectrum()
        if design is None:
            dg = np.ones_like(s)
        else:
            dg = np.interp(np.abs(s), self.L.s, self.L.dg)
        U_s = amp / np.sqrt(dg) / np.sqrt(dg)  # amp carries sqrt(dg) for the s integral; U dy = amp ds, dy = f dg ds
        yk = F_FAC * kern
        o = np.argsort(yk)
        U = np.interp(self.y, yk[o], U_s[o].real, left=0, right=0) + 1j * np.interp(self.y, yk[o], U_s[o].imag,
                                                                                     left=0, right=0)
        U[(self.y < yk.min()) | (self.y > yk.max())] = 0
        self.U_exit = U
        self.clip_at_S = clip_at_S
        self.k = 2 * np.pi / LAM
        k_yz = self.k * n_air  # projected wavenumber in air for x direction sine sx
        self.kz = np.sqrt((k_yz**2 - (2 * np.pi * np.fft.fftfreq(NY, DY)) ** 2).astype(complex))
        self.eta = self.r.eta
        self.Py_ref = self.r.pupil()  # RelayY reference
        self._scale = None

    def prop(self, U, z):
        """Angular-spectrum propagation by z (um); evanescent components are dropped (they would grow when z < 0)."""
        H = np.where(self.kz.imag == 0, np.exp(1j * self.kz.real * z), 0)
        return np.fft.ifft(np.fft.fft(U) * H)

    def field_S(self, z=None, a=None, mode="slot"):
        """Field at S. z: comb distance from S (um); a: comb half-aperture (slot) or edge position (um)."""
        z0 = self.z_ap
        if z is None or a is None:
            return self.prop(self.U_exit, -z0)
        Uc = self.prop(self.U_exit, z - z0)
        W = np.abs(self.y) < a if mode == "slot" else self.y < a
        return self.prop(Uc * W, -z)

    def pupil(self, US):
        """P(eta) = sum U_S(y) exp(i k NA eta y / D) dy over the central 4 mm, scaled to the RelayY normalisation."""
        m = np.abs(self.y) < 2000.0
        ph = np.exp(1j * self.k * NA / D * np.outer(self.eta, self.y[m]))
        P = ph @ US[m] * DY
        if self._scale is None:  # fix the constant against RelayY (exit field placed in S, as RelayY does)
            ref = self.U_exit
            Pr = ph @ ref[m] * DY
            self._scale = complex(np.vdot(Pr, self.Py_ref) / np.vdot(Pr, Pr))
            self.scale_check = float(np.max(np.abs(Pr * self._scale - self.Py_ref)) / np.max(np.abs(self.Py_ref)))
        return P * self._scale

    def metrics(self, Py, ya=np.arange(-12.0, 12.0001, 0.02)):
        I = self.r.atom_line(ya, Py=Py)
        i0 = np.argmax(I)
        w = I > 0.5 * I.max()
        cen = float(np.sum(ya[w] * I[w]) / np.sum(I[w]))
        return {"transmission": self.r.efficiency(Py=Py), "peak": float(I.max()), "centroid_nm": 1e3 * cen,
                "xt_3um": eps_at(ya, I, 3.0), "xt_4um": eps_at(ya, I, 4.0), "xt_8um": eps_at(ya, I, 8.0),
                "y_peak_um": float(ya[i0])}


def edge():
    cases = {}
    rows = []
    for apod_name, apod in (("hard", None), ("graded 0.8-0.9", raised_cosine(0.8, 0.9))):
        cy = CombY(apod=apod)
        ref_S = cy.metrics(cy.pupil(cy.prop(cy.U_exit, 0.0)))  # clip in S (RelayY convention)
        base = cy.metrics(cy.pupil(cy.field_S()))  # clip at the FAC exit, no comb
        relay = cy.metrics(cy.Py_ref)
        cases[apod_name] = {"relayY": relay, "fft_clip_in_S": ref_S, "no_comb": base, "scale_check": cy.scale_check}
        print(apod_name, "RelayY", {k: f"{v:.4g}" for k, v in relay.items()}, flush=True)
        print(apod_name, "FFT clip in S", {k: f"{v:.4g}" for k, v in ref_S.items()}, "scale err", f"{cy.scale_check:.1e}")
        print(apod_name, "no comb (clip at FAC exit)", {k: f"{v:.4g}" for k, v in base.items()}, flush=True)
        for mode in ("slot", "edge"):
            for z_mm in (1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 45.0, 60.0):
                for a in (100.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0):
                    m = cy.metrics(cy.pupil(cy.field_S(z_mm * 1e3, a, mode)))
                    r = {"apod": apod_name, "mode": mode, "z_mm": z_mm, "a_um": a,
                         "loss": 1 - m["transmission"] / base["transmission"], "peak_rel": m["peak"] / base["peak"],
                         "xt_3um": m["xt_3um"], "xt_4um": m["xt_4um"], "xt_8um": m["xt_8um"],
                         "row_shift_nm": m["centroid_nm"] - base["centroid_nm"]}
                    rows.append(r)
                    if a in (200.0, 250.0):
                        print(apod_name, mode, z_mm, a, "| loss %.2e peak %.4f xt4 %.2e xt3 %.2e shift %.2f nm" % (
                            r["loss"], r["peak_rel"], r["xt_4um"], r["xt_3um"], r["row_shift_nm"]), flush=True)
    # usable window: largest z and smallest a with loss <= 1% and xt_4um <= 2e-4 (hard edge, slot)
    ok = {}
    for apod_name in cases:
        for mode in ("slot", "edge"):
            for z_mm in sorted({r["z_mm"] for r in rows}):
                good = [r["a_um"] for r in rows if r["apod"] == apod_name and r["mode"] == mode and r["z_mm"] == z_mm
                        and r["loss"] <= 0.01 and r["xt_4um"] <= 2e-4]
                ok[f"{apod_name} | {mode} | z {z_mm:g} mm"] = min(good) if good else None
    # analytic Fresnel tail of the FAC edge at the comb: fraction of power beyond distance g from the geometric edge
    cy = CombY()
    I_exit = np.abs(cy.U_exit) ** 2
    edge_rel = float(I_exit[np.argmin(np.abs(cy.y - 0.99 * F_FAC * NA_FAC))] / I_exit.max())
    P_tot = float(I_exit.sum() * DY)
    # far-shadow Fresnel tail of a straight edge: intensity I_e lambda L / (4 pi^2 x^2) at distance x into the
    # shadow, i.e. a power I_e lambda L / (4 pi^2 g) per edge beyond g; two edges
    analytic = {f"z {z:g} mm, a {a:g} um": 2 * edge_rel * I_exit.max() * LAM * (z * 1e3 + F_FAC) /
                (4 * np.pi**2 * (a - F_FAC * NA_FAC)) / P_tot
                for z in (10.0, 30.0, 60.0) for a in (200.0, 250.0)}
    f, ax = fig(1, 2, h=3.0)
    for a, ls in ((200.0, "--"), (250.0, "-")):
        s = [r for r in rows if r["apod"] == "hard" and r["mode"] == "slot" and r["a_um"] == a]
        ax[0, 0].loglog([r["z_mm"] for r in s], [max(r["loss"], 1e-7) for r in s], ls, marker="o", ms=3,
                        label=f"loss, a = {a:g} um")
        ax[0, 1].semilogx([r["z_mm"] for r in s], [r["xt_4um"] for r in s], ls, marker="o", ms=3, label=f"a = {a:g} um")
    ax[0, 1].axhline(cases["hard"]["no_comb"]["xt_4um"], color="0.5", lw=0.7)
    ax[0, 0].set_xlabel("comb distance from S (mm)")
    ax[0, 0].set_ylabel("power clipped by the comb")
    ax[0, 1].set_xlabel("comb distance from S (mm)")
    ax[0, 1].set_ylabel("crosstalk at 4 um")
    ax[0, 0].legend(fontsize=7)
    ax[0, 1].legend(fontsize=7)
    save_result("m26_edge", script=SCRIPT, fig=f,
                inputs={"beta_um": BETA, "NA_FAC": NA_FAC, "f_FAC_um": F_FAC, "D": D, "row_pitch_um": DELTA,
                        "chip_pitch_at_S_um": D_CHIP, "grid_um": DY, "span_um": NY * DY,
                        "aperture_plane": "FAC exit, f_FAC before S (thin sine-condition lens, CHOICE)",
                        "window": "slot: |y| < a (interleaving comb, both edges at the same z); edge: y < a (one "
                                  "edge, block combining)"},
                outputs={"baselines": cases, "rows": rows, "usable_min_half_aperture_um": ok,
                         "analytic_tail_fraction": analytic, "edge_intensity_rel": edge_rel},
                notes="loss and peak_rel relative to the same chip without a comb (aperture at the FAC exit); "
                      "xt_*: largest atom-plane intensity within +/-0.5 um of that distance, relative to the peak "
                      "(1D, single tweezer, as m13 fill); row_shift_nm: centroid of the > 50% part of the profile.")
    return rows, ok


def _row_metrics(F, ls, Py_m, Py_n, deltas=(3.0, 4.0)):
    """m13 'crosstalk' metrics with given y pupils for chip m (row at 0) and chip m+1 (row at +delta)."""
    xt = np.arange(-60.0, 60.0 + 1e-9, 3.0)
    th = np.where(np.arange(len(xt)) % 2 == 0, 0.0, np.pi / 2)
    w = ls.weights(xt, phases=th)
    x = np.arange(-63, 63, 0.04)
    y = np.arange(-1.6, 9.0, 0.04)
    I = np.abs(F.field(w, x, y, Py=Py_m)) ** 2
    pk = np.mean([peak2d(I[:, np.abs(x - x0) < 0.6], x[np.abs(x - x0) < 0.6], y)[0] for x0 in xt])
    prof = I.max(axis=1) / pk
    out = {"row_peak_mean": float(pk),
           "eps_row": {str(d): float(prof[np.abs(np.abs(y) - d) <= 0.5].max()) for d in (3.0, 4.0, 6.0, 8.0)}}
    for delta in deltas:
        yy = np.arange(delta - 1.2, delta + 1.2, 0.02)
        xs_n = xt[15:26] + 1.5
        xx = np.arange(xs_n[0] - 1.5, xs_n[-1] + 1.5, 0.02)
        En = F.field(ls.weights(xs_n, phases=np.where(np.arange(len(xs_n)) % 2 == 0, 0.0, np.pi / 2)), xx, yy - delta,
                     Py=Py_n)
        Em = F.field(w, xx, yy, Py=Py_m)
        xs_m2 = xt.copy()
        xs_m2[20] += 0.25
        Em2 = F.field(ls.weights(xs_m2, phases=th), xx, yy, Py=Py_m)

        def tw(Iv):
            return np.array([peak2d(Iv[:, np.abs(xx - x0) < 0.6], xx[np.abs(xx - x0) < 0.6], yy) for x0 in xs_n])

        base = tw(np.abs(En) ** 2)
        rows = {}
        for mode in ("incoherent", "coherent"):
            dpk = dpos = dmk = dmp = 0.0
            for ph in ([0.0] if mode == "incoherent" else np.linspace(0, 2 * np.pi, 8, endpoint=False)):
                if mode == "incoherent":
                    Ia, Ib = np.abs(En) ** 2 + np.abs(Em) ** 2, np.abs(En) ** 2 + np.abs(Em2) ** 2
                else:
                    Ia, Ib = np.abs(En + np.exp(1j * ph) * Em) ** 2, np.abs(En + np.exp(1j * ph) * Em2) ** 2
                a, b = tw(Ia), tw(Ib)
                dpk = max(dpk, np.max(np.abs(a[:, 0] / base[:, 0] - 1)))
                dpos = max(dpos, np.max(np.hypot(a[:, 1] - base[:, 1], a[:, 2] - base[:, 2])))
                dmk = max(dmk, np.max(np.abs(b[:, 0] / a[:, 0] - 1)))
                dmp = max(dmp, np.max(np.hypot(b[:, 1] - a[:, 1], b[:, 2] - a[:, 2])))
            rows[mode] = {"static_peak_change": float(dpk), "static_pos_nm": float(1e3 * dpos),
                          "move_peak_change": float(dmk), "move_pos_nm": float(1e3 * dmp)}
        out[f"neighbour_{delta:g}um"] = rows
    return out


def twod():
    fp, ls, r, F = setup()
    cy = CombY()
    P0 = r.pupil()
    cases = {"no comb": (P0, P0)}
    for z_mm, a in ((23.0, 250.0), (23.0, 225.0), (15.0, 250.0), (30.0, 250.0)):
        P = cy.pupil(cy.field_S(z_mm * 1e3, a, "slot"))
        cases[f"slot z {z_mm:g} mm a {a:g} um"] = (P, P)  # chip m (slot) and m+1 (tooth) see the same window
    # block combining: the boundary between chip m (edge above it) and chip m+1 (edge below it, mirror image)
    Pe = cy.pupil(cy.field_S(23e3, 250.0, "edge"))
    cases["edge z 23 mm a 250 um"] = (Pe, Pe[::-1])
    res = {}
    w0, w60 = unit(ls.weights([0.0])), unit(ls.weights([60.0]))
    m00, m60 = metrics(F, w0, 0.0, Py=P0), metrics(F, w60, 60.0, Py=P0)
    for name, (Pm, Pn) in cases.items():
        s0, s60 = metrics(F, w0, 0.0, Py=Pm), metrics(F, w60, 60.0, Py=Pm)
        row = {"single_x0": {"peak_rel": s0["peak"] / m00["peak"], "y_nm": s0["y_nm"], "fwhm_y": s0["fwhm_y"],
                             "eps3": s0["eps3"], "eps4": s0["eps4"], "eps6": s0["eps6"]},
               "single_x60": {"peak_rel": s60["peak"] / m60["peak"], "y_nm": s60["y_nm"], "eps4": s60["eps4"]}}
        row.update(_row_metrics(F, ls, Pm, Pn))
        res[name] = row
        print(name, "| single x0: peak %.4f eps4 %.2e y %.1f nm | x60: peak %.4f eps4 %.2e | row eps3/4/8 %.2e %.2e %.2e" % (
            row["single_x0"]["peak_rel"], row["single_x0"]["eps4"], row["single_x0"]["y_nm"],
            row["single_x60"]["peak_rel"], row["single_x60"]["eps4"], row["eps_row"]["3.0"], row["eps_row"]["4.0"],
            row["eps_row"]["8.0"]), flush=True)
        for k in ("neighbour_3um", "neighbour_4um"):
            print("   ", k, {m: {kk: float("%.3g" % vv) for kk, vv in v.items()} for m, v in row[k].items()}, flush=True)
    # local-window sensitivity for comb roll / edge straightness (single tweezer at x = 0, 2D)
    sens = []
    for a in (205.0, 215.0, 225.0, 235.0, 250.0):
        P = cy.pupil(cy.field_S(23e3, a, "slot"))
        s0 = metrics(F, w0, 0.0, Py=P)
        sens.append({"a_um": a, "peak_rel": s0["peak"] / m00["peak"], "eps4": s0["eps4"], "eps3": s0["eps3"]})
        print("a", a, sens[-1], flush=True)
    # Sigma-space footprint: chief rays parallel beyond the slit plane (telecentric), X = f1 s; one tweezer's
    # geometric half-width at L beyond Sigma is L (A/2) / f1; four engines at 3.0 mm pitch with slit zones +/- f1 s_cut
    f1 = 3120.0
    A = fp.opa.A
    X_half_4eng = 1.5 * 3000.0 + f1 * fp.opa.s_cut
    geo = {"f1_um": f1, "emitter_aperture_um": A, "sigma_from_S_mm": (2 * f1 - 2 * F_FAC) * 1e-3,  # thin lenses
           "tweezer_half_footprint_um_at_17mm_beyond_sigma": 17e3 * (A / 2) / f1,
           "X_half_span_4_engines_um": X_half_4eng,
           "roll_mrad_for_10um": 10.0 / X_half_4eng * 1e3,
           "edge_change_across_one_tweezer_um_at_2mrad": 2e-3 * 2 * 17e3 * (A / 2) / f1}
    print(geo, flush=True)
    save_result("m26_twod", script=SCRIPT,
                inputs={"model": "m13 Focus2D (D1 BPM x response, circular NA 0.5 pupil), y pupil from CombY",
                        "row_m": "41 tweezers at 3 um, phases 0 / pi/2", "row_m+1": "11 tweezers offset 1.5 um",
                        "move": "one row-m tweezer moves 0.25 um", "coherent": "worst of 8 relative phases"},
                outputs={"cases": res, "local_window_sensitivity_z23mm": sens, "geometry": geo},
                notes="slot: both chips windowed by |y| < a at the same z (interleaving comb); edge: chip m windowed "
                      "by y < a, chip m+1 by the mirror image (block combining boundary). Roll and straightness use "
                      "the local window a(X) at each tweezer's chief-ray position X = f1 s (telecentric Sigma space).")
    return res


def _sigma_from_facet_um(d):
    """x slit plane from the facet: back focal distance of the traced x lens (results/m14_xlens.json, aplanat SiO2,
    measured from the emitter plane in air) plus the x-focus shift of the collimator glass, t (1 - 1/n) per section
    on axis (the collimator has no x power)."""
    import json
    zb = json.load(open(RESULTS / "m14_xlens.json"))["outputs"]["aplanat SiO2"]["back_focal_z_um"]
    if callable(getattr(d, "sigma_z", None)):  # folded reflective collimator: its own unfolded glass path
        sig = d.sigma_z(zb)
        return sig, zb, sig - zb
    zs = [0.0] + [sf.z_v for sf in d.surfaces]
    shift = sum((zs[j + 1] - zs[j]) * (1 - 1 / n) for j, n in enumerate(d.indices[:-1]) if n > 1)
    return zb + shift, zb, shift


def lens():
    """The comb with the traced collimator (on axis, sx = 0): baseline, scan, worst row in 2D, real geometry."""
    d = fac_adopted()
    cy = CombY(design=d)
    base = cy.metrics(cy.pupil(cy.field_S()))
    relay = cy.metrics(cy.Py_ref)
    thin = CombY()
    tbase = thin.metrics(thin.pupil(thin.field_S()))
    zS = _bfp(d)
    sig, zb_x, shift = _sigma_from_facet_um(d)
    z_c = (sig - zS) * 1e-3 + 0.5  # mm, comb start 0.5 mm beyond the slit plane (CHOICE)
    geo = {"S_from_facet_um": zS, "exit_plane_from_facet_um": d.z_out, "x_lens_back_focal_um": zb_x,
           "collimator_glass_x_shift_um": shift, "sigma_from_facet_um": sig, "sigma_from_S_mm": (sig - zS) * 1e-3,
           "z_c_mm": z_c, "worst_row_z_mm_32_rows": z_c + 16.0}
    print("geometry", geo, flush=True)
    print("baseline (lens)", {k: float("%.4g" % v) for k, v in base.items()}, "RelayY xt4 %.3g" % relay["xt_4um"], flush=True)
    rows = []
    for mode in ("slot", "edge"):
        for z_mm in (5.0, 10.0, 15.0, 20.0, round(z_c + 16.0, 2), 30.0, 45.0):
            for a in (200.0, 225.0, 250.0, 275.0):
                out = {"mode": mode, "z_mm": z_mm, "a_um": a}
                for name, c, b in (("lens", cy, base), ("thin", thin, tbase)):
                    m = c.metrics(c.pupil(c.field_S(z_mm * 1e3, a, mode)))
                    out[name] = {"loss": 1 - m["transmission"] / b["transmission"], "peak_rel": m["peak"] / b["peak"],
                                 "xt_4um": m["xt_4um"], "xt_3um": m["xt_3um"],
                                 "row_shift_nm": m["centroid_nm"] - b["centroid_nm"]}
                rows.append(out)
                if a in (225.0, 250.0):
                    print(mode, z_mm, a, "| lens loss %.2e peak %.4f xt4 %.2e | thin loss %.2e peak %.4f xt4 %.2e" % (
                        out["lens"]["loss"], out["lens"]["peak_rel"], out["lens"]["xt_4um"], out["thin"]["loss"],
                        out["thin"]["peak_rel"], out["thin"]["xt_4um"]), flush=True)
    # 2D at the worst row (32 rows comb, or 64 rows comb + PBS)
    fp, ls, _, _ = setup()
    F = Focus2D(fp, cy.r)
    zw = (z_c + 16.0) * 1e3
    P0 = cy.Py_ref
    cases = {"no comb": (P0, P0)}
    for a in (250.0, 225.0):
        P = cy.pupil(cy.field_S(zw, a, "slot"))
        cases[f"slot z {zw / 1e3:.1f} mm a {a:g} um"] = (P, P)
    w0 = unit(ls.weights([0.0]))
    m0 = metrics(F, w0, 0.0, Py=P0)
    twod_rows = {}
    for name, (Pm, Pn) in cases.items():
        s0 = metrics(F, w0, 0.0, Py=Pm)
        row = {"single_peak_rel": s0["peak"] / m0["peak"], "single_eps4": s0["eps4"], "single_eps3": s0["eps3"]}
        row.update(_row_metrics(F, ls, Pm, Pn))
        twod_rows[name] = row
        print("2D", name, "| peak %.4f eps4 %.2e | row eps4 %.2e |" % (row["single_peak_rel"], row["single_eps4"],
              row["eps_row"]["4.0"]), {m: {k: float("%.3g" % v) for k, v in vv.items()}
                                       for m, vv in row["neighbour_4um"].items()}, flush=True)
    save_result("m26_lens", script=SCRIPT,
                inputs={"collimator": "collimator.fac_design(n_SiO2, n2=1.8): Weierstrass block glued to the facet, "
                                      "10 um gap, 250 um aplanat; AR coated; facet radiating into fused silica",
                        "exit_plane": "d.z_out (60 um beyond the last vertex)", "S": "collimator.back_focal_plane",
                        "comb_start": "0.5 mm beyond the x slit plane (CHOICE)"},
                outputs={"baseline_lens": base, "baseline_relayY_lens": relay, "baseline_thin": tbase,
                         "geometry": geo, "rows": rows, "twod_worst_row": twod_rows})
    return rows


def _best_focus_at(F, w, xa, zs=np.arange(-3.0, 3.01, 0.25)):
    """Atom-plane defocus (um) that maximises the peak of a tweezer at xa (2D, window around xa)."""
    x = np.arange(xa - 1.0, xa + 1.0, 0.02)
    y = np.arange(-3.0, 3.0, 0.02)
    pk = [F.intensity(w, x, y, dz_a=z).max() for z in zs]
    i = int(np.argmax(pk))
    return float(zs[i]), float(pk[i])


def _profile_2d(F, w, xa):
    """y profile (max over |x - xa| < 1 um) at the atoms, for diagnosis of strongly aberrated cases."""
    x = np.arange(xa - 1.0, xa + 1.0, 0.04)
    y = np.arange(-8.0, 8.0, 0.04)
    I = F.intensity(w, x, y)
    return y, I.max(axis=1)


def skew(design=None, name="m26_skew", xas=None):
    """x dependence: the cylindrical collimator for tweezers at x_a (x direction sine sx = x_a / x_scale).
    design: a traced collimator (default: the adopted refractive one); name: result file."""
    d = fac_adopted() if design is None else design
    fp, ls, _, _ = setup()
    x_scale, s_cut = fp.opa.x_scale, fp.opa.s_cut
    ya = np.arange(-12.0, 12.0001, 0.01)
    sig, _, _ = _sigma_from_facet_um(d)
    zw = ((sig - _bfp(d)) * 1e-3 + 0.5 + 16.0) * 1e3
    rows = []
    dzs = np.arange(-8.0, 8.001, 0.25)
    peak_vs_dz = {}
    xas = (0.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 75.0, 90.0, s_cut * x_scale) if xas is None else xas
    for xa in xas:
        sx = xa / x_scale
        cy = CombY(design=d, sx=sx)
        r = cy.r
        ideal = RelayY(cy.fa, BETA, NA=NA, NA_FAC=NA_FAC)
        I_id = ideal.atom_line(ya)
        I = r.atom_line(ya)
        L = cy.L
        nm = cy.fa.n_m
        basis = (nm * np.sqrt(1 - (L.s / nm) ** 2) - nm) / LAM  # W of a facet displaced by +1 um (first order)
        dz_eq = float(-np.sum(L.W * basis) / np.sum(basis**2))
        resid = L.W + dz_eq * basis
        I_fix = r.atom_line(ya, dz=dz_eq)  # facet refocused for this tweezer alone
        base = cy.metrics(cy.pupil(cy.field_S()))
        comb = cy.metrics(cy.pupil(cy.field_S(zw, 250.0, "slot")))
        F = Focus2D(fp, r)
        Fi = Focus2D(fp, ideal)
        w = unit(ls.weights([xa]))
        m2 = metrics(F, w, xa)
        m2i = metrics(Fi, w, xa)
        zb, pb = _best_focus_at(F, w, xa)
        row = {"x_a_um": xa, "sx": sx, "peak_rel_ideal": float(I.max() / I_id.max()), "dz_equiv_um": dz_eq,
               "W_pv_waves": float(np.ptp(L.W)), "W_resid_after_focus_pv_waves": float(np.ptp(resid)),
               "peak_rel_ideal_refocused_facet": float(I_fix.max() / I_id.max()),
               "xt_3um": eps_at(ya, I, 3.0), "xt_4um": eps_at(ya, I, 4.0), "centroid_nm": base["centroid_nm"],
               "s_max": float(L.s_max), "g_edge_over_s": float(L.g[-1] / L.s[-1]),
               "twod_peak_rel_ideal": m2["peak"] / m2i["peak"], "twod_eps4": m2["eps4"],
               "twod_best_focus_um": zb, "twod_peak_at_best_focus_rel_ideal": pb / m2i["peak"],
               "comb_worst_row": {"loss": 1 - comb["transmission"] / base["transmission"],
                                  "peak_rel": comb["peak"] / base["peak"], "xt_4um": comb["xt_4um"]}}
        rows.append(row)
        peak_vs_dz[f"{xa:.1f}"] = [float(r.atom_line(ya, dz=dz).max() / I_id.max()) for dz in dzs]
        print("x_a %5.1f sx %.3f | peak %.3f (2D %.3f, best focus %.1f um -> %.3f) | dz_eq %+.2f um, W pv %.2f, resid %.3f "
              "| refocused %.3f | xt4 %.1e (2D %.1e) | comb: loss %.1e xt4 %.1e" % (
                  xa, sx, row["peak_rel_ideal"], row["twod_peak_rel_ideal"], zb, row["twod_peak_at_best_focus_rel_ideal"],
                  dz_eq, row["W_pv_waves"], row["W_resid_after_focus_pv_waves"], row["peak_rel_ideal_refocused_facet"],
                  row["xt_4um"], row["twod_eps4"], row["comb_worst_row"]["loss"], row["comb_worst_row"]["xt_4um"]),
              flush=True)
    sxs = np.array([r["sx"] for r in rows])
    dzq = np.array([r["dz_equiv_um"] for r in rows])
    c = float(np.sum(dzq * sxs**2) / np.sum(sxs**4))  # dz_equiv = c sx^2
    # one common facet offset for the whole row: maximise the worst peak over |x_a| <= 90 um
    keys = [k for k in peak_vs_dz if float(k) <= 90.0 + 1e-6]
    worst = [min(peak_vs_dz[k][j] for k in keys) for j in range(len(dzs))]
    jb = int(np.argmax(worst))
    usable = {}
    for crit in (0.97, 0.9):
        ok = [r["x_a_um"] for r in rows if r["peak_rel_ideal"] >= crit]
        usable[f"peak >= {crit}"] = max(ok)
    ok = [r["x_a_um"] for r in rows if r["xt_4um"] <= 2e-4]
    usable["xt_4um <= 2e-4"] = max(ok)
    n = N_SIO2
    thin_c = -F_FAC * (n + 1) / (2 * n)  # thin glass lens in air: delta f / f = -(n + 1) sx^2 / (2 n)
    need = {"tolerance_um": 0.5, "c_max_for_90um_um": 0.5 / (90.0 / x_scale) ** 2,
            "c_max_for_90um_with_balance_um": 1.0 / (90.0 / x_scale) ** 2}
    out = {"x_scale_um": x_scale, "s_cut": s_cut, "rows": rows, "dz_equiv_coefficient_um": c,
           "thin_lens_air_estimate_um": thin_c, "usable_x_a_um": usable,
           "balance": {"facet_offset_um": float(dzs[jb]), "worst_peak_rel_ideal": worst[jb],
                       "peak_vs_offset": {"dz_um": dzs.tolist(), **peak_vs_dz}},
           "requirement": need}
    print({k: v for k, v in out.items() if k not in ("rows", "balance")}, "balance", dzs[jb], worst[jb], flush=True)
    f, ax = fig(1, 2, h=3.0)
    ax[0, 0].plot([r["x_a_um"] for r in rows], [r["peak_rel_ideal"] for r in rows], "o-", ms=3, label="1D")
    ax[0, 0].plot([r["x_a_um"] for r in rows], [r["twod_peak_rel_ideal"] for r in rows], "s--", ms=3, label="2D")
    ax[0, 0].set_xlabel("tweezer position x_a (um)")
    ax[0, 0].set_ylabel("peak / ideal collimator")
    ax[0, 0].legend(fontsize=7)
    ax[0, 1].semilogy([r["x_a_um"] for r in rows], [r["xt_4um"] for r in rows], "o-", ms=3)
    ax[0, 1].axhline(2e-4, color="0.5", lw=0.7)
    ax[0, 1].set_xlabel("tweezer position x_a (um)")
    ax[0, 1].set_ylabel("crosstalk at 4 um")
    save_result(name, script=SCRIPT, fig=f,
                inputs={"collimator": "adopted (m14)" if design is None else d.name, "method": "cylindrical surfaces conserve k_x: the y-z problem is "
                        "the 2D trace with effective indices sqrt(n^2 - sx^2) (exact); facet spectrum FacetY with the "
                        "effective source index; AR coated; downstream y relay ideal for every sx",
                        "sx": "x_a / x_scale (D1 slit relay)", "comb": f"slot a = 250 um at z = {zw / 1e3:.1f} mm"},
                outputs=out,
                notes="peak_rel_ideal: y-line peak relative to an ideal sine-condition collimator with the same facet "
                      "spectrum and sx. dz_equiv: facet focus error (in the source medium) whose first-order wavefront "
                      "best matches the traced W(s). refocused: the same tweezer with the facet moved by dz_equiv.")
    return out


def skew_offset():
    """A common facet offset (block length) for the whole row: peak and crosstalk vs x_a (1D)."""
    d = fac_adopted()
    ya = np.arange(-12.0, 12.0001, 0.01)
    out = {}
    for off in (0.0, -0.5, -0.75, -1.0, -1.25, -1.5):
        rows = []
        for xa in (0.0, 15.0, 30.0, 40.0, 45.0, 50.0, 60.0):
            cy = CombY(design=d, sx=xa / 288.0)
            I = cy.r.atom_line(ya, dz=off)
            Ii = RelayY(cy.fa, BETA, NA=NA, NA_FAC=NA_FAC).atom_line(ya)
            rows.append({"x_a_um": xa, "peak_rel_ideal": float(I.max() / Ii.max()), "xt_4um": eps_at(ya, I, 4.0),
                         "xt_3um": eps_at(ya, I, 3.0)})
        out[f"{off:g}"] = rows
        print(off, " | ".join("x %g: %.3f, %.1e" % (r["x_a_um"], r["peak_rel_ideal"], r["xt_4um"]) for r in rows), flush=True)
    save_result("m26_skew_offset", script=SCRIPT,
                inputs={"collimator": "adopted (m14)", "x_scale_um": 288.0,
                        "offset": "facet moved along z by the same amount for every tweezer (first-order focus phase, "
                                  "RelayY dz, um in the source medium)"},
                outputs={"rows": out})
    return out


def z_c_real():
    """Nearest comb point from S (mm): 0.5 mm (CHOICE) beyond the x slit plane of the traced x lens and collimator."""
    d = fac_adopted()
    sig, _, _ = _sigma_from_facet_um(d)
    return (sig - _bfp(d)) * 1e-3 + 0.5


def z_max_rule(a):
    """Largest comb distance (mm) at which the in-band beam (half-angle NA/D at S) still fits the window a (um)."""
    return (a - F_FAC * NA_FAC) * D / NA * 1e-3


def layout():
    """First-order layouts (DERIVED) for 32 / 48 / 64 rows at 4 um, and the worst-row comb penalty (SIM)."""
    Z_C = z_c_real()
    F_y = design(DELTA, 16)["F_y_mm"]
    pupil_r = F_OBJ * NA  # mm
    cy = CombY()
    base = cy.metrics(cy.pupil(cy.field_S()))

    def worst(z_mm, a):
        m = cy.metrics(cy.pupil(cy.field_S(z_mm * 1e3, a)))
        return {"loss": 1 - m["transmission"] / base["transmission"], "peak_rel": m["peak"] / base["peak"],
                "xt_4um": m["xt_4um"], "xt_3um": m["xt_3um"]}

    out = {"F_y_mm": F_y, "pupil_radius_mm": pupil_r, "z_c_mm": Z_C, "baseline": base,
           "z_max_rule_mm": {f"a {a:g} um": z_max_rule(a) for a in (200.0, 225.0, 250.0)}, "layouts": []}
    for n in (32, 48, 64):
        d = design(DELTA, n)
        m_t, f_2 = d["y_telescope"]["m"], d["f_2_mm"]
        half_S = n * D_CHIP * 1e-3 / 2
        for name, n_st, pitch, comb_rows, n_pbs in (("single stack", 1, 0.5, 0, 0),
                                                     ("2 stacks, PBS", 2, 0.5, 0, 1),
                                                     ("2 stacks, comb", 2, 1.0, n, 0),
                                                     ("4 stacks, comb + PBS", 4, 1.0, n // 2, 1)):
            if n % n_st:
                continue
            row = {"rows": n, "layout": name, "stacks": n_st, "chips_per_stack": n // n_st, "chip_pitch_mm": pitch,
                   "stack_height_mm": n // n_st * pitch, "pbs": n_pbs, "combs": (2 if n_pbs else 1) if comb_rows else 0,
                   "field_half_at_f_a_mm": half_S, "field_angle_at_pupil_mrad": 1e3 * half_S / m_t / f_2,
                   "block_height_at_atoms_um": n * DELTA, "polarisation_insensitive_objective": bool(n_pbs)}
            if comb_rows:
                zr = (Z_C, Z_C + comb_rows * D_CHIP * 1e-3)
                row["comb_z_range_mm"] = zr
                row["worst_row_a250"] = worst(zr[1], 250.0)
                row["worst_row_a225"] = worst(zr[1], 225.0)
            out["layouts"].append(row)
            print(row, flush=True)
    # comb tilt about x: deflection 2 theta at S -> pupil shift F_y 2 theta. The m13 FAC decentre delta shifts the
    # facet image by delta / beta pupil radii, with peak 0.992 / 0.969 / 0.881 at delta = 0.05 / 0.1 / 0.2 um
    # (results/m13_tolerance.json)
    shift = np.array([0.0, 0.05, 0.1, 0.2]) / BETA
    peak = np.array([1.0, 0.992, 0.969, 0.881])
    tol = {}
    for loss in (0.01, 0.03):
        sh = float(np.interp(loss, 1 - peak, shift))
        theta = sh * pupil_r / (2 * F_y)
        tol[f"{loss:.0%} loss"] = {"pupil_shift_radii": sh, "comb_tilt_mrad": 1e3 * theta,
                                   "row_shift_nm_at_z_20mm": 1e3 * 2 * theta * 20e3 / D}
    out["comb_tilt_about_x"] = tol
    out["row_shift_nm_per_urad_at_z_20mm"] = 2e-6 * 20e3 * 1e3 / D
    out["path_equality"] = ("45 deg comb along z = z_c + y: stack B (facets in the plane y = -z_c, beams along +y) has "
                            "the path Z - z_c - y_B = Z to any plane z = Z for every row, equal to stack A's. Moving the "
                            "comb by dz along z changes B's path by dz (m13 piston: 100 um keeps the peak >= 0.992) and "
                            "shifts every window by dz in y (tolerance from the a scan)")
    out["io_per_chip"] = {"fast_actuators_per_engine": 383, "4 engines": 4 * 383, "8 engines": 8 * 383}
    save_result("m26_layout", script=SCRIPT,
                inputs={"design": "scripts/m13_stack.design (thin lenses, 25x y telescope, f_a 100 mm)",
                        "z_c_mm": Z_C, "m13_tolerance_peaks": peak.tolist(), "m13_decentres_um": [0, 0.05, 0.1, 0.2]},
                outputs=out,
                notes="worst_row: the row at the far end of the comb (largest z), slot half-aperture a; loss, peak_rel "
                      "relative to no comb. z_max_rule: (a - f_FAC NA_FAC) D / NA, the distance at which the light the "
                      "pupil accepts (half-angle NA/D at S) reaches the window edge.")
    return out


if __name__ == "__main__":
    parts = sys.argv[1:] or ["edge", "layout"]
    for p in parts:
        {"edge": edge, "layout": layout, "twod": twod, "lens": lens, "skew": skew, "skew_offset": skew_offset}[p]()
