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
  layout  first-order layouts for 32 / 48 / 64 rows (m13 design), comb z-range, path equality, tilt tolerances,
          chip pitch and I/O per chip.

Output: results/m26_edge.json (+ .png), results/m26_layout.json
"""
from __future__ import annotations

import sys

import numpy as np

from picsim.opa.stack import FacetY, RelayY, raised_cosine
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

sys.path.insert(0, str(RESULTS.parent / "scripts"))
from m13_stack import BETA, D_CHIP, F_OBJ, LAM, NA, NA_FAC, design, eps_at  # noqa: E402

SCRIPT = "scripts/m26_stack_combine.py"
DELTA = 4.0  # um, row pitch at the atoms
D = D_CHIP / DELTA
F_FAC = BETA * D / NA  # um
DY = 0.5  # um, grid at S
NY = 2**16  # 32.8 mm span


class CombY:
    """1D field at S of one chip, with an optional comb window at distance z from S."""

    def __init__(self, apod=None, clip_at_S=False):
        self.fa = FacetY(LAM, 1.2, 0.4)
        self.r = RelayY(self.fa, BETA, NA=NA, NA_FAC=NA_FAC, apod=apod)
        self.y = (np.arange(NY) - NY // 2) * DY
        s, a, _ = self.r.spectrum()
        sy = self.y / F_FAC
        U = np.interp(sy, s, a.real, left=0, right=0) + 1j * np.interp(sy, s, a.imag, left=0, right=0)
        U[np.abs(sy) >= NA_FAC] = 0
        self.U_exit = U  # sine-condition lens: collimated, amplitude a(y / f), clipped at the lens exit
        self.clip_at_S = clip_at_S  # True: RelayY convention (clip in S itself)
        self.k = 2 * np.pi / LAM
        self.kz = np.sqrt((self.k**2 - (2 * np.pi * np.fft.fftfreq(NY, DY)) ** 2).astype(complex))
        self.eta = self.r.eta
        self.Py_ref = self.r.pupil()  # RelayY reference
        self._scale = None

    def prop(self, U, z):
        """Angular-spectrum propagation by z (um); evanescent components are dropped (they would grow when z < 0)."""
        H = np.where(self.kz.imag == 0, np.exp(1j * self.kz.real * z), 0)
        return np.fft.ifft(np.fft.fft(U) * H)

    def field_S(self, z=None, a=None, mode="slot"):
        """Field at S. z: comb distance from S (um); a: comb half-aperture (slot) or edge position (um)."""
        z0 = 0.0 if self.clip_at_S else -F_FAC  # plane of the aperture relative to S
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
        if self._scale is None:  # fix the constant against RelayY with the clip in S
            ref = self.prop(self.U_exit, 0.0)
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
    analytic = {f"z {z:g} mm, a {a:g} um": 2 * edge_rel * I_exit.max() * LAM * (z * 1e3 + F_FAC) /
                (2 * np.pi**2 * (a - F_FAC * NA_FAC)) / P_tot
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


Z_C = 7.0  # mm, nearest comb point from S: beyond the x slit plane (2 f1 = 6.24 mm from the facets), CHOICE


def z_max_rule(a):
    """Largest comb distance (mm) at which the in-band beam (half-angle NA/D at S) still fits the window a (um)."""
    return (a - F_FAC * NA_FAC) * D / NA * 1e-3


def layout():
    """First-order layouts (DERIVED) for 32 / 48 / 64 rows at 4 um, and the worst-row comb penalty (SIM)."""
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
        {"edge": edge, "layout": layout}[p]()
