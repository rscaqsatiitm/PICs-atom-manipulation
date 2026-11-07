"""Module 24: light leakage between engines tiled along x with a shared slit plane (m22 tiling).

Geometry (1D, x across the row): engine i's emitter array (D1: 192 x 1.5 um, 1.2 um facets) sits in the front focal
plane of its own lenslet (focal length f1, clear aperture a = fill x Pe, centred on the engine); the slit plane
Sigma is the lenslets' common back focal plane. Engine i uses |X - X_i| < f1 s_cut of Sigma; its neighbour j = i + 1
uses |X - X_i - Pe| < f1 s_cut. Sigma is conjugate to the atoms, so the intensity that engine i puts into its
neighbour's slit zone appears as a ghost in the neighbour's segment.

Two paths:
  (a) through engine i's own lenslet: wave optics. The emitter near field is propagated over f1 (angular spectrum,
      exact), clipped by the lenslet aperture, given an ideal lens phase, and propagated over f1 to Sigma. Reported:
      peak intensity in the neighbour's slit zone relative to the tweezer spot peak, and the power fraction;
  (b) through neighbouring lenslets: exact ray trace (thin lenses, exact angles) of rays at all angles that miss
      lenslet i; reported: the closest landing point to the neighbour's slit zone.

Output: results/m24_tiling_wave.json (+ .png)
"""
from __future__ import annotations

import json

import numpy as np

from picsim import params as prm
from picsim.opa.focal_plane import FocalPlane, LinearSynth
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m24_tiling_wave.py"
LAM = prm.LAM_TRANSPORT.value  # um
F1 = 3120.0  # um (m13 x relay)
D1 = {"N_ports": 128, "N_arms": 192, "w_emit": 1.2, "p_out": 1.5, "w_port": 1.7, "w_arm": 1.7}


def emitter_field(fp, c, x):
    """Near field of the engine on grid x (um): sum of arm amplitudes c times the facet profile."""
    opa = fp.opa
    u = np.zeros(len(x), dtype=complex)
    e = opa.element
    for cn, xn in zip(c, opa.xn):
        m = np.abs(x - xn) < 5.0
        u[m] += cn * e(x[m] - xn)
    return u


def propagate(u, dx, z):
    k = 2 * np.pi / LAM
    fx = np.fft.fftfreq(len(u), dx)
    kz = np.sqrt((k**2 - (2 * np.pi * fx) ** 2).astype(complex))
    return np.fft.ifft(np.fft.fft(u) * np.exp(1j * kz * z))


def wave_case(fp, c, Pe, fill, dx=0.25, half=None):
    half = half or max(3 * Pe, 12000.0)
    n = int(2 ** np.ceil(np.log2(2 * half / dx)))
    x = (np.arange(n) - n // 2) * dx
    u0 = emitter_field(fp, c, x)
    k = 2 * np.pi / LAM
    u1 = propagate(u0, dx, F1)
    a = fill * Pe / 2
    u1 = np.where(np.abs(x) <= a, u1, 0) * np.exp(-1j * k * x**2 / (2 * F1))
    u2 = propagate(u1, dx, F1)
    I = np.abs(u2) ** 2
    s_cut = fp.opa.s_cut
    acc = fp.opa.A / (2 * F1)  # angular acceptance at Sigma that reaches the objective pupil (telecentric relay)

    def accepted(zone):
        """Field in a slit zone, keeping only directions inside the pupil acceptance."""
        E = np.where(zone, u2, 0)
        fx = np.fft.fftfreq(len(E), dx) * LAM  # direction sine
        return np.abs(np.fft.ifft(np.where(np.abs(fx) <= acc, np.fft.fft(E), 0))) ** 2

    own = np.abs(x) < F1 * s_cut
    nb = np.abs(x - Pe) < F1 * s_cut
    nb2 = np.abs(x + Pe) < F1 * s_cut
    Io, In, In2 = accepted(own), accepted(nb), accepted(nb2)
    p_own = I[own].sum()
    return {"Pe_um": Pe, "fill": fill,
            "peak_ratio_neighbour_raw": float(max(I[nb].max(), I[nb2].max()) / I[own].max()),
            "peak_ratio_neighbour_accepted": float(max(In.max(), In2.max()) / Io.max()),
            "power_ratio_neighbour_accepted": float(max(In.sum(), In2.sum()) / Io.sum()),
            "own_slit_power_fraction_of_aperture": float(p_own / max(I.sum(), 1e-300))}, (x, I)


def ray_case(Pe, fill, s_cut, acc, n_lens=4):
    """Rays from the engine centre at all propagating angles that miss lenslet i: landing points on Sigma through
    lenslets j = +-1 ... +-n_lens; distance to the nearest neighbour slit zone (negative = inside)."""
    th = np.linspace(-np.pi / 2 + 1e-4, np.pi / 2 - 1e-4, 200001)
    y = F1 * np.tan(th)  # height at the lenslet plane
    a = fill * Pe / 2
    worst = np.inf
    for j in range(-n_lens, n_lens + 1):
        if j == 0:
            continue
        Xj = j * Pe
        m = np.abs(y - Xj) <= a
        if not m.any():
            continue
        # thin lens j: direction changes by -atan((y - Xj)/F1) (exact for a ray through the lens's front focal point
        # structure is not assumed; this is the ideal-lens angle mapping tan(th') = tan(th) - (y - Xj)/F1)
        tp = np.tan(th[m]) - (y[m] - Xj) / F1
        X = y[m] + F1 * tp
        ok = np.abs(tp) <= acc  # only rays inside the pupil acceptance reach the atoms
        if not ok.any():
            continue
        for k in range(-n_lens - 2, n_lens + 3):
            if k == 0:
                continue
            d = np.abs(X[ok] - k * Pe) - F1 * s_cut
            worst = min(worst, float(d.min()))
    return worst


def main():
    T = np.load(RESULTS / "m10_T_opt_A.npz")["T"]
    fp = FocalPlane(A=T, **D1)
    ls = LinearSynth(fp)
    s_cut = fp.opa.s_cut
    rows, curves = [], {}
    for xa in (0.0, 60.0, 90.0):
        w = ls.weights([xa])
        c = w @ fp.A
        for Pe, fill in ((4430.0, 1.0), (4430.0, 0.8), (5000.0, 1.0), (6000.0, 1.0), (8000.0, 1.0), (3500.0, 1.0),
                         (3000.0, 1.0), (2600.0, 1.0), (2400.0, 1.0)):
            r, (x, I) = wave_case(fp, c, Pe, fill)
            r["tweezer_x_um"] = xa
            r["ray_min_distance_to_other_slit_um"] = ray_case(Pe, fill, s_cut, fp.opa.A / (2 * F1))
            rows.append(r)
            print({k: (f"{v:.3e}" if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
            if xa == 90.0 and fill == 1.0 and Pe in (4430.0, 3500.0):
                sel = slice(None, None, 64)
                curves[f"Pe {Pe}"] = [x[sel].tolist(), (I[sel] / I.max()).tolist()]
    f, ax = fig(1, 1, h=3.0)
    for k, (xx, yy) in curves.items():
        ax[0, 0].semilogy(np.array(xx) / 1e3, np.maximum(yy, 1e-14), lw=0.7, label=k)
    ax[0, 0].axvspan((4430 - F1 * s_cut) / 1e3, (4430 + F1 * s_cut) / 1e3, color="0.85")
    ax[0, 0].set_xlabel("slit plane X (mm)")
    ax[0, 0].set_ylabel("intensity / tweezer peak")
    ax[0, 0].legend(fontsize=7)
    save_result("m24_tiling_wave", script=SCRIPT, fig=f,
                inputs={"engine": "D1 (BPM matrix), LinearSynth tweezer at x_a", "f1_um": F1, "s_cut": s_cut,
                        "grid_um": 0.25, "lens": "ideal paraxial phase on the lenslet aperture",
                        "fill": "lenslet clear aperture / engine pitch"},
                outputs={"rows": rows, "limits": {"lenslet_half_aperture_for_slit_um": float(F1 * np.tan(np.arcsin(s_cut))),
                                                  "min_engine_pitch_um": float(2 * F1 * np.tan(np.arcsin(s_cut))),
                                                  "slit_zone_width_um": float(2 * F1 * s_cut)}},
                notes="peak_ratio_*: largest intensity in a neighbour's slit zone relative to the tweezer spot peak in Sigma "
                      "(Sigma is conjugate to the atoms); 'accepted' keeps only directions within the pupil acceptance "
                      "A/(2 f1). ray_min_distance: closest approach of an accepted ray to another engine's slit zone "
                      "(inf: no accepted ray through any neighbouring lenslet).")
    # own-slit throughput: the lenslet must pass |s| < s_cut, i.e. a half-aperture of f1 tan(asin(s_cut))
    out_min = {"lenslet_half_aperture_for_slit_um": float(F1 * np.tan(np.arcsin(s_cut))),
               "min_engine_pitch_um": float(2 * F1 * np.tan(np.arcsin(s_cut))),
               "slit_zone_width_um": float(2 * F1 * s_cut)}
    print(out_min, flush=True)
    return rows, out_min


if __name__ == "__main__":
    main()
