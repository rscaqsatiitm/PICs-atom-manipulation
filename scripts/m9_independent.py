"""Module 9: independently steerable tweezers (AOD-equivalent) on the SiN platform.

Why an AOD gives independent tweezers: (i) each RF tone shifts its light by its own frequency, so tweezers add
incoherently on average; (ii) the acoustic transit time across the aperture converts each tone's frequency
into its own angle (resolvable spots = RF bandwidth x transit time). Integrated routes studied here:

  E  sliding-window focal-plane array: a dense input-port array (one port per resolvable spot) on a star
     coupler; each tweezer is a Gaussian-weighted group of ports, moved by sliding the window. In the port
     (position) basis the tweezers occupy separate ports, so their fields add linearly and independently.
  F  AOD emulation: complex (IQ) modulation of every arm with one RF tone per tweezer (frequency tags).
  A' coherent complex synthesis on the arms without tags (same linear superposition, static).
  D  wavelength-addressed dispersive OPA: one laser line per tweezer, passive delay-line array.
  G  integrated surface-acoustic-wave deflector; T time multiplexing (order-of-magnitude bounds only).

Outputs: results/m9_E_window.json, m9_E_independence.json, m9_papr.json, m9_routes.json (+ figures)
"""
from __future__ import annotations

import json
import sys

import numpy as np
from scipy.optimize import minimize

from picsim import params as prm
from picsim.opa.array import OPA, channel_response, spot_metrics
from picsim.opa.star import StarCoupler
from picsim.plotting import fig
from picsim.provenance import RESULTS as Path_R
from picsim.provenance import P, save_result

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from m8_twenty import emitter  # noqa: E402

SCRIPT = "scripts/m9_independent.py"
LAM = prm.LAM_TRANSPORT.value
NA = prm.NA_OBJ.value
N = 128
P_A, W_ARM = 2.0, 1.7  # arm pitch and width at the FPR edge (CHOICE, as m8)
P_PORT, W_PORT = 2.0, 1.7  # dense input ports (CHOICE): pitch and width at the FPR edge
P_OUT = 1.5
SITE = prm.SITE_PITCH.value


# ------------------------------------------------------------------ E: model
class FocalPlaneOPA:
    """Dense-port star coupler + 128-arm OPA. Port k maps to atom position x_k = k * dx_port."""

    def __init__(self, lam=LAM, n=N):
        self.opa = OPA(N=n, pitch=P_OUT, lam=lam, element=emitter(lam, 0.6), NA=NA)
        n_s = StarCoupler(lam, prm.T_SIN.value, 2, P_A, 100.0, W_PORT, W_ARM).n_s
        # one port per resolvable direction: p_port = lambda_s R / (N p_a)  ->  R
        self.R = n * P_A * P_PORT / (lam / n_s)
        self.sc = StarCoupler(lam, prm.T_SIN.value, n, P_A, self.R, W_PORT, W_ARM)
        self.k = np.arange(n) - (n - 1) / 2
        self.A = np.array([self.sc.arm_amplitudes(kk * P_PORT) for kk in self.k])[::-1]  # (ports, arms)
        self.dx_port = lam / (2 * NA)  # = 2 x_range / N
        self.x = np.arange(-1.05 * self.opa.x_range, 1.05 * self.opa.x_range, 0.02)
        self.G = channel_response(self.opa, self.x)
        self.e_pow = float(np.sum(self.opa.element(self.opa.x) ** 2) * self.opa.dx)

    def weights(self, u, sigma, amp=1.0, phase=0.0):
        """Gaussian port weights centred at fractional port index u (port k at x = k dx_port), unit power."""
        w = np.exp(-((self.k - u) ** 2) / (2 * sigma**2))
        return amp * np.exp(1j * phase) * w / np.sqrt(np.sum(w**2))

    def intensity(self, port_amps):
        c = port_amps @ self.A
        return np.abs(c @ self.G) ** 2 / (self.opa.lam * self.opa.x_scale * self.e_pow)

    def peak_near(self, I, x0, win=0.6):
        m = np.abs(self.x - x0) < win
        return spot_metrics(self.x[m], I[m])


def window_study(fp: FocalPlaneOPA):
    rows = []
    for sigma in (0.8, 0.85, 0.9, 1.0, 1.2):
        us = np.linspace(0, 1, 21)
        xs, pks, fw, nb, dl = [], [], [], [], []
        for u in us:
            I = fp.intensity(fp.weights(u, sigma))
            s = spot_metrics(fp.x, I)
            xs.append(s["x_peak"]), pks.append(s["peak"]), fw.append(s["fwhm"]), dl.append(np.trapezoid(I, fp.x))
            nb.append(max(np.interp(s["x_peak"] + SITE, fp.x, I), np.interp(s["x_peak"] - SITE, fp.x, I)) / s["peak"])
        xs, pks = np.array(xs), np.array(pks)
        lin = np.polyfit(us, xs, 1)
        rows.append({"sigma_ports": sigma, "slope_um_per_port": float(lin[0]),
                     "nonlinearity_nm": float(np.max(np.abs(xs - np.polyval(lin, us))) * 1e3),
                     "peak_ripple": float(np.ptp(pks) / pks.mean()), "fwhm_um": float(np.mean(fw)),
                     "neighbour_max": float(np.max(nb)), "delivered": float(np.mean(dl)),
                     "x_of_u": xs.tolist(), "peak_of_u": (pks / pks.mean()).tolist()})
    # across the field (sigma = 0.9): delivered fraction and peak vs position
    field = []
    for u in np.arange(-60, 61, 5.0):
        I = fp.intensity(fp.weights(u, 0.9))
        s = spot_metrics(fp.x, I)
        field.append({"x_um": float(s["x_peak"]), "delivered": float(np.trapezoid(I, fp.x)), "peak": float(s["peak"])})
    return rows, field


def calibrated(fp: FocalPlaneOPA, sigma, step=0.125):
    """Lookup tables over the whole field: window centre u -> spot position x(u) and peak(u).

    Both depend deterministically on u (sub-port ripple, element-factor and capture roll-off), so a tweezer is
    placed at x by inverting x(u) and its depth is equalised by scaling its weights by sqrt(peak_ref/peak(u)).
    """
    us = np.arange(-62.0, 62.0 + step / 2, step)
    xs, pk = [], []
    for u in us:
        s = spot_metrics(fp.x, fp.intensity(fp.weights(u, sigma)))
        xs.append(s["x_peak"]), pk.append(s["peak"])
    xs, pk = np.array(xs), np.array(pk)
    ref = pk.max()

    def u_for(x_target):
        return float(np.interp(x_target, xs, us))

    def amp_for(u):
        return float(np.sqrt(ref / np.interp(u, us, pk)))

    return u_for, amp_for


def power_factors(cal, x_half=60.0):
    """Laser power needed to keep the depth equal over |x| <= x_half, relative to a tweezer at the centre.

    cal = (u_for, amp_for): u_for maps a position x (um) to the window centre u (ports), and amp_for(u)^2 is the
    power factor relative to the brightest window position. Returns (field mean, worse edge), both relative to x = 0.
    """
    u_for, amp_for = cal
    c = amp_for(u_for(0.0)) ** 2
    xs = np.arange(-x_half, x_half + 0.5, 1.0)
    f = np.array([amp_for(u_for(x)) ** 2 / c for x in xs])
    return float(f.mean()), float(max(f[0], f[-1]))


def independence(fp: FocalPlaneOPA, sigma=0.9, seed=0, cal=None):
    rng = np.random.default_rng(seed)
    u_for, amp_for = cal if cal is not None else calibrated(fp, sigma)
    # 20 tweezers at random positions >= 4 um apart within +/-60 um, random relative phases
    while True:
        xt = np.sort(rng.uniform(-60, 60, 20))
        if np.min(np.diff(xt)) >= SITE:
            break
    th = rng.uniform(0, 2 * np.pi, 20)

    def field(xs):
        return sum(fp.weights(u_for(x), sigma, amp_for(u_for(x)), t) for x, t in zip(xs, th))

    I0 = fp.intensity(field(xt))
    base = [fp.peak_near(I0, x) for x in xt]
    pos_err0 = max(abs(b["x_peak"] - x) for b, x in zip(base, xt))
    spread0 = (max(b["peak"] for b in base) - min(b["peak"] for b in base)) / np.mean([b["peak"] for b in base])
    moves = []
    j = int(np.argmax(np.diff(xt)))  # tweezer with the largest free gap to its right
    gap = xt[j + 1] - xt[j]
    for d in (0.1, 0.5, 1.0, 2.0, min(4.0, gap - 2.0)):
        xs = xt.copy()
        xs[j] += d
        I = fp.intensity(field(xs))
        mv = fp.peak_near(I, xs[j])
        others = [fp.peak_near(I, x) for i, x in enumerate(xs) if i != j]
        ob = [b for i, b in enumerate(base) if i != j]
        moves.append({"move_um": float(d), "moved_pos_err_nm": float((mv["x_peak"] - xs[j]) * 1e3),
                      "moved_peak_change": float(mv["peak"] / base[j]["peak"] - 1),
                      "others_pos_change_max_nm": float(max(abs(o["x_peak"] - b["x_peak"]) for o, b in zip(others, ob)) * 1e3),
                      "others_peak_change_max": float(max(abs(o["peak"] / b["peak"] - 1) for o, b in zip(others, ob)))})
    # switch one tweezer off (remove its weights): effect on the others and residual at its site
    xs_off = np.delete(xt, j)
    th_off = np.delete(th, j)
    I_off = fp.intensity(sum(fp.weights(u_for(x), sigma, amp_for(u_for(x)), t) for x, t in zip(xs_off, th_off)))
    oth = [fp.peak_near(I_off, x) for x in xs_off]
    ob = [b for i, b in enumerate(base) if i != j]
    off = {"residual_at_site": float(np.interp(xt[j], fp.x, I_off) / base[j]["peak"]),
           "others_peak_change_max": float(max(abs(o["peak"] / b["peak"] - 1) for o, b in zip(oth, ob)))}
    # two tweezers approaching each other, relative phase 0, pi/2, pi
    appr = []
    for dphi in (0.0, np.pi / 2, np.pi):
        for sep in (8.0, 4.0, 3.0, 2.5, 2.0, 1.5):
            w = fp.weights(u_for(-sep / 2), sigma, amp_for(u_for(-sep / 2)), 0.0) + \
                fp.weights(u_for(sep / 2), sigma, amp_for(u_for(sep / 2)), dphi)
            I = fp.intensity(w)
            a, b = fp.peak_near(I, -sep / 2, 0.4), fp.peak_near(I, sep / 2, 0.4)
            Is = fp.intensity(fp.weights(u_for(sep / 2), sigma, amp_for(u_for(sep / 2))))
            ref = fp.peak_near(Is, sep / 2)["peak"]
            appr.append({"dphi": float(dphi), "sep_um": sep, "pull_nm": float((b["x_peak"] - sep / 2) * 1e3),
                         "peak_change": float(b["peak"] / ref - 1)})
    return {"positions_um": xt.tolist(), "sigma": sigma, "static_pos_err_max_nm": float(pos_err0 * 1e3),
            "static_peak_spread": float(spread0), "moved_index": j, "moves": moves, "switch_off": off,
            "approach": appr}


def amplitude_noise(fp: FocalPlaneOPA, sigma=0.9, trials=200, seed=2):
    rng = np.random.default_rng(seed)
    out = []
    for sa, sp in ((0.003, 0.003), (0.01, 0.01), (0.03, 0.03)):
        w0 = fp.weights(0.3, sigma)
        s0 = spot_metrics(fp.x, fp.intensity(w0))
        dx, dI = [], []
        for _ in range(trials):
            w = w0 * (1 + sa * rng.standard_normal(N)) * np.exp(1j * sp * rng.standard_normal(N))
            s = fp.peak_near(fp.intensity(w), s0["x_peak"])
            dx.append(s["x_peak"] - s0["x_peak"]), dI.append(s["peak"] / s0["peak"] - 1)
        out.append({"sigma_amp": sa, "sigma_phase_rad": sp, "pos_rms_nm": float(np.std(dx) * 1e3),
                    "peak_rms": float(np.std(dI))})
    return out


# ------------------------------------------------------------------ F and A': peak-to-average of arm drive
def papr(K=20, n=N, seed=0, tagged=False):
    """Efficiency penalty of complex arm synthesis: mean |c|^2 / max |c|^2, phases of the tweezers optimised.

    c_n(t) = sum_k exp(i (n dphi_k + theta_k + 2 pi k t))   (t in units of the tag period if tagged)
    """
    rng = np.random.default_rng(seed)
    x = np.sort(rng.choice(np.arange(-60, 61, 4.0), K, replace=False))
    dphi = 2 * np.pi * x / (n * (LAM / (2 * NA)))  # phase step per arm for position x
    nn = np.arange(n)[:, None]
    ts = np.linspace(0, 1, 64, endpoint=False)[:, None, None] if tagged else np.zeros((1, 1, 1))
    kk = np.arange(K)[None, None, :]

    def peak_over_mean(theta, p=16):
        c = np.exp(1j * (nn[None] * dphi[None, None, :] + theta[None, None, :] + 2 * np.pi * kk * ts)).sum(-1)
        a = np.abs(c) ** 2 / K
        return np.mean(a**p) ** (1 / p), a.max()

    best = None
    for _ in range(6):
        th0 = rng.uniform(0, 2 * np.pi, K)
        res = minimize(lambda th: peak_over_mean(th)[0], th0, method="Powell", options={"maxiter": 3000, "xtol": 1e-3})
        pk = peak_over_mean(res.x)[1]
        best = pk if best is None else min(best, pk)
    return {"K": K, "tagged": tagged, "peak_over_mean": float(best), "efficiency": float(1 / best)}


# ------------------------------------------------------------------ D, G, T: numbers
def routes_numbers():
    c = 299792458.0
    ng = 2.14  # SIM m1_design_modes, 1064_0.4
    dndT = 2.36e-5  # SIM m1_design_modes
    field = 2 * N * LAM / (4 * NA)  # full steering width at the atoms (um)
    D = []
    for per_site_GHz in (5, 20, 40, 100):
        fsr = per_site_GHz * 1e9 * field / SITE  # Hz, frequency span that maps onto the whole field
        tau0 = 1 / fsr
        dL = c * tau0 / ng  # m, extra length per arm
        L_max = (N - 1) * dL
        drift = (2 * np.pi / (LAM * 1e-6)) * dndT * L_max / (N - 1)  # rad/K per arm (uniform heating)
        dx_per_K = drift * LAM / (2 * np.pi * P_OUT) * (N * P_OUT / (2 * NA))  # um/K at the atoms
        D.append({"GHz_per_site": per_site_GHz, "FSR_THz": fsr / 1e12, "dL_per_arm_mm": dL * 1e3,
                  "longest_arm_extra_mm": L_max * 1e3, "total_extra_length_m": dL * N * (N - 1) / 2,
                  "shift_um_per_K": dx_per_K, "tuning_nm_full_field": fsr * (LAM * 1e-6) ** 2 / c * 1e9})
    f_obj = prm.WD.value * 1e3  # um, focal length ~ working distance (DERIVED order of magnitude)
    z_R = np.pi * prm.W0_TWEEZER.value**2 / LAM
    dlam_max = LAM * 1.0 / f_obj  # um, for |focal shift| <= 1 um of a purely diffractive lens (f ~ 1/lambda)
    chrom = {"diffractive_focal_shift_um_per_nm": f_obj / (LAM * 1e3), "Rayleigh_um": z_R,
             "max_detuning_pm_for_1um_focal_shift": dlam_max * 1e6,
             "max_detuning_GHz": c * dlam_max * 1e-6 / (LAM * 1e-6) ** 2 / 1e9}
    v_saw = 5.0e3  # m/s, RECALLED order of magnitude for AlN-on-oxide surface waves
    G = [{"aperture_mm": a, "RF_bandwidth_MHz_for_128_spots": N * v_saw / (a * 1e-3) / 1e6} for a in (1, 2, 5)]
    f_trap = 84e3
    T = {"cycle_rate_MHz_min": 10 * f_trap / 1e6, "pattern_rate_MHz_for_20": 20 * 10 * f_trap / 1e6,
         "slot_ns": 1e9 / (20 * 10 * f_trap), "peak_intensity_factor": 20}
    return {"D_dispersive": D, "D_chromatic_constraint": chrom, "G_saw": G, "T_time_multiplex": T}


if __name__ == "__main__":
    parts = sys.argv[1:] or ["E", "papr", "routes"]
    if "E" in parts:
        fp = FocalPlaneOPA()
        print(f"R = {fp.R:.0f} um, port pitch {P_PORT} um -> {fp.dx_port:.4f} um at the atoms", flush=True)
        rows, field = window_study(fp)
        for r in rows:
            print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items() if k not in ("x_of_u", "peak_of_u")}, flush=True)
        print("field", [(round(f["x_um"], 1), round(f["delivered"], 3)) for f in field], flush=True)
        f, ax = fig(1, 2, h=3.0)
        for r in rows:
            us = np.linspace(0, 1, 21)
            ax[0, 0].plot(us, (np.array(r["x_of_u"]) - np.polyval(np.polyfit(us, r["x_of_u"], 1), us)) * 1e3,
                          label=f"sigma {r['sigma_ports']}")
            ax[0, 1].plot(us, r["peak_of_u"], label=f"sigma {r['sigma_ports']}")
        ax[0, 0].set_xlabel("window centre (fraction of port pitch)")
        ax[0, 0].set_ylabel("position nonlinearity (nm)")
        ax[0, 1].set_xlabel("window centre (fraction of port pitch)")
        ax[0, 1].set_ylabel("peak / mean")
        ax[0, 0].legend(fontsize=6)
        save_result("m9_E_window", script=SCRIPT, fig=f,
                    inputs={"N": N, "port_pitch": "CHOICE 2.0 um", "port_width": "CHOICE 1.7 um",
                            "arm_pitch": "CHOICE 2.0 um", "emitter_pitch": "CHOICE 1.5 um", "NA": prm.NA_OBJ,
                            "R": P(fp.R, "um", "DERIVED", "R = N p_a p_port / lambda_s (one port per resolvable spot)"),
                            "model": "picsim.opa.star (2D paraxial FPR, isolated port and arm modes) + relay"},
                    outputs={"R_um": fp.R, "dx_port_um": fp.dx_port, "rows": rows, "field": field},
                    notes="Uncalibrated sliding Gaussian window; ripple and nonlinearity are deterministic in the "
                          "window position and are removed by a lookup table (see m9_E_independence).")
        cal = calibrated(fp, 0.9)
        ind = independence(fp, cal=cal)
        # power cost of equalising depth across the field: mean of 1/peak relative to the centre
        ind["power_factor_mean_over_field"], ind["power_factor_edge"] = power_factors(cal)
        print(json.dumps({k: v for k, v in ind.items() if k != "positions_um"}, indent=0, default=float), flush=True)
        noise = amplitude_noise(fp)
        print("noise", noise, flush=True)
        save_result("m9_E_independence", script=SCRIPT,
                    inputs={"sigma_ports": 0.9, "tweezers": 20, "placement": "random, >= 4 um apart, +/-60 um",
                            "relative_phases": "random", "calibration": "lookup of position and peak vs window position"},
                    outputs={"independence": ind, "amplitude_phase_noise": noise})
    if "power" in parts:
        from picsim import transport as tr
        fp = FocalPlaneOPA()
        m2 = json.loads((Path_R / "m2_mmi.json").read_text())["outputs"]
        m2c = json.loads((Path_R / "m2_couplers.json").read_text())["outputs"]
        ex22 = [r for r in m2["rows"] if r["lam"] == "1064" and r["W"] == 3.0][0]["excess_2x2_dB"]
        ind = json.loads((Path_R / "m9_E_independence.json").read_text())["outputs"]["independence"]
        chip_db = m2c["best"]["1064_0.4_lensed_2.5"]["dB"] - 2 * 7 * ex22 - 0.3 + 10 * np.log10(0.7)
        s0 = spot_metrics(fp.x, fp.intensity(fp.weights(0.0, 0.9)))
        w0y = prm.W0_TWEEZER.value
        peak2d = s0["peak"] * np.sqrt(2 / np.pi) / w0y * 10 ** (chip_db / 10)  # 1/um^2 per W at the fibre
        I_g = 2 / (np.pi * w0y**2)
        P_centre = prm.DEPTH_TRANSPORT.value / (tr.trap_depth_K(1.0, w0y, prm.ALPHA_CS_1064.value) * peak2d / I_g)
        out = {"chip_items_dB": {"facet": m2c["best"]["1064_0.4_lensed_2.5"]["dB"], "tree 1x128": -2 * 7 * ex22,
                                 "propagation": -0.3, "meta-objective 0.7": 10 * np.log10(0.7)},
               "centre_fwhm_um": s0["fwhm"], "P_fibre_per_tweezer_centre_mW": P_centre * 1e3,
               "P_fibre_per_tweezer_field_mean_mW": P_centre * ind["power_factor_mean_over_field"] * 1e3,
               "P_fibre_per_tweezer_edge_mW": P_centre * ind["power_factor_edge"] * 1e3,
               "P_fibre_20_tweezers_W_field_mean": 20 * P_centre * ind["power_factor_mean_over_field"]}
        print(json.dumps(out, indent=1), flush=True)
        save_result("m9_E_power", script=SCRIPT,
                    inputs={"depth": prm.DEPTH_TRANSPORT, "alpha": prm.ALPHA_CS_1064, "w0_y": prm.W0_TWEEZER,
                            "sigma_ports": 0.9, "meta_objective": "CHOICE 0.7", "tree": "7-stage tunable tree, 3 um 2x2 "
                            "MMIs (SIM m2_mmi)"},
                    outputs=out,
                    notes="Depth equalised over the field by the tunable tree; mean factor from m9_E_independence.")
    if "papr" in parts:
        pr = [papr(tagged=False), papr(tagged=True)]
        print("papr", pr, flush=True)
        save_result("m9_papr", script=SCRIPT,
                    inputs={"K": 20, "N": N, "positions": "20 random 4 um sites within +/-60 um",
                            "optimisation": "tweezer phases (Powell on a soft-max of |c|^2), best of 6 starts"},
                    outputs={"rows": pr},
                    notes="Complex synthesis on the arms must scale every arm to the largest |c_n|; the efficiency is "
                          "mean/max. With frequency tags (route F) the maximum is over arms and time.")
    if "routes" in parts:
        rn = routes_numbers()
        print(json.dumps(rn, indent=1, default=float), flush=True)
        save_result("m9_routes", script=SCRIPT,
                    inputs={"n_g": P(2.14, "-", "SIM", "m1_design_modes 1064_0.4"),
                            "dneff_dT": P(2.36e-5, "1/K", "SIM", "m1_design_modes 1064_0.4"),
                            "f_objective": P(prm.WD.value, "mm", "DERIVED", "focal length taken as the 15 mm working "
                                             "distance (order of magnitude)"),
                            "v_SAW": P(5.0e3, "m/s", "RECALLED", "AlN surface-acoustic-wave velocity, order of "
                                       "magnitude; verify for the AlN/SiN/SiO2 stack"),
                            "f_trap": P(84e3, "Hz", "DERIVED", "radial, 1.2 mK, w0 = 1.04 um (m3_update_rate)")},
                    outputs=rn)
