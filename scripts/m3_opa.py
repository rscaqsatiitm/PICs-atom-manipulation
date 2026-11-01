"""Module 3: OPA steering engine at 1064.4 nm, end-fire facet array vs surface grating antennas.

Parts (results/):
  m3_field_tradeoff      channels N and pupil fill vs steering range, spot size, sites addressable, power/tweezer
  m3_pitch_efficiency    slit (central-zone) efficiency and field roll-off vs emitter pitch; minimum pitch
                         set by inter-channel coupling for each emitter type
  m3_spot                atom-plane spot profile along the steering axis; Gaussian channel apodisation option
  m3_multibeam           K simultaneous tweezers: stray light at untargeted sites, uniformity, efficiency
  m3_noise               position jitter and intensity change vs channel phase error (MC vs analytic); bits;
                         PZT hysteresis budget
  m3_update_rate         motional excitation of a moved atom vs OPA update rate (staircase trajectory)
  m3_loss_budget         emitted-power chain from chip input to the atoms vs the proposal's ~10 dB
"""
from __future__ import annotations

import json
import sys

import numpy as np

from picsim import materials as mat
from picsim import modes as M
from picsim import params as prm
from picsim import transport as tr
from picsim.opa import multibeam as mb
from picsim.opa import noise
from picsim.opa.array import OPA, channel_response, spot_metrics
from picsim.parallel import pmap
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m3_opa.py"
LAM = prm.LAM_TRANSPORT.value
NA = prm.NA_OBJ.value
PITCH_SITE = prm.SITE_PITCH.value


# ------------------------------------------------------------------ emitter near field from the channel mode
def emitter_profile(width, t=prm.T_SIN.value):
    """Lateral near-field amplitude |Ex|(x) of the TE0 channel mode at mid-height, as an interpolating function."""
    ms = M.solve(M.XSection(cores=[(0, width)], t=t, x_left=-5, x_right=5), LAM, num_modes=2)
    te = M.fundamental(ms)
    x = np.linspace(-4.5, 4.5, 901)
    ex = np.abs(M.field_on_grid(te.mode, x, np.array([t / 2]), 0))[0]
    ex = ex / ex.max()
    f = lambda xx: np.interp(xx, x, ex, left=0.0, right=0.0)
    # 1/e^2 intensity half-width for reporting
    w_e = 0.5 * np.ptp(x[ex**2 > np.exp(-2)])
    return f, float(w_e), float(te.neff.real)


def design_width(key="1064_0.4"):
    d = json.loads((RESULTS / "m1_width_sweep.json").read_text())["outputs"]["design"]
    return d[key]["w_design"]


# ------------------------------------------------------------------ A. field / resolution trade-off
def field_tradeoff(elem):
    rows = []
    for N in (32, 64, 128, 256):
        for fill in (1.0, 0.75, 0.5):
            o = OPA(N=N, pitch=2.0, lam=LAM, element=elem, NA=NA, fill=fill)
            ph = o.steer_phases(0.0)
            xa = np.linspace(-3, 3, 1201)
            sm = spot_metrics(xa, o.atom_intensity(ph, xa))
            # usable range: main lobe keeps >= 50% of its on-axis peak (slit clipping + element factor roll-off)
            s_try = np.linspace(0, o.s_cut, 60)
            pk0 = sm["peak"]
            usable = 0.0
            for s0 in s_try:
                x0 = s0 * o.x_scale
                xx = np.linspace(x0 - 2, x0 + 2, 161)
                if o.atom_intensity(o.steer_phases(s0), xx).max() >= 0.5 * pk0:
                    usable = x0
            # power per tweezer for equal peak depth scales with 1/(peak intensity per W) along x
            rows.append({"N": N, "fill": fill, "aperture_um": o.A, "x_half_range_um": o.x_range,
                         "x_half_usable_50pct_um": usable, "fwhm_x_um": sm["fwhm"], "w_1e2_x_um": sm["w_1e2"],
                         "sites_in_usable_row": int(2 * np.floor(usable / PITCH_SITE) + 1),
                         "fraction_of_pm100um": usable / prm.FIELD_HALF.value,
                         "peak_per_W_rel": sm["peak"]})
    p_ref = [r for r in rows if r["N"] == 64 and r["fill"] == 1.0][0]["peak_per_W_rel"]
    for r in rows:
        r["power_for_equal_depth_rel"] = p_ref / r["peak_per_W_rel"]
    f, ax = fig(1, 2, h=3.0)
    for fill, c in zip((1.0, 0.75, 0.5), ("C0", "C1", "C2")):
        r = [x for x in rows if x["fill"] == fill]
        ax[0, 0].loglog([x["N"] for x in r], [x["x_half_usable_50pct_um"] for x in r], "-o", color=c, ms=3,
                        label=f"fill {fill}")
        ax[0, 1].semilogx([x["N"] for x in r], [x["sites_in_usable_row"] for x in r], "-o", color=c, ms=3)
    ax[0, 0].axhline(prm.FIELD_HALF.value, color="k", lw=0.8)
    ax[0, 0].set_xlabel("channels N")
    ax[0, 0].set_ylabel("usable half-range at atoms (um)")
    ax[0, 0].legend()
    ax[0, 1].axhline(prm.K_TWEEZERS.value, color="k", lw=0.8)
    ax[0, 1].set_xlabel("channels N")
    ax[0, 1].set_ylabel(f"sites on one row ({PITCH_SITE:.0f} um pitch)")
    save_result("m3_field_tradeoff", script=SCRIPT, fig=f,
                inputs={"lam": prm.LAM_TRANSPORT, "NA": prm.NA_OBJ, "site_pitch": prm.SITE_PITCH,
                        "field_target": prm.FIELD_HALF, "pitch_um": 2.0,
                        "usable_criterion": "steered main-lobe peak >= 50% of on-axis peak"},
                outputs={"rows": rows},
                notes="1D scalar Fourier optics; aplanatic objective (pupil radius F NA). The number of resolvable "
                      "positions equals N; the range at the atoms is N lambda/(4 NA fill).")
    return rows


# ------------------------------------------------------------------ B. pitch: efficiency vs crosstalk
def pitch_efficiency(elem, w_e):
    cpl = json.loads((RESULTS / "m1_coupling.json").read_text())["outputs"]["rows"]
    c = sorted([r for r in cpl if r["key"] == "1064_0.4"], key=lambda r: r["pitch"])
    pitch_c = np.array([r["pitch"] for r in c])
    # beyond ~1e5 um the supermode splitting (< 1e-6) is at the FEM floor and non-monotonic: treat the
    # largest value reached as a lower bound and keep L_c monotonic in pitch
    Lc = np.maximum.accumulate(np.minimum(np.array([r["Lc_um"] for r in c]), 1e5))

    def Lc_at(p):
        return float(np.exp(np.interp(p, pitch_c, np.log(Lc))))

    rows = []
    for p in (1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0):
        o = OPA(N=64, pitch=p, lam=LAM, element=elem, NA=NA)
        eff0 = o.slit_efficiency(o.steer_phases(0.0))
        s_edge = 0.8 * o.s_cut
        eff_e = o.slit_efficiency(o.steer_phases(s_edge))
        xx = np.linspace(-2, 2, 161)
        pk0 = o.atom_intensity(o.steer_phases(0.0), xx).max()
        pke = o.atom_intensity(o.steer_phases(s_edge), xx + s_edge * o.x_scale).max()
        rows.append({"pitch_um": p, "slit_eff_center": eff0, "slit_eff_at_0.8_range": eff_e,
                     "peak_rel_at_0.8_range": pke / pk0, "Lc_um": Lc_at(p) if p >= pitch_c[0] else None,
                     "L_parallel_max_for_-20dB_um": 0.064 * Lc_at(p) if p >= pitch_c[0] else None})
    # parallel length each emitter type needs
    L_endfire = 50.0  # um, straight section at the facet after the fan-in (CHOICE)
    lengths = {"endfire": L_endfire}
    raw = RESULTS / "m3_grating_fdtd_raw.json"
    if raw.exists():
        g = json.loads(raw.read_text())
        for e in sorted({r["etch_um"] for r in g}):
            best = max([r for r in g if r["etch_um"] == e], key=lambda r: r["directionality_c"])
            lengths[f"grating_etch{e*1e3:.0f}nm_box{best['t_box_um']:.0f}um"] = float(
                np.log(10) / (2 * best["alpha_per_um"]))  # length radiating 90% of the power
    out = {"rows": rows, "w_e_1e2_um": w_e, "parallel_lengths_um": lengths, "min_pitch_um": {}}
    p_fine = np.linspace(pitch_c[0], pitch_c[-1], 400)
    Lc_fine = np.exp(np.interp(p_fine, pitch_c, np.log(Lc)))
    for name, Lp in lengths.items():
        # power coupling over Lp below -20 dB: sin^2(pi Lp / (2 Lc)) < 0.01  ->  Lp <= 0.064 Lc
        ok = p_fine[0.064 * Lc_fine >= Lp]
        out["min_pitch_um"][name] = float(ok.min()) if len(ok) else f">{p_fine[-1]:.2f}"
    f, ax = fig(1, 2, h=3.0)
    P = [r["pitch_um"] for r in rows]
    ax[0, 0].plot(P, [r["slit_eff_center"] for r in rows], "-o", ms=3, label="on axis")
    ax[0, 0].plot(P, [r["slit_eff_at_0.8_range"] for r in rows], "-s", ms=3, label="at 0.8 of range")
    ax[0, 0].plot(P, [r["peak_rel_at_0.8_range"] for r in rows], "--", label="peak at 0.8 range / on axis")
    ax[0, 0].set_xlabel("emitter pitch (um)")
    ax[0, 0].set_ylabel("fraction")
    ax[0, 0].legend(fontsize=7)
    ax[0, 1].semilogy(pitch_c, 0.064 * Lc, "-o", ms=3, label="max parallel length for -20 dB")
    for (name, Lp), c in zip(lengths.items(), ("C1", "C2", "C3", "C4")):
        ax[0, 1].axhline(Lp, color=c, lw=1, label=name)
    ax[0, 1].set_xlabel("pitch (um)")
    ax[0, 1].set_ylabel("length (um)")
    ax[0, 1].legend(fontsize=7)
    save_result("m3_pitch_efficiency", script=SCRIPT, fig=f,
                inputs={"N": 64, "element": "TE0 channel mode lateral |Ex| at 1064.4 nm, t = 400 nm",
                        "coupling": "results/m1_coupling.json (1064_0.4)", "L_endfire_um": "CHOICE 50 um",
                        "grating": "results/m3_grating_fdtd_raw.json"},
                outputs=out)
    return out


# ------------------------------------------------------------------ C. spot shape and apodisation
def spot(elem):
    out = {}
    f, ax = fig(1, 1, w=4.5, h=3.0)
    xa = np.linspace(-6, 6, 2401)
    for name, apod in (("uniform", None), ("gauss_edge_1e2", "g")):
        o = OPA(N=64, pitch=2.0, lam=LAM, element=elem, NA=NA)
        if apod == "g":
            o.apod = np.exp(-(o.xn / (o.A / 2)) ** 2)  # 1/e^2 intensity at the aperture edge
        ph = o.steer_phases(0.0)
        I = o.atom_intensity(ph, xa)
        sm = spot_metrics(xa, I)
        # first sidelobe: largest local maximum outside the main lobe
        main = np.abs(xa - sm["x_peak"]) < 1.2 * sm["fwhm"]
        side = I[~main].max() / sm["peak"]
        out[name] = {**{k: float(v) for k, v in sm.items()}, "sidelobe_rel": float(side),
                     "transmission": o.transmission(ph)}
        ax[0, 0].semilogy(xa, I / sm["peak"], label=f"{name}: FWHM {sm['fwhm']:.2f} um, sidelobe {side:.3f}")
    # equivalent Gaussian of the proposal (w0 = 1.04 um): FWHM = w0 sqrt(2 ln 2)
    out["proposal_gaussian_fwhm_um"] = prm.W0_TWEEZER.value * np.sqrt(2 * np.log(2))
    ax[0, 0].set_ylim(1e-4, 1.5)
    ax[0, 0].set_xlabel("x at atoms (um)")
    ax[0, 0].set_ylabel("I / I_peak")
    ax[0, 0].legend(fontsize=6)
    save_result("m3_spot", script=SCRIPT, fig=f,
                inputs={"N": 64, "pitch_um": 2.0, "NA": prm.NA_OBJ, "w0_proposal": prm.W0_TWEEZER},
                outputs=out,
                notes="Uniform channel amplitudes give a sinc^2-like profile along the steering axis with ~4-5% "
                      "first sidelobes; edge apodisation trades spot width for sidelobes.")
    return out


# ------------------------------------------------------------------ D. multi-tweezer synthesis
def _mb_job(args):
    N, K, seed, w_dsgn = args
    elem, _, _ = EMITTER
    o = OPA(N=N, pitch=2.0, lam=LAM, element=elem, NA=NA)
    usable = 0.9 * o.x_range
    sites = np.arange(np.ceil(-usable / PITCH_SITE), np.floor(usable / PITCH_SITE) + 1) * PITCH_SITE
    if K > len(sites):
        return {"N": N, "K": K, "seed": seed, "feasible": False, "n_sites": len(sites)}
    rng = np.random.default_rng(1000 * N + 10 * K + seed)
    xt = np.sort(rng.choice(sites, K, replace=False))
    res = {"N": N, "K": K, "seed": seed, "feasible": True, "n_sites": len(sites)}
    for tag, kw in (("gsw", dict(beta=0, gamma=0, eta=0)), ("refined", dict(beta=10, gamma=20, eta=20))):
        ph = mb.synthesise(o, xt, PITCH_SITE, seed=seed, **kw)
        m = mb.evaluate(o, ph, xt, PITCH_SITE)
        res.update({f"{tag}_{k}": v for k, v in m.items()})
    return res


def multibeam():
    jobs = []
    for N, Ks in ((32, (1, 2, 4, 6, 8)), (64, (1, 2, 4, 8, 12, 16, 17)), (128, (1, 4, 8, 16, 21, 24, 32))):
        for K in Ks:
            for seed in range(4):
                jobs.append((N, K, seed, None))
    rows = pmap(_mb_job, jobs)
    summ = []
    f, ax = fig(1, 3, w=9, h=3.0)
    for N, c in zip((32, 64, 128), ("C0", "C1", "C2")):
        Ks = sorted({r["K"] for r in rows if r["N"] == N and r["feasible"]})
        for tag, ls in (("gsw", "--"), ("refined", "-")):
            st = [np.max([r[f"{tag}_stray_sites"] for r in rows if r["N"] == N and r["K"] == K and r["feasible"]]) for K in Ks]
            un = [np.max([r[f"{tag}_uniformity"] for r in rows if r["N"] == N and r["K"] == K and r["feasible"]]) for K in Ks]
            ef = [np.mean([r[f"{tag}_eff_spots"] for r in rows if r["N"] == N and r["K"] == K and r["feasible"]]) for K in Ks]
            ax[0, 0].plot(Ks, st, ls + "o", color=c, ms=3, label=f"N={N} {tag}")
            ax[0, 1].plot(Ks, un, ls + "o", color=c, ms=3)
            ax[0, 2].plot(Ks, ef, ls + "o", color=c, ms=3)
            for K, a, b, e in zip(Ks, st, un, ef):
                summ.append({"N": N, "K": K, "method": tag, "max_stray_sites": a, "max_uniformity": b, "mean_eff_spots": e})
    ax[0, 0].set_ylabel("max stray at untargeted sites / tweezer")
    ax[0, 1].set_ylabel("max (Imax-Imin)/(Imax+Imin)")
    ax[0, 2].set_ylabel("power in tweezers / emitted")
    for a in ax[0]:
        a.set_xlabel("simultaneous tweezers K")
    ax[0, 0].legend(fontsize=6)
    save_result("m3_multibeam", script=SCRIPT, fig=f,
                inputs={"pitch_um": 2.0, "NA": prm.NA_OBJ, "site_pitch": prm.SITE_PITCH,
                        "targets": "K distinct lattice sites drawn at random within 0.9 of the steering half-range; "
                                   "4 draws per (N, K)",
                        "methods": "gsw: weighted Gerchberg-Saxton (Di Leonardo et al., Opt. Express 15, 1913 (2007), "
                                   "RECALLED); refined: + L-BFGS on brightness, uniformity, centring and dark "
                                   "surroundings (beta=10, gamma=20, eta=20, q=2)"},
                outputs={"rows": rows, "summary": summ},
                notes="One row of the lattice along the steering axis. K is limited by the number of sites in range.")
    return summ


# ------------------------------------------------------------------ E. phase noise, quantisation, hysteresis
def phase_noise(elem):
    rows = []
    for N in (32, 64, 128):
        o = OPA(N=N, pitch=2.0, lam=LAM, element=elem, NA=NA)
        for sig in (0.01, 0.03, 0.1, 0.3):
            dx, dI = noise.monte_carlo(o, 0.3 * o.s_cut, sig, trials=150)
            rows.append({"N": N, "sigma_rad": sig, "dx_rms_um": float(dx.std()),
                         "dx_analytic_um": float(noise.analytic_position_rms(sig, N, LAM, NA)),
                         "dI_mean": float(dI.mean()), "dI_rms": float(dI.std())})
    # required phase accuracy for 10 nm rms position noise (analytic, verified by MC above)
    req = {N: float(prm.POS_NOISE.value / noise.analytic_position_rms(1.0, N, LAM, NA)) for N in (32, 64, 128)}
    bits = {b: float(noise.quantisation_sigma(b)) for b in (6, 8, 10, 12)}
    # PZT hysteresis: fractional error h of the applied phase; for phases spread over [0, 2 pi) the rms error
    # is h * 2 pi / sqrt(3) (uniform distribution)  -> h_max = sigma_req / (2 pi / sqrt(3))
    hyst = {N: req[N] / (2 * np.pi / np.sqrt(3)) for N in req}
    f, ax = fig(1, 1, w=4.5, h=3.0)
    for N in (32, 64, 128):
        r = [x for x in rows if x["N"] == N]
        ax[0, 0].loglog([x["sigma_rad"] for x in r], [x["dx_rms_um"] * 1e3 for x in r], "o", ms=4, label=f"MC N={N}")
        ax[0, 0].loglog([x["sigma_rad"] for x in r], [x["dx_analytic_um"] * 1e3 for x in r], "-", lw=0.8)
    ax[0, 0].axhline(prm.POS_NOISE.value * 1e3, color="k", lw=0.8)
    ax[0, 0].set_xlabel("channel phase error rms (rad)")
    ax[0, 0].set_ylabel("tweezer position rms (nm)")
    ax[0, 0].legend(fontsize=7)
    save_result("m3_noise", script=SCRIPT, fig=f,
                inputs={"pos_noise_target": prm.POS_NOISE, "NA": prm.NA_OBJ, "steer": "0.3 of half-range"},
                outputs={"rows": rows, "sigma_required_rad": req, "quantisation_sigma_rad": bits,
                         "pzt_hysteresis_fraction_max": hyst},
                notes="Static (frozen) channel errors. For a time-varying error the same mapping gives the "
                      "position noise spectrum, which sets heating (see m7_budget).")
    return req, hyst


# ------------------------------------------------------------------ F. update rate and transport heating
def update_rate():
    rows = []
    U_static = 0.3e-3
    for U in (U_static, prm.DEPTH_TRANSPORT.value):
        f_tr = tr.radial_freq_Hz(U, prm.W0_TWEEZER.value)
        for T_move in (50e-6, 100e-6, 200e-6, 400e-6):
            for d in (PITCH_SITE * 1e-6, 4 * PITCH_SITE * 1e-6):
                n_s = tr.n_bar_smooth(T_move, d, f_tr)
                for fu in (50e3, 100e3, 200e3, 300e3, 500e3, 1e6, 2e6):
                    for tau, tname in ((0.0, "ideal"), (1 / (2 * np.pi * prm.PZT_BW.value), "pzt_0.4MHz")):
                        rows.append({"U_mK": U * 1e3, "f_trap_kHz": f_tr / 1e3, "T_move_us": T_move * 1e6,
                                     "d_um": d * 1e6, "f_update_kHz": fu / 1e3, "actuator": tname,
                                     "n_bar": float(tr.n_bar_staircase(T_move, d, fu, f_tr, tau)),
                                     "n_bar_smooth": float(n_s)})
    f, ax = fig(1, 2, h=3.0)
    for j, U in enumerate((0.3, prm.DEPTH_TRANSPORT.value * 1e3)):
        for T_move, c in zip((50, 100, 200, 400), ("C0", "C1", "C2", "C3")):
            r = [x for x in rows if np.isclose(x["U_mK"], U) and x["T_move_us"] == T_move and x["d_um"] == PITCH_SITE
                 and x["actuator"] == "pzt_0.4MHz"]
            ax[0, j].loglog([x["f_update_kHz"] for x in r], np.maximum([x["n_bar"] for x in r], 1e-6), "-o", color=c,
                            ms=3, label=f"T = {T_move} us")
        ax[0, j].axvline(100, color="k", lw=0.8)
        ax[0, j].set_title(f"U = {U:.2g} mK, f_r = {r[0]['f_trap_kHz']:.0f} kHz; d = {PITCH_SITE:.0f} um")
        ax[0, j].set_xlabel("OPA update rate (kHz)")
    ax[0, 0].set_ylabel("added motional quanta n")
    ax[0, 0].legend(fontsize=7)
    save_result("m3_update_rate", script=SCRIPT, fig=f,
                inputs={"w0": prm.W0_TWEEZER, "alpha_Cs_1064": prm.ALPHA_CS_1064, "depth_transport": prm.DEPTH_TRANSPORT,
                        "U_static_mK": "CHOICE 0.3 (lower end of proposal range)", "pzt_bw": prm.PZT_BW,
                        "trajectory": "minimum-jerk, sample-and-hold at f_update, first-order actuator"},
                outputs={"rows": rows},
                notes="Harmonic approximation, radial axis, atom initially in the ground state; the staircase "
                      "aliases the move spectrum to |f_update - f_trap|.")
    return rows


# ------------------------------------------------------------------ G. loss budget to the atoms
def loss_budget(pe, mbs, grating):
    """Power chain from fibre to atoms (dB). Uses simulated values where available."""
    m2 = RESULTS / "m2_couplers.json"
    m2m = RESULTS / "m2_mmi.json"
    facet = json.loads(m2.read_text())["outputs"]["best_lensed_1064_dB"] if m2.exists() else -1.5
    mmi_ex = json.loads(m2m.read_text())["outputs"]["excess_loss_1064_dB"] if m2m.exists() else 0.15
    row = [r for r in pe["rows"] if r["pitch_um"] == 2.0][0]
    k21 = [s for s in mbs if s["N"] == 128 and s["K"] == 21 and s["method"] == "refined"]
    eff_mb = k21[0]["mean_eff_spots"] if k21 else np.nan
    eff_mb_slit = eff_mb  # eff_spots already includes slit and pupil
    items = {
        "fibre-to-chip edge coupler (SIM, m2)": facet,
        "6-stage MMI tree excess (SIM, m2)": 6 * (-abs(mmi_ex)),
        "propagation 1 cm at 0.3 dB/cm (CHOICE)": -0.3,
        "emitters + slit + pupil + K=21 synthesis, end-fire (SIM)": 10 * np.log10(eff_mb_slit),
        "meta-objective transmission 0.7 (CHOICE, Table 1 analogue)": 10 * np.log10(0.7),
    }
    if grating:
        items_g = dict(items)
        items_g["grating directionality (SIM, m3 FDTD)"] = 10 * np.log10(grating)
    total_e = sum(items.values())
    out = {"endfire_items_dB": items, "endfire_total_dB": total_e,
           "proposal_assumption_dB": -prm.LOSS_CHIP_OPA.value,
           "per_tweezer_on_chip_mW_needed_endfire": float(prm.P_TWEEZER_ATOMS.value * 1e3 / 10 ** (total_e / 10))}
    if grating:
        out["grating_items_dB"] = items_g
        out["grating_total_dB"] = sum(items_g.values())
    save_result("m3_loss_budget", script=SCRIPT, inputs={"P_tweezer_atoms": prm.P_TWEEZER_ATOMS,
                                                          "loss_assumed": prm.LOSS_CHIP_OPA}, outputs=out,
                notes="Per-tweezer share of the emitted power for K = 21 at N = 128, 2 um pitch. The 6.6 mW per "
                      "tweezer target uses the proposal's w0 = 1.04 um; the OPA spot differs (see m3_spot).")
    return out


# ------------------------------------------------------------------ H. grating-antenna summary (FDTD)
def grating_summary():
    """Tabulate the 2D FDTD grating results and derive emission-angle wavelength sensitivity.

    Grating equation (first order, into air): sin(theta) = n_avg - lambda/Lambda, hence
    d sin(theta)/d lambda = (n_avg - n_g)/lambda - 1/Lambda (n_g: group index of the grating section,
    approximated by the duty-weighted slab group index).
    """
    from picsim import slab

    raw = json.loads((RESULTS / "m3_grating_fdtd_raw.json").read_text())
    n_sin, n_ox = float(mat.n_sin(LAM)), float(mat.n_sio2(LAM))
    rows = []
    for r in sorted(raw, key=lambda r: (r["etch_um"], r["t_box_um"])):
        e = r["etch_um"]
        ng = []
        for tt in (prm.T_SIN.value, prm.T_SIN.value - e):
            npl = slab.slab_neff(tt, LAM + 0.001, float(mat.n_sin(LAM + 0.001)), float(mat.n_sio2(LAM + 0.001)))
            nmi = slab.slab_neff(tt, LAM - 0.001, float(mat.n_sin(LAM - 0.001)), float(mat.n_sio2(LAM - 0.001)))
            n0 = slab.slab_neff(tt, LAM, n_sin, n_ox)
            ng.append(n0 - LAM * (npl - nmi) / 0.002)
        n_g = 0.5 * (ng[0] + ng[1])
        ds = (r["n_avg_slab"] - n_g) / LAM - 1 / r["period_um"]  # per um
        th = np.radians(r["theta_air_deg"])
        dth_deg_per_nm = float(np.degrees(ds / np.cos(th)) * 1e-3)
        lam = np.array(r["lam_um"])
        d = np.array(r["up"]) / (np.array(r["up"]) + np.array(r["down"]))
        rows.append({"etch_um": e, "t_box_um": r["t_box_um"], "period_um": r["period_um"],
                     "theta_air_deg": r["theta_air_deg"], "alpha_per_um": r["alpha_per_um"],
                     "L_90pct_um": float(np.log(10) / (2 * r["alpha_per_um"])),
                     "directionality": r["directionality_c"], "reflection": r["refl_c"],
                     "directionality_min_pm15nm": float(d.min()), "directionality_max_pm15nm": float(d.max()),
                     "dtheta_dlambda_deg_per_nm": dth_deg_per_nm})
    f, ax = fig(1, 2, h=3.0)
    for e, c in zip(sorted({r["etch_um"] for r in rows}), ("C0", "C1", "C2")):
        rr = [r for r in rows if r["etch_um"] == e]
        ax[0, 0].plot([r["t_box_um"] for r in rr], [r["directionality"] for r in rr], "-o", color=c, ms=3,
                      label=f"etch {e*1e3:.0f} nm")
        ax[0, 1].semilogy([r["t_box_um"] for r in rr], [r["L_90pct_um"] for r in rr], "-o", color=c, ms=3)
    ax[0, 0].set_xlabel("buried oxide (um)")
    ax[0, 0].set_ylabel("up / (up + down)")
    ax[0, 0].legend(fontsize=7)
    ax[0, 1].set_xlabel("buried oxide (um)")
    ax[0, 1].set_ylabel("length for 90% emission (um)")
    save_result("m3_grating", script=SCRIPT, fig=f,
                inputs={"fdtd": "results/m3_grating_fdtd_raw.json (scripts/m3_grating_fdtd.py, MEEP 1.35, 2D, "
                                "40 px/um, Si substrate, 2 um top oxide, duty 0.5, nominal 10 deg)"},
                outputs={"rows": rows},
                notes="2D model of a uniform grating. Directionality oscillates with buried-oxide thickness through "
                      "the Si-substrate reflection; a bottom reflector or a two-layer grating is needed for > 0.8.")
    return rows


EMITTER = None

if __name__ == "__main__":
    parts = sys.argv[1:] or ["grating", "field", "pitch", "spot", "multibeam", "noise", "update", "budget"]
    if "grating" in parts:
        for r in grating_summary():
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
        if parts == ["grating"]:
            sys.exit(0)
    w = design_width()
    EMITTER = emitter_profile(w)
    elem, w_e, neff = EMITTER
    print(f"emitter width {w} um, 1/e^2 half-width {w_e:.3f} um, neff {neff:.4f}", flush=True)
    if "field" in parts:
        for r in field_tradeoff(elem):
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    pe = None
    if "pitch" in parts:
        pe = pitch_efficiency(elem, w_e)
        print(json.dumps({k: v for k, v in pe.items() if k != "rows"}, default=float), flush=True)
        for r in pe["rows"]:
            print(r, flush=True)
    if "spot" in parts:
        print(json.dumps(spot(elem), indent=0), flush=True)
    mbs = None
    if "multibeam" in parts:
        mbs = multibeam()
        for s in mbs:
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()}, flush=True)
    if "noise" in parts:
        print(phase_noise(elem), flush=True)
    if "update" in parts:
        update_rate()
        print("update-rate done", flush=True)
    if "budget" in parts:
        if pe is None:
            pe = json.loads((RESULTS / "m3_pitch_efficiency.json").read_text())["outputs"]
        if mbs is None:
            mbs = json.loads((RESULTS / "m3_multibeam.json").read_text())["outputs"]["summary"]
        raw = RESULTS / "m3_grating_fdtd_raw.json"
        gdir = None
        if raw.exists():
            g = json.loads(raw.read_text())
            gdir = max(r["directionality_c"] for r in g)
        print(json.dumps(loss_budget(pe, mbs, gdir), indent=1), flush=True)
