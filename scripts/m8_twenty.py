"""Module 8: twenty steerable, individually switchable tweezers at 1064.4 nm and 894.6 nm.

Two architectures through one full-NA aperture from one laser:
  (A) shared OPA with phase-only multi-spot synthesis; a tweezer is switched off by re-synthesis;
  (C) multi-port star-coupler OPA: one input port per tweezer, each with its own MZI switch, a free-propagation
      region (FPR) feeding N arms and emitters; a common phase ramp on the arms moves all lit tweezers
      together (parallel move), the lit ports select which sites carry tweezers.

Outputs (results/):
  m8_star_design      FPR radius, port pitch, capture vs port, per-port spot, neighbour crosstalk
  m8_star_multi       all/partial ports lit with random relative phases: peak spread and stray at dark sites
  m8_star_steer       collective shift: position tracking and peak change of all lit tweezers
  m8_resynth          architecture A: residual at a switched-off site, disturbance of the others, phase path
  m8_power            per-tweezer and total on-chip power at 1064.4 nm; addressing power at 894.6 nm
  m8_summary          comparison table and channel/driver counts
"""
from __future__ import annotations

import json
import sys

import numpy as np

from picsim import modes as M
from picsim import params as prm
from picsim import transport as tr
from picsim.opa import multibeam as mb
from picsim.opa.array import OPA, channel_response, spot_metrics
from picsim.opa.star import StarCoupler, design_radius
from picsim.plotting import fig
from picsim.provenance import RESULTS, P, save_result

SCRIPT = "scripts/m8_twenty.py"
K = 20
SITE = prm.SITE_PITCH.value
NA = prm.NA_OBJ.value
T_FILM = prm.T_SIN.value
P_A = 2.0  # arm pitch at the FPR edge (CHOICE)
W_ARM = 1.7  # arm width at the FPR edge, 0.3 um gap (CHOICE)
P_OUT = 1.5  # emitter pitch (CHOICE; >= 1.46 um end-fire minimum from m3_pitch_efficiency)
CASES = {"1064": {"lam": prm.LAM_TRANSPORT.value, "w": 0.6}, "895": {"lam": prm.LAM_ADDRESS.value, "w": 0.5}}
TARGETS = (np.arange(K) - (K - 1) / 2) * SITE  # 20 sites, -38 ... +38 um


def emitter(lam, width):
    ms = M.solve(M.XSection(cores=[(0, width)], t=T_FILM, x_left=-5, x_right=5), lam, num_modes=2)
    te = M.fundamental(ms)
    x = np.linspace(-4.5, 4.5, 901)
    ex = np.abs(M.field_on_grid(te.mode, x, np.array([T_FILM / 2]), 0))[0]
    ex = ex / ex.max()
    return lambda xx: np.interp(xx, x, ex, left=0.0, right=0.0)


# Final configuration per wavelength (chosen from m8_apod and the port-gap check):
#   1064.4 nm: envelope exp(-6) at the outer arms, 1.4 um input tapers -> port pitch 1.89 um, gap 0.49 um
#   894.6 nm:  envelope exp(-6), 1.0 um input tapers -> port pitch 1.64 um, gap 0.64 um (crosstalk ~3e-6)
# The edge level is set from the 1/e^2 width of the port far field; for wide (non-Gaussian) taper modes the
# tails differ, so both were chosen from explicit scans (session log, documented in docs/twenty_tweezers.md).
FINAL = {"1064": {"edge": np.exp(-6), "w_in": 1.4}, "895": {"edge": np.exp(-6), "w_in": 1.0}}


def build(case, N, edge=np.exp(-2), w_in=None):
    lam, w = CASES[case]["lam"], CASES[case]["w"]
    w_in = w if w_in is None else w_in
    R = design_radius(lam, T_FILM, N, P_A, w_in, W_ARM, edge=edge)
    sc = StarCoupler(lam, T_FILM, N, P_A, R, w_in, W_ARM)
    opa = OPA(N=N, pitch=P_OUT, lam=lam, element=emitter(lam, w), NA=NA)
    # port offsets that place beams on the 4 um sites: s = n_s d p_a / (R p_out)
    d = TARGETS / opa.x_scale * R * P_OUT / (sc.n_s * P_A)
    A = np.array([sc.arm_amplitudes(dj) for dj in d])  # (K, N) complex arm amplitudes per port
    # sign convention: the FPR phase -k_s d y/R steers to -s; flip the port order so that port j -> TARGETS[j]
    A = A[::-1]
    return lam, sc, opa, d, A


def port_field(opa, amps, x):
    """Atom-plane intensity (1/um, normalised to emitted power) for complex channel amplitudes."""
    a = np.abs(amps)
    keep = opa.apod
    opa.apod = a / a.max()
    I = opa.atom_intensity(np.angle(amps), x) * 1.0
    p_rel = np.sum(a**2) / np.sum((a / a.max()) ** 2)  # undo the apod normalisation (power scales with |a|^2)
    opa.apod = keep
    return I, p_rel


def design(case, N, edge=np.exp(-2), w_in=None):
    lam, sc, opa, d, A = build(case, N, edge, w_in)
    x = np.arange(-opa.x_range * 1.05, opa.x_range * 1.05, 0.02)
    G = channel_response(opa, x)  # (N, len(x)), unit channel amplitudes
    # emitted power for given channel amplitudes: element powers add (non-overlapping emitters)
    e_pow = float(np.sum(np.abs(opa.element(opa.x)) ** 2) * opa.dx)
    rows = []
    for j in range(K):
        c = A[j]
        I = np.abs(c @ G) ** 2 / (opa.lam * opa.x_scale * e_pow)  # per unit port input power
        sm = spot_metrics(x, I)
        nb = [np.interp(TARGETS[j] + k * SITE, x, I) / sm["peak"] for k in (-2, -1, 1, 2)]
        rows.append({"port": j, "d_um": float(d[::-1][j]), "x_target": float(TARGETS[j]),
                     "x_peak": float(sm["x_peak"]), "fwhm": float(sm["fwhm"]), "peak_per_W_in": float(sm["peak"]),
                     "capture": float(np.sum(np.abs(c) ** 2)), "delivered_frac": float(np.trapezoid(I, x)),
                     "nb_max": float(max(nb)), "nnb_max": float(max(nb[0], nb[3]))})
    out = {"lam": lam, "N": N, "R_um": sc.R, "n_s": sc.n_s, "port_pitch_um": float(abs(d[1] - d[0])),
           "port_span_um": float(np.ptp(d)), "x_range_um": opa.x_range, "envelope_edge": sc.envelope_edge(),
           "rows": rows}
    return out, (lam, sc, opa, d, A, x, G, e_pow)


def multi(ctx, n_draw=40, seed=0):
    lam, sc, opa, d, A, x, G, e_pow = ctx
    rng = np.random.default_rng(seed)
    # equalise single-port peaks with port powers (router setting)
    pk = np.array([spot_metrics(x, np.abs(A[j] @ G) ** 2)["peak"] for j in range(K)])
    amp = np.sqrt(pk.max() / pk)
    res = []
    for k_on in (K, 10, 5):
        spread, stray = [], []
        for _ in range(n_draw):
            on = np.sort(rng.choice(K, k_on, replace=False))
            ph = rng.uniform(0, 2 * np.pi, K)
            c = (amp[on] * np.exp(1j * ph[on])) @ A[on]
            I = np.abs(c @ G) ** 2
            peaks = np.array([I[np.abs(x - TARGETS[j]) < 0.5].max() for j in on])
            ref = pk.max()  # equalised single-port peak
            spread.append((peaks.max() - peaks.min()) / (peaks.max() + peaks.min()))
            off = np.setdiff1d(np.arange(K), on)
            if len(off):
                stray.append(max(np.interp(TARGETS[j], x, I) for j in off) / peaks.mean())
            else:
                stray.append(0.0)
            # peak change relative to single port (interference with other ports' tails)
        res.append({"k_on": k_on, "peak_spread_max": float(np.max(spread)), "peak_spread_median": float(np.median(spread)),
                    "stray_dark_sites_max": float(np.max(stray)), "stray_dark_sites_median": float(np.median(stray))})
    return res, amp


def steer(ctx, amp):
    lam, sc, opa, d, A, x, G, e_pow = ctx
    rows = []
    c0 = amp @ A  # all ports lit, equal relative phase
    I0 = np.abs(c0 @ G) ** 2
    pk0 = np.array([I0[np.abs(x - t) < 0.5].max() for t in TARGETS])
    for shift in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 24.0):
        ramp = np.exp(1j * opa.steer_phases(shift / opa.x_scale))
        I = np.abs((c0 * ramp) @ G) ** 2
        xs, pks = [], []
        for t in TARGETS:
            win = np.abs(x - (t + shift)) < 0.6
            if not np.any(win):
                xs.append(np.nan)
                pks.append(0.0)
                continue
            sm = spot_metrics(x[win], I[win])
            xs.append(sm["x_peak"] - t)
            pks.append(sm["peak"])
        xs, pks = np.array(xs), np.array(pks)
        rows.append({"shift_um": shift, "shift_error_max_nm": float(np.nanmax(np.abs(xs - shift)) * 1e3),
                     "peak_rel_min": float(np.min(pks / pk0)), "peak_rel_max": float(np.max(pks / pk0))})
    return rows


def resynth(N=128, seed=0):
    """Architecture A: K = 20 by synthesis; switch one tweezer off by re-synthesis (warm start)."""
    lam = prm.LAM_TRANSPORT.value
    opa = OPA(N=N, pitch=P_OUT, lam=lam, element=emitter(lam, CASES["1064"]["w"]), NA=NA)
    xt = TARGETS
    ph20 = mb.synthesise(opa, xt, SITE, seed=seed)
    x = np.arange(-opa.x_range, opa.x_range, opa.spot_pitch / 20)
    out = []
    for i_off in (0, 9):
        keep = np.delete(xt, i_off)
        pts = np.arange(-opa.x_range, opa.x_range, opa.spot_pitch / 2)
        dark = pts[np.min(np.abs(pts[:, None] - keep[None, :]), axis=1) > 1.5 * opa.spot_pitch]
        ph19 = mb.refine(opa, ph20, keep, dark, beta=10.0, gamma=20.0, eta=20.0, q=2)
        I20 = opa.atom_intensity(ph20, x)
        I19 = opa.atom_intensity(ph19, x)
        pk = lambda I, ts: np.array([I[np.abs(x - t) < 0.5].max() for t in ts])
        p20, p19 = pk(I20, keep), pk(I19, keep)
        resid = float(np.interp(xt[i_off], x, I19) / p19.mean())
        # linear phase path from ph20 to ph19 (wrapped differences): disturbance of kept tweezers
        dph = np.angle(np.exp(1j * (ph19 - ph20)))
        dev = []
        for f in np.linspace(0, 1, 11):
            If = opa.atom_intensity(ph20 + f * dph, x)
            pf = pk(If, keep)
            dev.append(float(np.max(np.abs(pf / p20 - 1))))
        out.append({"switched_off_site_um": float(xt[i_off]), "residual_rel": resid,
                    "kept_peak_change_max": float(np.max(np.abs(p19 / p20 - 1))),
                    "path_kept_deviation_max": float(np.max(dev)), "max_phase_step_rad": float(np.max(np.abs(dph))),
                    "stray_sites_after": mb.evaluate(opa, ph19, keep, SITE)["stray_sites"]})
    return out


def power(designs):
    """On-chip power for 20 tweezers at 1064.4 nm (1.2 mK) and addressing power at 894.6 nm."""
    m2 = json.loads((RESULTS / "m2_mmi.json").read_text())["outputs"]
    m2c = json.loads((RESULTS / "m2_couplers.json").read_text())["outputs"]
    ex22 = {lk: [r for r in m2["rows"] if r["lam"] == lk and r["W"] == 3.0][0]["excess_2x2_dB"] for lk in ("1064", "895")}
    facet = {"1064": m2c["best"]["1064_0.4_lensed_2.5"]["dB"], "895": m2c["best"]["895_0.4_lensed_2.5"]["dB"]}
    out = {}
    w0y = prm.W0_TWEEZER.value
    for lk, dsg in designs.items():
        rows = dsg["rows"]
        # 1D peak per unit port power x Gaussian y-profile (cylindrical relay fills NA in y) -> per W at the port
        peak2d = np.array([r["peak_per_W_in"] for r in rows]) * np.sqrt(2 / np.pi) / w0y  # 1/um^2 per W
        router = -2 * 5 * ex22[lk]  # 5-stage tunable-splitter tree (1x32), two 3 um 2x2 MMIs per stage
        switch = -2 * ex22[lk] * (2 if lk == "895" else 1)  # 895: two cascaded MZIs; 1064: one trimmed MZI
        chip_db = facet[lk] + router + switch - 0.3  # + 1 cm at 0.3 dB/cm (CHOICE)
        obj = 10 * np.log10(0.7)  # meta-objective transmission (CHOICE)
        eff = 10 ** ((chip_db + obj) / 10) * peak2d  # peak intensity at atoms per W at the fibre (1/um^2)
        res = {"chip_items_dB": {"facet": facet[lk], "router 1x32": router, "port switch": switch, "propagation": -0.3,
                                 "meta-objective 0.7": obj},
               "port_capture_min": float(min(r["capture"] for r in rows)),
               "port_capture_max": float(max(r["capture"] for r in rows)),
               "delivered_frac_min": float(min(r["delivered_frac"] for r in rows)),
               "delivered_frac_max": float(max(r["delivered_frac"] for r in rows))}
        if lk == "1064":
            I_g = 2 / (np.pi * w0y**2)
            U_per_W_gauss = tr.trap_depth_K(1.0, w0y, prm.ALPHA_CS_1064.value)  # K per W for the Gaussian
            P_fibre = prm.DEPTH_TRANSPORT.value / (U_per_W_gauss * eff / I_g)  # W at fibre per tweezer
            res.update({"P_fibre_per_tweezer_mW_min": float(P_fibre.min() * 1e3),
                        "P_fibre_per_tweezer_mW_max": float(P_fibre.max() * 1e3),
                        "P_fibre_total_20_W": float(P_fibre.sum()),
                        "end_to_end_dB_best_port": float(10 * np.log10(eff.max() / I_g)),
                        "end_to_end_dB_worst_port": float(10 * np.log10(eff.min() / I_g))})
        else:
            p_site = max(prm.P_SITE_ADDR.value)
            I_g = 2 / (np.pi * w0y**2)
            P_fibre = p_site * I_g / eff
            res.update({"P_fibre_per_site_uW_for_20uW_equiv_max": float(P_fibre.max() * 1e6),
                        "P_fibre_total_20_mW": float(P_fibre.sum() * 1e3)})
        out[lk] = res
    return out


def apod_sweep():
    """Envelope edge level vs neighbour crosstalk, spot width, delivery and multi-port peak spread (N = 128)."""
    rows = []
    for case in ("1064", "895"):
        for e_pow in (2, 3, 4, 5, 6, 8):
            dsg, ctx = design(case, 128, edge=np.exp(-e_pow))
            ml, _ = multi(ctx, n_draw=30)
            r = dsg["rows"]
            rows.append({"case": case, "edge": f"exp(-{e_pow})", "R_um": dsg["R_um"], "port_pitch_um": dsg["port_pitch_um"],
                         "fwhm_um": float(np.mean([x["fwhm"] for x in r])),
                         "nb_max": float(max(x["nb_max"] for x in r)), "nnb_max": float(max(x["nnb_max"] for x in r)),
                         "delivered_min": float(min(x["delivered_frac"] for x in r)),
                         "delivered_max": float(max(x["delivered_frac"] for x in r)),
                         "spread20_max": ml[0]["peak_spread_max"], "stray10_max": ml[1]["stray_dark_sites_max"]})
            print("apod", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    return rows


def all_lit_phases(n_trials=50, seed=0):
    """Static position error of the 20 spots with all ports lit (depth-equalised), for equal and for random
    relative port phases. A spot is displaced by the tails of the other ports' fields, so the error depends on the
    relative port phases, which are not stabilised in the star-coupler design."""
    from picsim.opa.array import spot_metrics
    out = {}
    for case in ("1064", "895"):
        cfg = FINAL[case]
        _, ctx = design(case, 128, edge=cfg["edge"], w_in=cfg["w_in"])
        lam, sc, opa, d, A, x, G, e_pow = ctx
        pk = np.array([spot_metrics(x, np.abs(A[j] @ G) ** 2)["peak"] for j in range(K)])
        amp = np.sqrt(pk.max() / pk)

        def errs(ph):
            I = np.abs((amp * np.exp(1j * ph)) @ A @ G) ** 2
            e = []
            for t in TARGETS:
                win = np.abs(x - t) < 0.6
                e.append(spot_metrics(x[win], I[win])["x_peak"] - t)
            return np.array(e) * 1e3

        single = []
        for j, t in enumerate(TARGETS):
            I = np.abs(A[j] @ G) ** 2
            win = np.abs(x - t) < 0.6
            single.append(spot_metrics(x[win], I[win])["x_peak"] - t)
        rng = np.random.default_rng(seed)
        rand = np.array([errs(rng.uniform(0, 2 * np.pi, K)) for _ in range(n_trials)])
        out[case] = {"single_port_max_nm": float(np.max(np.abs(single)) * 1e3),
                     "all_lit_equal_phases_max_nm": float(np.max(np.abs(errs(np.zeros(K))))),
                     "all_lit_random_phases_max_nm": float(np.max(np.abs(rand))),
                     "all_lit_random_phases_rms_nm": float(np.sqrt(np.mean(rand**2))),
                     "all_lit_random_phases_outer_rms_nm": float(np.sqrt(np.mean(rand[:, [0, -1]] ** 2))),
                     "all_lit_random_phases_inner_max_nm": float(np.max(np.abs(rand[:, 2:-2])))}
        print(case, out[case], flush=True)
    save_result("m8_all_lit", script=SCRIPT,
                inputs={"design": "FINAL (128 arms, e^-6 envelope)", "trials": n_trials, "seed": seed,
                        "window": "+/-0.6 um around each target"},
                outputs=out)
    return out


def summary(designs, multis, steers, rs, pw):
    tab = {
        "A_shared_synthesis_1064": {
            "steering": "independent per tweezer (computed phases)",
            "on_off": "re-synthesis, PZT settle ~2 us; residual at switched-off site "
                      f"{max(r['residual_rel'] for r in rs):.1e}; kept tweezers change by "
                      f"<= {max(r['kept_peak_change_max'] for r in rs):.1%}",
            "extinction_limit": "set by synthesis stray (~1e-3 to 1e-2); not adequate for 50 dB gating",
            "channels": "N = 128 PZT phase shifters", "switches": "none"},
        "C_star_coupler": {
            "steering": "common to all lit tweezers (parallel move); site selection by lit ports",
            "on_off": "per-port MZI: thermo-optic at 1064.4 nm (adiabatic ramps), cascaded AlN at 894.6 nm (ns)",
            "extinction_limit": "switch ER (>= 50 dB cascaded) and neighbours' optical tails",
            "channels": "N arm phase shifters + 20 port switches + 31 router MZIs per wavelength", "switches": "20"},
    }
    return tab


if __name__ == "__main__":
    parts = sys.argv[1:] or ["design", "resynth"]
    if "alllit" in parts:
        all_lit_phases()
    designs, multis, steers = {}, {}, {}
    if "design" in parts:
        for case in ("1064", "895"):
            for N in (96, 128):
                dsg, ctx = design(case, N)
                key = f"{case}_N{N}"
                ok = abs(dsg["rows"][0]["x_target"]) < dsg["x_range_um"] * 0.9
                ml, amp = multi(ctx)
                st = steer(ctx, amp)
                designs[key], multis[key], steers[key] = dsg, ml, st
                r = dsg["rows"]
                print(key, f"R={dsg['R_um']:.0f} um, port pitch {dsg['port_pitch_um']:.2f} um, range +/-{dsg['x_range_um']:.0f}",
                      f"capture {min(x['capture'] for x in r):.3f}-{max(x['capture'] for x in r):.3f}",
                      f"delivered {min(x['delivered_frac'] for x in r):.3f}-{max(x['delivered_frac'] for x in r):.3f}",
                      f"fwhm {np.mean([x['fwhm'] for x in r]):.3f}",
                      f"pos err max {max(abs(x['x_peak']-x['x_target']) for x in r)*1e3:.1f} nm",
                      f"nb {max(x['nb_max'] for x in r):.2e} nnb {max(x['nnb_max'] for x in r):.2e}", flush=True)
                print("   multi", ml, flush=True)
                print("   steer", st, flush=True)
        f, ax = fig(1, 2, h=3.0)
        for j, case in enumerate(("1064", "895")):
            for N, c in ((96, "C0"), (128, "C1")):
                r = designs[f"{case}_N{N}"]["rows"]
                ax[0, j].plot([x["x_target"] for x in r], [x["delivered_frac"] for x in r], "-o", color=c, ms=3,
                              label=f"N={N} delivered")
                ax[0, j].semilogy([x["x_target"] for x in r], [max(x["nb_max"], 1e-8) for x in r], "--s", color=c, ms=3,
                                  label=f"N={N} neighbour/peak")
            ax[0, j].axhline(1e-4, color="k", lw=0.6)
            ax[0, j].set_title(f"{CASES[case]['lam']*1e3:.1f} nm, 20 ports")
            ax[0, j].set_xlabel("tweezer site (um)")
            ax[0, j].legend(fontsize=6)
        save_result("m8_star_design", script=SCRIPT, fig=f,
                    inputs={"K": K, "site_pitch": prm.SITE_PITCH, "NA": prm.NA_OBJ,
                            "arm_pitch_FPR": "CHOICE 2.0 um", "arm_width_FPR": "CHOICE 1.7 um (0.3 um gap)",
                            "emitter_pitch": "CHOICE 1.5 um (>= 1.46 um end-fire minimum, m3_pitch_efficiency)",
                            "FPR_radius": "set so the arm-power envelope is 1/e^2 at the outer arms",
                            "model": "picsim.opa.star (2D paraxial FPR) + picsim.opa.array relay"},
                    outputs={"designs": designs},
                    notes="delivered_frac: fraction of port input power reaching the atom plane through FPR, arms, "
                          "emitters, slit and pupil (meta-objective transmission not included).")
        save_result("m8_star_multi", script=SCRIPT,
                    inputs={"relative_port_phases": "uniformly random (unstabilised), 40 draws",
                            "port_powers": "equalised single-port peaks"},
                    outputs={"multi": multis})
        save_result("m8_star_steer", script=SCRIPT, inputs={"all_ports_lit": True}, outputs={"steer": steers})
        pw = power({k.split("_")[0]: v for k, v in designs.items() if k.endswith("N128")})
        print(json.dumps(pw, indent=1), flush=True)
        save_result("m8_power", script=SCRIPT,
                    inputs={"depth": prm.DEPTH_TRANSPORT, "alpha": prm.ALPHA_CS_1064, "w0_y": prm.W0_TWEEZER,
                            "P_site_addr": prm.P_SITE_ADDR, "meta_objective": "CHOICE 0.7",
                            "router": "1x32 tunable-splitter tree, 3 um 2x2 MMIs (SIM m2_mmi)",
                            "facet": "SIM m2_couplers, lensed fibre"},
                    outputs=pw,
                    notes="Peak intensity from the simulated port spot along x and a Gaussian (w0 = 1.04 um) along y. "
                          "The tunable router sends power only to lit ports, so per-tweezer power is independent of the "
                          "number lit.")
    rs = None
    if "resynth" in parts:
        rs = resynth()
        for r in rs:
            print("resynth", r, flush=True)
        save_result("m8_resynth", script=SCRIPT,
                    inputs={"N": 128, "pitch": "CHOICE 1.5 um", "K": K, "method": "synthesise (GSW + refine), then "
                            "refine with the switched-off site among the dark points, warm-started"},
                    outputs={"rows": rs})
    if designs and rs:
        save_result("m8_summary", script=SCRIPT, inputs={}, outputs=summary(designs, multis, steers, rs, pw))
    if "final" in parts:
        fd, fm, fs = {}, {}, {}
        for case, cfg in FINAL.items():
            dsg, ctx = design(case, 128, edge=cfg["edge"], w_in=cfg["w_in"])
            ml, amp = multi(ctx, n_draw=60)
            st = steer(ctx, amp)
            fd[case], fm[case], fs[case] = dsg, ml, st
            r = dsg["rows"]
            print(case, f"R={dsg['R_um']:.0f} port pitch {dsg['port_pitch_um']:.2f} capture "
                  f"{min(x['capture'] for x in r):.3f}-{max(x['capture'] for x in r):.3f} delivered "
                  f"{min(x['delivered_frac'] for x in r):.3f}-{max(x['delivered_frac'] for x in r):.3f} fwhm "
                  f"{np.mean([x['fwhm'] for x in r]):.3f} pos err {max(abs(x['x_peak']-x['x_target']) for x in r)*1e3:.1f} nm "
                  f"nb {max(x['nb_max'] for x in r):.2e}", flush=True)
            print("  multi", ml, flush=True)
            print("  steer", st, flush=True)
        pw = power(fd)
        print(json.dumps(pw, indent=1), flush=True)
        # switching requirements
        nb895 = max(x["nb_max"] for x in fd["895"]["rows"])
        sw = {
            "895_neighbour_crosstalk": nb895,
            "895_neighbour_error_one_beam": float((np.pi**2 / 4) * nb895),
            "895_switch_ER_dB_for_idle_error_2.5e-5": 50.0,
            "895_switch_ER_dB_where_leak_equals_neighbour_tail": float(-10 * np.log10(nb895)),
            "895_switch": "two cascaded AlN MZIs per port (>= 50 dB in all trials, m4_extinction); ZVD drive",
            "1064_switch_ER_dB_for_residual_1uK_at_1.2mK": float(10 * np.log10(1.2e-3 / 1e-6)),
            "1064_switch": "one trimmed thermo-optic MZI per port (rise ~34 us, m5_heaters); ramps >= 100 us",
        }
        counts = {
            "1064": {"ports": 20, "router_MZIs_1x32": 31, "port_switches": 20, "arm_phase_shifters_PZT": 128,
                     "emitters": 128, "FPR_radius_um": fd["1064"]["R_um"]},
            "895": {"ports": 20, "router_MZIs_1x32": 31, "port_switches_AlN": 40, "arm_phase_shifters": 128,
                    "emitters": 128, "FPR_radius_um": fd["895"]["R_um"]},
        }
        save_result("m8_final", script=SCRIPT,
                    inputs={"FINAL": {k: {"edge": float(v["edge"]), "w_in": v["w_in"]} for k, v in FINAL.items()},
                            "N": 128, "p_out": P_OUT, "p_a": P_A, "w_arm": W_ARM, "K": K, "site_pitch": prm.SITE_PITCH},
                    outputs={"designs": fd, "multi": fm, "steer": fs, "power": pw, "switching": sw, "counts": counts},
                    notes="Adjacent 1.4 um port tapers: coupling length 1.33 mm at 1064.4 nm (0.49 um gap); at 894.6 nm "
                          "1.4 um tapers with a 0.45 um gap give 3.25 mm, and the final 1.0 um tapers have a wider "
                          "0.64 um gap (femwell, run in session).")
        print(json.dumps(sw, indent=1), flush=True)
    if "scan" in parts:
        rows = []
        for case, edges, widths in (("1064", (3, 4, 5, 6), (1.2, 1.4)), ("895", (6, 8), (1.0, 1.2, 1.4))):
            for e in edges:
                for w in widths:
                    dsg, ctx = design(case, 128, edge=np.exp(-e), w_in=w)
                    ml, _ = multi(ctx, n_draw=20)
                    r = dsg["rows"]
                    rows.append({"case": case, "edge_nominal": f"exp(-{e})", "w_in": w, "R_um": dsg["R_um"],
                                 "edge_realised": dsg["envelope_edge"], "port_pitch_um": dsg["port_pitch_um"],
                                 "gap_um": dsg["port_pitch_um"] - w, "fwhm_um": float(np.mean([x["fwhm"] for x in r])),
                                 "nb_max": float(max(x["nb_max"] for x in r)),
                                 "peak_per_W_mean": float(np.mean([x["peak_per_W_in"] for x in r])),
                                 "spread20_max": ml[0]["peak_spread_max"]})
                    print("scan", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
        save_result("m8_taper_scan", script=SCRIPT,
                    inputs={"N": 128, "edge": "nominal arm-power envelope at the outer arm",
                            "w_in": "input taper width at the FPR (um)"},
                    outputs={"rows": rows},
                    notes="Basis for the FINAL choice: 1064.4 nm exp(-6), 1.4 um; 894.6 nm exp(-6), 1.0 um.")
    if "apod" in parts:
        rows = apod_sweep()
        f, ax = fig(1, 2, h=3.0)
        for case, c in (("1064", "C0"), ("895", "C1")):
            r = [x for x in rows if x["case"] == case]
            ax[0, 0].semilogy([x["fwhm_um"] for x in r], [x["nb_max"] for x in r], "-o", color=c, ms=3, label=case)
            ax[0, 1].plot([x["fwhm_um"] for x in r], [x["spread20_max"] for x in r], "-o", color=c, ms=3, label=case)
        ax[0, 0].axhline(1e-4, color="k", lw=0.6)
        ax[0, 0].set_xlabel("spot FWHM along x (um)")
        ax[0, 0].set_ylabel("intensity at 4 um neighbour / peak")
        ax[0, 1].set_xlabel("spot FWHM along x (um)")
        ax[0, 1].set_ylabel("peak spread, 20 lit, random phases")
        ax[0, 0].legend()
        save_result("m8_apod", script=SCRIPT, fig=f,
                    inputs={"N": 128, "edge_levels": "arm-power envelope at the outer arm, exp(-2) ... exp(-8)"},
                    outputs={"rows": rows},
                    notes="Stronger apodisation narrows the port far field (smaller FPR radius and port pitch), "
                          "widens the spot and suppresses the tails that make lit ports interfere.")
