"""Module 16: conversion of the decorrelating phase noise into intensity noise (PM-to-AM) along the laser chain.

The seed carries Gaussian phase noise (3 rad rms, correlation rate B = 2 GHz; docs/row_coherence.md). A pure
phase modulation leaves the intensity constant; any spectral filtering converts part of it to intensity noise.
Each element below is applied in the frequency domain to a long record of the modulated field, and the relative
intensity noise (RIN) density is evaluated at 2 nu_trap = 168 kHz, where it drives parametric heating
Gamma = pi^2 nu^2 S_RIN(2 nu) (O'Hara et al., arXiv:physics/0003049; nu = 84 kHz radial at 1.2 mK).

Elements (parameters are CHOICE ranges, not measured values):
  gain_slope   amplifier gain tilt g1 (dB/nm) across the modulated spectrum
  bandpass     Gaussian band-pass (ASE filter, fibre Bragg grating) of FWHM BW, centred or offset
  multipath    a delayed copy of amplitude r = sqrt(R1 R2) and delay tau (two reflections; fibre or on chip)
  pm_fibre     polarisation crosstalk (extinction PER) in PM fibre of length L, then a polariser (the chip TE
               guide): a copy of amplitude sqrt(PER) delayed by the fibre's differential group delay
  dispersion   group-velocity dispersion beta2 L of the delivery fibre
  chain        a combination of representative values

Output: results/m16_pm_am.json (+ .png)
"""
from __future__ import annotations

import sys

import numpy as np

from picsim.plotting import fig
from picsim.provenance import save_result

SCRIPT = "scripts/m16_pm_am.py"
LAM = 1.064435e-6
C0 = 2.998e8
NU_TRAP = 84e3
S_MAX = 0.1 / (np.pi**2 * NU_TRAP**2)  # RIN density for Gamma = 0.1 /s (CHOICE criterion)


def phase_noise(rng, n, fs, B, sigma):
    """Gaussian phase noise, rms sigma, Gaussian correlation of width 1/(2B) (as in m15)."""
    w = rng.standard_normal(n)
    f = np.fft.rfftfreq(n, 1 / fs)
    # Gaussian impulse response of sigma_t = 1/(2B)  ->  transfer exp(-(2 pi f sigma_t)^2 / 2)
    H = np.exp(-0.5 * (2 * np.pi * f / (2 * B)) ** 2)
    x = np.fft.irfft(np.fft.rfft(w) * H, n)
    return x / x.std() * sigma


def rin_density(I, fs, f_lo=100e3, f_hi=300e3):
    """One-sided RIN density (1/Hz) averaged over [f_lo, f_hi], and total RIN rms."""
    d = I / I.mean() - 1
    F = np.fft.rfft(d)
    f = np.fft.rfftfreq(len(d), 1 / fs)
    S = 2 * np.abs(F) ** 2 / (fs * len(d))
    m = (f >= f_lo) & (f <= f_hi)
    return float(S[m].mean()), float(d.std()), f, S


def run(H_of_nu, n=2**21, fs=16e9, B=2e9, sigma=3.0, reps=3, seed=1):
    rng = np.random.default_rng(seed)
    nu = np.fft.fftfreq(n, 1 / fs)
    Hn = H_of_nu(nu)
    out = []
    for _ in range(reps):
        E = np.exp(1j * phase_noise(rng, n, fs, B, sigma))
        Eo = np.fft.ifft(np.fft.fft(E) * Hn)
        s, rms, _, _ = rin_density(np.abs(Eo) ** 2, fs)
        out.append((s, rms))
    s = float(np.mean([o[0] for o in out]))
    return {"S_2nu_per_Hz": s, "S_dB_per_Hz": float(10 * np.log10(max(s, 1e-300))), "RIN_rms": float(np.mean([o[1] for o in out])),
            "Gamma_per_s": float(np.pi**2 * NU_TRAP**2 * s), "meets_0.1_per_s": bool(s <= S_MAX)}


def dl(nu):
    return -LAM**2 * nu / C0 * 1e9  # wavelength offset in nm for an optical frequency offset nu


def H_gain(g1):
    return lambda nu: 10 ** (g1 * dl(nu) / 20)


def H_band(bw_pm, off_pm=0.0):
    sig = bw_pm * 1e-3 / 2.3548
    return lambda nu: np.exp(-0.25 * ((dl(nu) - off_pm * 1e-3) / sig) ** 2)  # amplitude: intensity FWHM = bw


def H_multi(r, tau, theta=0.0):
    return lambda nu: 1 + r * np.exp(1j * theta) * np.exp(-2j * np.pi * nu * tau)


def H_disp(beta2_ps2_km, L_m):
    b = beta2_ps2_km * 1e-24 / 1e3 * L_m
    return lambda nu: np.exp(0.5j * b * (2 * np.pi * nu) ** 2)


def main():
    out = {"criterion": {"Gamma_max_per_s": 0.1, "S_max_per_Hz": S_MAX, "S_max_dB_per_Hz": float(10 * np.log10(S_MAX))},
           "modulation": {"sigma_rad": 3.0, "B_GHz": 2.0}}
    out["check_no_filter"] = run(lambda nu: np.ones_like(nu), reps=1)
    rows = []

    def add(kind, label, H, **kw):
        r = run(H, **kw)
        rows.append({"element": kind, "case": label, **r})
        print(kind, label, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)

    for g in (0.05, 0.5, 5.0):
        add("gain_slope", f"{g} dB/nm", H_gain(g))
    for bw, off in ((100.0, 0.0), (100.0, 25.0), (30.0, 0.0), (30.0, 7.5), (1000.0, 250.0)):
        add("bandpass", f"FWHM {bw} pm, offset {off} pm", H_band(bw, off))
    for R in (1e-2, 1e-3, 1e-4):
        for tau in (1e-11, 1e-10, 3e-10, 1e-9, 1e-8, 5e-8):
            r_amp = R  # r = sqrt(R1 R2) with R1 = R2 = R
            worst = None
            for th in (0.0, np.pi / 2):
                rr = run(H_multi(r_amp, tau, th), reps=2)
                if worst is None or rr["S_2nu_per_Hz"] > worst["S_2nu_per_Hz"]:
                    worst = rr
            rows.append({"element": "multipath", "case": f"R1 = R2 = {R:g}, tau = {tau*1e9:g} ns", **worst})
            print("multipath", R, tau, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in worst.items()}, flush=True)
    for per_db in (-20.0, -25.0, -30.0):
        for L in (1.0, 5.0, 20.0):
            tau = 1.2e-12 * L  # PM fibre DGD ~1.2 ps/m at 1064 nm (beat length ~3 mm, RECALLED)
            worst = None
            for th in (0.0, np.pi / 2):
                rr = run(H_multi(10 ** (per_db / 20), tau, th), reps=2)
                if worst is None or rr["S_2nu_per_Hz"] > worst["S_2nu_per_Hz"]:
                    worst = rr
            rows.append({"element": "pm_fibre", "case": f"PER {per_db} dB, L = {L} m", **worst})
            print("pm_fibre", per_db, L, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in worst.items()}, flush=True)
    for L in (5.0, 100.0, 1000.0):
        add("dispersion", f"beta2 20 ps^2/km, L = {L} m", H_disp(20.0, L))

    # representative chain: amplifier tilt, ASE band-pass, PM fibre, two connector reflections, on-chip etalon
    def chain(nu, th1=np.pi / 2, th2=np.pi / 2):
        return (H_gain(0.5)(nu) * H_band(1000.0, 250.0)(nu) * H_multi(10 ** (-20 / 20), 6e-12, np.pi / 2)(nu)
                * H_multi(1e-3, 2.5e-8, th1)(nu) * H_multi(np.sqrt(1e-3 * 0.023), 1.9e-10, th2)(nu) * H_disp(20.0, 5.0)(nu))
    worst = None
    for th1 in (0.0, np.pi / 2):
        for th2 in (0.0, np.pi / 2):
            rr = run(lambda nu: chain(nu, th1, th2), reps=2)
            if worst is None or rr["S_2nu_per_Hz"] > worst["S_2nu_per_Hz"]:
                worst = rr
    rows.append({"element": "chain", "case": "0.5 dB/nm tilt; 1 nm ASE filter offset 250 pm; PM fibre PER -20 dB, 5 m; "
                 "two -30 dB reflections 5 m apart; on-chip etalon (input facet -30 dB, glued emitter facet 2.3%, "
                 "190 ps); 5 m dispersion", **worst})
    print("chain", worst, flush=True)
    # B dependence for the worst multipath delay class
    bdep = []
    for B in (0.5e9, 2e9, 8e9):
        for tau in (1e-10, 1e-9):
            rr = run(H_multi(1e-2, tau, np.pi / 2), B=B, fs=max(16e9, 8 * B), reps=2)
            bdep.append({"B_GHz": B / 1e9, "tau_ns": tau * 1e9, "r": 1e-2, **rr})
            print("B", B, tau, rr, flush=True)
    out["rows"] = rows
    out["B_dependence_multipath_r0.01"] = bdep
    # figure: multipath S vs tau
    f, ax = fig(1, 1, h=3.0)
    for R in (1e-2, 1e-3, 1e-4):
        sel = [r for r in rows if r["element"] == "multipath" and r["case"].startswith(f"R1 = R2 = {R:g},")]
        taus = [float(r["case"].split("tau = ")[1].split(" ns")[0]) for r in sel]
        ax[0, 0].loglog(taus, [max(r["S_2nu_per_Hz"], 1e-22) for r in sel], "-o", ms=3, label=f"R1 = R2 = {R:g}")
    ax[0, 0].axhline(S_MAX, color="k", lw=0.8, ls="--", label="Gamma = 0.1 /s")
    ax[0, 0].set_xlabel("delay of the parasitic copy (ns)")
    ax[0, 0].set_ylabel("RIN density at 168 kHz (1/Hz)")
    ax[0, 0].legend(fontsize=7)
    save_result("m16_pm_am", script=SCRIPT, fig=f,
                inputs={"record": "2^21 samples at 16 GS/s (131 us), 2-3 realisations", "band": "100-300 kHz",
                        "phase_noise": "Gaussian, 3 rad rms, B = 2 GHz", "nu_trap_Hz": NU_TRAP,
                        "element_values": "CHOICE ranges; PM-fibre DGD 1.2 ps/m and beta2 20 ps^2/km RECALLED"},
                outputs=out,
                notes="Multipath and PM-fibre values are the worse of relative phases 0 and pi/2 of the copy.")
    return out


if __name__ == "__main__":
    main()
