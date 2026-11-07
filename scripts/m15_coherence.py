"""Module 15: independent rows without per-chip acousto-optic frequency offsets.

The rows of the stacked relay interfere where one chip's tails reach the next row (module 13). Options that
need no per-chip frequency shifter:
  apod   a passive Gaussian y apodiser in the y-telescope plane conjugate to the objective pupil (common to all
         chips): crosstalk floor vs power cost (1D scan)
  coh    2D coherent perturbation of the neighbouring rows (worst of 8 relative phases), static and with a moving
         tweezer, for: baseline; apodisers; alternating linear polarisation (vector Richards-Wolf focusing,
         x / y on alternate chips); polarisation + apodiser
  prbs   common pseudo-random phase modulation of the seed before the 1:16 split, plus fibre delays: residual
         coherence (m-sequence, modulation-depth error, finite bandwidth, random multi-level phase), modulation
         spectrum, and chromatic smearing of the tweezers by the broadened spectrum
All with the adopted collimator (Weierstrass + aplanat, glued, AR coated), beta = 0.7 um unless stated, and the
D1 chip response (BPM matrix).

Outputs: results/m15_*.json (+ .png)
"""
from __future__ import annotations

import json
import sys

import numpy as np

from picsim import materials as mat
from picsim import params as prm
from picsim.opa import collimator as C
from picsim.opa.focal_plane import FocalPlane, LinearSynth
from picsim.opa.stack import FacetY, Focus2D, RelayY, fwhm, peak2d
from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m15_coherence.py"
LAM = prm.LAM_TRANSPORT.value
NA = prm.NA_OBJ.value
N_SIO2 = float(mat.n_sio2(LAM))
D1 = {"N_ports": 128, "N_arms": 192, "w_emit": 1.2, "p_out": 1.5, "w_port": 1.7, "w_arm": 1.7}
DIST = (4.0, 8.0, 12.0, 16.0)
_S = {}


def relay(beta=0.7):
    if "lens" not in _S:
        _S["lens"] = C.fac_design(N_SIO2).lensmap(LAM, 0.9, coated=True)
        _S["facet"] = FacetY(LAM, 1.2, 0.4, n_m=N_SIO2)
    return RelayY(_S["facet"], beta, NA=NA, NA_FAC=0.9, lens=_S["lens"])


def gauss_apod(a):
    return (lambda eta: np.ones_like(eta)) if a is None else (lambda eta: np.exp(-(eta / a) ** 2))


def eps_at(y, I, d, half=0.5):
    m = np.abs(np.abs(y) - d) <= half
    return float(I[m].max() / I.max())


# ------------------------------------------------------------------ apodiser scan (1D)
def apod():
    y = np.arange(-20, 20, 0.005)
    base = relay(0.7)
    I0 = base.atom_line(y)
    pk0 = I0.max()
    rows = []
    for a in (None, 1.0, 0.8, 0.65, 0.55, 0.45):
        for beta in (0.5, 0.6, 0.7, 0.8, 1.0):
            r = relay(beta)
            Py = r.pupil() * gauss_apod(a)(r.eta)
            I = r.atom_line(y, Py)
            row = {"a": a, "beta_um": beta, "transmission": r.efficiency(Py), "fwhm_um": fwhm(y, I),
                   "peak_rel": float(I.max() / pk0), "eps": {str(d): eps_at(y, I, d) for d in (3.0, 4.0, 6.0, 8.0, 12.0, 16.0)}}
            rows.append(row)
            print(a, beta, round(row["transmission"], 3), round(row["fwhm_um"], 3), round(row["peak_rel"], 3),
                  {k: f"{v:.1e}" for k, v in row["eps"].items()}, flush=True)
    save_result("m15_apod", script=SCRIPT,
                inputs={"apodiser": "amplitude exp(-(eta/a)^2) in a plane conjugate to the pupil (y only)",
                        "collimator": "Weierstrass + aplanat, AR", "reference": "beta 0.7 um, no apodiser"},
                outputs={"rows": rows},
                notes="1D y line; peak_rel is per unit guided power relative to the reference (power per tweezer "
                      "scales as its inverse).")
    return rows


# ------------------------------------------------------------------ 2D coherent metrics
def setup2d():
    if "fp" not in _S:
        T = np.load(RESULTS / "m10_T_opt_A.npz")["T"]
        _S["fp"] = FocalPlane(A=T, **D1)
        _S["ls"] = LinearSynth(_S["fp"])
    return _S["fp"], _S["ls"]


def coherent_metrics(F, Py, pol_m, pol_n, dist=DIST):
    """Row m (41 tweezers at 3 um over +/-60 um) and 11 tweezers of row n = m + 1 (offset 1.5 um in x) at row
    separation d. pol_*: None (scalar), 'x' or 'y' (vector). Returns, per d, the worst over 8 relative phases of
    the static and moving (row-m tweezer +0.25 um) changes of the row-n tweezers, plus the incoherent values."""
    fp, ls = setup2d()
    xt = np.arange(-60.0, 60.0 + 1e-9, 3.0)
    th = np.where(np.arange(len(xt)) % 2 == 0, 0.0, np.pi / 2)
    xs_n = xt[15:26] + 1.5
    xx = np.arange(xs_n[0] - 1.5, xs_n[-1] + 1.5, 0.02)
    y0 = np.arange(-1.2, 1.2, 0.02)
    xm2 = xt.copy()
    xm2[20] += 0.25

    def fld(w, yy, pol):
        if pol is None:
            return F.field(w, xx, yy, Py=Py)[None]
        return F.field_vector(w, xx, yy, pol=pol, Py=Py)

    En = fld(ls.weights(xs_n, phases=np.where(np.arange(len(xs_n)) % 2 == 0, 0.0, np.pi / 2)), y0, pol_n)
    yall = np.concatenate([y0 + d for d in dist])
    Em_all = fld(ls.weights(xt, phases=th), yall, pol_m)
    Em2_all = fld(ls.weights(xm2, phases=th), yall, pol_m)

    def tweezers(I):
        out = []
        for x0 in xs_n:
            m = np.abs(xx - x0) < 0.6
            out.append(peak2d(I[:, m], xx[m], y0))
        return np.array(out)

    In = np.sum(np.abs(En) ** 2, axis=0)
    base = tweezers(In)
    res = {}
    for i, d in enumerate(dist):
        sl = slice(i * len(y0), (i + 1) * len(y0))
        Em, Em2 = Em_all[:, sl], Em2_all[:, sl]
        r = {}
        for mode in ("incoherent", "coherent"):
            dpk = dpos = mpk = mpos = 0.0
            for ph in ([0.0] if mode == "incoherent" else np.linspace(0, 2 * np.pi, 8, endpoint=False)):
                if mode == "incoherent":
                    Ia = In + np.sum(np.abs(Em) ** 2, axis=0)
                    Ib = In + np.sum(np.abs(Em2) ** 2, axis=0)
                else:
                    Ia = np.sum(np.abs(En + np.exp(1j * ph) * Em) ** 2, axis=0)
                    Ib = np.sum(np.abs(En + np.exp(1j * ph) * Em2) ** 2, axis=0)
                a, b = tweezers(Ia), tweezers(Ib)
                dpk = max(dpk, np.max(np.abs(a[:, 0] / base[:, 0] - 1)))
                dpos = max(dpos, np.max(np.hypot(a[:, 1] - base[:, 1], a[:, 2] - base[:, 2])))
                mpk = max(mpk, np.max(np.abs(b[:, 0] / a[:, 0] - 1)))
                mpos = max(mpos, np.max(np.hypot(b[:, 1] - a[:, 1], b[:, 2] - a[:, 2])))
            r[mode] = {"static_peak": float(dpk), "static_pos_nm": float(dpos * 1e3), "move_peak": float(mpk),
                       "move_pos_nm": float(mpos * 1e3)}
        res[str(d)] = r
    return res, float(base[:, 0].mean())


def coh():
    fp, ls = setup2d()
    out = {}
    configs = [("baseline, scalar", 0.7, None, None, None),
               ("baseline, vector, same polarisation", 0.7, None, "x", "x"),
               ("alternating polarisation (x / y), vector", 0.7, None, "x", "y"),
               ("apodiser a=0.8, beta 0.7", 0.7, 0.8, None, None),
               ("apodiser a=0.65, beta 0.7", 0.7, 0.65, None, None),
               ("alternating polarisation + apodiser a=0.8, beta 0.7", 0.7, 0.8, "x", "y")]
    pk_ref = None
    for name, beta, a, pm, pn in configs:
        r = relay(beta)
        F = Focus2D(fp, r)
        Py = r.pupil() * gauss_apod(a)(r.eta)
        res, pk = coherent_metrics(F, Py, pm, pn)
        if pk_ref is None:
            pk_ref = pk
        out[name] = {"beta_um": beta, "apod_a": a, "pol": [pm, pn], "row_peak_rel": pk / pk_ref, "by_distance": res}
        print(name, round(pk / pk_ref, 3), json.dumps({d: {m: {k: round(v, 5) for k, v in vv.items()}
                                                           for m, vv in x.items()} for d, x in res.items()}), flush=True)
    save_result("m15_coherent", script=SCRIPT,
                inputs={"row_m": "41 tweezers, 3 um, phases 0/pi/2", "row_n": "11 tweezers offset 1.5 um",
                        "move": "row-m tweezer +0.25 um", "phases": "worst of 8 relative chip phases",
                        "collimator": "Weierstrass + aplanat, AR", "vector": "Richards-Wolf, aplanatic objective, "
                        "intensity summed over Ex, Ey, Ez (scalar light shift)"},
                outputs=out,
                notes="With alternating polarisation, rows 4 and 12 um apart are orthogonal, 8 and 16 um parallel.")
    return out


# ------------------------------------------------------------------ PRBS phase modulation
def mseq(deg):
    taps = {7: (7, 6), 9: (9, 5), 10: (10, 7), 11: (11, 9), 13: (13, 12, 10, 9)}[deg]
    state = [1] * deg
    out = []
    for _ in range(2**deg - 1):
        out.append(state[-1])
        fb = 0
        for t in taps:
            fb ^= state[t - 1]
        state = [fb] + state[:-1]
    return np.array(out)


def coherence_residual(phase, k):
    """|<exp(i (phi(t) - phi(t - k samples)))>| over one period: residual mutual coherence of two chips whose
    feeds differ by k samples of delay."""
    return float(np.abs(np.mean(np.exp(1j * (phase - np.roll(phase, k))))))


def smooth(ph, os_, rise):
    """Band-limit a phase waveform (Gaussian impulse response with 10-90% rise of `rise` bits)."""
    if not rise:
        return ph
    sig = rise * os_ / 2.563
    k = np.exp(-0.5 * (np.arange(-int(4 * sig) - 1, int(4 * sig) + 2) / sig) ** 2)
    pad = len(k)
    return np.convolve(np.concatenate([ph[-pad:], ph, ph[:pad]]), k / k.sum(), mode="same")[pad:-pad]


def zadoff_chu(N, u=1):
    k = np.arange(N)
    return np.pi * u * k * (k + 1) / N  # N odd


def unwrap_levels(ph):
    """Take the shortest phase path between successive levels (what a phase modulator driven continuously does)."""
    d = np.angle(np.exp(1j * np.diff(ph)))
    return np.concatenate([[ph[0]], ph[0] + np.cumsum(d)])


def gauss_noise_phase(rng, n, os_, sigma):
    """Gaussian phase noise of rms sigma (rad), correlation time ~1 bit (os_ samples)."""
    w = rng.standard_normal(n * os_ + 8 * os_)
    k = np.exp(-0.5 * (np.arange(-3 * os_, 3 * os_ + 1) / (os_ / 2)) ** 2)
    x = np.convolve(w, k, mode="same")[4 * os_:4 * os_ + n * os_]
    return x / x.std() * sigma


def prbs():
    """Common phase modulation of the seed + per-chip fibre delays: residual mutual coherence of two chips."""
    rng = np.random.default_rng(11)
    out = {}
    os_ = 16
    rows = []
    for name, maker in (
            ("binary 0/pi m-sequence, N=1023", lambda: np.repeat(mseq(10) * np.pi, os_).astype(float)),
            ("Zadoff-Chu, N=1021", lambda: np.repeat(unwrap_levels(zadoff_chu(1021)), os_)),
            ("Gaussian phase noise, 2 rad rms", lambda: gauss_noise_phase(rng, 2**15, os_, 2.0)),
            ("Gaussian phase noise, 3 rad rms", lambda: gauss_noise_phase(rng, 2**15, os_, 3.0)),
            ("Gaussian phase noise, 4 rad rms", lambda: gauss_noise_phase(rng, 2**15, os_, 4.0))):
        for rise in (0.0, 0.3, 0.6):
            ph = smooth(maker(), os_, rise) if "noise" not in name else maker()
            res_int = max(coherence_residual(ph, d * os_) for d in range(2, 16))
            res_any = max(coherence_residual(ph, int(d * os_)) for d in np.arange(2.0, 16.0, 0.37))
            carrier = float(np.abs(np.mean(np.exp(1j * ph))))
            rows.append({"scheme": name, "rise_bits": rise, "residual_integer_delays": res_int,
                         "residual_any_delay": res_any, "carrier": carrier})
            print(rows[-1], flush=True)
            if "noise" in name:
                break
    # depth error on the binary sequence
    for e in (0.05, 0.10):
        ph = np.repeat(mseq(10) * np.pi * (1 + e), os_).astype(float)
        rows.append({"scheme": f"binary 0/pi m-sequence, depth +{int(e*100)}%", "rise_bits": 0.0,
                     "residual_integer_delays": max(coherence_residual(ph, d * os_) for d in range(2, 16)),
                     "residual_any_delay": None, "carrier": float(np.abs(np.mean(np.exp(1j * ph))))})
    out["schemes"] = rows
    # optical spectrum of Gaussian phase noise (3 rad, correlation ~1 bit): rms width in units of the bit rate
    ph = gauss_noise_phase(rng, 2**15, os_, 3.0)
    S = np.abs(np.fft.fftshift(np.fft.fft(np.exp(1j * ph)))) ** 2
    fq = np.fft.fftshift(np.fft.fftfreq(len(ph), d=1 / os_))  # in units of B
    S /= S.sum()
    out["noise3_spectrum_rms_width_over_B"] = float(np.sqrt(np.sum(S * fq**2)))
    # low-frequency density of the interference term (heating at 2 nu_trap) for the 3 rad noise
    cross = np.cos(ph - np.roll(ph, 8 * os_))
    Sc = np.abs(np.fft.rfft(cross - cross.mean())) ** 2 / len(cross)  # per sample^2; one-sided density x 2 / f_s
    out["noise3_cross_term_lowfreq_density_x_B"] = float(2 * np.mean(Sc[1:50]) / os_)
    out["noise3_spectrum_fq"] = fq[::64].tolist()
    out["noise3_spectrum_S"] = S[::64].tolist()
    # chromatic smearing: axial scan of a single tweezer, averaged over the sinc^2 spectrum
    fp, ls = setup2d()
    r = relay(0.7)
    F = Focus2D(fp, r)
    w = ls.weights([0.0])
    w = w / np.linalg.norm(w)
    zs = np.arange(-4.0, 4.01, 0.05)
    x = np.arange(-0.6, 0.6, 0.02)
    yv = np.arange(-0.8, 0.8, 0.02)
    Iz = np.array([F.intensity(w, x, yv, dz_a=z).max() for z in zs])
    dzdl = 16.0  # um of focal shift per nm (diffractive objective, f proportional to 1/lambda; DERIVED, f = 17.3 mm)
    chrom = []
    nu0 = 2.998e8 / (LAM * 1e-6)
    wB = out["noise3_spectrum_rms_width_over_B"]
    for B in (0.5e9, 2e9, 10e9):
        f = np.linspace(-6 * wB * B, 6 * wB * B, 2401)
        S = np.exp(-0.5 * (f / (wB * B)) ** 2)  # Gaussian optical spectrum of the phase-noise-modulated field
        S /= S.sum()
        dl_nm = -f / nu0 * LAM * 1e3
        dz = dl_nm * dzdl
        Iavg = np.array([np.sum(S * np.interp(z - dz, zs, Iz, left=0, right=0)) for z in zs])
        x60 = 60.0 * np.sqrt(np.sum(S * (dl_nm / (LAM * 1e3)) ** 2)) * 1e3  # rms lateral smear at 60 um, nm
        chrom.append({"B_GBd": B / 1e9, "optical_rms_width_GHz": wB * B / 1e9, "peak_rel": float(Iavg.max() / Iz.max()),
                      "axial_rms_um": float(np.sqrt(np.sum(S * dz**2))), "lateral_rms_at_60um_nm": float(x60)})
        print(chrom[-1], flush=True)
    out["chromatic"] = chrom
    save_result("m15_prbs", script=SCRIPT,
                inputs={"modulation": "common phase modulation of the seed before the 1:16 split; per-chip fibre "
                        "delays", "samples_per_bit": 16, "B": "symbol rate (binary, Zadoff-Chu) or inverse correlation "
                        "time (Gaussian noise)", "objective_chromatic": "16 um/nm axial (diffractive), lateral x "
                        "proportional to lambda", "fibre_group_index": 1.45},
                outputs=out,
                notes="Residual coherence multiplies the coherent crosstalk terms of m15_coherent. Delays from 2 to 15 "
                      "bits (integer, or arbitrary in steps of 0.37 bit). Cross-term low-frequency density in units of 1/B.")
    return out


# ------------------------------------------------------------------ summary
def summary():
    """Perturbation of a tweezer by all other rows within +/-16 um (both sides), per option.

    Each row pair contributes its worst-phase (coherent) or incoherent value from m15_coherent; with drifting
    relative phases the rms over phases is worst / sqrt(2), and independent pairs add in quadrature.
    """
    co = json.loads((RESULTS / "m15_coherent.json").read_text())["outputs"]
    pr = json.loads((RESULTS / "m15_prbs.json").read_text())["outputs"]
    noise = [r for r in pr["schemes"] if r["scheme"] == "Gaussian phase noise, 3 rad rms"][0]
    resid = noise["residual_any_delay"]  # finite-record upper bound; ensemble value exp(-9) = 1.2e-4

    def total(pairs, key_pk, key_pos, scale=1.0, coherent=True):
        f = 1 / np.sqrt(2) if coherent else 1.0
        pk = np.sqrt(2 * sum((p[key_pk] * scale * f) ** 2 for p in pairs))
        pos = np.sqrt(2 * sum((p[key_pos] * scale * f) ** 2 for p in pairs))
        return float(pk), float(pos)

    def pairs(cfg, mode, dmap=None):
        out = []
        for d in ("4.0", "8.0", "12.0", "16.0"):
            src = co[dmap[d]] if dmap else co[cfg]
            out.append(src["by_distance"][d][mode])
        return out

    same = "baseline, vector, same polarisation"
    alt = "alternating polarisation (x / y), vector"
    rows = []

    def add(name, prs, coherent=True, scale=1.0, power=1.0, active=""):
        s_pk, s_pos = total(prs, "static_peak", "static_pos_nm", scale, coherent)
        m_pk, m_pos = total(prs, "move_peak", "move_pos_nm", scale, coherent)
        rows.append({"option": name, "active_parts": active, "power_factor": power, "static_depth_rms": s_pk,
                     "static_pos_rms_nm": s_pos, "move_depth_rms": m_pk, "move_pos_rms_nm": m_pos})

    add("all chips coherent (no measure)", pairs(same, "coherent"))
    add("per-chip AOM offsets (incoherent)", pairs(same, "incoherent"), coherent=False, active="16 AOMs at 2.6-3.4 W")
    # residual coherence r scales only the interference part: incoherent + r (worst-phase coherent - incoherent)
    inc, coh = pairs(same, "incoherent"), pairs(same, "coherent")
    mixed = [{k: abs(i[k] + resid * (c[k] - i[k])) for k in i} for i, c in zip(inc, coh)]
    add("common phase noise (3 rad) + fibre delays", mixed, coherent=False,
        active="1 EOM on the seed, noise source, <= 5 m extra fibre")
    add("alternating polarisation", pairs(None, "coherent", {"4.0": alt, "8.0": same, "12.0": alt, "16.0": same}),
        active="none (half-wave plate strips)")
    apo = "apodiser a=0.65, beta 0.7"
    add("y apodiser a=0.65", pairs(apo, "coherent"), power=1 / co[apo]["row_peak_rel"], active="none (graded mirror)")
    pa = "alternating polarisation + apodiser a=0.8, beta 0.7"
    a8 = "apodiser a=0.8, beta 0.7"
    add("polarisation + apodiser a=0.8", pairs(None, "coherent", {"4.0": pa, "8.0": a8, "12.0": pa, "16.0": a8}),
        power=1 / co[pa]["row_peak_rel"], active="none")
    for r in rows:
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    save_result("m15_summary", script=SCRIPT,
                inputs={"rows": "both sides, separations 4, 8, 12, 16 um", "phase": "rms over drifting relative phases = "
                        "worst / sqrt(2), pairs in quadrature", "noise_residual": resid,
                        "apodiser_parallel_rows": "same-polarisation values of the apodised scalar case"},
                outputs={"rows": rows},
                notes="Criteria used elsewhere: others' depth change < 3e-3 and position < 5 nm when a tweezer moves; "
                      "position noise budget 10 nm rms (proposal).")
    return rows


if __name__ == "__main__":
    parts = sys.argv[1:] or ["apod", "coh", "prbs", "summary"]
    for p in parts:
        {"apod": apod, "coh": coh, "prbs": prbs, "summary": summary}[p]()
