"""Module 31: addressing at 894.6 nm. The PIC switch router against crossed AODs and a free-space modulator (Zhang et al.).

Parts:
  patterns   arbitrary on/off patterns on an n x n block: crossed AODs light row-set x column-set rectangles, so a
             pattern needs as many rounds as a partition of its ones into rectangles (lower bound: matrix rank over
             the reals; upper bound: greedy partition, distinct rows, distinct columns). The PIC router (one switch
             per site) and a free-space modulator (Zhang et al.) need one round. Time per pattern from the access / switching times
  leakage    idle-site light: PIC (m4 cascade, measured device), free-space modulator (Zhang et al.), crossed AODs with k tones
             (third-order intermodulation spots at 2 f_i - f_j land on sites of the same line; level swept)
  frequency  per-site frequency tags of an AOD (tone spacing = bandwidth / resolvable spots) and the phase they
             accumulate in coherent (Raman) addressing; the PIC and the DMD keep one frequency

Outputs: results/m31_patterns.json, m31_leakage.json, m31_frequency.json (+ m31_patterns.png)
"""
from __future__ import annotations

import json
import sys

import numpy as np

from picsim.plotting import fig
from picsim.provenance import RESULTS, save_result

SCRIPT = "scripts/m31_addressing.py"
# timing (SOURCE unless stated)
T_AOD_ACCESS = (4.5e-6, 7.5e-6)  # s, 1064 nm TeO2 AODs, 3 / 5 mm aperture (sintec.sg/ao/497.html, read 2026-11-14)
T_PIC_RISE = 26e-9  # s, 10-90% rise of the AlN-on-SiN switch, arXiv:2508.09920 (abstract read 2026-11-15)
T_PIC_SETTLE = 0.19e-6  # s, ZVD-shaped settling with a 5% error in the mechanical mode (m5, assumed f0 and Q)
T_DMD = 1 / 43e3  # s, update period at 43 kHz, arXiv:2310.08539 (abstract read 2026-11-15)
PULSES = (0.1e-6, 1e-6, 10e-6)  # s, gate pulse durations (CHOICE)


# ---------------------------------------------------------------------------------------------- patterns
def rect_bounds(M):
    """Rounds of a crossed AOD for 0/1 matrix M (exact partition of the ones into combinatorial rectangles).
    Returns (lower, upper): lower = real rank; upper = min(distinct nonzero rows, distinct nonzero columns, greedy)."""
    M = np.asarray(M, dtype=int)
    if M.sum() == 0:
        return 0, 0
    lower = int(np.linalg.matrix_rank(M.astype(float)))
    rows = {tuple(r) for r in M if r.any()}
    cols = {tuple(c) for c in M.T if c.any()}
    # greedy: repeatedly take the rectangle (rows sharing a column pattern subset) covering most uncovered ones
    R = M.copy()
    g = 0
    while R.any():
        best, best_cells = 0, None
        for i in np.flatnonzero(R.any(axis=1)):
            cset = R[i].astype(bool)
            rset = np.all(R[:, cset] == 1, axis=1)  # rows whose uncovered ones include all of cset
            n = int(rset.sum() * cset.sum())
            if n > best:
                best, best_cells = n, (rset, cset)
        rset, cset = best_cells
        R[np.ix_(rset, cset)] = 0
        g += 1
    return lower, min(len(rows), len(cols), g)


def patterns():
    rng = np.random.default_rng(31)
    rows = []
    for n in (4, 8, 16, 24, 40):
        for rho in (0.05, 0.1, 0.25, 0.5):
            lo, up = [], []
            for _ in range(40 if n <= 16 else 12):
                M = rng.random((n, n)) < rho
                a, b = rect_bounds(M)
                lo.append(a)
                up.append(b)
            r = {"n": n, "density": rho, "rounds_lower_mean": float(np.mean(lo)), "rounds_upper_mean": float(np.mean(up)),
                 "rounds_upper_max": int(np.max(up))}
            for tp in PULSES:
                aod = [r["rounds_lower_mean"] * (T_AOD_ACCESS[0] + tp), r["rounds_upper_mean"] * (T_AOD_ACCESS[1] + tp)]
                r[f"time_us_pulse_{tp * 1e6:g}us"] = {"AOD": [t * 1e6 for t in aod],
                                                      "PIC": (T_PIC_RISE + T_PIC_SETTLE + tp) * 1e6,
                                                      "free_space_modulator": (T_DMD + tp) * 1e6}
            rows.append(r)
            print("n %2d density %.2f | AOD rounds %.1f to %.1f (max %d) | 1 us pulse: AOD %.0f-%.0f us, PIC %.2f us, DMD %.1f us"
                  % (n, rho, r["rounds_lower_mean"], r["rounds_upper_mean"], r["rounds_upper_max"],
                     *r["time_us_pulse_1us"]["AOD"], r["time_us_pulse_1us"]["PIC"], r["time_us_pulse_1us"]["free_space_modulator"]),
                  flush=True)
    f, ax = fig(1, 1, h=3.0)
    for rho in (0.05, 0.1, 0.25, 0.5):
        rr = [r for r in rows if r["density"] == rho]
        ax[0, 0].fill_between([r["n"] for r in rr], [r["rounds_lower_mean"] for r in rr],
                              [r["rounds_upper_mean"] for r in rr], alpha=0.35, label=f"density {rho:g}")
    ax[0, 0].axhline(1, color="k", lw=0.8)
    ax[0, 0].set_xlabel("block size n (n x n sites)")
    ax[0, 0].set_ylabel("rounds per arbitrary pattern")
    ax[0, 0].legend(fontsize=6)
    save_result("m31_patterns", script=SCRIPT, fig=f,
                inputs={"patterns": "random, independent sites", "aod_access_s": T_AOD_ACCESS, "pic_rise_s": T_PIC_RISE,
                        "pic_settle_s": T_PIC_SETTLE, "dmd_period_s": T_DMD, "pulses_s": PULSES},
                outputs={"rows": rows},
                notes="Crossed AODs illuminate row-set x column-set rectangles; an arbitrary pattern lit once per site "
                      "needs an exact partition of its ones into rectangles. Bands: real rank (lower) to the best of "
                      "greedy, distinct rows and distinct columns (upper). PIC router and free-space modulator: one round. Black "
                      "line: one round. The figure shows the band between the bounds.")
    return rows


# ---------------------------------------------------------------------------------------------- leakage
def imd3_hits(sites):
    """Idle sites of a line hit by third-order products 2 f_i - f_j of tones at the lit sites (equal tone spacing per
    site), and the number of distinct products."""
    s = sorted(sites)
    prod = {2 * a - b for a in s for b in s if a != b}
    return prod - set(s)


def leakage():
    rng = np.random.default_rng(7)
    m4 = json.loads((RESULTS / "m4_switches.json").read_text())["outputs"] if (RESULTS / "m4_switches.json").exists() else {}
    line = []
    for n in (16, 40):
        for k in (2, 4, 8, 16):
            if k >= n:
                continue
            hits = []
            for _ in range(200):
                lit = rng.choice(n, size=k, replace=False)
                h = [x for x in imd3_hits(lit) if 0 <= x < n]
                hits.append(len(h))
            line.append({"sites_per_line": n, "tones": k, "idle_sites_hit_mean": float(np.mean(hits)),
                         "idle_fraction_hit": float(np.mean(hits) / (n - k))})
            print("AOD line of %d sites, %2d tones: %.1f idle sites receive 2f1-f2 products (%.0f%% of idle sites)" % (
                n, k, line[-1]["idle_sites_hit_mean"], 100 * line[-1]["idle_fraction_hit"]), flush=True)
    table = {
        "PIC router + cascaded switch (simulated)": {"idle_extinction_dB": ">= 88 (5th percentile, m4), capped at 70 "
                                                     "by the assumed stray-light floor", "source": "SIM (m4_switches)"},
        "PIC switch, measured (8 channels, 795 nm)": {"extinction_dB": 71.4, "crosstalk_on_chip_dB": -68.0,
                                                      "crosstalk_after_free_space_dB": -50.8, "rise_ns": 26,
                                                      "source": "arXiv:2508.09920 abstract, read 2026-11-15"},
        "free-space modulator, Zhang et al. (>10,000 spots; device not named in the abstract)": {"contrast_dB": 46, "crosstalk_nearest_dB": -44, "rate_kHz": 43,
                                              "source": "arXiv:2310.08539 abstract, read 2026-11-15"},
        "crossed AOD": {"intermodulation_dBc": "[VALUE?] not found in the sources read; swept -30 to -60",
                        "rounds": "m31_patterns", "access_us": "4.5 to 7.5 (1064 nm TeO2 listing)"}}
    save_result("m31_leakage", script=SCRIPT,
                inputs={"imd": "tones on an equally spaced grid (one tone per site), third-order products 2 f_i - f_j",
                        "m4": "results/m4_switches.json"},
                outputs={"imd3_line": line, "table": table, "m4_available": bool(m4)},
                notes="With equally spaced tones every third-order product lands on a site of the same line, so the "
                      "idle-site leakage of a multi-tone AOD equals its intermodulation level at those sites.")
    return line


# ---------------------------------------------------------------------------------------------- frequency
def tone_spacing_per_site(bw=20e6, sites_per_field=None):
    """RF offset between neighbouring 4 um sites: the band spread over the sites of one field (module 28: 40)."""
    if sites_per_field is None:
        sys.path.insert(0, str(RESULTS.parent / "scripts"))
        from m28_aod_baseline import AOD_FIELD
        sites_per_field = AOD_FIELD
    return bw / sites_per_field

def frequency():
    """(Corrected: an earlier version divided the band by the 150 resolvable spots and called it the offset per
    site; the 150 spots span a 40-site field at 4 um pitch (module 28), so the offset per site is band / 40.)"""
    sys.path.insert(0, str(RESULTS.parent / "scripts"))
    from m28_aod_baseline import AOD_FIELD
    bw, spots = 20e6, 150  # Hz, 25-45 MHz operating band and 150 spots (45035-5-6.5DEG-1.06, SOURCE as above)
    df = tone_spacing_per_site(bw, AOD_FIELD)  # per 4 um site: 150 spots x 1.06 um = 160 um = 40 sites
    df_spot = bw / spots
    rows = []
    for d_sites in (1, 10, 40):
        for tp in PULSES:
            rows.append({"site_separation": d_sites, "pulse_s": tp, "delta_f_Hz": d_sites * df,
                         "phase_rad": 2 * np.pi * d_sites * df * tp})
    for r in rows:
        print("sites apart %2d, pulse %5.1f us: frequency offset %.2f MHz, phase %.2f rad" % (
            r["site_separation"], r["pulse_s"] * 1e6, r["delta_f_Hz"] / 1e6, r["phase_rad"]), flush=True)
    save_result("m31_frequency", script=SCRIPT,
                inputs={"bandwidth_Hz": bw, "resolvable_spots": spots, "sites_per_field": AOD_FIELD,
                        "source": "45035-5-6.5DEG-1.06 (sintec listing); 40 sites per field from module 28"},
                outputs={"tone_spacing_per_site_Hz": df, "tone_spacing_per_resolvable_spot_Hz": df_spot, "rows": rows},
                notes="A crossed AOD shifts each addressed site by its own RF frequency. For a light-shift (Z) gate this "
                      "is harmless; for coherent two-photon (Raman) addressing against a global beam the relative phase "
                      "walks as 2 pi delta_f t and must be tracked per site. The PIC router and a DMD keep one optical "
                      "frequency for all sites.")
    return rows


if __name__ == "__main__":
    for part in sys.argv[1:] or ["patterns", "leakage", "frequency"]:
        globals()[part]()
