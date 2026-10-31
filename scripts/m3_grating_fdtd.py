"""Module 3b: 2D FDTD (MEEP) of a partially etched SiN grating antenna at 1064.4 nm.

Run with the MEEP interpreter (see environment-meep.yml):
    PYTHONPATH=. <meep-env>/bin/python scripts/m3_grating_fdtd.py [etch_um ...]

Geometry (x along propagation, y vertical, TE chip mode = Ez in 2D):
  Si substrate (PML inside) | BOX SiO2 t_box | SiN t (teeth etched by depth e from the top, duty 0.5) |
  top SiO2 t_top | air (PML)
The grating period is chosen from the 2D slab indices for a nominal air emission angle THETA_DEG
(first order, sin(theta) = n_avg - lambda/Lambda), and the actual angle is measured from the spatial DFT of
the field above the chip.

Outputs per case: emission angle in air, decay constant alpha (1/um, field), fractions of the guided input
power emitted up (into air), down (into the Si substrate), transmitted and reflected; up/(up+down) is the
directionality, independent of grating length for a uniform grating. Spectra over +/-15 nm give dtheta/dlambda
and the directionality bandwidth.

Limits: 2D (infinitely wide emitter); the lateral (x-steering axis) element factor of a ~1 um wide emitter
is taken from the channel mode in scripts/m3_opa.py. Uniform grating (exponential near field) only.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

import meep as mp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from picsim import materials as mat  # noqa: E402
from picsim import params as prm  # noqa: E402
from picsim import slab  # noqa: E402

LAM = prm.LAM_TRANSPORT.value
T = prm.T_SIN.value
T_TOP = 2.0
THETA_DEG = 10.0
RES = int(os.environ.get("MEEP_RES", 40))  # pixels per um
DPML = 1.0


def run_case(etch, t_box, L_grat=float(os.environ.get("GRAT_L", 30.0)), theta_deg=THETA_DEG, duty=0.5):
    n_sin, n_ox, n_si = float(mat.n_sin(LAM)), float(mat.n_sio2(LAM)), float(prm.N_SI_1064.value)
    # nominal period from duty-weighted slab index (first-order phase matching to air)
    n_full = slab.slab_neff(T, LAM, n_sin, n_ox)
    n_etch = slab.slab_neff(T - etch, LAM, n_sin, n_ox)
    n_avg = duty * n_full + (1 - duty) * n_etch
    period = LAM / (n_avg - np.sin(np.radians(theta_deg)))
    n_per = int(round(L_grat / period))
    L = n_per * period

    t_si = 1.0
    pad_x = 4.0
    sx = L + 2 * pad_x + 2 * DPML
    y_bot = -t_box - t_si - DPML
    y_top = T + T_TOP + 1.5 + DPML
    sy = y_top - y_bot
    yc = 0.5 * (y_top + y_bot)
    Y = lambda y: y - yc  # meep cell is centred at 0
    x0 = -L / 2

    SiN, Ox, Si = mp.Medium(index=n_sin), mp.Medium(index=n_ox), mp.Medium(index=n_si)
    geom = [
        mp.Block(size=mp.Vector3(mp.inf, t_si + DPML), center=mp.Vector3(0, Y(-t_box - (t_si + DPML) / 2)), material=Si),
        mp.Block(size=mp.Vector3(mp.inf, t_box + T + T_TOP), center=mp.Vector3(0, Y((-t_box + T + T_TOP) / 2)),
                 material=Ox),
        mp.Block(size=mp.Vector3(mp.inf, T), center=mp.Vector3(0, Y(T / 2)), material=SiN),
    ]
    # etched trenches (filled by top oxide), (1 - duty) of each period
    w_tr = (1 - duty) * period
    for i in range(n_per):
        xc = x0 + i * period + duty * period + w_tr / 2
        geom.append(mp.Block(size=mp.Vector3(w_tr, etch), center=mp.Vector3(xc, Y(T - etch / 2)), material=Ox))

    fcen = 1 / LAM
    df = fcen * 0.03
    src_x = -sx / 2 + DPML + 1.0
    sources = [mp.EigenModeSource(mp.GaussianSource(fcen, fwidth=df), center=mp.Vector3(src_x, Y(T / 2)),
                                  size=mp.Vector3(0, 3.0), eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True)]
    sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), geometry=geom, sources=sources, resolution=RES,
                        boundary_layers=[mp.PML(DPML)], default_material=mp.Medium(index=1.0))
    nf = 31
    yw = Y(T / 2)
    mon = {
        "in": sim.add_flux(fcen, df, nf, mp.FluxRegion(center=mp.Vector3(src_x + 0.6, yw), size=mp.Vector3(0, 3.0))),
        "trans": sim.add_flux(fcen, df, nf, mp.FluxRegion(center=mp.Vector3(sx / 2 - DPML - 0.5, yw),
                                                          size=mp.Vector3(0, 3.0))),
        "up": sim.add_flux(fcen, df, nf, mp.FluxRegion(center=mp.Vector3(0, Y(T + T_TOP + 0.75)),
                                                       size=mp.Vector3(sx - 2 * DPML, 0))),
        "down": sim.add_flux(fcen, df, nf, mp.FluxRegion(center=mp.Vector3(0, Y(-t_box - 0.3)),
                                                         size=mp.Vector3(sx - 2 * DPML, 0), weight=-1)),
    }
    dft = sim.add_dft_fields([mp.Ez], fcen, 0, 1, center=mp.Vector3(0, Y(T + T_TOP + 0.75)),
                             size=mp.Vector3(sx - 2 * DPML, 0))
    # normalisation run would be cleaner; the input monitor sits before the grating and records
    # incident minus reflected flux, so we also record reflection via the straight-guide reference below
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(sx / 2 - DPML - 0.5, yw), 1e-6))
    freqs = np.array(mp.get_flux_freqs(mon["in"]))
    F = {k: np.array(mp.get_fluxes(v)) for k, v in mon.items()}
    ez = sim.get_dft_array(dft, mp.Ez, 0)
    sim.reset_meep()

    # reference: straight waveguide, same source, gives incident power
    sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), geometry=geom[:3], sources=sources, resolution=RES,
                        boundary_layers=[mp.PML(DPML)], default_material=mp.Medium(index=1.0))
    ref = sim.add_flux(fcen, df, nf, mp.FluxRegion(center=mp.Vector3(src_x + 0.6, yw), size=mp.Vector3(0, 3.0)))
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(sx / 2 - DPML - 0.5, yw), 1e-6))
    P0 = np.array(mp.get_fluxes(ref))
    sim.reset_meep()

    up, down, tr = F["up"] / P0, F["down"] / P0, F["trans"] / P0
    refl = 1 - F["in"] / P0
    # emission angle from spatial spectrum of the near field above the chip
    xs = np.linspace(-(sx - 2 * DPML) / 2, (sx - 2 * DPML) / 2, len(ez))
    win = np.abs(xs) < L / 2
    nfft = 1 << 16
    spec = np.abs(np.fft.fftshift(np.fft.fft(ez * win, nfft))) ** 2
    kx = np.fft.fftshift(np.fft.fftfreq(nfft, d=xs[1] - xs[0])) * 2 * np.pi
    s = kx / (2 * np.pi / LAM)
    prop = np.abs(s) < 1
    s_peak = float(s[prop][np.argmax(spec[prop])])
    # field decay constant: T = exp(-2 alpha L)
    ic = nf // 2
    alpha = float(-np.log(max(tr[ic], 1e-12)) / (2 * L))
    lam_arr = 1 / freqs
    # angle vs wavelength from phase matching with measured neff-shift neglected: use measured s at centre and
    # the grating equation slope d s / d lambda = -(1/Lambda) + d n_avg/d lambda (~ -1/Lambda - ng-ish term)
    return {
        "etch_um": etch, "t_box_um": t_box, "period_um": period, "n_periods": n_per, "L_um": L,
        "n_avg_slab": n_avg, "theta_air_deg": float(np.degrees(np.arcsin(s_peak))),
        "alpha_per_um": alpha, "L_1e_power_um": 1 / (2 * alpha) if alpha > 0 else float("inf"),
        "lam_um": lam_arr.tolist(), "up": up.tolist(), "down": down.tolist(), "trans": tr.tolist(),
        "refl": refl.tolist(),
        "up_c": float(up[ic]), "down_c": float(down[ic]), "trans_c": float(tr[ic]), "refl_c": float(refl[ic]),
        "directionality_c": float(up[ic] / (up[ic] + down[ic])),
    }


if __name__ == "__main__":
    etches = [float(a) for a in sys.argv[1:]] or [0.05, 0.1, 0.2]
    boxes = [2.0, 3.0, 4.0]
    out = []
    for e in etches:
        for tb in boxes:
            r = run_case(e, tb)
            print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                              if not isinstance(v, list)}), flush=True)
            out.append(r)
    dest = REPO / "results" / "m3_grating_fdtd_raw.json"
    prev = json.loads(dest.read_text()) if dest.exists() else []
    keep = [p for p in prev if (p["etch_um"], p["t_box_um"]) not in {(o["etch_um"], o["t_box_um"]) for o in out}]
    dest.write_text(json.dumps(keep + out, indent=1))
