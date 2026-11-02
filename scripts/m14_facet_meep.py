"""Module 14b: 2D FDTD (MEEP) of the emitter facet in its vertical cross-section, into air or into glass.

Plane (z along the guide -> MEEP x; y vertical -> MEEP y); chip TE field (along the film, out of this plane) ->
Ez (ODD_Z). The 400 nm SiN core is replaced by a slab whose TE index matches the 3D femwell mode
(n_eff = 1.7995, lateral confinement folded into the core index), in SiO2, with the Si substrate below a 4 um
BOX. The guide ends at x = 0 in a medium of index n_m (1: air; 1.4496: index-matched glass block).

Outputs (results/m14_facet_meep_<medium>.json): reflected and transmitted power relative to the incident mode,
the angular power spectrum P(s) of the transmitted field (s = n_m sin theta) from the field on a line just
beyond the facet, the fraction within |s| < 0.9, and the Kirchhoff estimate for the same slab mode.
Run:  PYTHONPATH=. <meep-env>/bin/python scripts/m14_facet_meep.py air|glass [resolution]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import meep as mp
import numpy as np
from scipy.optimize import brentq

LAM = 1.064435
N_OX = 1.4496257710158857
N_SI = 3.55
N_EFF_3D = 1.7994707727722363  # femwell, results/m13_facet.json
T = 0.4


def slab_te_neff(nc, ncl, t, lam):
    k0 = 2 * np.pi / lam
    f = lambda ne: np.tan(k0 * np.sqrt(nc**2 - ne**2) * t / 2) - np.sqrt(ne**2 - ncl**2) / np.sqrt(nc**2 - ne**2)
    hi = min(nc - 1e-9, ncl + (nc - ncl) * 0.999)
    grid = np.linspace(ncl + 1e-6, hi, 4000)
    v = f(grid)
    for i in range(len(grid) - 1, 0, -1):
        if np.isfinite(v[i]) and np.isfinite(v[i - 1]) and v[i] * v[i - 1] < 0 and abs(v[i]) < 50:
            return brentq(f, grid[i - 1], grid[i])
    raise ValueError


def slab_mode(nc, ncl, t, lam, ne, y):
    k0 = 2 * np.pi / lam
    kap, gam = k0 * np.sqrt(nc**2 - ne**2), k0 * np.sqrt(ne**2 - ncl**2)
    inside = np.abs(y) <= t / 2
    return np.where(inside, np.cos(kap * y), np.cos(kap * t / 2) * np.exp(-gam * (np.abs(y) - t / 2)))


def main(medium="air", res=40):
    n_m = 1.0 if medium == "air" else N_OX
    nc = brentq(lambda n: slab_te_neff(n, N_OX, T, LAM) - N_EFF_3D, N_EFF_3D + 1e-3, 2.2)
    ne = slab_te_neff(nc, N_OX, T, LAM)
    dpml = 1.5
    sx, sy = 16.0 + 2 * dpml, 26.0 + 2 * dpml
    cell = mp.Vector3(sx, sy)
    xL, xR = -sx / 2, sx / 2
    x_facet = xL + dpml + 7.0
    geom = [
        mp.Block(mp.Vector3(x_facet - xL, sy), center=mp.Vector3((xL + x_facet) / 2, 0), material=mp.Medium(index=N_OX)),
        mp.Block(mp.Vector3(x_facet - xL, sy / 2 - 4.0 - T / 2), center=mp.Vector3((xL + x_facet) / 2, (-sy / 2 + (-4.0 - T / 2)) / 2),
                 material=mp.Medium(index=N_SI)),
        mp.Block(mp.Vector3(x_facet - xL, T), center=mp.Vector3((xL + x_facet) / 2, 0), material=mp.Medium(index=nc)),
    ]
    fcen = 1 / LAM
    src_x = xL + dpml + 1.5
    sources = [mp.EigenModeSource(mp.GaussianSource(fcen, fwidth=0.1 * fcen), center=mp.Vector3(src_x, 0),
                                  size=mp.Vector3(0, 6.0), eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True)]
    sim = mp.Simulation(cell_size=cell, geometry=geom, sources=sources, boundary_layers=[mp.PML(dpml)],
                        resolution=res, default_material=mp.Medium(index=n_m), dimensions=2)
    mon_x = src_x + 2.0
    mode_mon = sim.add_mode_monitor(fcen, 0, 1, mp.ModeRegion(center=mp.Vector3(mon_x, 0), size=mp.Vector3(0, 6.0)))
    line_x = x_facet + 0.5
    ny = int(round((sy - 2 * dpml) * res))
    dft = sim.add_dft_fields([mp.Ez, mp.Hy], fcen, 0, 1, center=mp.Vector3(line_x, 0),
                             size=mp.Vector3(0, sy - 2 * dpml))
    tflux = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(line_x, 0), size=mp.Vector3(0, sy - 2 * dpml)))
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(line_x + 2, 0), 1e-7))
    coef = sim.get_eigenmode_coefficients(mode_mon, [1], eig_parity=mp.ODD_Z).alpha[0, 0]
    P_in, P_back = abs(coef[0]) ** 2, abs(coef[1]) ** 2
    P_t = mp.get_fluxes(tflux)[0]
    Ez = np.array(sim.get_dft_array(dft, mp.Ez, 0))
    y = np.linspace(-(sy - 2 * dpml) / 2, (sy - 2 * dpml) / 2, len(Ez))
    # angular spectrum of the transmitted field (homogeneous medium n_m beyond the facet)
    nfft = 2**15
    dy = y[1] - y[0]
    pad = np.zeros(nfft, dtype=complex)
    pad[:len(Ez)] = Ez
    s = np.fft.fftshift(np.fft.fftfreq(nfft, d=dy)) * LAM
    A = np.fft.fftshift(np.fft.fft(pad))
    prop = np.abs(s) < n_m
    dens = np.where(prop, np.abs(A) ** 2 * np.sqrt(np.clip(1 - (s / n_m) ** 2, 0, None)), 0)
    frac09 = float(dens[np.abs(s) < 0.9].sum() / dens.sum())
    s1 = s[prop & (s >= 0)]
    d1 = dens[prop & (s >= 0)] / dens.max()
    na_1e2 = float(np.interp(np.exp(-2), d1[::-1], s1[::-1]))
    # Kirchhoff estimate for the same slab mode
    yk = np.arange(-6, 6, 0.002)
    Yk = slab_mode(nc, N_OX, T, LAM, ne, yk)
    padk = np.zeros(nfft)
    padk[:len(Yk)] = Yk
    sk = np.fft.fftshift(np.fft.fftfreq(nfft, d=0.002)) * LAM
    Ak = np.abs(np.fft.fftshift(np.fft.fft(padk))) ** 2
    pk = np.abs(sk) < n_m
    frac_prop = float((Ak[pk] * np.sqrt(1 - (sk[pk] / n_m) ** 2)).sum() / Ak.sum())
    T_kirch = 4 * ne * n_m / (ne + n_m) ** 2 * frac_prop
    densk = np.where(pk, Ak * np.sqrt(np.clip(1 - (sk / n_m) ** 2, 0, None)), 0)
    frac09_k = float(densk[np.abs(sk) < 0.9].sum() / densk.sum())
    out = {"medium": medium, "n_m": n_m, "resolution": res, "core_index_2d": nc, "neff_2d": ne,
           "R": P_back / P_in, "T": P_t / P_in, "T_plus_R": (P_t + P_back) / P_in,
           "fraction_within_s0.9": frac09, "T_within_s0.9": P_t / P_in * frac09, "NA_1e2": na_1e2,
           "kirchhoff": {"T": T_kirch, "fraction_within_s0.9": frac09_k, "T_within_s0.9": T_kirch * frac09_k},
           "spectrum_s": s1[::8].tolist(), "spectrum_P": d1[::8].tolist()}
    print(json.dumps({k: v for k, v in out.items() if not k.startswith("spectrum")}, indent=1), flush=True)
    Path("results").joinpath(f"m14_facet_meep_{medium}.json").write_text(json.dumps(
        {"script": "scripts/m14_facet_meep.py", "inputs": {"medium": medium, "resolution": res, "slab": "SiN core "
         "index fitted to n_eff 1.7995, t 0.4 um, SiO2, Si below 4 um BOX", "polarisation": "Ez (chip TE)"},
         "outputs": out}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "air", int(sys.argv[2]) if len(sys.argv) > 2 else 40)
