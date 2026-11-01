"""Port-to-arm transfer matrix of a star coupler by wide-angle FD-BPM (picsim.bpm) on an FPRLayout.

Each port is launched in the fundamental lateral mode of its narrow section, tilted along the port axis; all
ports are propagated together (one right-hand side each). At the measurement plane the field is projected onto
the tilted fundamental mode of each narrow arm. T[k, n] is the amplitude in arm n for unit power in port k.
"""
from __future__ import annotations

import numpy as np

from ..bpm import BPM2D, coverage, lateral_profile, tm_overlap, tm_power
from ..slab import fd_modes
from .fpr_layout import FPRLayout


def narrow_mode(width, lam, n_core, n_clad, half=5.0, dx=0.005):
    xi = np.arange(-half - width / 2, half + width / 2, dx)
    ne, v = fd_modes(xi, lateral_profile(xi, [(0.0, width)], n_core, n_clad), lam, pol="TM", num=1)
    m = v[:, 0] * np.sign(v[len(xi) // 2, 0])
    return xi, m, float(ne[0])


def tilted(x, xc, alpha, xi, m, beta):
    """Fundamental mode of a guide centred at xc whose axis makes angle alpha with +z, sampled on x."""
    prof = np.interp((x - xc) * np.cos(alpha), xi, m, left=0.0, right=0.0)
    return prof * np.exp(-1j * beta * np.sin(alpha) * (x - xc))


def transfer_matrix(lay: FPRLayout, lam, n_core, n_clad, dx=0.04, dz=0.25, ports=None, margin=40.0,
                    absorber=30.0, k_abs=0.05, progress=False):
    ports = np.arange(lay.N_ports) if ports is None else np.asarray(ports)
    X = max(lay.x_slab, np.abs(lay.Ax).max() + lay.L_arm * np.abs(np.sin(lay.phi)).max(),
            np.abs(lay.Px).max()) + margin
    x = np.arange(-X, X, dx)
    xi, m, ne = narrow_mode(lay.w_narrow, lam, n_core, n_clad)
    beta = 2 * np.pi / lam * ne
    bpm = BPM2D(x, lam, n_core, absorber=absorber, k_abs=k_abs)
    z = lay.z_launch
    H = np.empty((len(x), len(ports)), dtype=complex)
    for j, k in enumerate(ports):
        xc, al = lay.port_axis(k, z)
        H[:, j] = tilted(x, xc, al, xi, m, beta)

    def nprof(zz):
        a, b = lay.intervals(zz)
        f = coverage(x, a, b) if len(a) else np.zeros_like(x)
        return np.sqrt(1.0 / (f / n_core**2 + (1 - f) / n_clad**2))

    n0 = nprof(z)
    p_in = tm_power(H, n0)
    z_end = lay.z_measure
    nsteps = int(np.ceil((z_end - z) / dz))
    for i in range(nsteps):
        zm = z + (i + 0.5) * dz
        H = bpm.step(H, nprof(zm), dz)
        if progress and i % 500 == 0:
            print(f"  BPM step {i}/{nsteps}", flush=True)
    z = z + nsteps * dz
    n1 = nprof(z)
    T = np.empty((len(ports), lay.N_arms), dtype=complex)
    for nn in range(lay.N_arms):
        xc, al = lay.arm_axis(nn, z)
        T[:, nn] = tm_overlap(H, tilted(x, xc, al, xi, m, beta), n1) / np.sqrt(p_in)
    return T, {"x": x, "p_out_window": tm_power(H, n1) / p_in, "z_end": z, "neff_narrow": ne}
