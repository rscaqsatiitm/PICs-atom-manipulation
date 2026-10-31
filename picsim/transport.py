"""Motional excitation of an atom in a harmonic tweezer whose centre follows a prescribed trajectory.

For a harmonic trap of angular frequency omega whose centre moves along x0(t), an atom starting in the ground
state ends in a coherent state with mean phonon number (exact for a harmonic potential)

    n_bar = (m / (2 hbar omega)) |int xdd_0(t) exp(-i omega t) dt|^2 ,

where xdd_0 is the acceleration of the trap centre. Anharmonicity and depth changes during the move are
neglected; this is valid while the induced oscillation amplitude stays well below the waist.

Trajectories: a minimum-jerk profile x(t) = d (10 u^3 - 15 u^4 + 6 u^5), u = t/T, optionally sampled and held at
the OPA update rate f_u (staircase), then low-pass filtered by the phase-shifter response (first order, tau).
"""
from __future__ import annotations

import numpy as np
from scipy.constants import hbar, k as kB

from . import params as prm

AU_POL = 1.64877727436e-41  # C m^2 / V, atomic unit of polarisability (CODATA)
EPS0 = 8.8541878128e-12
C = 299792458.0


def trap_depth_K(P, w0_um, alpha_au):
    """Peak depth U0/kB (K) of a Gaussian tweezer of power P (W) and waist w0 (um)."""
    I0 = 2 * P / (np.pi * (w0_um * 1e-6) ** 2)
    return alpha_au * AU_POL * I0 / (2 * EPS0 * C) / kB


def power_for_depth(U_K, w0_um, alpha_au):
    return U_K / trap_depth_K(1.0, w0_um, alpha_au)


def radial_freq_Hz(U_K, w0_um, mass=None):
    m = mass or prm.CS_MASS.value
    return np.sqrt(4 * U_K * kB / (m * (w0_um * 1e-6) ** 2)) / (2 * np.pi)


def min_jerk(t, T, d):
    u = np.clip(t / T, 0, 1)
    return d * (10 * u**3 - 15 * u**4 + 6 * u**5)


def trajectory(T, d, f_update=None, tau=0.0, dt=None, t_settle=None):
    """Trap-centre trajectory (m) on a uniform time grid (s): min-jerk, optional sample-and-hold, low pass."""
    dt = dt or min(T / 4000, 1e-8)
    t_settle = t_settle if t_settle is not None else max(5 * tau, 0.2 * T)
    t = np.arange(0, T + t_settle, dt)
    if f_update:
        t_hold = np.floor(t * f_update) / f_update
        x = min_jerk(t_hold, T, d)
    else:
        x = min_jerk(t, T, d)
    if tau > 0:
        a = dt / (tau + dt)
        y = np.empty_like(x)
        y[0] = x[0]
        for i in range(1, len(x)):
            y[i] = y[i - 1] + a * (x[i] - y[i - 1])
        x = y
    return t, x


def n_bar(t, x, f_trap, mass=None):
    """Mean phonon number added by the trajectory, from the Fourier component of the velocity.

    int xdd exp(-i w t) dt = [xd exp(-i w t)] + i w int xd exp(-i w t) dt; xd vanishes at both ends.
    """
    m = mass or prm.CS_MASS.value
    w = 2 * np.pi * f_trap
    v = np.gradient(x, t)
    A = 1j * w * np.trapezoid(v * np.exp(-1j * w * t), t)
    return m * np.abs(A) ** 2 / (2 * hbar * w)


def n_bar_staircase(T, d, f_update, f_trap, tau=0.0, mass=None):
    """Mean phonon number for a sample-and-hold trajectory, using the exact step decomposition.

    A staircase of steps dx_j at times t_j, each filtered by a first-order response (tau), has velocity
    v(t) = sum_j dx_j exp(-(t - t_j)/tau)/tau, hence its Fourier component at omega is
    sum_j dx_j exp(-i w t_j) / (1 + i w tau). This avoids resolving the steps on a time grid.
    """
    m = mass or prm.CS_MASS.value
    w = 2 * np.pi * f_trap
    tj = np.arange(0, T + 1 / f_update, 1 / f_update)
    xj = min_jerk(tj, T, d)
    dx = np.diff(np.concatenate([[0.0], xj]))
    F = np.sum(dx * np.exp(-1j * w * tj)) / (1 + 1j * w * tau)
    A = 1j * w * F
    return m * np.abs(A) ** 2 / (2 * hbar * w)


def n_bar_smooth(T, d, f_trap, mass=None, n=20001):
    t = np.linspace(0, T, n)
    return n_bar(t, min_jerk(t, T, d), f_trap, mass)


def heating_rate_position_noise(S_x, f_trap, mass=None):
    """Mean energy heating rate (quanta/s) from trap-position noise with one-sided PSD S_x (m^2/Hz) at f_trap.

    Savard, O'Hara and Thomas, PRA 56, R1095 (1997) (RECALLED): <E_dot> = (pi/2) m w^4 S_x(w) (w = 2 pi f_trap),
    hence quanta/s = (pi/2) m w^3 S_x / hbar.
    """
    m = mass or prm.CS_MASS.value
    w = 2 * np.pi * f_trap
    return np.pi / 2 * m * w**3 * S_x / hbar


def heating_rate_intensity_noise(S_eps, f_trap):
    """Parametric heating rate constant Gamma (1/s) from fractional intensity noise PSD S_eps (1/Hz) at 2 f_trap.

    Savard et al. 1997 (RECALLED): Gamma = pi^2 f_trap^2 S_eps(2 f_trap); <E>(t) = <E>(0) exp(Gamma t).
    """
    return np.pi**2 * f_trap**2 * S_eps
