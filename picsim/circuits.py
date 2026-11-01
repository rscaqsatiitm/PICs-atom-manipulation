"""Circuit models of MZI switches, cascades and routing trees (2x2 transfer matrices, numpy).

Conventions: a 2x2 coupler with power cross-coupling kappa and excess transmission tau has
    C(kappa) = sqrt(tau) [[sqrt(1-kappa), i sqrt(kappa)], [i sqrt(kappa), sqrt(1-kappa)]].
An MZI is C(k2) . diag(a1 e^{i phi}, a2) . C(k1) with arm amplitude transmissions a1, a2. Input in port 0,
output taken from the cross port (index 1): at phi = 0 and kappa = 1/2 the cross port is fully ON.

A general (unbalanced) MZI cannot null its output; the minimum transmission over phi is
    T_min = (|t1 t2| a1 - |k1 k2| a2)^2 ... in amplitude form, computed numerically below.
"""
from __future__ import annotations

import numpy as np


def coupler(kappa, tau=1.0):
    k = np.sqrt(kappa)
    t = np.sqrt(1 - kappa)
    return np.sqrt(tau) * np.array([[t, 1j * k], [1j * k, t]])


def mzi(phi, k1=0.5, k2=0.5, a1=1.0, a2=1.0, tau=1.0):
    return coupler(k2, tau) @ np.diag([a1 * np.exp(1j * phi), a2]) @ coupler(k1, tau)


def mzi_cross_T(phi, k1=0.5, k2=0.5, a1=1.0, a2=1.0):
    """Cross-port power transmission vs phase (vectorised over phi)."""
    phi = np.asarray(phi, float)
    # cross amplitude: [C2 D C1]_{10} = i(sqrt(k2 (1-k1)) a1 e^{i phi} + sqrt((1-k2) k1) a2) ... (derived)
    A = np.sqrt(k2) * np.sqrt(1 - k1) * a1 * np.exp(1j * phi) + np.sqrt(1 - k2) * np.sqrt(k1) * a2
    return np.abs(A) ** 2


def mzi_on_off(k1, k2, a1=1.0, a2=1.0, dphi_off=0.0):
    """ON (max) and OFF (at the ideal null phase pi + dphi_off) cross transmissions."""
    # max at phi = 0, min at phi = pi for real positive coefficients
    T_on = mzi_cross_T(0.0, k1, k2, a1, a2)
    T_off = mzi_cross_T(np.pi + dphi_off, k1, k2, a1, a2)
    return float(T_on), float(T_off)


def extinction_db(T_on, T_off, floor=0.0):
    return 10 * np.log10(T_on / (T_off + floor))


def mc_extinction(n, sigma_k, sigma_phi, sigma_a=0.0, stages=1, tunable=False, floor=0.0, seed=0,
                  k_mean=0.5):
    """Monte Carlo extinction (dB) of a switch channel.

    stages: number of cascaded MZIs, all switched together (ON: all at phi = 0; OFF: all at null).
    tunable: if True, each MZI's input splitter is itself a thermo-optically trimmed MZI, set to the exact
             null condition (compensating coupler and arm-loss imbalance); the residual is set by the trim phase
             error sigma_phi (see trimmed_kappa).
    floor: incoherent stray-light floor relative to the input (e.g. 1e-7 for -70 dB).
    """
    rng = np.random.default_rng(seed)
    out = np.empty(n)
    for i in range(n):
        T_on, T_off = 1.0, 1.0
        for _ in range(stages):
            k1 = np.clip(k_mean + sigma_k * rng.standard_normal(), 1e-6, 1 - 1e-6)
            k2 = np.clip(k_mean + sigma_k * rng.standard_normal(), 1e-6, 1 - 1e-6)
            a1 = 1 + sigma_a * rng.standard_normal()
            if tunable:
                # trim the input splitter to the exact null condition sqrt(k2 (1-k1)) a1 = sqrt((1-k2) k1),
                # which compensates both coupler and arm-loss imbalance; residual set by the trim phase error
                k_null = k2 * a1**2 / (k2 * a1**2 + (1 - k2))
                k1 = trimmed_kappa(k_null, sigma_phi * rng.standard_normal())
            on, off = mzi_on_off(k1, k2, a1, 1.0, sigma_phi * rng.standard_normal())
            T_on *= on
            T_off *= off
        out[i] = extinction_db(T_on, T_off, floor)
    return out


def trimmed_kappa(k_target, dphi):
    """Splitting ratio of an MZI-based tunable splitter set to k_target with a phase error dphi.

    A balanced MZI splitter has kappa = sin^2(theta/2); an error dphi in theta shifts kappa by
    sin(theta) dphi / 2 to first order.
    """
    theta = 2 * np.arcsin(np.sqrt(k_target))
    return float(np.sin((theta + dphi) / 2) ** 2)


# ------------------------------------------------------------------ actuator dynamics
def second_order_step(t, f0, Q, t_steps=(0.0,), amps=(1.0,)):
    """Response of a unity-gain second-order system (resonance f0, quality factor Q) to a sum of steps."""
    w0 = 2 * np.pi * f0
    z = 1 / (2 * Q)
    wd = w0 * np.sqrt(max(1 - z**2, 1e-12))
    y = np.zeros_like(t)
    for ts, A in zip(t_steps, amps):
        tt = np.clip(t - ts, 0, None)
        resp = 1 - np.exp(-z * w0 * tt) * (np.cos(wd * tt) + z / np.sqrt(1 - z**2) * np.sin(wd * tt))
        y += A * np.where(t >= ts, resp, 0.0)
    return y


def zv_shaper(f0, Q):
    """Zero-vibration input shaper (two impulses) for a second-order mode: times and amplitudes."""
    z = 1 / (2 * Q)
    K = np.exp(-z * np.pi / np.sqrt(1 - z**2))
    wd = 2 * np.pi * f0 * np.sqrt(1 - z**2)
    return (0.0, np.pi / wd), (1 / (1 + K), K / (1 + K))


def zvd_shaper(f0, Q):
    """Zero-vibration-derivative shaper (three impulses): robust to resonance-frequency error."""
    z = 1 / (2 * Q)
    K = np.exp(-z * np.pi / np.sqrt(1 - z**2))
    wd = 2 * np.pi * f0 * np.sqrt(1 - z**2)
    s = 1 + 2 * K + K**2
    return (0.0, np.pi / wd, 2 * np.pi / wd), (1 / s, 2 * K / s, K**2 / s)


def switching_metrics(f0, Q, shaped=False, f0_true=None, ER_target_db=50.0, t_max=5e-6, n=200001):
    """Optical 10-90% rise (OFF->ON) and time to reach ER_target after ON->OFF, for a balanced MZI.

    The phase follows the actuator step response; OFF->ON: phi from pi to 0; ON->OFF: phi from 0 to pi.
    Cross-port intensity T = cos^2(phi/2).
    """
    f0_true = f0_true or f0
    t = np.linspace(0, t_max, n)
    if shaped == "zvd":
        ts, A = zvd_shaper(f0, Q)
    elif shaped:
        ts, A = zv_shaper(f0, Q)
    else:
        ts, A = (0.0,), (1.0,)
    y = second_order_step(t, f0_true, Q, ts, A)  # 0 -> 1
    # OFF -> ON: phi = pi (1 - y)
    T_on = np.cos(np.pi * (1 - y) / 2) ** 2
    i10 = np.argmax(T_on >= 0.1)
    i90 = np.argmax(T_on >= 0.9)
    rise = t[i90] - t[i10]
    # ON -> OFF: phi = pi y; extinction relative to ON
    T_off = np.cos(np.pi * y / 2) ** 2
    lvl = 10 ** (-ER_target_db / 10)
    above = np.where(T_off > lvl)[0]
    t_settle = t[above[-1]] if len(above) and above[-1] < n - 1 else (np.nan if len(above) else 0.0)
    return {"rise_10_90_s": float(rise), "t_to_ER_s": float(t_settle)}


# ------------------------------------------------------------------ routing tree
def router_leakage(n_stages, sigma_phi, n=2000, seed=0, target=0):
    """Binary tree of tunable MZI splitters routing all light to output 'target'.

    Each node is a tunable splitter set to send everything down the target branch; a phase error sigma_phi
    leaks a fraction ~ (sigma_phi/2)^2. Returns per-output power fractions (n, 2**n_stages).
    """
    rng = np.random.default_rng(seed)
    N = 2**n_stages
    P = np.zeros((n, N))
    bits = [(target >> (n_stages - 1 - s)) & 1 for s in range(n_stages)]
    for i in range(n):
        p = np.ones(1)
        for s in range(n_stages):
            # each node at this stage: fraction to child 1 is kappa
            kap = np.empty(2**s)
            want = 1.0 - 1e-12 if bits[s] == 1 else 1e-12
            for j in range(2**s):
                kap[j] = trimmed_kappa(want, sigma_phi * rng.standard_normal())
            new = np.empty(2 ** (s + 1))
            new[0::2] = p * (1 - kap)
            new[1::2] = p * kap
            p = new
        P[i] = p
    return P
