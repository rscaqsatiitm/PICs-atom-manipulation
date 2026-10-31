"""Phase-only synthesis of K simultaneous tweezers with an N-channel OPA (weighted Gerchberg-Saxton).

The OPA has N complex amplitudes with fixed moduli (set by the splitter tree) and free phases. The field at
target m (sin(theta) = s_m) is E_m = sum_n a_n exp(i phi_n) H_nm with H_nm = exp(-i k x_n s_m) E_el(s_m).
The weighted GS iteration (Di Leonardo et al., Opt. Express 15, 1913 (2007), RECALLED) equalises |E_m|.
"""
from __future__ import annotations

import numpy as np

from .array import OPA, channel_response


def element_factor_at(opa: OPA, s):
    return np.interp(s, opa.s, np.abs(opa._E_el))


def gsw(opa: OPA, x_targets, iters=200, seed=0, target_weights=None, H=None):
    """Return channel phases producing spots at atom-plane positions x_targets with equal (or weighted)
    intensities. H defaults to the exact channel response through slit and pupil."""
    rng = np.random.default_rng(seed)
    x_t = np.asarray(x_targets, float)
    H = channel_response(opa, x_t) if H is None else H
    a = opa.apod
    tw = np.ones(len(x_t)) if target_weights is None else np.asarray(target_weights, float)
    w = np.ones(len(x_t))
    phi = rng.uniform(0, 2 * np.pi, opa.N)
    for _ in range(iters):
        E = (a * np.exp(1j * phi)) @ H
        amp = np.abs(E) / np.sqrt(tw)
        w = w * np.mean(amp) / np.maximum(amp, 1e-30)
        back = H.conj() @ (w * np.sqrt(tw) * E / np.maximum(np.abs(E), 1e-30))
        phi = np.angle(back)
    return phi


def evaluate(opa: OPA, phases, x_targets, site_pitch, x_eval=None, r_excl=None):
    """Metrics of a multi-tweezer pattern at the atom plane (steering axis).

    Returns a dict:
      eff_spots      fraction of emitted power within +/- r_excl of the targets
      transmission   fraction of emitted power through slit and pupil
      uniformity     (max - min) / (max + min) of target peak intensities
      stray_any      brightest intensity farther than r_excl from all targets / mean target peak
      stray_sites    brightest intensity at untargeted lattice sites (multiples of site_pitch) / mean target peak
    """
    x_targets = np.asarray(x_targets, float)
    r_excl = r_excl if r_excl is not None else 0.5 * site_pitch
    if x_eval is None:
        x_eval = np.arange(-1.1 * opa.x_range, 1.1 * opa.x_range, opa.spot_pitch / 20)
    I = opa.atom_intensity(phases, x_eval)
    dx = x_eval[1] - x_eval[0]
    near = np.min(np.abs(x_eval[:, None] - x_targets[None, :]), axis=1) <= r_excl
    pk = []
    for xt in x_targets:
        win = np.abs(x_eval - xt) <= 0.5 * opa.spot_pitch
        pk.append(I[win].max())
    pk = np.array(pk)
    sites = np.arange(np.ceil(-opa.x_range / site_pitch), np.floor(opa.x_range / site_pitch) + 1) * site_pitch
    untargeted = sites[np.min(np.abs(sites[:, None] - x_targets[None, :]), axis=1) > 0.5 * site_pitch]
    I_sites = opa.atom_intensity(phases, untargeted) if len(untargeted) else np.array([0.0])
    return {
        "eff_spots": float(np.sum(I[near]) * dx),
        "transmission": float(opa.transmission(phases)),
        "uniformity": float((pk.max() - pk.min()) / (pk.max() + pk.min())),
        "peak_mean": float(pk.mean()),
        "stray_any": float(I[~near].max() / pk.mean()) if np.any(~near) else 0.0,
        "stray_sites": float(I_sites.max() / pk.mean()),
    }


def refine(opa: OPA, phases0, x_targets, x_dark, beta=10.0, gamma=20.0, eta=20.0, q=2, iters=600):
    """Gradient refinement of channel phases on the exact relay response.

    Cost on sampled intensities (targets t, target +/- delta, dark points d):
        C = -log mu + gamma (<I_t^2>/mu^2 - 1) + beta sum_d (I_d/mu)^q + eta sum_t ((I_t+ - I_t-)/I_t)^2,
    mu = <I_t>. The terms reward brightness, equal target intensities, dark surroundings and spots centred on
    their targets (delta = 1/4 resolution cell). Gradient by the chain rule,
        dC/dphi_n = 2 Re(i c_n sum_p H_np conj(E_p) dC/dI_p).
    """
    from scipy.optimize import minimize

    xt = np.asarray(x_targets, float)
    xd = np.asarray(x_dark, float)
    K, D = len(xt), len(xd)
    dlt = 0.25 * opa.spot_pitch
    H = channel_response(opa, np.concatenate([xt, xt + dlt, xt - dlt, xd]))
    a = opa.apod

    def cost(phi):
        c = a * np.exp(1j * phi)
        E = c @ H
        I = np.abs(E) ** 2
        It, Ip, Im, Id = I[:K], I[K:2 * K], I[2 * K:3 * K], I[3 * K:]
        mu = It.mean()
        m2 = np.mean(It**2)
        r = Id / mu
        Dm = (Ip - Im) / It
        C = -np.log(mu) + gamma * (m2 / mu**2 - 1) + beta * np.sum(r**q) + eta * np.sum(Dm**2)
        w = np.zeros_like(I)
        w[:K] = (-1 / (K * mu) + gamma * (2 * It / (K * mu**2) - 2 * m2 / (K * mu**3))
                 - beta * q * np.sum(r**q) / (mu * K) - 2 * eta * Dm**2 / It)
        w[K:2 * K] = 2 * eta * Dm / It
        w[2 * K:3 * K] = -2 * eta * Dm / It
        if D:
            w[3 * K:] = beta * q * r ** (q - 1) / mu
        g = 2 * np.real(1j * c * (H @ (np.conj(E) * w)))
        return C, g

    res = minimize(cost, phases0, jac=True, method="L-BFGS-B", options={"maxiter": iters})
    return res.x


def synthesise(opa: OPA, x_targets, site_pitch, beta=10.0, gamma=20.0, eta=20.0, seed=0, mode="dense", q=2,
               r_core=None):
    """GSW initial guess followed by gradient refinement. Returns phases.

    mode='sites': darken only untargeted lattice sites.
    mode='dense': darken a grid (spacing lambda/(4 NA F_fill)) everywhere farther than r_core from any target
                  (default r_core = 1.5 resolution cells, which keeps each spot's main lobe free).
    """
    x_targets = np.asarray(x_targets, float)
    if mode == "sites":
        pts = np.arange(np.ceil(-opa.x_range / site_pitch), np.floor(opa.x_range / site_pitch) + 1) * site_pitch
        excl = 0.5 * site_pitch
    else:
        pts = np.arange(-opa.x_range, opa.x_range, opa.spot_pitch / 2)
        excl = r_core if r_core is not None else 1.5 * opa.spot_pitch
    dark = pts[np.min(np.abs(pts[:, None] - x_targets[None, :]), axis=1) > excl]
    phi = gsw(opa, x_targets, seed=seed)
    if beta > 0 or gamma > 0 or eta > 0:
        phi = refine(opa, phi, x_targets, dark, beta=beta, gamma=gamma, eta=eta, q=q if mode == "dense" else 1)
    return phi
