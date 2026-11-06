"""Module 19: quantum check of the low-heating transport points of module 18.

1D split-step Fourier propagation of a Cs atom starting in the ground state of its trap:
  y   handoff axis: static site (1075 nm Gaussian, w0 0.85 um) plus the moving tweezer's y profile (from the 3D grid
      of m18 at x = z = 0), gravity along y. Sequence as in m18: ramp up 100 us, 2 um to the lane in 100 us, wait
      400 us (the x move), back in 100 us, ramp down 100 us, hold 100 us; 200 kHz sample-and-hold.
  x   move axis: the moving tweezer's x profile alone (in the lane the static sites are 1.6e-5 of their depth away),
      12 um or 40 um minimum-jerk in 400 us, 200 kHz sample-and-hold.
The final state is projected on the eigenstates of the final trap (static site for y, moving trap for x):
ground-state probability and mean energy gain. The same 1D problem is solved classically from a Wigner (Gaussian)
sample of the initial ground state, to test the classical model of m18 at these heating levels.

Output: results/m19_quantum.json
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.constants import hbar, k as kB, g as G0
from scipy.linalg import eigh_tridiagonal

from picsim import params as prm
from picsim.provenance import RESULTS, save_result
from picsim.transport3d import min_jerk_u

SCRIPT = "scripts/m19_quantum.py"
M = prm.CS_MASS.value
W0_S = 0.85e-6
F_U = 200e3


def profiles():
    d = np.load(RESULTS / "m18_trap.npz")
    g, x, y, z = d["g"].astype(float), d["x"], d["y"], d["z"]
    i, j, k = np.argmin(abs(x)), np.argmin(abs(y)), np.argmin(abs(z))
    return (x, g[:, j, k]), (y, g[i, :, k])


def interp(prof, q):
    a, b = prof
    return np.interp(q, a, b, left=0.0, right=0.0)


def ground_state(V, dq, n_states=61):
    n = len(V)
    t = hbar**2 / (2 * M * dq**2)
    w, v = eigh_tridiagonal(V + 2 * t, -t * np.ones(n - 1), select="i", select_range=(0, n_states - 1))
    return w, v / np.sqrt(dq)


def schedule_y(t, Um, T_ramp=100e-6, T_lane=100e-6, T_wait=400e-6, T_hold=100e-6, lane=2e-6, staircase=True):
    if staircase:
        t = np.floor(t * F_U) / F_U
    b = np.cumsum([0, T_ramp, T_lane, T_wait, T_lane, T_ramp, T_hold])
    if t < b[1]:
        return 0.0, Um * np.sin(0.5 * np.pi * t / T_ramp) ** 2
    if t < b[2]:
        return lane * min_jerk_u((t - b[1]) / T_lane), Um
    if t < b[3]:
        return lane, Um
    if t < b[4]:
        return lane * (1 - min_jerk_u((t - b[3]) / T_lane)), Um
    if t < b[5]:
        return 0.0, Um * np.cos(0.5 * np.pi * (t - b[4]) / T_ramp) ** 2
    return 0.0, 0.0


def potential_y(q, t, Us, Um, py, T_ramp=100e-6, T_lane=100e-6):
    c, dep = schedule_y(t, Um, T_ramp=T_ramp, T_lane=T_lane)
    return -Us * np.exp(-2 * q**2 / W0_S**2) - dep * interp(py, q - c) + M * G0 * q


def run_quantum(Vt, q, t_end, psi0, dt):
    dq = q[1] - q[0]
    k = 2 * np.pi * np.fft.fftfreq(len(q), dq)
    Kh = np.exp(-1j * hbar * k**2 / (2 * M) * dt)
    psi = psi0.astype(complex)
    t = 0.0
    while t < t_end:
        V = Vt(q, t + 0.5 * dt)
        psi *= np.exp(-0.5j * V * dt / hbar)
        psi = np.fft.ifft(Kh * np.fft.fft(psi))
        psi *= np.exp(-0.5j * V * dt / hbar)
        t += dt
    return psi


def run_classical(Vt, q, t_end, x0, v0, dt):
    x, v = x0.copy(), v0.copy()
    dq = q[1] - q[0]

    def F(x, t):
        Vg = Vt(q, t)
        return -np.interp(x, q, np.gradient(Vg, dq))

    a = F(x, 0.0) / M
    t = 0.0
    while t < t_end:
        v += 0.5 * dt * a
        x += dt * v
        t += dt
        a = F(x, t) / M
        v += 0.5 * dt * a
    return x, v


def analyse(psi, Vf, q, w_f, v_f, E0_i):
    dq = q[1] - q[0]
    c = v_f.T @ psi * dq
    P = np.abs(c) ** 2
    E = np.sum(P * w_f) / max(P.sum(), 1e-300)
    return {"P_ground": float(P[0]), "P_bound_states_used": float(P.sum()),
            "E_gain_uK": float((E - w_f[0]) / kB * 1e6), "E0_final_minus_initial_uK": float((w_f[0] - E0_i) / kB * 1e6)}


def case_y(Us_mK, Um_mK, py, n_cl=2000, T_ramp=100e-6, T_lane=100e-6):
    Us, Um = Us_mK * 1e-3 * kB, Um_mK * 1e-3 * kB
    q = np.arange(-3e-6, 6e-6, 4e-9)
    V0 = -Us * np.exp(-2 * q**2 / W0_S**2) + M * G0 * q
    w0, v0 = ground_state(V0, q[1] - q[0])
    t_end = 2 * T_ramp + 2 * T_lane + 400e-6 + 100e-6
    Vt = lambda qq, t: potential_y(qq, t, Us, Um, py, T_ramp, T_lane)
    psi = run_quantum(Vt, q, t_end, v0[:, 0], 5e-9)
    qres = analyse(psi, V0, q, w0, v0, w0[0])
    # classical from a Wigner sample of the ground state (harmonic approximation of the static site)
    om = np.sqrt(4 * Us / (M * W0_S**2))
    rng = np.random.default_rng(0)
    x0 = rng.normal(0, np.sqrt(hbar / (2 * M * om)), n_cl)
    vv = rng.normal(0, np.sqrt(hbar * om / (2 * M)), n_cl)
    E_i = 0.5 * M * vv**2 + np.interp(x0, q, V0)
    x1, v1 = run_classical(Vt, q, t_end, x0, vv, 50e-9)
    E_f = 0.5 * M * v1**2 + np.interp(x1, q, V0)
    bound = E_f < np.interp(3e-6, q, V0)
    return {"axis": "y (handoff + lane)", "Us_mK": Us_mK, "Um_mK": Um_mK, "T_ramp_us": T_ramp * 1e6,
            "T_lane_us": T_lane * 1e6, "hbar_omega_uK": float(hbar * om / kB * 1e6),
            "quantum": qres, "classical_wigner": {"E_gain_uK": float((E_f[bound].mean() - E_i[bound].mean()) / kB * 1e6),
                                                  "lost_fraction": float(1 - bound.mean())}}


def case_x(Um_mK, d_um, px, T=400e-6, n_cl=2000):
    Um = Um_mK * 1e-3 * kB
    q = np.arange(-3e-6, (d_um + 3) * 1e-6, 4e-9)
    traj = lambda t: d_um * 1e-6 * min_jerk_u((np.floor(t * F_U) / F_U) / T)
    Vt = lambda qq, t: -Um * interp(px, qq - traj(t))
    V0 = Vt(q, 0.0)
    w0, v0 = ground_state(V0, q[1] - q[0])
    t_end = T * 1.25
    psi = run_quantum(Vt, q, t_end, v0[:, 0], 5e-9)
    Vf = Vt(q, t_end)
    wf, vf = ground_state(Vf, q[1] - q[0])
    qres = analyse(psi, Vf, q, wf, vf, w0[0])
    curv = np.gradient(np.gradient(V0, q), q)[np.argmin(np.abs(q))]
    om = np.sqrt(curv / M)
    rng = np.random.default_rng(0)
    x0 = rng.normal(0, np.sqrt(hbar / (2 * M * om)), n_cl)
    vv = rng.normal(0, np.sqrt(hbar * om / (2 * M)), n_cl)
    E_i = 0.5 * M * vv**2 + np.interp(x0, q, V0)
    x1, v1 = run_classical(Vt, q, t_end, x0, vv, 50e-9)
    E_f = 0.5 * M * v1**2 + np.interp(x1, q, Vf)
    return {"axis": "x (move)", "Um_mK": Um_mK, "d_um": d_um, "hbar_omega_uK": float(hbar * om / kB * 1e6),
            "quantum": qres, "classical_wigner": {"E_gain_uK": float((E_f.mean() - E_i.mean()) / kB * 1e6)}}


def profile_z():
    d = np.load(RESULTS / "m18_trap.npz")
    g, x, y, z = d["g"].astype(float), d["x"], d["y"], d["z"]
    return z, g[np.argmin(abs(x)), np.argmin(abs(y)), :]


def case_z(Us_mK, Um_mK, pz, dz_um, T_ramp, T_lane, n_cl=500):
    """Axial motion during handoff: static site (Gaussian beam on axis, focus at dz) weighted by the lateral
    separation of the atom (at the moving-trap centre) from the static site, plus the moving trap's axial profile."""
    Us, Um = Us_mK * 1e-3 * kB, Um_mK * 1e-3 * kB
    zR = np.pi * W0_S**2 / (prm.LAM_STATIC.value * 1e-6)
    q = np.arange(-12e-6, 12e-6, 4e-9)
    dz = dz_um * 1e-6

    def Vt(qq, t):
        c, dep = schedule_y(t, Um, T_ramp=T_ramp, T_lane=T_lane)
        return -Us * np.exp(-2 * c**2 / W0_S**2) / (1 + ((qq - dz) / zR) ** 2) - dep * interp(pz, qq)

    V0 = Vt(q, 0.0)
    w0, v0 = ground_state(V0, q[1] - q[0], n_states=200)
    t_end = 2 * T_ramp + 2 * T_lane + 400e-6 + 100e-6
    psi = run_quantum(Vt, q, t_end, v0[:, 0], 5e-9)
    qres = analyse(psi, V0, q, w0, v0, w0[0])
    om = np.sqrt(2 * Us / (M * zR**2))
    rng = np.random.default_rng(0)
    x0 = dz + rng.normal(0, np.sqrt(hbar / (2 * M * om)), n_cl)
    vv = rng.normal(0, np.sqrt(hbar * om / (2 * M)), n_cl)
    E_i = 0.5 * M * vv**2 + np.interp(x0, q, V0)
    x1, v1 = run_classical(Vt, q, t_end, x0, vv, 50e-9)
    E_f = 0.5 * M * v1**2 + np.interp(x1, q, V0)
    return {"axis": "z (handoff, axial)", "Us_mK": Us_mK, "Um_mK": Um_mK, "focus_offset_um": dz_um,
            "T_ramp_us": T_ramp * 1e6, "T_lane_us": T_lane * 1e6, "hbar_omega_uK": float(hbar * om / kB * 1e6),
            "quantum": qres, "classical_wigner": {"E_gain_uK": float((E_f.mean() - E_i.mean()) / kB * 1e6)}}


def axial():
    pz = profile_z()
    rows = []
    for Us, Um in ((0.1, 0.3), (0.1, 0.5)):
        for Tr, Tl in ((100e-6, 100e-6), (200e-6, 400e-6)):
            for dz in (0.0, 0.5, 1.0, 2.0):
                r = case_z(Us, Um, pz, dz, Tr, Tl)
                rows.append(r)
                print(Us, Um, Tr, Tl, dz, "Q dE %.3f P0 %.3f bound %.4f | Cl %.3f" % (
                    r["quantum"]["E_gain_uK"], r["quantum"]["P_ground"], r["quantum"]["P_bound_states_used"],
                    r["classical_wigner"]["E_gain_uK"]), flush=True)
    save_result("m19_axial", script=SCRIPT,
                inputs={"model": "1D along the optical axis; the static site's depth follows the lateral separation of "
                        "the moving-trap centre (atom assumed at that centre)", "grid": "4 nm, dt 5 ns",
                        "profiles": "moving: m18_trap.npz axial cut; static: Gaussian beam, w0 0.85 um, 1075 nm"},
                outputs={"rows": rows},
                notes="P_bound_states_used: weight in the 200 lowest static states (projection completeness).")
    return rows


def main():
    px, py = profiles()
    rows = []
    for Us, Um in ((0.1, 0.3), (0.1, 0.5), (0.2, 0.6), (0.2, 1.0)):
        r = case_y(Us, Um, py)
        rows.append(r)
        print(r, flush=True)
    for Us, Um in ((0.1, 0.3), (0.1, 0.5), (0.2, 0.6)):
        for Tr, Tl in ((100e-6, 200e-6), (200e-6, 200e-6), (200e-6, 400e-6)):
            r = case_y(Us, Um, py, n_cl=500, T_ramp=Tr, T_lane=Tl)
            rows.append(r)
            print(r, flush=True)
    for Um in (0.3, 0.5):
        for d in (12.0, 40.0):
            r = case_x(Um, d, px)
            rows.append(r)
            print(r, flush=True)
    save_result("m19_quantum", script=SCRIPT,
                inputs={"grid": "4 nm, split-step dt 5 ns", "profiles": "results/m18_trap.npz (x and y cuts through focus)",
                        "static_w0_um": W0_S * 1e6, "update": "200 kHz sample-and-hold", "classical": "2000-atom Wigner "
                        "sample of the initial ground state (harmonic approximation), velocity Verlet 50 ns"},
                outputs={"rows": rows},
                notes="E_gain: mean energy above the final ground state. 1D only: axial motion and the coupling between "
                      "axes are not included.")
    return rows


if __name__ == "__main__":
    {"all": main, "axial": axial}[sys.argv[1] if len(sys.argv) > 1 else "all"]()
