"""Thermo-optic phase shifters: 2D steady-state and transient heat conduction (scikit-fem, P1 elements).

Cross-section (x lateral, y vertical, um; SiN core bottom at y = 0):
    Si substrate (depth d_si, bottom held at ambient) | BOX t_box | SiN core | top oxide t_top | metal heater
Sides and top are adiabatic. Optional air trenches through top oxide and BOX beside the waveguide, and an
optional undercut (air) in the Si below the waveguide. The heater dissipates P' (W/m) uniformly.

Phase per unit heater power, with dneff/dT from the mode solver (uniform heating of core and cladding):
    phi / P = k0 dneff/dT * dT_core / P'   ->   P_pi = pi / (k0 dneff/dT dT_core/P')
which is independent of heater length for L >> cross-section.
"""
from __future__ import annotations

import contextlib
import io
from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
from shapely.geometry import box
from skfem import Basis, BilinearForm, ElementTriP0, ElementTriP1, LinearForm, asm, condense, solve
from skfem.helpers import dot
from skfem.io.meshio import from_meshio

from . import params as prm


@dataclass
class HeaterXS:
    t: float = 0.4
    w_core: float = 0.6
    t_box: float = 4.0
    t_top: float = 2.0
    w_h: float = 2.0
    t_h: float = 0.1
    d_si: float = 100.0
    half_width: float = 150.0
    trench: float | None = None  # centre offset of air trenches (um); None disables
    w_trench: float = 4.0
    undercut: float | None = None  # half-width of an air undercut below the BOX (um), depth 10 um
    neighbours: tuple = ()  # x positions of neighbouring waveguides (for crosstalk sampling only)
    source: str = "heater"  # "heater": metal strip dissipates; "core": optical absorption in the core

    def polygons(self):
        y_top = self.t + self.t_top
        P = OrderedDict()
        P["heater"] = box(-self.w_h / 2, y_top, self.w_h / 2, y_top + self.t_h)
        P["core"] = box(-self.w_core / 2, 0, self.w_core / 2, self.t)
        if self.trench is not None:
            for s, nm in ((-1, "trL"), (1, "trR")):
                xc = s * self.trench
                P[nm] = box(xc - self.w_trench / 2, -self.t_box, xc + self.w_trench / 2, y_top)
        if self.undercut is not None:
            P["uc"] = box(-self.undercut, -self.t_box - 10.0, self.undercut, -self.t_box)
        P["oxide"] = box(-self.half_width, -self.t_box, self.half_width, y_top)
        P["si"] = box(-self.half_width, -self.t_box - self.d_si, self.half_width, -self.t_box)
        return P

    def mesh(self):
        from femwell.mesh import mesh_from_OrderedDict

        res = {"heater": {"resolution": 0.05, "distance": 3.0}, "core": {"resolution": 0.05, "distance": 3.0}}
        for nm in ("trL", "trR", "uc"):
            res[nm] = {"resolution": 0.3, "distance": 3.0}
        with contextlib.redirect_stdout(io.StringIO()):
            m = mesh_from_OrderedDict(self.polygons(), res, default_resolution_max=8.0)
        return from_meshio(m)


K_OF = {"heater": prm.K_METAL.value, "core": prm.K_SIN.value, "oxide": prm.K_SIO2.value, "si": prm.K_SI.value,
        "trL": prm.K_AIR.value, "trR": prm.K_AIR.value, "uc": prm.K_AIR.value}
C_OF = {"heater": 2.5e6, "core": prm.RHO_CP_SIN.value, "oxide": prm.RHO_CP_SIO2.value, "si": prm.RHO_CP_SI.value,
        "trL": 1.2e3, "trR": 1.2e3, "uc": 1.2e3}


def _assemble(xs: HeaterXS):
    mesh = xs.mesh()
    b0 = Basis(mesh, ElementTriP0())
    kap = b0.zeros()
    cap = b0.zeros()
    for nm, els in mesh.subdomains.items():
        if nm not in K_OF:
            continue
        kap[els] = K_OF[nm]
        cap[els] = C_OF[nm]
    b = b0.with_element(ElementTriP1())

    @BilinearForm
    def cond(u, v, w):
        return w["k"] * dot(u.grad, v.grad)

    @BilinearForm
    def mass(u, v, w):
        return w["c"] * u * v

    @LinearForm
    def load(v, w):
        return w["q"] * v

    # lengths in um: scale conductivity W/(m K) -> W/(um K) and volumetric heat J/(m^3 K) -> J/(um^3 K)
    K = asm(cond, b, k=b0.interpolate(kap * 1e-6))
    Mm = asm(mass, b, c=b0.interpolate(cap * 1e-18))
    q = b0.zeros()
    if xs.source == "core":
        q[mesh.subdomains["core"]] = 1.0 / (xs.w_core * xs.t)
    else:
        q[mesh.subdomains["heater"]] = 1.0 / (xs.w_h * xs.t_h)  # 1 W/um per unit length -> W/um^3
    F = asm(load, b, q=b0.interpolate(q))
    y_bot = -xs.t_box - xs.d_si
    bottom = mesh.facets_satisfying(lambda x: np.isclose(x[1], y_bot))
    D = b.get_dofs(bottom)
    return mesh, b, K, Mm, F, D


def steady(xs: HeaterXS, x_probe=None):
    """Temperature rise per unit heater power per unit length, dT/P' in K/(W/m), at the core and probes."""
    mesh, b, K, Mm, F, D = _assemble(xs)
    T = solve(*condense(K, F, D=D))  # K per (W/um)
    probes = np.array([[0.0], [xs.t / 2]])
    if x_probe is not None:
        probes = np.hstack([probes, np.vstack([np.asarray(x_probe, float), np.full(len(x_probe), xs.t / 2)])])
    vals = b.probes(probes) @ T
    return vals * 1e-6, (mesh, b, T * 1e-6)  # K/(W/m)


def transient(xs: HeaterXS, t_end=2e-3, n=400):
    """Core temperature after a heater power step (implicit Euler, log-spaced steps). Returns t, dT/dT_ss."""
    mesh, b, K, Mm, F, D = _assemble(xs)
    T_ss = solve(*condense(K, F, D=D))
    pr = b.probes(np.array([[0.0], [xs.t / 2]]))
    ts = np.concatenate([[0.0], np.geomspace(1e-8, t_end, n)])
    T = b.zeros()
    out = [0.0]
    for i in range(1, len(ts)):
        dt = ts[i] - ts[i - 1]
        A = Mm / dt + K
        rhs = Mm @ T / dt + F
        T = solve(*condense(A, rhs, D=D))
        out.append(float((pr @ T)[0]))
    out = np.array(out) / float((pr @ T_ss)[0])
    return ts, out


def rise_time(t, y):
    i10 = np.argmax(y >= 0.1)
    i90 = np.argmax(y >= 0.9)
    return float(t[i90] - t[i10])
