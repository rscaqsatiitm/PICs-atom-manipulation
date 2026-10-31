"""Vector FEM mode solving of SiN channel waveguides (femwell, Nedelec order 2).

Coordinates in micrometres. The core sits on the buried oxide at y = 0 and is centred at x = 0
(or at the given offsets for multi-core cross-sections).
"""
from __future__ import annotations

import contextlib
import io
import os
from collections import OrderedDict
from dataclasses import dataclass, field

import numpy as np
from shapely.geometry import box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from . import materials as mat


@dataclass
class XSection:
    """Rectangular SiN cores in an oxide cladding.

    cores: list of (x_centre, width) in um, all of thickness ``t``.
    x_left/x_right: extent of the computational window (um).
    absorber: width (um) of a lossy layer at the right edge (for bend radiation); 0 disables.
    substrate: if True, a Si substrate below the BOX is added with a lossy bottom layer (leakage).
    """

    cores: list
    t: float
    x_left: float = -4.0
    x_right: float = 4.0
    t_box: float = 4.0
    t_top: float = 3.0
    clad_bottom: float = 3.0  # oxide below core included in window when substrate is False
    absorber: float = 0.0
    substrate: bool = False
    t_si: float = 1.5
    res_core: float = 0.03
    res_max: float = 0.4
    dT_core: float = 0.0
    dT_clad: float = 0.0
    extra: dict = field(default_factory=dict)

    def polygons(self):
        y_bot = -(self.t_box if self.substrate else self.clad_bottom)
        polys = OrderedDict()
        for i, (xc, w) in enumerate(self.cores):
            polys[f"core{i}"] = box(xc - w / 2, 0.0, xc + w / 2, self.t)
        if self.absorber > 0:
            polys["absorber"] = box(self.x_right - self.absorber, y_bot, self.x_right, self.t + self.t_top)
        if self.substrate:
            polys["si_abs"] = box(self.x_left, y_bot - self.t_si - 1.0, self.x_right, y_bot - self.t_si)
            polys["si"] = box(self.x_left, y_bot - self.t_si, self.x_right, y_bot)
        polys["clad"] = box(self.x_left, y_bot, self.x_right, self.t + self.t_top)
        return polys

    def resolutions(self):
        r = {f"core{i}": {"resolution": self.res_core, "distance": 1.5} for i in range(len(self.cores))}
        return r


def _mesh(xs: XSection):
    from femwell.mesh import mesh_from_OrderedDict

    with contextlib.redirect_stdout(io.StringIO()):
        m = mesh_from_OrderedDict(xs.polygons(), xs.resolutions(), default_resolution_max=xs.res_max)
    return from_meshio(m)


@dataclass
class ModeResult:
    neff: complex
    te_fraction: float
    gamma_core: float  # fraction of Poynting power in all cores
    aeff: float  # um^2, effective area of transverse E
    mode: object

    @property
    def loss_db_per_cm(self):
        lam = self.mode.wavelength  # um
        return 20 / np.log(10) * 2 * np.pi / lam * np.imag(self.neff) * 1e4


def solve(xs: XSection, lam: float, num_modes: int = 2, radius: float = np.inf, n_guess=None,
          order: int = 2, absorber_k: float = 0.5):
    """Return modes sorted by Re(neff), descending."""
    from femwell.maxwell.waveguide import compute_modes

    mesh = _mesh(xs)
    basis0 = Basis(mesh, ElementTriP0())
    eps = basis0.zeros(dtype=complex)
    n_clad = mat.n_sio2(lam, xs.dT_clad)
    eps[:] = n_clad**2
    centroids = mesh.p[:, mesh.t].mean(axis=1)  # (2, nelements); P0 dof index == element index
    y_bot = -(xs.t_box if xs.substrate else xs.clad_bottom)
    for name in mesh.subdomains:
        els = mesh.subdomains[name]
        if name.startswith("core"):
            eps[els] = mat.n_sin(lam, xs.dT_core) ** 2
        elif name == "absorber":
            # quadratic loss ramp suppresses reflection at the absorber edge
            s_ = (centroids[0, els] - (xs.x_right - xs.absorber)) / xs.absorber
            eps[els] = (n_clad + 1j * absorber_k * np.clip(s_, 0, 1) ** 2) ** 2
        elif name == "si":
            eps[els] = mat.n_si(lam) ** 2
        elif name == "si_abs":
            s_ = ((y_bot - xs.t_si) - centroids[1, els]) / 1.0
            eps[els] = (mat.n_si(lam) + 0.5j * np.clip(s_, 0, 1) ** 2) ** 2
    with contextlib.redirect_stdout(io.StringIO()):
        modes = compute_modes(basis0, eps, wavelength=lam, num_modes=num_modes, order=order,
                              radius=radius, n_guess=n_guess)
    out = []
    cores = [k for k in mesh.subdomains if k.startswith("core")]
    for m in modes:
        ptot = np.real(m.calculate_power())
        pcore = sum(np.real(m.calculate_power(elements=c)) for c in cores)
        out.append(ModeResult(neff=complex(m.n_eff), te_fraction=float(np.real(m.te_fraction)),
                              gamma_core=float(pcore / ptot), aeff=float(m.calculate_effective_area()),
                              mode=m))
    out.sort(key=lambda r: -r.neff.real)
    return out


def guided(modes, lam, margin=1e-3):
    """Modes with Re(neff) above the cladding index (plus margin)."""
    nc = mat.n_sio2(lam)
    return [m for m in modes if m.neff.real > nc + margin]


def fundamental(modes, pol="TE"):
    for m in modes:
        if (pol == "TE" and m.te_fraction > 0.5) or (pol == "TM" and m.te_fraction <= 0.5):
            return m
    return None


def field_on_grid(mode, x, y, component=0):
    """Interpolate a transverse E component (0: Ex, 1: Ey) of a femwell mode on a grid (complex)."""
    from skfem import ElementTriP1, ElementVector

    (et, et_basis), _ = mode.basis.split(mode.E)
    X, Y = np.meshgrid(x, y, indexing="xy")
    pts = np.vstack([X.ravel(), Y.ravel()])
    # project real and imaginary parts onto a P1 vector field for point evaluation
    b_p1 = et_basis.with_element(ElementVector(ElementTriP1()))
    out = np.zeros(pts.shape[1], dtype=complex)
    for part, fac in ((np.real, 1.0), (np.imag, 1j)):
        e_p1 = b_p1.project(et_basis.interpolate(part(et)))
        out += fac * b_p1.interpolator(e_p1)(pts)[component]
    return out.reshape(X.shape)


def quiet_gmsh():
    os.environ.setdefault("OMP_NUM_THREADS", "1")
