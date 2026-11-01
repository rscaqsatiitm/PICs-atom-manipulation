"""Multi-port star-coupler (free-propagation region, FPR) feeding an optical phased array.

Each tweezer has its own input port. Port j sits at lateral offset d_j on the input arc of a confocal FPR of
radius R; the N arms sit at y_n = (n - (N-1)/2) p_a on the output arc. In the paraxial Fraunhofer limit
(R much larger than the apertures) the slab field on the output arc is

    E_j(y) = sqrt(k_s / (2 pi R)) U(k_s y / R) exp(-i k_s d_j y / R),     k_s = 2 pi n_s / lambda,

with U the Fourier transform of the input-port mode u(xi); this normalisation conserves power. Arm n captures
the amplitude a_jn = int E_j(y) e(y - y_n) dy, with e the (normalised) fundamental mode of an isolated arm at
the FPR edge. Hence (i) every port illuminates all N arms with the same envelope |U(k_s y_n/R)|, which acts
as a Gaussian apodisation of the OPA, and (ii) port j imposes a linear phase k_s d_j p_a / R per arm, which
steers its beam to s_j = n_s d_j p_a / (R p_out) after the emitters (pitch p_out). A common phase ramp on
the arms shifts all beams together.

Approximations: 2D effective-index slab; isolated-arm modes (no inter-arm coupling at the FPR edge); flat
arcs (Rowland aberrations ignored); power not captured by arm fundamental modes is lost.
Lengths in micrometres.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..mmi import Lateral


def _mode(lat: Lateral, width, half=4.0):
    x = np.arange(-half, half + lat.dx / 2, lat.dx)
    _, v = lat.modes(x, [(0.0, width)])
    m = v[:, 0]
    m = m * np.sign(m[len(m) // 2])
    return x, m / np.sqrt(np.sum(np.abs(m) ** 2) * lat.dx)


@dataclass
class StarCoupler:
    lam: float
    t: float
    N: int
    p_a: float  # arm pitch at the FPR edge (um)
    R: float  # FPR radius (um)
    w_in: float  # input-port waveguide width at the FPR (um)
    w_arm: float  # arm width at the FPR edge (um)
    dy: float = 0.02
    lat: Lateral = field(init=False)

    def __post_init__(self):
        self.lat = Lateral(self.lam, self.t)
        self.n_s = self.lat.n_s
        self.k_s = 2 * np.pi * self.n_s / self.lam
        self.xi, self.u = _mode(self.lat, self.w_in)
        self.eta_x, self.e = _mode(self.lat, self.w_arm)
        self.y_n = (np.arange(self.N) - (self.N - 1) / 2) * self.p_a
        Y = self.N * self.p_a / 2 + 4 * self.p_a
        self.y = np.arange(-Y, Y + self.dy / 2, self.dy)
        # far field of the port mode on the output arc (without the port-offset phase)
        q = self.k_s * self.y / self.R
        U = np.exp(-1j * np.outer(q, self.xi)) @ self.u * self.lat.dx
        self.E0 = np.sqrt(self.k_s / (2 * np.pi * self.R)) * U

    def field(self, d):
        return self.E0 * np.exp(-1j * self.k_s * d * self.y / self.R)

    def arm_amplitudes(self, d):
        """Complex amplitudes coupled into the N arms from a port at offset d (unit input power)."""
        E = self.field(d)
        e_on_y = np.interp(self.y[:, None] - self.y_n[None, :], self.eta_x, self.e, left=0.0, right=0.0)
        return (E[:, None] * e_on_y).sum(axis=0) * self.dy

    def capture(self, d):
        return float(np.sum(np.abs(self.arm_amplitudes(d)) ** 2))

    def envelope_edge(self):
        """Arm-power envelope at the outermost arm relative to the centre (1/e^2 = 0.135)."""
        a = np.abs(self.arm_amplitudes(0.0)) ** 2
        return float(a[0] / a[self.N // 2])

    def s_of_port(self, d, p_out):
        """Beam direction sin(theta) after the emitters for a port at offset d."""
        return self.n_s * d * self.p_a / (self.R * p_out)


def _q_halfwidth(sc: StarCoupler):
    """1/e^2 half-width of the port far field in transverse wavenumber q (rad/um), window independent."""
    q = np.linspace(-6 * np.pi / sc.w_in, 6 * np.pi / sc.w_in, 4001)
    U = np.exp(-1j * np.outer(q, sc.xi)) @ sc.u * sc.lat.dx
    I = np.abs(U) ** 2
    return float(0.5 * np.ptp(q[I >= I.max() * np.exp(-2)]))


def far_field_halfwidth(sc: StarCoupler):
    """1/e^2 half-width (um) of the port far field on the output arc, y = q R / k_s."""
    return _q_halfwidth(sc) * sc.R / sc.k_s


def design_radius(lam, t, N, p_a, w_in, w_arm, edge=np.exp(-2)):
    """Choose R so that a Gaussian-equivalent arm-power envelope falls to 'edge' at the outermost arm.

    The 1/e^2 far-field half-width scales as y_e = q_e R / k_s; for I(y) = exp(-2 y^2 / y_e^2) the outer arm
    at y = (N-1) p_a / 2 sees 'edge' when y_e = y_outer sqrt(-2 / ln(edge)). For non-Gaussian port modes
    (wide tapers) the tails differ, so the realised edge level should be checked (StarCoupler.envelope_edge).
    """
    sc = StarCoupler(lam, t, 2, p_a, 100.0, w_in, w_arm)
    target = (N - 1) / 2 * p_a * np.sqrt(-2 / np.log(edge))
    return target * sc.k_s / _q_halfwidth(sc)
