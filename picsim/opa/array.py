"""1D optical phased array: emitter geometry, near field and Fourier-filtered relay to the atom plane.

Optical model (scalar, 1D along the steering axis x; the cross axis is treated separately):

  chip near field  u(x) = sum_n a_n exp(i phi_n) e(x - x_n)
  Fourier plane    U(s) = FT[u](k0 s), s = sin(theta); a slit keeps |s| < s_cut (default lambda/(2p))
  pupil            the filtered near field is imaged onto the meta-objective pupil with magnification M;
                   the pupil (diameter D = 2 F NA, Abbe sine condition) truncates it to |x| <= A/(2 F_fill)
  atom plane       E(x_a) = FT[pupil field] evaluated at s = x_a * 2 NA F_fill / A

Here A = N p is the OPA aperture and F_fill = M A / D the pupil fill factor (1: OPA image just fills the pupil).
Hence the atom-plane field is a scaled, slit-filtered copy of the OPA far field, and the steering range at the
atoms is |x_a| < s_cut A / (2 NA F_fill) = N lambda / (4 NA F_fill) for s_cut = lambda/(2p).

Lengths in micrometres.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def gaussian_element(w_e):
    """Emitter near-field amplitude: Gaussian with 1/e^2 intensity radius w_e (um)."""
    return lambda x: np.exp(-(x / w_e) ** 2)


@dataclass
class OPA:
    N: int
    pitch: float  # um
    lam: float  # um
    element: callable  # near-field amplitude profile e(x) of one emitter
    NA: float = 0.5
    fill: float = 1.0  # pupil fill factor F_fill
    s_cut: float | None = None  # Fourier slit half-width in s; default lambda/(2 p)
    apod: np.ndarray | None = None  # channel amplitudes (default uniform)
    dx: float = 0.02  # near-field grid step (um)

    def __post_init__(self):
        if self.s_cut is None:
            self.s_cut = min(self.lam / (2 * self.pitch), 1.0)
        self.xn = (np.arange(self.N) - (self.N - 1) / 2) * self.pitch
        if self.apod is None:
            self.apod = np.ones(self.N)
        self.A = self.N * self.pitch
        L = self.A + 12 * self.pitch
        # FFT grid: pad so that the s-sampling resolves the beam (d s << lambda / A)
        nfft = int(2 ** np.ceil(np.log2(16 * L / self.dx)))
        self.x = (np.arange(nfft) - nfft // 2) * self.dx
        self.s = np.fft.fftshift(np.fft.fftfreq(nfft, d=self.dx)) * self.lam  # s = lambda * f_x
        e = self.element(self.x)
        # continuous FT of one emitter, E_el(s) = int e(x) exp(-i k s x) dx, on the s-grid
        self._E_el = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(e))) * self.dx
        self._ds = self.s[1] - self.s[0]

    # ----------------------------------------------------------------- basics
    @property
    def x_scale(self):
        """Atom-plane position per unit s: x_a = s * A / (2 NA F_fill)."""
        return self.A / (2 * self.NA * self.fill)

    @property
    def spot_pitch(self):
        """Diffraction-limited resolution cell at the atoms, lambda/(2 NA F_fill)."""
        return self.lam / (2 * self.NA * self.fill)

    @property
    def x_range(self):
        """Half-range of steering at the atoms set by the slit, s_cut * x_scale."""
        return self.s_cut * self.x_scale

    def steer_phases(self, s0):
        """Linear phase ramp that steers the main lobe to sin(theta) = s0."""
        return 2 * np.pi / self.lam * self.xn * s0

    def x_to_s(self, x_a):
        return np.asarray(x_a) / self.x_scale

    def array_factor(self, phases, s=None):
        s = self.s if s is None else s
        k = 2 * np.pi / self.lam
        c = self.apod * np.exp(1j * np.asarray(phases))
        return np.exp(-1j * k * np.outer(s, self.xn)) @ c

    def far_field(self, phases):
        """Far-field amplitude U(s) on self.s (element factor x array factor).

        Only the propagating band |s| < 1 is computed; evanescent components never reach the slit or the
        atoms, and the emitted power is defined on |s| < 1.
        """
        if not hasattr(self, "_prop"):
            self._prop = np.nonzero(np.abs(self.s) < 1)[0]
        U = np.zeros_like(self._E_el)
        idx = self._prop
        U[idx] = self._E_el[idx] * self.array_factor(phases, self.s[idx])
        return U

    # ----------------------------------------------------------------- relay
    def _powers(self, phases):
        """Return emitted power, filtered near field u_t (after slit and pupil) and power fractions.

        Powers in units of int |u|^2 dx; by Parseval int |U|^2 ds = lambda int |u|^2 dx.
        """
        U = self.far_field(phases)
        prop = np.abs(self.s) < 1
        p_em = np.sum(np.abs(U[prop]) ** 2) * self._ds / self.lam
        Uf = np.where(np.abs(self.s) < self.s_cut, U, 0)
        p_slit = np.sum(np.abs(Uf) ** 2) * self._ds / self.lam
        u = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(Uf))) / self.dx
        u_t = np.where(np.abs(self.x) <= self.A / (2 * self.fill), u, 0)
        p_pupil = np.sum(np.abs(u_t) ** 2) * self.dx
        return p_em, u_t, p_slit / p_em, p_pupil / p_em

    def slit_efficiency(self, phases):
        """Fraction of radiated power passing the Fourier slit |s| < s_cut."""
        return self._powers(phases)[2]

    def transmission(self, phases):
        """Fraction of radiated power passing slit and pupil."""
        return self._powers(phases)[3]

    def atom_intensity(self, phases, x_a):
        """Atom-plane intensity (1/um) along the steering axis, normalised to emitted power.

        int I dx_a equals the fraction of emitted power transmitted through slit and pupil.
        """
        p_em, u_t, _, _ = self._powers(phases)
        nz = np.nonzero(u_t)[0]
        k = 2 * np.pi / self.lam
        s_eval = np.atleast_1d(self.x_to_s(x_a))
        E = np.empty(len(s_eval), dtype=complex)
        for i in range(0, len(s_eval), 128):  # chunked direct transform (memory)
            E[i:i + 128] = np.exp(-1j * k * np.outer(s_eval[i:i + 128], self.x[nz])) @ u_t[nz] * self.dx
        # int |E|^2 ds = lambda int |u_t|^2 dx and ds = dx_a / x_scale
        return np.abs(E) ** 2 / (self.lam * self.x_scale * p_em)


def channel_response(opa: "OPA", x_a):
    """Complex matrix G (N, M): atom-plane field at x_a from channel n with unit amplitude, through slit and
    pupil. For channel amplitudes c_n = apod_n exp(i phi_n), the intensity is |c @ G|^2 / (lambda x_scale p_em),
    identical to OPA.atom_intensity."""
    k = 2 * np.pi / opa.lam
    slit = np.abs(opa.s) < opa.s_cut
    pupil = np.abs(opa.x) <= opa.A / (2 * opa.fill)
    nz = np.nonzero(pupil)[0]
    D = np.exp(-1j * k * np.outer(opa.x[nz], np.atleast_1d(opa.x_to_s(x_a)))) * opa.dx  # (nz, M)
    G = np.empty((opa.N, D.shape[1]), dtype=complex)
    sl = np.nonzero(slit)[0]
    for n, xn in enumerate(opa.xn):
        U = np.zeros_like(opa._E_el)
        U[sl] = opa._E_el[sl] * np.exp(-1j * k * opa.s[sl] * xn)
        u = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(U))) / opa.dx
        G[n] = u[nz] @ D
    return G


def spot_metrics(x, I):
    """Peak position (parabolic interpolation), 1/e^2 half-width and FWHM of the main peak of I(x)."""
    i = int(np.argmax(I))
    if 0 < i < len(I) - 1:
        y0, y1, y2 = np.log(I[i - 1] + 1e-300), np.log(I[i] + 1e-300), np.log(I[i + 1] + 1e-300)
        d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    else:
        d = 0.0
    dx = x[1] - x[0]
    xp = x[i] + d * dx
    pk = I[i]

    def width(level):
        j = i
        while j > 0 and I[j] > level * pk:
            j -= 1
        k = i
        while k < len(I) - 1 and I[k] > level * pk:
            k += 1
        # linear interpolation at both edges
        xl = x[j] + (level * pk - I[j]) / (I[j + 1] - I[j] + 1e-300) * dx
        xr = x[k - 1] + (I[k - 1] - level * pk) / (I[k - 1] - I[k] + 1e-300) * dx
        return xr - xl

    return {"x_peak": xp, "peak": pk, "fwhm": width(0.5), "w_1e2": width(np.exp(-2)) / 2}
