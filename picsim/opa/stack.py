"""Stacked-chip anamorphic relay: the cross axis (y) of one chip, and the 2D focus through a circular pupil.

Geometry (lengths in um unless stated). M chips are stacked in y at pitch d_c, facets flush in one plane. Each
chip has its own fast-axis collimator (FAC, a cylindrical lens of focal length f_FAC with its front focal plane on
the facets). The FAC back focal planes form the plane S. A common y-relay of effective focal length F_y Fourier
transforms S onto the objective pupil, and the objective (focal length f_obj, numerical aperture NA, Abbe sine
condition) transforms the pupil onto the atoms. Hence, in y:

    facet --(FAC: Fourier)--> S --(F_y: Fourier)--> pupil --(objective: Fourier)--> atoms

  * S is conjugate to the atoms with magnification 1/D, D = F_y / f_obj; chip m lands at y_a = m d_c / D, so the row
    pitch is Delta = d_c / D.
  * The pupil is conjugate to the facet with magnification F_y / f_FAC. The pupil radius f_obj NA, referred back to
    the facet, is beta = NA f_FAC / D. beta is the design parameter that sets the y fill of the pupil.
  * The FAC aperture (numerical aperture NA_FAC) clips the facet's angular spectrum at |s| < NA_FAC, where s is the
    sine of the emission angle in air.

The x axis is the slit-filtered relay of picsim.opa.array (emitter plane imaged onto the pupil). Both axes share
the circular pupil, xi^2 + eta^2 <= 1, with xi = x/(A/2) and eta = y'/beta normalised pupil coordinates. The
atom-plane field is

    E(x_a, y_a) = sum_xi sum_eta u_x(xi) P_y(eta) circ(xi, eta) exp(-i k NA (xi x_a + eta y_a)),

scalar, without the objective's apodisation factor (as in the x model).

Per-chip errors enter the y pupil field as phases on the angular spectrum: FAC decentre delta (facet relative to
the FAC axis) shifts the facet image in the pupil by delta/beta; FAC focus error dz multiplies by
exp(i k dz cos theta); FAC aberration multiplies by exp(i 2 pi W(s/NA_FAC)) with W in waves.
"""
from __future__ import annotations

import numpy as np

from .. import modes as M
from .array import channel_pupil

_FACET_CACHE = {}


def raised_cosine(s1, s2):
    """Amplitude transmission 1 for |s| < s1, cosine roll-off to 0 at |s| = s2 (graded FAC aperture)."""
    def T(s):
        u = np.clip((np.abs(s) - s1) / (s2 - s1), 0, 1)
        return np.cos(0.5 * np.pi * u)
    return T


def facet_profile(lam, width, t, dy=0.005):
    """Vertical facet field of the TE0 emitter mode, integrated over x (the theta_x = 0 cut of the far field).

    Returns (y, Y) with y measured from the film centre and Y normalised to max |Y| = 1.
    """
    key = (round(lam, 6), round(width, 4), round(t, 4), dy)
    if key not in _FACET_CACHE:
        ms = M.solve(M.XSection(cores=[(0, width)], t=t, x_left=-4 - width / 2, x_right=4 + width / 2), lam,
                     num_modes=2)
        te = M.fundamental(ms)
        x = np.arange(-3.5 - width / 2, 3.5 + width / 2, 0.01)
        y = np.arange(-2.9, t + 2.9, dy)
        F = M.field_on_grid(te.mode, x, y, 0)  # Ex on (y, x)
        Y = F.sum(axis=1) * 0.01
        Y = Y * np.exp(-1j * np.angle(Y[np.argmax(np.abs(Y))]))
        _FACET_CACHE[key] = (y - t / 2, Y / np.abs(Y).max(), float(te.neff.real))
    return _FACET_CACHE[key]


class FacetY:
    """Angular spectrum A(s) in air of a vertical facet profile (femwell mode, or Gaussian of 1/e^2 radius w)."""

    def __init__(self, lam, width=1.2, t=0.4, gaussian_w=None, dy=0.005, nfft=2**16):
        self.lam = lam
        if gaussian_w is not None:
            y = (np.arange(-1200, 1201)) * dy
            Y = np.exp(-(y / gaussian_w) ** 2).astype(complex)
            self.neff = None
        else:
            y, Y, self.neff = facet_profile(lam, width, t, dy)
        self.y, self.Y = y, Y
        pad = np.zeros(nfft, dtype=complex)
        i0 = nfft // 2 - np.argmin(np.abs(y))
        pad[i0:i0 + len(y)] = Y
        s = np.fft.fftshift(np.fft.fftfreq(nfft, d=dy)) * lam
        A = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(pad))) * dy
        keep = np.abs(s) < 1
        self.s, self.A = s[keep], A[keep]
        self.ds = s[1] - s[0]
        self.cos = np.sqrt(1 - self.s**2)
        dens = np.abs(self.A) ** 2 * self.cos  # radiated power per unit s
        self.P_rad = float(np.sum(dens) * self.ds)
        self._dens = dens

    def fraction_within(self, NA_FAC):
        return float(np.sum(self._dens[np.abs(self.s) < NA_FAC]) * self.ds / self.P_rad)

    @property
    def NA_1e2(self):
        """Sine of the 1/e^2 half-angle of the far-field intensity |A|^2 cos(theta)."""
        d = self._dens / self._dens.max()
        pos = self.s >= 0
        return float(np.interp(np.exp(-2), d[pos][::-1], self.s[pos][::-1]))

    @property
    def w_equiv(self):
        """1/e^2 radius of the Gaussian with the same far-field divergence, lambda / (pi NA_1e2)."""
        return self.lam / (np.pi * self.NA_1e2)


class RelayY:
    """Cross-axis (y) relay of one chip: facet -> FAC -> S -> pupil -> atoms."""

    def __init__(self, facet: FacetY, beta, NA=0.5, NA_FAC=0.8, n_eta=801, apod=None):
        self.f, self.beta, self.NA, self.NA_FAC = facet, beta, NA, NA_FAC
        self.apod = apod  # amplitude transmission T(s) of the FAC aperture (graded edge); None: hard edge
        self.lam = facet.lam
        self.k = 2 * np.pi / self.lam
        self.eta = np.linspace(-1, 1, n_eta)
        self.deta = self.eta[1] - self.eta[0]
        # all radiated power, referred to the pupil coordinate eta: int |P_y|^2 d eta = lam P_rad / beta
        self.P_full = self.lam * facet.P_rad / beta

    def spectrum(self, dz=0.0, delta=0.0, W=None):
        """Angular spectrum after the FAC, sqrt(cos) power factor, aperture and phase errors included."""
        f = self.f
        a = f.A * np.sqrt(f.cos)
        if self.apod is not None:
            a = a * self.apod(f.s)
        m = np.abs(f.s) < self.NA_FAC
        ph = self.k * dz * f.cos - self.k * f.s * delta
        if W is not None:
            ph = ph + 2 * np.pi * W(f.s / self.NA_FAC)
        return f.s[m], a[m] * np.exp(1j * ph[m])

    def pupil(self, eta=None, **err):
        """Pupil field P_y(eta), eta in units of the pupil radius (facet image, magnified)."""
        eta = self.eta if eta is None else np.asarray(eta)
        s, a = self.spectrum(**err)
        yp = self.beta * eta
        out = np.empty(yp.shape, dtype=complex)
        flat = yp.ravel()
        o = out.ravel()
        for i in range(0, len(flat), 4096):
            o[i:i + 4096] = np.exp(1j * self.k * np.outer(flat[i:i + 4096], s)) @ a * self.f.ds
        return out

    @property
    def norm(self):
        """Factor that makes |E_y|^2 * norm a line intensity (1/um) per unit radiated facet power."""
        return self.NA / (self.lam * self.P_full)

    def atom_line(self, y_a, Py=None, dz_a=0.0, **err):
        """1D y profile at the atoms (no circle; the x pupil is not involved). Integral = transmitted fraction."""
        Py = self.pupil(**err) if Py is None else Py
        ph = np.exp(1j * self.k * dz_a * np.sqrt(1 - (self.NA * self.eta) ** 2))
        E = np.exp(-1j * self.k * self.NA * np.outer(np.atleast_1d(y_a), self.eta)) @ (Py * ph) * self.deta
        return np.abs(E) ** 2 * self.norm

    def efficiency(self, Py=None, **err):
        """Fraction of the radiated facet power inside |eta| <= 1 (FAC aperture and pupil, no circle)."""
        Py = self.pupil(**err) if Py is None else Py
        return float(np.sum(np.abs(Py) ** 2) * self.deta / self.P_full)


class Focus2D:
    """2D atom-plane field of one chip: x response of a FocalPlane (D1) and the y relay, circular pupil."""

    def __init__(self, fp, relay: RelayY, stride=5, dz_chip=0.0):
        self.fp, self.r = fp, relay
        opa = fp.opa
        xs, u = channel_pupil(opa, dz_chip)
        self.xs, self.u = xs[::stride], u[:, ::stride]
        self.dxp = opa.dx * stride
        self.xi = self.xs / (opa.A / 2)
        self.norm_x = 1.0 / (opa.lam * opa.x_scale * fp.e_pow)
        self.k = 2 * np.pi / opa.lam
        self.opa = opa

    def pupil_x(self, port_w):
        return (port_w @ self.fp.A) @ self.u

    def field(self, port_w, x_a, y_a, Py=None, circle=True, dz_a=0.0, shear=None):
        """Complex field on the grid (y_a, x_a); |E|^2 is the intensity (1/um^2) per unit port power and per unit
        radiated facet power fraction. shear(xi) gives an eta offset per pupil column (FAC roll or smile)."""
        ux = self.pupil_x(port_w)
        r = self.r
        if shear is None:
            Py = r.pupil() if Py is None else Py
            P = ux[:, None] * Py[None, :]
        else:
            P = ux[:, None] * self._sheared(shear)
        rho2 = self.xi[:, None] ** 2 + r.eta[None, :] ** 2
        if circle:
            P = np.where(rho2 <= 1.0, P, 0)
        if dz_a:
            P = P * np.exp(1j * self.k * dz_a * np.sqrt(np.clip(1 - r.NA**2 * rho2, 0, None)))
        Ex = np.exp(-1j * self.k * np.outer(self.opa.x_to_s(np.atleast_1d(x_a)), self.xs)) * self.dxp
        Ey = np.exp(-1j * self.k * r.NA * np.outer(r.eta, np.atleast_1d(y_a))) * r.deta
        E = Ex @ (P @ Ey) if Ey.shape[1] < Ex.shape[0] else (Ex @ P) @ Ey
        return E.T * np.sqrt(self.norm_x * r.norm)

    def _sheared(self, shear):
        r = self.r
        if not hasattr(self, "_Py_wide"):
            self._eta_w = np.linspace(-3, 3, 6 * (len(r.eta) - 1) // 2 + 1)
            self._Py_wide = r.pupil(self._eta_w)
        off = np.asarray(shear(self.xi))
        E = self._eta_w
        re = np.interp(r.eta[None, :] - off[:, None], E, self._Py_wide.real)
        im = np.interp(r.eta[None, :] - off[:, None], E, self._Py_wide.imag)
        return re + 1j * im

    def intensity(self, port_w, x_a, y_a, **kw):
        return np.abs(self.field(port_w, x_a, y_a, **kw)) ** 2


def peak2d(I, x, y):
    """Peak value and sub-grid position (parabolic in x and y) of a 2D intensity on the grid (y, x)."""
    j, i = np.unravel_index(np.argmax(I), I.shape)

    def par(a, b, c):
        a, b, c = np.log(a + 1e-300), np.log(b + 1e-300), np.log(c + 1e-300)
        den = a - 2 * b + c
        return 0.5 * (a - c) / den if den != 0 else 0.0

    dx = par(I[j, i - 1], I[j, i], I[j, i + 1]) if 0 < i < I.shape[1] - 1 else 0.0
    dy = par(I[j - 1, i], I[j, i], I[j + 1, i]) if 0 < j < I.shape[0] - 1 else 0.0
    return float(I[j, i]), float(x[i] + dx * (x[1] - x[0])), float(y[j] + dy * (y[1] - y[0]))


def fwhm(z, I):
    i = int(np.argmax(I))
    h = I[i] / 2
    lo = i
    while lo > 0 and I[lo] > h:
        lo -= 1
    hi = i
    while hi < len(I) - 1 and I[hi] > h:
        hi += 1
    zl = np.interp(h, [I[lo], I[lo + 1]], [z[lo], z[lo + 1]])
    zr = np.interp(h, [I[hi], I[hi - 1]], [z[hi], z[hi - 1]])
    return float(zr - zl)
