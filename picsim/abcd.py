"""Paraxial ray-transfer (ABCD) helpers, 2x2 and augmented 3x3 (for decentred lenses and prisms). Units are free
but consistent (mm in the relay designs)."""
from __future__ import annotations

import numpy as np


def abcd_free(d):
    return np.array([[1.0, d], [0.0, 1.0]])


def abcd_lens(f):
    return np.array([[1.0, 0.0], [-1.0 / f, 1.0]])


def chain(*ms):
    """Product of matrices in the order the ray meets them (first element acts first)."""
    out = np.eye(ms[0].shape[0])
    for m in ms:
        out = m @ out
    return out


def aug_free(d):
    return np.array([[1.0, d, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])


def aug_lens(f, h=0.0):
    """Thin lens of focal length f centred at transverse position h: theta' = theta - (y - h)/f."""
    return np.array([[1.0, 0.0, 0.0], [-1.0 / f, 1.0, h / f], [0.0, 0.0, 1.0]])


def trace(M, y, theta):
    """Apply an augmented 3x3 matrix to rays (y, theta)."""
    v = np.vstack([np.atleast_1d(y), np.atleast_1d(theta), np.ones(np.size(y))])
    out = M @ v
    return out[0], out[1]
