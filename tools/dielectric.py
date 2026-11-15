"""Causal wideband dielectric model (Djordjevic-Sarkar).

    eps(s) = eps_inf + k * log10((w2 + s) / (w1 + s)),   w = 2*pi*10^m

Between f1 = 10^m1 and f2 = 10^m2 the loss tangent is nearly constant and
eps' falls by k per decade (the Kramers-Kronig consequence of the loss);
below f1 and above f2 the material is lossless and eps' is flat. The function
is analytic in s, so it can be used directly in the circuit model and in the
numerical inverse Laplace transform.

Fitting from one datasheet point (eps', tan d at f_spec):
  * er_hf given: the loss band ends where eps' has fallen to er_hf (f2 solved);
  * er_hf = 0:   the loss band spans [f1, f_max] and eps_inf is solved.
"""

import math

import numpy as np
from scipy.optimize import brentq


class Dielectric:
    def __init__(self, eps_inf, k=0.0, m1=-2.0, m2=10.0):
        self.eps_inf, self.k, self.m1, self.m2 = eps_inf, k, m1, m2
        self.w1, self.w2 = 2 * math.pi * 10 ** m1, 2 * math.pi * 10 ** m2

    def eps(self, s):
        s = np.asarray(s, complex)
        if self.k == 0:
            return np.full(s.shape, self.eps_inf, complex)
        return self.eps_inf + self.k * np.log10((self.w2 + s) / (self.w1 + s))

    def at(self, f):
        """Complex permittivity eps' - j eps'' at frequency f."""
        return self.eps(2j * math.pi * np.asarray(f, float))

    def er(self, f):
        return self.at(f).real

    def tan_d(self, f):
        e = self.at(f)
        return -e.imag / e.real

    @property
    def f_band(self):
        return 10 ** self.m1, 10 ** self.m2

    @classmethod
    def fit(cls, er, tan_d, f_spec, er_hf=0.0, m1=-2.0, f_max=1e10):
        if tan_d <= 0:
            return cls(er)
        w0 = 2 * math.pi * f_spec

        def L(m2):
            return np.log10((2 * math.pi * 10 ** m2 + 1j * w0) / (2 * math.pi * 10 ** m1 + 1j * w0))

        if not er_hf:
            m2 = math.log10(f_max)
            k = -er * tan_d / L(m2).imag
            return cls(er - k * L(m2).real, k, m1, m2)

        def mismatch(m2):
            k = -er * tan_d / L(m2).imag
            return er_hf + k * L(m2).real - er

        lo, hi = math.log10(f_spec) - 1.0, 14.0
        if mismatch(lo) * mismatch(hi) > 0:
            raise ValueError(f"cannot fit eps'={er}, tan d={tan_d} at {f_spec} Hz with a high-frequency "
                             f"permittivity of {er_hf}: the loss is too small for that much dispersion")
        m2 = brentq(mismatch, lo, hi, xtol=1e-9)
        k = -er * tan_d / L(m2).imag
        return cls(er_hf, k, m1, m2)

    def describe(self, freqs=(50, 1e3, 1e4, 1e5, 1e6, 1e8)):
        return "  ".join(f"{_ff(f)}: {self.er(f):.2f}/{self.tan_d(f):.3f}" for f in freqs)


def _ff(f):
    for div, unit in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if f >= div:
            return f"{f/div:g}{unit}"
    return f"{f:g}"


def materials(cfg):
    """{'air': ..., 'resin': ..., 'FR4': ...} Dielectric models from the [materials] config."""
    mc = cfg["materials"]
    out = {"air": Dielectric(mc.get("air_er", 1.0))}
    for key, name in (("resin", "resin"), ("fr4", "FR4")):
        out[name] = Dielectric.fit(mc[f"{key}_er"], mc.get(f"{key}_tan_d", 0.0),
                                   mc.get(f"{key}_spec_freq", 1e6), mc.get(f"{key}_er_hf", 0.0),
                                   m1=math.log10(mc.get(f"{key}_loss_f_lo", 0.01)))
    return out


def solve_permittivities(cfg):
    """Real permittivity of each material at the field-solver frequency."""
    f = cfg["materials"].get("solve_freq", 1e6)
    return {name: float(m.er(f)) for name, m in materials(cfg).items()}
