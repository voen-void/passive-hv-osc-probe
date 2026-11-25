"""Causal wideband dielectric models.

Loss and dispersion are linked (Kramers-Kronig): a material with a high loss
tangent has a permittivity that falls with frequency. Both models here are
analytic functions of s, so they go straight into the circuit model and into
the numerical inverse Laplace transform.

Djordjevic-Sarkar (class Dielectric), fitted to one datasheet point:

    eps(s) = eps_inf + k * log10((w2 + s) / (w1 + s)),   w = 2*pi*10^m

  Between f1 = 10^m1 and f2 = 10^m2 the loss tangent is nearly constant and
  eps' falls by k per decade; below f1 and above f2 the material is lossless.
    * er_hf given: the loss band ends where eps' has fallen to er_hf;
    * er_hf = 0:   the loss band spans [f1, f_max] and eps_inf is solved.

Cole-Cole (class ColeCole), fitted to a measured sweep of a sample capacitor:

    eps(s) = eps_inf + d_eps / (1 + (s*tau)^beta)

  One broad relaxation: beta = 1 is a single Debye pole, small beta spreads
  the loss over many decades (the cast polyurethane resin has beta ~ 0.24).
  See measured_model() for how the sweep is turned into a permittivity.
"""

import csv
import functools
import math
import os

import numpy as np
from scipy.optimize import brentq, least_squares

EPS0 = 8.8541878128e-12


class _Model:
    """Shared evaluation helpers; subclasses provide eps(s)."""

    dispersive = True
    note = ""

    def at(self, f):
        """Complex permittivity eps' - j eps'' at frequency f."""
        return self.eps(2j * math.pi * np.asarray(f, float))

    def er(self, f):
        return self.at(f).real

    def tan_d(self, f):
        e = self.at(f)
        return -e.imag / e.real

    def describe(self, freqs=(50, 1e3, 1e4, 1e5, 1e6, 1e8)):
        return "  ".join(f"{_ff(f)}: {self.er(f):.2f}/{self.tan_d(f):.3f}" for f in freqs)


class Dielectric(_Model):
    def __init__(self, eps_inf, k=0.0, m1=-2.0, m2=10.0):
        self.eps_inf, self.k, self.m1, self.m2 = eps_inf, k, m1, m2
        self.w1, self.w2 = 2 * math.pi * 10 ** m1, 2 * math.pi * 10 ** m2
        self.dispersive = bool(k)

    def eps(self, s):
        s = np.asarray(s, complex)
        if self.k == 0:
            return np.full(s.shape, self.eps_inf, complex)
        return self.eps_inf + self.k * np.log10((self.w2 + s) / (self.w1 + s))

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


class ColeCole(_Model):
    def __init__(self, eps_inf, d_eps, tau, beta):
        self.eps_inf, self.d_eps, self.tau, self.beta = eps_inf, d_eps, tau, beta

    def eps(self, s):
        s = np.asarray(s, complex)
        return self.eps_inf + self.d_eps / (1 + (s * self.tau) ** self.beta)

    @property
    def f_peak(self):
        return 1 / (2 * math.pi * self.tau)


def _ff(f):
    for div, unit in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if f >= div:
            return f"{f/div:g}{unit}"
    return f"{f:g}"


# ----------------------------------------------------------------------------
# model from a measured sample

def read_sweep(path):
    """Corrected sweep of a sample capacitor, as written by resin_fit.py:
    columns f_Hz, C_pF, tan_d, sigma_C_pF, sigma_tan_d, use. Lines starting
    with # are comments. Returns the rows with use = 1 as arrays."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(l for l in fh if not l.startswith("#"))]
    rows = [r for r in rows if r.get("use", "1").strip() not in ("0", "")]
    if len(rows) < 6:
        raise ValueError(f"{path}: fewer than 6 usable points")
    col = lambda k, d=None: np.array([float(r[k]) if r.get(k) not in (None, "") else d for r in rows], float)
    f, c, d = col("f_Hz"), col("C_pF"), col("tan_d")
    return f, c, d, col("sigma_C_pF", 0.0), col("sigma_tan_d", 0.0)


def without_gap(c, d, gap):
    """Take a lossless series capacitance `gap` [pF] (air between electrode
    and sample) off a measured C [pF] and loss factor."""
    if not gap or not math.isfinite(gap):
        return c, d
    y = 1 / (1 / (c * (1 - 1j * d)) - 1 / gap)
    return y.real, -y.imag / y.real


def fit_cole_cole(f, c, d, sc=None, sd=None):
    """Weighted least-squares Cole-Cole fit to C(f) [pF] and tan d(f).
    Returns ((c_inf, d_c, tau, beta), reduced chi-square)."""
    w = 2j * np.pi * np.asarray(f, float)
    sc = np.maximum(sc if sc is not None else 0.0, 1e-3 * c)
    sd = np.maximum(sd if sd is not None else 0.0, 2e-4)

    def model(p):
        return p[0] + p[1] / (1 + (w * 10 ** p[2]) ** p[3])

    def resid(p):
        y = model(p)
        return np.concatenate(((y.real - c) / sc, (-y.imag / y.real - d) / sd))

    c_hi, span = float(np.min(c)), float(np.max(c) - np.min(c))
    best = None
    for lt in (-1.0, -3.0, -5.0):                # a few starts: the loss peak may lie outside the data
        r = least_squares(resid, [0.9 * c_hi, 2 * span + 1e-3, lt, 0.4],
                          bounds=([0.0, 0.0, -9.0, 0.03], [c_hi, 50 * (span + c_hi), 4.0, 1.0]), x_scale="jac")
        if best is None or r.cost < best.cost:
            best = r
    p = best.x
    chi2 = 2 * best.cost / max(1, 2 * len(c) - 4)
    return (float(p[0]), float(p[1]), float(10 ** p[2]), float(p[3])), float(chi2)


def measured_model(path, er, f_spec, tan_d=None, gap="none", area=None, thickness=None):
    """Cole-Cole permittivity from a sweep of a sample capacitor.

    The sweep gives C(f) and tan d(f) of the sample between its electrodes.
    Their shape over frequency is the material's; two things it does not give
    are settled by the arguments:

      scale   the factor between capacitance and permittivity.
              area [mm^2] and thickness [mm] of the sample: eps = C d / (eps0 A).
              Without them the datasheet point is used: eps' = er at f_spec.
              A thickness given as (min, max) is an uneven sample: too
              uncertain to set the scale, so the datasheet point still does,
              and the geometry only serves as a check (.fit["t_eff_mm"], the
              flat-sample thickness the scale stands for, against the range).
      gap     air between electrode and sample is a lossless capacitance in
              series; it lowers C, tan d and the dispersion by the same factor.
              "none"       the electrodes touch the sample everywhere;
              a number     the gap capacitance [pF];
              "datasheet"  the gap that brings the fitted tan d at f_spec up
                           to the datasheet's `tan_d` (none if it is already
                           there).

    Returns a ColeCole model; .fit holds the fit's details."""
    f, c, d, sc, sd = read_sweep(path)

    def fitted(gap_pf):
        cs, ds = without_gap(c, d, gap_pf)
        k = cs / c                                   # the sigmas scale with the values
        p, chi2 = fit_cole_cole(f, cs, ds, sc * k, sd * k)
        return p, chi2

    def tan_at_spec(p):
        y = p[0] + p[1] / (1 + (2j * math.pi * f_spec * p[2]) ** p[3])
        return -y.imag / y.real

    gap_pf = None
    if isinstance(gap, str) and gap.lower() == "datasheet":
        if not tan_d:
            raise ValueError("gap = 'datasheet' needs the datasheet loss factor")
        if tan_at_spec(fitted(None)[0]) < tan_d:
            lo, hi = math.log(1.02 * float(np.max(c))), math.log(1e6)
            g = lambda lg: tan_at_spec(fitted(math.exp(lg))[0]) - tan_d
            if g(lo) > 0:
                gap_pf = math.exp(brentq(g, lo, hi, xtol=1e-4))
    elif not isinstance(gap, str) and gap:
        gap_pf = float(gap)

    t_range = None
    if isinstance(thickness, (list, tuple)):
        t_range = (float(min(thickness)), float(max(thickness)))
        thickness = None

    (c_inf, d_c, tau, beta), chi2 = fitted(gap_pf)
    c_spec = (c_inf + d_c / (1 + (2j * math.pi * f_spec * tau) ** beta)).real
    if area and thickness:
        scale = thickness * 1e-3 / (EPS0 * area * 1e-6) * 1e-12       # eps per pF
        how = f"sample geometry {area:g} mm2 x {thickness:g} mm"
    else:
        scale = er / c_spec
        how = f"eps' = {er:g} at {_ff(f_spec)}Hz (datasheet)"
    m = ColeCole(c_inf * scale, d_c * scale, tau, beta)
    cs, ds = without_gap(c, d, gap_pf)
    y = m.at(f)
    m.fit = {
        "path": path, "points": int(len(f)), "f_min": float(f.min()), "f_max": float(f.max()),
        "chi2": chi2, "gap_pF": gap_pf, "scale_from": how, "pF_per_eps": 1 / scale,
        "rms_c_pct": float(100 * np.sqrt(np.mean((y.real / (cs * scale) - 1) ** 2))),
        "rms_tan_d": float(np.sqrt(np.mean((-y.imag / y.real - ds) ** 2))),
    }
    if area:
        # C = eps0 eps A / t: the thickness a flat sample of this electrode area
        # would have at this scale. For an uneven one it is A / integral(dA / t),
        # which thin regions dominate; fringing at the electrode edge is left out.
        m.fit["area_mm2"] = float(area)
        m.fit["t_eff_mm"] = EPS0 * area * 1e-6 * scale * 1e12 * 1e3
        if gap_pf:
            m.fit["gap_um"] = EPS0 * area * 1e-6 / (gap_pf * 1e-12) * 1e6      # mean air gap of that capacitance
        if t_range:
            lo, hi = t_range
            m.fit["t_range_mm"] = t_range
            # effective thickness if every thickness between lo and hi takes an equal share of the area
            m.fit["t_even_mm"] = (hi - lo) / math.log(hi / lo) if hi > lo > 0 else lo
    m.note = (f"Cole-Cole fit to {os.path.basename(path)} ({len(f)} points, {_ff(f.min())}Hz-{_ff(f.max())}Hz): "
              f"eps_inf {m.eps_inf:.2f}, d_eps {m.d_eps:.2f}, loss peak {m.f_peak:.3g} Hz, beta {beta:.3f}; "
              f"scale from {how}; electrode gap " + (f"{gap_pf:.0f} pF in series" if gap_pf else "none")
              + f". Outside {_ff(f.min())}Hz-{_ff(f.max())}Hz it is an extrapolation.")
    return m


@functools.lru_cache(maxsize=8)
def _measured_cached(path, mtime, er, f_spec, tan_d, gap, area, thickness):
    return measured_model(path, er, f_spec, tan_d, gap, area, thickness)


def materials(cfg):
    """{'air': ..., 'resin': ..., 'FR4': ...} models from the [materials] config."""
    mc = cfg["materials"]
    out = {"air": Dielectric(mc.get("air_er", 1.0))}
    for key, name in (("resin", "resin"), ("fr4", "FR4")):
        data = mc.get(f"{key}_data")
        if data:
            base = cfg.get("_base") or os.path.dirname(os.path.abspath(__file__))
            path = os.path.normpath(os.path.join(base, data))
            thick = mc.get(f"{key}_sample_thickness", 0.0)
            out[name] = _measured_cached(path, os.path.getmtime(path), mc[f"{key}_er"],
                                         mc.get(f"{key}_spec_freq", 50.0), mc.get(f"{key}_tan_d", 0.0),
                                         mc.get(f"{key}_electrode_gap", "none"),
                                         mc.get(f"{key}_sample_area", 0.0),
                                         tuple(thick) if isinstance(thick, list) else thick)
        else:
            out[name] = Dielectric.fit(mc[f"{key}_er"], mc.get(f"{key}_tan_d", 0.0),
                                       mc.get(f"{key}_spec_freq", 1e6), mc.get(f"{key}_er_hf", 0.0),
                                       m1=math.log10(mc.get(f"{key}_loss_f_lo", 0.01)))
    return out


def solve_permittivities(cfg):
    """Real permittivity of each material at the field-solver frequency
    ([materials] solve_er = {resin = ...} overrides single materials)."""
    f = cfg["materials"].get("solve_freq", 1e6)
    er = {name: float(m.er(f)) for name, m in materials(cfg).items()}
    er.update({k: float(v) for k, v in cfg["materials"].get("solve_er", {}).items()})
    return er
