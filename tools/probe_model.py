"""Linear circuit model of the probe and its analyses.

The circuit is built from
  * discrete parts read from the PCB (R, C, trimmer, pot),
  * the parasitic capacitance matrix from extract_caps.py,
  * capacitor ESL,
  * an optional series (back-termination) resistor,
  * a lossy coaxial cable (exact transmission-line two-port),
  * the oscilloscope input (R || C, optional bandwidth limit).

It is solved by nodal analysis at arbitrary complex frequencies s, which
gives the frequency response, the input impedance and (through a numerical
inverse Laplace transform) step responses.
"""

import json
import math
from dataclasses import dataclass, field

import numpy as np

import kicad_board
from kicad_board import parse_si

C0 = 299792458.0


# ----------------------------------------------------------------------------
# generic nodal-analysis network

class Network:
    def __init__(self):
        self.nodes = []          # node names, ground is "0"
        self.idx = {}
        self.elems = []          # (kind, a, b, params)

    def node(self, name):
        if name == "0":
            return -1
        if name not in self.idx:
            self.idx[name] = len(self.nodes)
            self.nodes.append(name)
        return self.idx[name]

    def add(self, kind, a, b, **p):
        self.node(a), self.node(b)
        self.elems.append((kind, a, b, p))

    # element admittances --------------------------------------------------
    @staticmethod
    def admittance(kind, p, s):
        if kind == "R":
            return np.full(s.shape, 1.0 / p["r"], complex)
        if kind == "C":
            z = 1.0 / (s * p["c"])
            if p.get("esl"):
                z = z + s * p["esl"]
            if p.get("esr"):
                z = z + p["esr"]
            return 1.0 / z
        raise ValueError(kind)

    @staticmethod
    def tline(p, s):
        """Series and shunt admittance of the pi-equivalent of a lossy line."""
        z0, vf, l = p["z0"], p["vf"], p["length"]
        Lp = z0 / (vf * C0)
        Cp = 1.0 / (z0 * vf * C0)
        w0 = 2 * math.pi * 1e8
        r_ac = 2 * z0 * p["loss_100MHz"] / 8.686
        K = r_ac / math.sqrt(w0 / 2)
        Zp = p["r_dc"] + K * np.sqrt(s) + s * Lp
        Yp = s * Cp
        gl = np.sqrt(Zp * Yp) * l
        zc = np.sqrt(Zp / Yp)
        return 1.0 / (zc * np.sinh(gl)), np.tanh(gl / 2) / zc

    def matrix(self, s):
        n = len(self.nodes)
        Y = np.zeros((len(s), n, n), complex)

        def stamp(a, b, y):
            ia, ib = self.idx.get(a, -1), self.idx.get(b, -1)
            if ia >= 0:
                Y[:, ia, ia] += y
            if ib >= 0:
                Y[:, ib, ib] += y
            if ia >= 0 and ib >= 0:
                Y[:, ia, ib] -= y
                Y[:, ib, ia] -= y

        for kind, a, b, p in self.elems:
            if kind == "TL":
                ys, yp = self.tline(p, s)
                stamp(a, b, ys)
                stamp(a, "0", yp)
                stamp(b, "0", yp)
            else:
                stamp(a, b, self.admittance(kind, p, s))
        return Y

    def solve(self, s, fixed, chunk=2048):
        """Node voltages for fixed node voltages `fixed` {name: V}.

        Returns (V [len(s), n], I_fixed {name: current into the network}).
        """
        s = np.atleast_1d(np.asarray(s, complex))
        fi = [self.idx[k] for k in fixed]
        fv = np.array([fixed[k] for k in fixed], complex)
        ui = [i for i in range(len(self.nodes)) if i not in fi]
        V = np.zeros((len(s), len(self.nodes)), complex)
        I = np.zeros((len(s), len(fi)), complex)
        for c0 in range(0, len(s), chunk):
            sl = slice(c0, c0 + chunk)
            Y = self.matrix(s[sl])
            Yuu = Y[:, ui][:, :, ui]
            Yuf = Y[:, ui][:, :, fi]
            rhs = -(Yuf @ fv)
            Vu = np.linalg.solve(Yuu, rhs[..., None])[..., 0]
            Vs = V[sl]          # view
            Vs[:, ui] = Vu
            Vs[:, fi] = fv
            Yfu = Y[:, fi][:, :, ui]
            Yff = Y[:, fi][:, :, fi]
            I[sl] = Yff @ fv + (Yfu @ Vu[..., None])[..., 0]
        return V, {k: I[:, n] for n, k in enumerate(fixed)}


# ----------------------------------------------------------------------------
# probe model

@dataclass
class Params:
    cable_type: str = "RG58"
    length: float = 1.5
    series_r: float = 0.0
    r_scope: float = 1e6
    c_scope: float = 15e-12
    scope_bw: float = 0.0
    source_r: float = 0.0          # DUT source impedance
    tip_r: float = 0.0             # damping resistor in series with the HV tip
    seg_r: float = 0.0             # damping resistor in series with every divider capacitor
    trim: float = None             # trimmer capacitance [F], None = config
    extra_c: float = 0.0           # extra ideal capacitance OUT-GND [F]
    top_caps: list = None          # per-segment discrete capacitance [F] (None = PCB values)
    pot: float = None              # pot resistance [ohm], None = config / auto
    parasitics: bool = True
    esl: bool = True


def params_from_config(cfg, **over):
    cab, sc = cfg["cable"], cfg["scope"]
    p = Params(cable_type=cab["type"], length=cab["length"], series_r=cab["series_r"],
               r_scope=sc["r_in"], c_scope=sc["c_in"], scope_bw=sc.get("bandwidth", 0.0),
               source_r=cfg["circuit"].get("source_resistance", 0.0),
               tip_r=cfg["circuit"].get("tip_r", 0.0), seg_r=cfg["circuit"].get("segment_r", 0.0))
    for k, v in over.items():
        if v is not None:
            setattr(p, k, v)
    return p


class Probe:
    """Probe circuit built from the PCB + parasitics + parameters."""

    def __init__(self, cfg, board, capdata, params):
        self.cfg, self.board, self.capdata, self.p = cfg, board, capdata, params
        cc = cfg["circuit"]
        self.hv, self.out = cc["input_net"], cc["output_net"]
        self.labels = kicad_board.net_labels(board, self.hv, self.out)
        self.ladder = kicad_board.ladder_order(board, self.hv, self.out)
        self.nseg = len(self.ladder) - 1
        self._classify_parts()
        if params.pot is None:
            setting = cc.get("parts", {}).get(self.pot_ref, "auto") if self.pot_ref else None
            if setting == "auto" or setting is None:
                self.p.pot = self.solve_pot(cc["target_ratio"]) if self.pot_ref else 0.0
            else:
                self.p.pot = parse_si(str(setting))
        self.net = self._build()

    # --- parts ------------------------------------------------------------
    def _classify_parts(self):
        cc = self.cfg["circuit"]
        seg_of = {frozenset(pair): k for k, pair in enumerate(zip(self.ladder, self.ladder[1:]))}
        self.top = {}        # segment -> [cap refs]
        self.top_r = {}      # segment -> [res refs]
        self.bottom_c, self.bottom_r = [], []
        self.trim_ref, self.trim_range, self.pot_ref = None, None, None
        for p in self.board.parts.values():
            nets = frozenset(p.nets)
            if p.ref.startswith("RV"):
                self.pot_ref = p.ref
                continue
            if "-" in p.value and p.ref.startswith("C"):
                lo, hi = (parse_si(v) for v in p.value.split("-"))
                self.trim_ref, self.trim_range = p.ref, (lo, hi)
                continue
            if nets in seg_of:
                d = self.top if p.ref.startswith("C") else self.top_r
                d.setdefault(seg_of[nets], []).append(p.ref)
            elif p.ref.startswith("C"):
                self.bottom_c.append(p.ref)
            else:
                self.bottom_r.append(p.ref)
        over = cc.get("parts", {})
        if self.trim_ref and self.p.trim is None:
            t = over.get(self.trim_ref, "mid")
            self.p.trim = sum(self.trim_range) / 2 if t == "mid" else parse_si(str(t))

    def part_value(self, ref):
        over = self.cfg["circuit"].get("parts", {})
        if ref in over and ref not in (self.trim_ref, self.pot_ref):
            return parse_si(str(over[ref]))
        return parse_si(self.board.parts[ref].value)

    def esl_of(self, ref):
        if not self.p.esl:
            return 0.0
        fp = self.board.parts[ref].footprint
        for k, v in self.cfg["circuit"].get("esl", {}).items():
            if k in fp:
                return float(v)
        return 0.0

    def top_cap_values(self):
        """Discrete capacitance per segment [F] (sum of parallel parts)."""
        if self.p.top_caps is not None:
            v = self.p.top_caps
            return list(v) if np.ndim(v) else [float(v)] * self.nseg
        return [sum(self.part_value(r) for r in self.top.get(k, [])) for k in range(self.nseg)]

    def bottom_fixed_c(self):
        return sum(self.part_value(r) for r in self.bottom_c)

    # --- network ----------------------------------------------------------
    def _build(self, dc_only=False):
        p, L = self.p, self.labels
        n = Network()
        g = lambda net: "0" if net == "GND" else L[net]

        # divider segments
        topc = self.top_cap_values()
        for k, (a, b) in enumerate(zip(self.ladder, self.ladder[1:])):
            for r in self.top_r.get(k, []):
                n.add("R", g(a), g(b), r=self.part_value(r), ref=r)
            refs = self.top.get(k, [])
            esl = self.esl_of(refs[0]) if refs else 0.0
            if topc[k] > 0:
                n.add("C", g(a), g(b), c=topc[k], esl=esl, esr=p.seg_r, ref="+".join(refs) or f"Cseg{k+1}")
        # bottom network
        for r in self.bottom_r:
            a, b = self.board.parts[r].nets
            n.add("R", g(a), g(b), r=self.part_value(r), ref=r)
        for c in self.bottom_c:
            a, b = self.board.parts[c].nets
            n.add("C", g(a), g(b), c=self.part_value(c), esl=self.esl_of(c), ref=c)
        if self.trim_ref and p.trim:
            a, b = self.board.parts[self.trim_ref].nets
            n.add("C", g(a), g(b), c=p.trim, esl=self.esl_of(self.trim_ref), ref=self.trim_ref)
        if p.extra_c:
            # behaves like one more SMD capacitor (same ESL) when positive
            esl = self.cfg["circuit"].get("extra_esl", 1.0e-9) if p.extra_c > 0 and p.esl else 0.0
            n.add("C", "OUT", "0", c=p.extra_c, esl=esl, ref="Cextra")
        if self.pot_ref:
            a, b = self.board.parts[self.pot_ref].nets
            n.add("R", g(a), g(b), r=max(p.pot, 1e-3), ref=self.pot_ref)

        # parasitics from the field solver
        if p.parasitics and self.capdata is not None:
            names = self.capdata["names"]
            M = np.array(self.capdata["maxwell"])
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    c = -M[i, j]
                    if c > 1e-17 and names[i] in L and names[j] in L:
                        n.add("C", g(names[i]), g(names[j]), c=c, ref=f"Cp_{L[names[i]]}_{L[names[j]]}")

        # source
        n.node("HV")
        if p.source_r + p.tip_r > 0:
            n.add("R", "SRC", "HV", r=p.source_r + p.tip_r, ref="Rsrc_tip")
            self.fixed = {"SRC": 1.0}
        else:
            self.fixed = {"HV": 1.0}

        # output: series R -> cable -> scope
        node = "OUT"
        if p.series_r > 0:
            n.add("R", "OUT", "CBL", r=p.series_r, ref="Rseries")
            node = "CBL"
        if p.length > 0:
            cab = dict(self.cfg["cable"]["types"][p.cable_type])
            cab["length"] = p.length
            n.add("TL", node, "SCOPE", **cab)
        else:
            n.add("R", node, "SCOPE", r=1e-3)
        n.add("R", "SCOPE", "0", r=p.r_scope, ref="Rscope")
        n.add("C", "SCOPE", "0", c=p.c_scope, ref="Cscope")
        return n

    # --- analyses ---------------------------------------------------------
    def solve(self, s):
        return self.net.solve(s, self.fixed)

    def H(self, f=None, s=None):
        """Transfer V(scope)/V(source)."""
        if s is None:
            s = 2j * np.pi * np.asarray(f, float)
        V, _ = self.solve(s)
        h = V[:, self.net.idx["SCOPE"]]
        if self.p.scope_bw > 0:
            h = h / (1 + s / (2 * np.pi * self.p.scope_bw))
        return h

    def H0(self):
        return self.H(f=[1e-3])[0].real

    def ratio(self):
        return 1.0 / self.H0()

    def input_admittance(self, f):
        s = 2j * np.pi * np.asarray(f, float)
        V, I = self.solve(s)
        k = list(self.fixed)[0]
        return I[k] / self.fixed[k]

    def node_voltages(self, f):
        """Phasor voltage of each ladder net at frequency f (source = 1 V)."""
        V, _ = self.solve(np.array([2j * np.pi * f]))
        out = {}
        for net in self.ladder:
            name = self.labels[net]
            out[net] = V[0, self.net.idx[name]] if name in self.net.idx else 0.0
        return out

    def solve_pot(self, target):
        """Pot resistance giving the target DC ratio (bisection)."""
        lo, hi = 0.0, self.cfg["circuit"].get("pot_max", 1e6)
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            self.p.pot = mid
            self.net = self._build()
            r = self.ratio()
            if r > target:
                lo = mid   # more resistance at the bottom -> smaller ratio
            else:
                hi = mid
        if hi >= self.cfg["circuit"].get("pot_max", 1e6) - 1:
            print(f"WARNING: target ratio {target} not reachable with the pot")
        return 0.5 * (lo + hi)

    def step(self, t_max, n=2**14, rise=0.0):
        """Unit step response (normalised to the DC gain) on [0, t_max].

        Numerical Laplace transform with damping and a Hanning window
        (Wedepohl/Wilcox method), which handles both the fast transmission
        line effects and the slow RC tails without wrap-around.

        rise > 0 gives the source a Gaussian edge with that 10-90 % rise time,
        centred at t0 = 2*rise. This also keeps out the GHz content where the
        lumped model is not meaningful.
        """
        T = 2.0 * t_max
        dw = 2 * np.pi / T
        c = 2 * math.log(n) / T
        w = np.arange(n) * dw
        s = c + 1j * w
        F = self.H(s=s) / s
        if rise > 0:
            sigma = rise / 2.563
            F = F * np.exp(0.5 * (sigma * s) ** 2 - s * 2 * rise)
        win = 0.5 * (1 + np.cos(np.pi * np.arange(n) / n))
        a = F * win
        a[0] *= 0.5
        t = np.arange(n) * (T / n)
        f = (2.0 / T) * np.exp(c * t) * np.real(n * np.fft.ifft(a))
        keep = t <= t_max
        return t[keep], f[keep] / self.H0()


# ----------------------------------------------------------------------------
# metrics

F_VALID = 1e9   # the lumped quasi-static model is not meaningful far beyond this


def response_metrics(probe, f=None):
    if f is None:
        f = np.logspace(0, np.log10(F_VALID), 1200)
    h = probe.H(f) / probe.H0()
    mag = np.abs(h)
    db = 20 * np.log10(mag)
    def first_below(level):
        i = np.nonzero(mag < level)[0]
        return f[i[0]] if len(i) else np.inf

    def first_outside(tol_db):
        i = np.nonzero(np.abs(db) > tol_db)[0]
        return f[i[0]] if len(i) else np.inf
    f3 = first_below(10 ** (-3 / 20))
    band = f < f3
    lf = (f >= probe.cfg["optimize"]["f_lo"]) & (f <= probe.cfg["optimize"]["f_hi"])
    Yin_1M = probe.input_admittance([1e6])[0]
    Yin_10k = probe.input_admittance([1e4])[0]
    return {
        "ratio_dc": probe.ratio(),
        "f_3dB": f3,
        "f_1dB": first_outside(1.0),
        "peaking_dB": float(db[band].max()) if band.any() else 0.0,
        "lf_error_pct": float(np.max(np.abs(mag[lf] - 1)) * 100),
        "hf_peak_dB": float(db[f >= 1e7].max()),
        "c_in_10k_pF": float(Yin_10k.imag / (2 * np.pi * 1e4) * 1e12),
        "c_in_1M_pF": float(Yin_1M.imag / (2 * np.pi * 1e6) * 1e12),
    }


def step_metrics(probe, t_max=1e-6, rise=None):
    """Rise time and overshoot for a source step with a finite rise time."""
    if rise is None:
        rise = probe.cfg.get("analysis", {}).get("source_rise", 1e-9)
    t, y = probe.step(t_max, n=2**14, rise=rise)
    final = 1.0
    i10 = np.argmax(y >= 0.1 * final)
    i90 = np.argmax(y >= 0.9 * final)
    tr = t[i90] - t[i10] if y.max() >= 0.9 else np.inf
    return {"rise_ns": tr * 1e9, "overshoot_pct": max(0.0, (y.max() - final) * 100)}


# ----------------------------------------------------------------------------
# compensation

def compensation_error(probe, f_lo, f_hi, n=60):
    f = np.logspace(math.log10(f_lo), math.log10(f_hi), n)
    h = probe.H(f) / probe.H0()
    return float(np.sqrt(np.mean((np.abs(h) - 1) ** 2)))


def optimize_bottom(probe, f_lo=None, f_hi=None, lo=-2e-9, hi=20e-9):
    """Find the extra OUT-GND capacitance (added to the current trimmer
    setting) that makes the response flattest over [f_lo, f_hi].
    Returns the optimum total trimmer-equivalent capacitance [F]."""
    from scipy.optimize import minimize_scalar
    o = probe.cfg["optimize"]
    f_lo = f_lo or o["f_lo"]
    f_hi = f_hi or o["f_hi"]
    base = probe.p.extra_c
    fixed_bottom = probe.bottom_fixed_c() + (probe.p.trim or 0)
    lo = max(lo, -fixed_bottom + 1e-12)

    def err(c):
        probe.p.extra_c = c
        probe.net = probe._build()
        return compensation_error(probe, f_lo, f_hi)

    r = minimize_scalar(err, bounds=(lo, hi), method="bounded", options={"xatol": 1e-14})
    probe.p.extra_c = r.x
    probe.net = probe._build()
    return base + r.x, r.fun


E_SERIES = {
    "E6": [1.0, 1.5, 2.2, 3.3, 4.7, 6.8],
    "E12": [1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2],
    "E24": [1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9,
            4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1],
}


def series_values(name, lo=1e-12, hi=1e-7):
    out = []
    dec = 1e-12
    while dec < hi:
        out += [v * dec for v in E_SERIES[name] if lo <= v * dec <= hi]
        dec *= 10
    return out


def suggest_fixed(target, series="E12", max_parts=2):
    """Combinations of up to max_parts standard values summing close to target."""
    vals = series_values(series, 1e-12, max(target * 1.2, 2e-12))
    best = []
    for a in vals:
        best.append(((a,), a))
    if max_parts >= 2:
        for i, a in enumerate(vals):
            for b in vals[i:]:
                best.append(((b, a), a + b))
    best.sort(key=lambda t: abs(t[1] - target))
    return best[:5]


def graded_caps(probe, c_bottom_segment):
    """Discrete segment capacitances that make the capacitive (HF) voltage
    distribution along the ladder linear, given the parasitic matrix.

    Works from the bottom segment (fixed to c_bottom_segment) upward: the
    current in segment k must supply segment k+1 plus all parasitic currents
    leaving node k when the node voltages are exactly linear.
    """
    names = probe.capdata["names"]
    M = np.array(probe.capdata["maxwell"])
    ix = {n: i for i, n in enumerate(names)}
    lad = probe.ladder
    nseg = len(lad) - 1
    V = {net: 1 - k / nseg for k, net in enumerate(lad)}

    def cp(a, b):
        return -M[ix[a], ix[b]] if a in ix and b in ix else 0.0

    seg_par = [cp(a, b) for a, b in zip(lad, lad[1:])]
    total = [0.0] * nseg
    total[-1] = c_bottom_segment + seg_par[-1]
    for k in range(nseg - 1, 0, -1):     # node lad[k] sits between segment k-1 and k
        node = lad[k]
        leak = 0.0
        for other in names:
            if other in (node, lad[k - 1], lad[k + 1]):
                continue
            leak += cp(node, other) * (V[node] - V.get(other, 0.0))
        total[k - 1] = total[k] + nseg * leak
    return [t - sp for t, sp in zip(total, seg_par)]


def load_capdata(path):
    with open(path) as f:
        return json.load(f)
