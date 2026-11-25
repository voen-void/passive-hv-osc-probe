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

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")   # batched small solves: ~80x faster single-threaded
import json
import math
from dataclasses import dataclass, field

import numpy as np

import dielectric
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
        self.materials = {}      # material -> function(s) = eps(s) / eps at the field-solve frequency

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
            if p.get("l"):                       # resistor with series inductance
                return 1.0 / (max(p["r"], 1e-6) + s * p["l"])
            return np.full(s.shape, 1.0 / max(p["r"], 1e-6), complex)
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
            if kind == "CMAT":
                # nodal capacitance matrices of the parasitics, one per material;
                # each follows that material's (dispersive, lossy) permittivity
                ix = [self.idx[nm] for nm in p["nodes"]]
                Yp = 0
                for m, Cm in p["mats"].items():
                    r = self.materials[m](s) if m in self.materials else 1.0
                    Yp = Yp + (s * r)[:, None, None] * Cm[None, :, :]
                Y[np.ix_(range(len(s)), ix, ix)] += Yp
            elif kind == "TL":
                ys, yp = self.tline(p, s)
                stamp(a, b, ys)
                stamp(a, "0", yp)
                stamp(b, "0", yp)
            else:
                stamp(a, b, self.admittance(kind, p, s))
        return Y

    def solve(self, s, fixed, chunk=2048, inject=None):
        """Node voltages for fixed node voltages `fixed` {name: V}, plus
        optional current sources `inject` {name: current into that node}.

        Returns (V [len(s), n], I_fixed {name: current into the network}).
        """
        s = np.atleast_1d(np.asarray(s, complex))
        fi = [self.idx[k] for k in fixed]
        fv = np.array([fixed[k] for k in fixed], complex)
        ui = [i for i in range(len(self.nodes)) if i not in fi]
        pos = {i: p for p, i in enumerate(ui)}
        V = np.zeros((len(s), len(self.nodes)), complex)
        I = np.zeros((len(s), len(fi)), complex)
        for c0 in range(0, len(s), chunk):
            sl = slice(c0, c0 + chunk)
            Y = self.matrix(s[sl])
            Yuu = Y[:, ui][:, :, ui]
            Yuf = Y[:, ui][:, :, fi]
            rhs = -(Yuf @ fv)
            for name, cur in (inject or {}).items():
                rhs[:, pos[self.idx[name]]] += cur
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
    series_r: float = None         # back-termination [ohm]; None = as on the PCB (R29)
    r_scope: float = 1e6
    c_scope: float = 15e-12
    scope_bw: float = 0.0
    source_r: float = 0.0          # DUT source impedance
    tip_r: float = 0.0             # damping resistor in series with the HV tip
    seg_r: float = None            # damping resistor in series with every divider capacitor [ohm];
                                   # None = the resistors on the PCB (none there: no damping)
    trim: float = None             # trimmer capacitance [F], None = config
    extra_c: float = 0.0           # extra ideal capacitance OUT-GND [F]
    top_caps: list = None          # per-segment discrete capacitance [F] (None = PCB values)
    pot: float = None              # pot resistance [ohm], None = config / auto
    parasitics: bool = True
    esl: bool = True
    esr: bool = True
    dielectric_loss: bool = True


def params_from_config(cfg, **over):
    cab, sc = cfg["cable"], cfg["scope"]
    sr = cab.get("series_r", "pcb")
    seg = cfg["circuit"].get("segment_r", "pcb")
    p = Params(cable_type=cab["type"], length=cab["length"], series_r=None if sr == "pcb" else float(sr),
               r_scope=sc["r_in"], c_scope=sc["c_in"], scope_bw=sc.get("bandwidth", 0.0),
               source_r=cfg["circuit"].get("source_resistance", 0.0),
               tip_r=cfg["circuit"].get("tip_r", 0.0), seg_r=None if seg == "pcb" else float(seg))
    for k, v in over.items():
        if v is not None:
            setattr(p, k, v)
    return p


class Probe:
    """Probe circuit built from the PCB + parasitics + parameters."""

    def __init__(self, cfg, board, capdata, params):
        self.cfg, self.board, self.capdata, self.p = cfg, board, capdata, params
        cc = cfg["circuit"]
        self.hv = kicad_board.resolve_net(board, cc["input_net"])
        self.out = kicad_board.resolve_net(board, cc["output_net"])
        self.labels = kicad_board.net_labels(board, self.hv, self.out)
        self.ladder = kicad_board.ladder_order(board, self.hv)
        # damping resistors on the board, in series with the divider capacitors
        self.damping, self.merged = kicad_board.segment_damping(board, self.ladder)
        self.bottom = self.ladder[-1]
        self.bot = self.labels[self.bottom]      # node label of the divider bottom
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
        self.seg_res = {}    # segment -> [damping resistor refs, in series with the capacitor]
        self.bottom_c, self.bottom_r, self.other = [], [], []
        self.trim_ref, self.trim_range, self.pot_ref, self.term_ref = None, None, None, None
        damp_r = {r: k for k, brs in self.damping.items() for r, _, _ in brs}
        damp_c = {c: k for k, brs in self.damping.items() for _, c, _ in brs}
        for p in self.board.parts.values():
            nets = frozenset(p.nets)
            if p.ref.startswith("RV"):
                self.pot_ref = p.ref
                continue
            if "-" in p.value and p.ref.startswith("C"):
                lo, hi = (parse_si(v) for v in p.value.split("-"))
                self.trim_ref, self.trim_range = p.ref, (lo, hi)
                continue
            if p.ref in damp_c:
                self.top.setdefault(damp_c[p.ref], []).append(p.ref)
            elif p.ref in damp_r:
                self.seg_res.setdefault(damp_r[p.ref], []).append(p.ref)
            elif nets in seg_of:
                d = self.top if p.ref.startswith("C") else self.top_r
                d.setdefault(seg_of[nets], []).append(p.ref)
            elif p.ref.startswith("C"):
                self.bottom_c.append(p.ref)
            elif p.ref.startswith("R") and len(p.nets) == 2:
                self.bottom_r.append(p.ref)
                if nets == frozenset((self.bottom, self.out)):
                    self.term_ref = p.ref       # back-termination resistor (R29)
            else:
                self.other.append(p.ref)        # diodes etc.: modelled by their capacitance
        self.divider_refs = {r for d in (self.top, self.top_r, self.seg_res) for refs in d.values() for r in refs}
        over = cc.get("parts", {})
        if self.trim_ref and self.p.trim is None:
            t = over.get(self.trim_ref, "mid")
            self.p.trim = sum(self.trim_range) / 2 if t == "mid" else parse_si(str(t))

    def part_value(self, ref):
        if ref == self.term_ref and self.p.series_r is not None:
            return self.p.series_r
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

    def esr_of(self, ref):
        if not self.p.esr:
            return 0.0
        fp = self.board.parts[ref].footprint
        for k, v in self.cfg["circuit"].get("esr", {}).items():
            if k in fp:
                return float(v)
        return 0.0

    def other_capacitance(self, ref):
        v = self.cfg["circuit"].get("part_capacitance", {}).get(ref)
        return parse_si(str(v)) if v is not None else None

    def series_r_value(self):
        if self.term_ref:
            return self.part_value(self.term_ref)
        return self.p.series_r or 0.0

    def dielectric_models(self):
        """(models, eps used by the field solve) if dispersion/loss is active, else None."""
        if not self.p.dielectric_loss or "partial" not in (self.capdata or {}):
            return None
        models = dielectric.materials(self.cfg)
        meta = self.capdata.get("meta", {})
        eps_solve = meta.get("eps_solve")
        if eps_solve is None and "resin_er" in meta:     # files from older versions
            eps_solve = {"air": 1.0, "resin": meta["resin_er"], "FR4": meta["fr4_er"]}
        eps_solve = eps_solve or dielectric.solve_permittivities(self.cfg)
        return models, eps_solve

    def top_cap_values(self):
        """Discrete capacitance per segment [F] (sum of parallel parts)."""
        if self.p.top_caps is not None:
            v = self.p.top_caps
            return list(v) if np.ndim(v) else [float(v)] * self.nseg
        return [sum(self.part_value(r) for r in self.top.get(k, [])) for k in range(self.nseg)]

    def bottom_fixed_c(self):
        return sum(self.part_value(r) for r in self.bottom_c)

    def seg_r_values(self):
        """Damping resistance in series with each segment's capacitor [ohm]:
        the board's resistors, or p.seg_r for every segment if that is set.
        Several damped branches in one segment are taken in parallel."""
        if self.p.seg_r is not None:
            return [float(self.p.seg_r)] * self.nseg
        out = []
        for k in range(self.nseg):
            g = sum(1.0 / max(self.part_value(r), 1e-6) for r in self.seg_res.get(k, []))
            out.append(1.0 / g if g else 0.0)
        return out

    def _node_x(self, net):
        """Mean x [mm] of the SMD pads of a ladder net (with the nets merged into it)."""
        nets = {net} | {m for m, to in self.merged.items() if to == net}
        xs = [x for ref, fp, layer, n, poly in self.board.pad_shapes if n in nets and ref in self.divider_refs
              for x, _ in poly]
        return float(np.mean(xs)) if xs else None

    def seg_l_values(self):
        """Series inductance of each divider segment [H], or None to use the
        capacitors' ESL from [circuit.esl].

        [circuit] segment_l = "auto": the segment current runs one pitch of the
        divider with its return on the tube wall, so it sees the inductance of
        a coaxial line of that length, (mu0 / 2 pi) * pitch * ln(r_tube / r),
        with r = segment_radius, the equivalent radius of the parts and pours
        that carry it. For the 10.9 mm pitch in the 27.6 mm tube that is about
        6 nH; the 1 nH datasheet ESL of an 1808 capacitor is its own 4.5 mm
        over a return a fraction of a millimetre away. A number sets every
        segment; "parts" keeps the datasheet ESL."""
        v = self.cfg["circuit"].get("segment_l", "parts")
        if v == "parts" or not self.p.esl:
            return None
        if v != "auto":
            return [float(v)] * self.nseg
        r_eq = self.cfg["circuit"].get("segment_radius", 0.8)
        r_tube = self.board.tube["r_in"]
        xs = [self._node_x(n) for n in self.ladder]
        out = []
        for a, b in zip(xs, xs[1:]):
            if a is None or b is None or r_tube <= r_eq:
                return None
            out.append(2e-7 * abs(a - b) * 1e-3 * math.log(r_tube / r_eq))
        return out

    def tip_l_value(self):
        """Inductance in series with the HV tip [H]: the board's own strip from
        the tip to the first divider parts ([circuit] tip_l; "auto" = partial
        inductance of a strip of that length and tip_strip_width) plus
        source_l, whatever connects the probe to the circuit under test."""
        cc = self.cfg["circuit"]
        if not (self.p.esl and self.p.parasitics):
            return 0.0
        v = cc.get("tip_l", 0.0)
        if v == "auto":
            x_div = self._node_x(self.hv)
            xs = [x for (net, _), cu in self.board.copper.items() if net == self.hv for poly in cu.polys for x, _ in poly]
            v = 0.0
            if x_div is not None and xs:
                length = max(abs(max(xs) - x_div), abs(min(xs) - x_div))        # [mm] to the far end of the HV copper
                w = cc.get("tip_strip_width", 1.0)
                if length > w:
                    v = 2e-10 * length * (math.log(2 * length / w) + 0.5 + 0.2235 * w / length)
        return float(v) + float(cc.get("source_l", 0.0))

    def seg_r_text(self):
        v = sorted(set(self.seg_r_values()))
        txt = f"{v[0]:g}" if len(v) == 1 else f"{v[0]:g}..{v[-1]:g}"
        return txt + (" ohm (PCB)" if self.p.seg_r is None and self.seg_res else " ohm")

    # --- network ----------------------------------------------------------
    def _build(self, dc_only=False):
        p, L = self.p, self.labels
        n = Network()
        g = lambda net: "0" if net == "GND" else L[net]

        # divider segments
        topc = self.top_cap_values()
        seg_r = self.seg_r_values()
        seg_l = self.seg_l_values() if p.parasitics else None
        for k, (a, b) in enumerate(zip(self.ladder, self.ladder[1:])):
            for r in self.top_r.get(k, []):
                n.add("R", g(a), g(b), r=self.part_value(r), ref=r)
            refs = self.top.get(k, [])
            esl = self.esl_of(refs[0]) if refs else 0.0
            if seg_l is not None:
                esl = max(esl, seg_l[k])         # the geometric value contains the part's own
            if topc[k] > 0:
                esr = (self.esr_of(refs[0]) if refs else 0.0) + seg_r[k]
                n.add("C", g(a), g(b), c=topc[k], esl=esl, esr=esr, ref="+".join(refs) or f"Cseg{k+1}")
        # bottom network
        for r in self.bottom_r:
            a, b = self.board.parts[r].nets
            n.add("R", g(a), g(b), r=self.part_value(r), ref=r)
        for c in self.bottom_c:
            a, b = self.board.parts[c].nets
            n.add("C", g(a), g(b), c=self.part_value(c), esl=self.esl_of(c), esr=self.esr_of(c), ref=c)
        if self.trim_ref and p.trim:
            a, b = self.board.parts[self.trim_ref].nets
            n.add("C", g(a), g(b), c=p.trim, esl=self.esl_of(self.trim_ref), esr=self.esr_of(self.trim_ref),
                  ref=self.trim_ref)
        for o in self.other:
            c = self.other_capacitance(o)
            if c:
                a, b = self.board.parts[o].nets[:2]
                n.add("C", g(a), g(b), c=c, ref=o)
        if p.extra_c:
            # behaves like one more SMD capacitor (same ESL/ESR) when positive
            cc = self.cfg["circuit"]
            pos = p.extra_c > 0
            n.add("C", self.bot, "0", c=p.extra_c, ref="Cextra",
                  esl=cc.get("extra_esl", 1.0e-9) if pos and p.esl else 0.0,
                  esr=cc.get("extra_esr", 0.05) if pos and p.esr else 0.0)
        if self.pot_ref:
            a, b = self.board.parts[self.pot_ref].nets
            n.add("R", g(a), g(b), r=max(p.pot, 1e-3), ref=self.pot_ref)

        # parasitics from the field solver
        if p.parasitics and self.capdata is not None:
            # The Maxwell matrix without its GND row/column is the nodal
            # capacitance matrix of the parasitics (all pair and ground caps).
            names = self.capdata["names"]
            missing = [L[n] for n in self.ladder if n not in names]
            if missing and not getattr(self, "_warned_missing", False):
                self._warned_missing = True
                print(f"WARNING: the capacitance file has no entry for {len(missing)} divider net(s) "
                      f"({', '.join(missing[:4])}{' ...' if len(missing) > 4 else ''}): it is from another "
                      f"board. Run extract_caps.py again.")
            keep = [i for i, nm in enumerate(names) if nm != "GND" and nm in L]
            nodes = [g(names[i]) for i in keep]
            for nd in nodes:
                n.node(nd)
            sub = lambda A: np.asarray(A)[np.ix_(keep, keep)]
            dm = self.dielectric_models()
            if dm:
                models, eps_solve = dm
                mats = {m: sub(v) for m, v in self.capdata["partial"].items()}
                n.materials = {m: (lambda s, d=models[m], e=eps_solve[m]: d.eps(s) / e) for m in mats}
            else:
                mats = {"static": sub(self.capdata["maxwell"])}
            n.add("CMAT", nodes[0], nodes[0], nodes=nodes, mats=mats,
                  maxwell=sub(self.capdata["maxwell"]), ref="Cp")

        # source
        n.node("HV")
        tip_l = self.tip_l_value()
        if p.source_r + p.tip_r > 0 or tip_l > 0:
            n.add("R", "SRC", "HV", r=p.source_r + p.tip_r, l=tip_l, ref="Rsrc_tip")
            self.fixed = {"SRC": 1.0}
        else:
            self.fixed = {"HV": 1.0}

        # output: series R -> cable -> scope
        node = "OUT"
        if not self.term_ref and p.series_r:
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

def optimize_bottom(probe, f_lo=None, f_hi=None, lo=-2e-9, hi=20e-9, n=120, tiebreak=True):
    """Find the extra BOT-GND capacitance (on top of the PCB's fixed parts and
    the current trimmer setting) with the smallest peak gain error
    max |H/H(DC)| - 1 over [f_lo, f_hi] (n log-spaced points, the reported
    'LF error'). The peak error is quasi-convex in the capacitance, so a
    bounded scalar search finds the minimax optimum.

    Tie-break: the peak is often set by a point the bottom capacitance cannot
    move (the low-frequency dip from the resin dispersion), so it is nearly
    flat over a range of capacitances. With tiebreak, the capacitance with the
    smallest RMS error over the same points is taken from those within
    [optimize] peak_tolerance of the minimum peak. The grading and rounding
    searches compare designs by the strict minimax error (tiebreak=False).

    The extra part is a shunt admittance y at BOT, so the network is solved
    only twice (source on; unit current into BOT with the source off) and the
    response for any y follows exactly from superposition:
        V_scope(y) = V_scope - Z_scope,bot * y * V_bot / (1 + Z_bot,bot * y)
    Returns (extra capacitance [F], peak error there); probe.p.extra_c is set to it."""
    from scipy.optimize import brentq, minimize_scalar
    o = probe.cfg["optimize"]
    f_lo = f_lo or o["f_lo"]
    f_hi = f_hi or o["f_hi"]
    fixed_bottom = probe.bottom_fixed_c() + (probe.p.trim or 0)
    lo = max(lo, -fixed_bottom + 1e-12)

    probe.p.extra_c = 0.0
    net = probe._build()
    f = np.concatenate(([1e-3], np.logspace(math.log10(f_lo), math.log10(f_hi), n)))   # [0] for H(DC), as H0()
    s = 2j * np.pi * f
    V0, _ = net.solve(s, probe.fixed)
    Z, _ = net.solve(s, {k: 0.0 for k in probe.fixed}, inject={probe.bot: 1.0})
    ib, isc = net.idx[probe.bot], net.idx["SCOPE"]
    v_bot, v_scope, z_bb, z_sb = V0[:, ib], V0[:, isc], Z[:, ib], Z[:, isc]
    bw = 1 / (1 + s / (2 * np.pi * probe.p.scope_bw)) if probe.p.scope_bw > 0 else 1.0
    # same element as _build() adds for extra_c: SMD part (ESL/ESR) when positive, ideal when negative
    cc = probe.cfg["circuit"]
    esl = cc.get("extra_esl", 1.0e-9) if probe.p.esl else 0.0
    esr = cc.get("extra_esr", 0.05) if probe.p.esr else 0.0

    def gain_err(c):
        y = 1 / (1 / (s * c) + s * esl + esr) if c > 0 else s * c
        h = (v_scope - z_sb * y * v_bot / (1 + z_bb * y)) * bw
        return np.abs(h[1:] / h[0].real) - 1

    def err(c):
        return float(np.max(np.abs(gain_err(c))))

    r = minimize_scalar(err, bounds=(lo, hi), method="bounded", options={"xatol": 1e-14})
    x = r.x
    tol = o.get("peak_tolerance", 0.01) if tiebreak else 0.0
    if tol > 0 and r.fun > 0:
        # capacitances within tolerance form one interval, as the peak is quasi-convex
        limit = (1 + tol) * r.fun

        def edge(end):
            if err(end) <= limit:
                return end
            return brentq(lambda c: err(c) - limit, min(r.x, end), max(r.x, end), xtol=1e-16)

        a, b = edge(lo), edge(hi)
        if b > a:
            rms = lambda c: float(np.sqrt(np.mean(gain_err(c) ** 2)))
            x = minimize_scalar(rms, bounds=(a, b), method="bounded", options={"xatol": 1e-14}).x
    probe.p.extra_c = x
    probe.net = probe._build()
    return x, err(x)


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


def cap_limit(probe):
    """Largest divider capacitor that can be bought: [optimize] cap_max, a
    value like "270p" or in farad. None = no limit."""
    v = probe.cfg.get("optimize", {}).get("cap_max")
    if not v:
        return None
    return kicad_board.parse_si(v) if isinstance(v, str) else float(v)


def graded_caps(probe, c_bottom_segment, kappa=1.0):
    """Discrete segment capacitances that make the capacitive (HF) voltage
    distribution along the ladder linear, given the parasitic matrix.

    Works from the bottom segment (fixed to c_bottom_segment) upward: the
    current in segment k must supply segment k+1 plus all parasitic currents
    leaving node k when the node voltages are exactly linear.

    kappa scales the parasitic currents: 1 uses the capacitances of the field
    solve; with a dispersive resin the effective value over the compensation
    band differs, which optimize_grading() accounts for.

    With [optimize] cap_max the values are cut off there. The segments at the
    HV end then get less than the grading asks for; a stronger grading (larger
    kappa) makes up for it further down, which optimize_grading() finds.
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
        total[k - 1] = total[k] + nseg * kappa * leak
    caps = [t - sp for t, sp in zip(total, seg_par)]
    lim = cap_limit(probe)
    return [min(c, lim) for c in caps] if lim else caps


def round_to_series(caps, series, c_max=None):
    vals = np.array(series_values(series, 1e-12, 1e-8))
    if c_max:
        vals = vals[vals <= c_max * (1 + 1e-9)]
    return [float(vals[np.argmin(np.abs(np.log(vals / c)))]) for c in caps]


def refine_rounding(probe, caps, series, max_passes=4, min_gain=None):
    """Round each capacitor down or up to `series` for the smallest peak gain
    error, instead of independently to the nearest value.

    Starts from nearest rounding and switches one capacitor at a time to its
    other neighbour (the series value on the other side of the unrounded
    value) while that lowers the error after re-optimising the bottom
    capacitance; repeated until no switch helps. Every value stays within
    one series step of the unrounded grading, so the fast-edge distribution
    stays close to it. Returns (caps, peak error).

    A switch is taken only if it lowers the error by more than the fraction
    `min_gain` (default: [optimize] peak_tolerance, at least 0.1 %). Many
    lists are equally good, and gains below that are not real: taking every
    gain of 0.1 % made the list change in two places when the parasitics
    moved by 0.3 %.
    """
    if min_gain is None:
        min_gain = max(probe.cfg.get("optimize", {}).get("peak_tolerance", 0.01), 1e-3)
    vals = np.array(series_values(series, 1e-12, 1e-8))
    lim = cap_limit(probe)
    if lim:
        vals = vals[vals <= lim * (1 + 1e-9)]

    def bracket(c):
        lo, hi = vals[vals <= c * (1 + 1e-9)], vals[vals >= c * (1 - 1e-9)]
        return {float(lo.max()) if len(lo) else float(vals[0]), float(hi.min()) if len(hi) else float(vals[-1])}

    def err(cs):
        probe.p.top_caps = list(cs)
        probe.p.extra_c = 0.0
        probe.net = probe._build()
        return optimize_bottom(probe, tiebreak=False)[1]

    options = [bracket(c) for c in caps]
    cur = round_to_series(caps, series, lim)
    best = err(cur)
    for _ in range(max_passes):
        improved = False
        for k, opts in enumerate(options):
            for v in opts - {cur[k]}:
                trial = cur[:k] + [v] + cur[k + 1:]
                e = err(trial)
                if e < best * (1 - min_gain):
                    cur, best, improved = trial, e, True
        if not improved:
            break
    return cur, best


def optimize_grading(probe, c_bottom_segment, series=None, kappa=None, kappa_range=(0.8, 1.8), tolerance=1.1):
    """Grade the divider for the flattest response over the compensation band.

    The grading shape comes from graded_caps(); its strength kappa sets how
    strong the parasitics are assumed to be (1 = their value at the field-solve
    frequency; for a static dielectric 1 is optimal). With a dispersive resin
    the parasitics are larger in the kHz band, so a stronger grading flattens
    the low-frequency response, but it also makes the fast-edge voltage
    distribution uneven and costs bandwidth.

    kappa is scanned; the peak gain error over [f_lo, f_hi] (bottom capacitance
    optimised each time) has a broad minimum, and the weakest grading within
    `tolerance` x the minimum is used. Pass kappa to fix it. Values are then
    rounded to `series` (each up or down, see refine_rounding) and the bottom
    capacitance re-optimised.
    Returns (caps, kappa, extra bottom C, peak error).
    """
    def peak_err(k):
        probe.p.top_caps = graded_caps(probe, c_bottom_segment, k)
        probe.p.extra_c = 0.0
        probe.net = probe._build()
        return optimize_bottom(probe, tiebreak=False)[1]

    if kappa is None:
        # with a largest value the best grading is stronger than without
        k_hi = kappa_range[1] + (0.8 if cap_limit(probe) else 0.0)
        ks = np.arange(kappa_range[0], k_hi + 1e-9, 0.05)
        errs = np.array([peak_err(k) for k in ks])
        kappa = float(ks[np.nonzero(errs <= tolerance * errs.min())[0][0]])
        probe.grading_scan = list(zip(ks.tolist(), errs.tolist()))
    caps = graded_caps(probe, c_bottom_segment, kappa)
    if series:
        caps, _ = refine_rounding(probe, caps, series)
    probe.p.top_caps = caps
    probe.p.extra_c = 0.0
    probe.net = probe._build()
    extra, e = optimize_bottom(probe)
    return caps, kappa, extra, e


def file_hash(path):
    import hashlib
    with open(path, "rb") as f:
        return hashlib.sha1(f.read()).hexdigest()


def load_capdata(path, pcb=None):
    """Capacitance matrices from extract_caps.py. With `pcb` given, says so if
    that board file is not the one the field solve read."""
    with open(path) as f:
        cap = json.load(f)
    solved = cap.get("meta", {}).get("pcb_sha1")
    if pcb and solved and solved != file_hash(pcb):
        print(f"NOTE: {os.path.basename(pcb)} has changed since the field solve of {cap['meta'].get('date', '?')} "
              f"({os.path.basename(path)}). If copper, outlines or part positions changed, run extract_caps.py again.")
    return cap
