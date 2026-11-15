"""Circuit-level simulation of the probe: response, compensation, design sweeps
and LTspice export. Uses the parasitic capacitances from extract_caps.py.

    python probe_sim.py info
    python probe_sim.py response [--optimize] [--length 2] [--series-r 47] [--top-cap 100p | --graded 47p]
    python probe_sim.py optimize            # compensation capacitor per cable length
    python probe_sim.py sweep               # top-cap / back-termination design study
    python probe_sim.py damping --graded 47p --round E12   # damping resistor study for one design
    python probe_sim.py spice [--optimize] [--run]   # LTspice netlist (+ cross-check)

Common options override probe.toml for one run.
"""

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")   # batched small solves: ~80x faster single-threaded
import argparse
import csv
import re
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import kicad_board
import probe_model as pm
from extract_caps import load_config
from kicad_board import parse_si

# categorical slots, fixed order (blue, orange, aqua, yellow, magenta, green, violet, red)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "lines.linewidth": 1.6,
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb",
    "axes.facecolor": "#fcfcfb", "legend.frameon": False, "font.size": 9,
})


def fmt_c(c):
    if c is None:
        return "-"
    a = abs(c)
    if a >= 1e-9:
        return f"{c*1e9:.3g} nF"
    return f"{c*1e12:.3g} pF"


def fmt_f(f):
    if not np.isfinite(f):
        return ">3 GHz"
    for div, unit in ((1e9, "GHz"), (1e6, "MHz"), (1e3, "kHz")):
        if f >= div:
            return f"{f/div:.3g} {unit}"
    return f"{f:.3g} Hz"


# ----------------------------------------------------------------------------

class Ctx:
    def __init__(self, a):
        self.cfg = load_config(a.config)
        self.res = self.cfg["files"]["results"]
        self.board = kicad_board.load_board(self.cfg["files"]["pcb"])
        capfile = a.capacitance or self.cfg["files"]["capacitance"]
        self.capdata = pm.load_capdata(capfile) if os.path.exists(capfile) else None
        if self.capdata is None:
            print(f"NOTE: {capfile} not found - run extract_caps.py first. Simulating without parasitics.")
        self.a = a

    def params(self, **over):
        a = self.a
        o = dict(length=a.length, series_r=a.series_r, cable_type=a.cable, trim=a.trim,
                 c_scope=a.c_scope, parasitics=not a.no_parasitics, tip_r=a.tip_r, seg_r=a.seg_r,
                 scope_bw=a.scope_bw)
        if a.no_loss:
            o["dielectric_loss"] = False
        if a.no_esr:
            o["esr"] = False
        if a.top_cap:
            o["top_caps"] = parse_si(a.top_cap)
        o.update(over)
        return pm.params_from_config(self.cfg, **o)

    def probe(self, **over):
        p = pm.Probe(self.cfg, self.board, self.capdata, self.params(**over))
        if getattr(self.a, "graded", None) and "top_caps" not in over:
            self.apply_grading(p, parse_si(self.a.graded), self.a.round, kappa=self.a.kappa)
        return p

    @staticmethod
    def apply_grading(p, base, series=None, fit=True, kappa=None):
        """Grade the divider capacitors; with fit=True the grading strength is
        chosen for the flattest response (see probe_model.optimize_grading)."""
        if fit or kappa:
            caps, kappa, _, _ = pm.optimize_grading(p, base, series, kappa=kappa)
            p.grading_kappa = kappa
            return caps
        caps = pm.graded_caps(p, base)
        if series:
            caps = pm.round_to_series(caps, series)
        p.p.top_caps = caps
        p.net = p._build()
        return caps


def optimized(p, quiet=False):
    """Optimise the bottom capacitance in place; returns (total C to add vs PCB, error)."""
    extra, err = pm.optimize_bottom(p)
    return extra, err


def recommend(p, extra):
    """Turn the optimum extra capacitance into fixed parts + trimmer setting."""
    trim_lo, trim_hi = p.trim_range or (0, 0)
    current_fixed = p.bottom_fixed_c()
    needed_total = current_fixed + (p.p.trim or 0) + extra       # fixed + trimmer
    trim_mid = 0.5 * (trim_lo + trim_hi)
    target_fixed = needed_total - trim_mid
    lines = [f"  required fixed+trimmer capacitance: {fmt_c(needed_total)} "
             f"(PCB now: {fmt_c(current_fixed)} fixed + trimmer {fmt_c(trim_lo)}..{fmt_c(trim_hi)})"]
    if needed_total < trim_lo:
        lines.append("  -> less than the trimmer minimum: remove the fixed bottom capacitors and/or "
                     "reduce cable capacitance; the trimmer alone is too large.")
        return lines, needed_total
    if trim_lo <= needed_total - current_fixed <= trim_hi:
        lines.append(f"  -> reachable with the current fixed parts, trimmer at {fmt_c(needed_total - current_fixed)}")
    if target_fixed > 0:
        o = p.cfg["optimize"]
        for parts, tot in pm.suggest_fixed(target_fixed, o["series"], o["max_fixed_parts"])[:3]:
            t = needed_total - tot
            ok = "ok" if trim_lo <= t <= trim_hi else "out of trimmer range"
            lines.append(f"  -> fixed {' + '.join(fmt_c(v) for v in parts):<22} trimmer at {fmt_c(t)} ({ok})")
    return lines, needed_total


# ----------------------------------------------------------------------------
# commands

def cmd_info(ctx):
    p = ctx.probe()
    L = p.labels
    print(f"divider: {p.nseg} segments, HV net {p.hv} -> output {p.out}")
    topc = p.top_cap_values()
    r_top = sum(p.part_value(r) for k in p.top_r for r in p.top_r[k])
    c_series = 1 / sum(1 / c for c in topc)
    print(f"top: R = {r_top/1e6:.3g} MOhm, discrete C in series = {fmt_c(c_series)}")
    print(f"bottom: R parts {', '.join(p.bottom_r)}, pot {p.pot_ref} = {p.p.pot/1e3:.1f} kOhm, "
          f"C parts {', '.join(p.bottom_c)} = {fmt_c(p.bottom_fixed_c())}, trimmer {p.trim_ref} {fmt_c(p.p.trim)}")
    cab = ctx.cfg["cable"]["types"][p.p.cable_type]
    c_cab = p.p.length / (cab["z0"] * cab["vf"] * pm.C0)
    term = f"{p.term_ref} = {p.series_r_value():g} Ohm" if p.term_ref else f"{p.series_r_value():g} Ohm (added)"
    print(f"back termination: {term}")
    for o in p.other:
        c = p.other_capacitance(o)
        print(f"{o} ({p.board.parts[o].value}): " + (f"modelled as {fmt_c(c)}" if c else
              "NOT MODELLED - set its capacitance in [circuit.part_capacitance]"))
    dm = p.dielectric_models()
    if dm is None:
        print("dielectrics: static (no dispersion/loss), parasitics as solved")
    else:
        models, eps_solve = dm
        print("dielectrics (eps' / tan d), parasitics scaled per material from the field solve at "
              + ", ".join(f"{k} {v:.2f}" for k, v in eps_solve.items()) + ":")
        for k, d in models.items():
            if d.k:
                print(f"  {k:<6} {d.describe()}")
    print(f"cable: {p.p.cable_type} {p.p.length} m ({fmt_c(c_cab)}); "
          f"scope {p.p.r_scope/1e6:g} MOhm || {fmt_c(p.p.c_scope)}")
    print(f"DC ratio: {p.ratio():.2f}")
    if ctx.capdata:
        names = ctx.capdata["names"]
        M = np.array(ctx.capdata["maxwell"])
        gi = names.index("GND")
        cg = sum(-M[names.index(n), gi] for n in p.ladder[1:-1])
        print(f"parasitic: sum of divider node capacitance to the shield = {fmt_c(cg)} "
              f"(vs {fmt_c(c_series)} series) -> alpha = sqrt(Cg/Cs) = {np.sqrt(cg/c_series):.2f}")
    ideal = ctx.probe(parasitics=False, top_caps=p.top_cap_values())
    e_i, _ = optimized(ideal)
    need_ideal = ideal.bottom_fixed_c() + ideal.p.trim + e_i
    print(f"ideal divider (no parasitics) needs fixed+trimmer = {fmt_c(need_ideal)}")
    m = pm.response_metrics(p)
    print(f"as built: LF error {m['lf_error_pct']:.1f} %, C_in {m['c_in_1M_pF']:.2f} pF @1 MHz")
    if ctx.capdata:
        e, err = optimized(p)
        m = pm.response_metrics(p)
        print(f"with parasitics, best bottom C leaves LF error {m['lf_error_pct']:.1f} %")
        for l in recommend(p, e)[0]:
            print(l)


def _plot_response(ctx, variants, path, title):
    f = np.logspace(0, np.log10(pm.F_VALID), 900)
    fig, axs = plt.subplots(3, 2, figsize=(13, 12))
    ax_m, ax_p, ax_tf, ax_ts, ax_v, ax_z = axs.ravel()
    t_fast = ctx.a.t_fast
    for n, (name, p) in enumerate(variants):
        col = SERIES[n % len(SERIES)]
        h = p.H(f) / p.H0()
        ax_m.semilogx(f, 20 * np.log10(np.abs(h)), color=col, label=name)
        k = f <= 1e7
        ax_p.semilogx(f[k], (np.abs(h[k]) - 1) * 100, color=col, label=name)
        t, y = p.step(t_fast, n=2**13, rise=ctx.a.rise)
        ax_tf.plot(t * 1e9, y, color=col, label=name)
        t, y = p.step(20e-3, n=2**16)
        k = t > 0.5e-6
        ax_ts.semilogx(t[k] * 1e6, y[k], color=col, label=name)
        Yin = p.input_admittance(f)
        ax_z.loglog(f, 1 / np.abs(Yin), color=col, label=name)
    ax_m.set(xlabel="frequency [Hz]", ylabel="|H / H(DC)| [dB]", title="Frequency response (normalised to DC ratio)")
    ax_m.set_ylim(max(ax_m.get_ylim()[0], -30), min(ax_m.get_ylim()[1], 30))
    ax_p.set(xlabel="frequency [Hz]", ylabel="gain error [%]",
             title="Compensation flatness up to 10 MHz (0 % = perfect)")
    ax_p.set_ylim(max(ax_p.get_ylim()[0], -60), min(ax_p.get_ylim()[1], 60))
    ax_tf.set(xlabel="time [ns]", ylabel="output / DC final value",
              title=f"Step response, fast (source rise {ctx.a.rise*1e9:g} ns)")
    ax_ts.set(xlabel="time [µs]", ylabel="output / DC final value", title="Step response, slow (compensation)")
    ax_z.set(xlabel="frequency [Hz]", ylabel="|Z_in| [Ohm]", title="Probe input impedance (DUT loading)")

    # voltage distribution of the first variant
    p = variants[0][1]
    seg = np.arange(1, p.nseg + 1)
    for n, (f0, lab) in enumerate(((1e-3, "DC"), (ctx.cfg["stress"]["hf_frequency"], "fast edge"))):
        v = p.node_voltages(f0)
        vv = np.array([v[net].real for net in p.ladder])
        share = (vv[:-1] - vv[1:]) / (vv[0] - vv[-1]) * 100
        ax_v.bar(seg + (n - 0.5) * 0.4, share, width=0.38, color=SERIES[n], label=lab)
    ax_v.axhline(100 / p.nseg, color=INK2, lw=0.8, ls="--")
    ax_v.set(xlabel="segment (1 = HV end)", ylabel="share of total voltage [%]",
             title=f"Voltage per segment: {variants[0][0]}")
    ax_v.set_xticks(seg)
    for ax in axs.ravel():
        ax.legend(fontsize=8)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def print_metrics(name, p, extra=None):
    m = pm.response_metrics(p)
    s = pm.step_metrics(p)
    print(f"{name}")
    print(f"  DC ratio {m['ratio_dc']:.1f}   C_in {m['c_in_1M_pF']:.2f} pF @1 MHz   "
          f"LF error {m['lf_error_pct']:.2f} %   f(-3 dB) {fmt_f(m['f_3dB'])}   f(+/-1 dB) {fmt_f(m['f_1dB'])}   "
          f"peaking {m['peaking_dB']:.1f} dB   HF peak {m['hf_peak_dB']:.1f} dB   rise {s['rise_ns']:.2f} ns   overshoot {s['overshoot_pct']:.0f} %"
          + (f"   bottom C change {fmt_c(extra)}" if extra is not None else ""))
    return m, s


def cmd_response(ctx):
    p = ctx.probe()
    variants = []
    extra = None
    if ctx.a.optimize:
        extra, _ = optimized(p)
    if ctx.a.graded or ctx.a.top_cap:
        refs = [r for k in sorted(p.top) for r in p.top[k]]
        print("divider capacitors: " + ", ".join(f"{r} {fmt_c(c)}" for r, c in zip(refs, p.top_cap_values())))
        if hasattr(p, "grading_kappa"):
            print(f"grading strength kappa = {p.grading_kappa:.3f}")
    variants.append(("this design" + (" (bottom C optimised)" if extra is not None else ""), p))
    print_metrics(variants[0][0], p, extra)
    if ctx.capdata and not ctx.a.no_parasitics:
        q = ctx.probe(parasitics=False, top_caps=p.top_cap_values())   # same capacitors
        if ctx.a.optimize:
            optimized(q)
        variants.append(("same capacitors, no parasitics", q))
        print_metrics(variants[1][0], q)
    path = os.path.join(ctx.res, ctx.a.name + ".png")
    _plot_response(ctx, variants, path,
                   f"{p.p.cable_type} {p.p.length} m, back termination {p.series_r_value():g} Ω, "
                   f"scope {fmt_c(p.p.c_scope)}")
    print(f"wrote {path}")


def cmd_optimize(ctx):
    rows = []
    for length in ctx.cfg["sweep"]["cable_lengths"]:
        p = ctx.probe(length=length)
        extra, err = optimized(p)
        print(f"\ncable {p.p.cable_type} {length} m, back termination {p.series_r_value():g} Ohm")
        lines, total = recommend(p, extra)
        for l in lines:
            print(l)
        m, s = print_metrics("  result", p)
        rows.append([length, total * 1e12, m["lf_error_pct"], m["f_3dB"], m["f_1dB"], m["peaking_dB"],
                     s["rise_ns"], s["overshoot_pct"], m["c_in_1M_pF"]])
    path = os.path.join(ctx.res, ctx.a.name + ".csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cable_m", "bottom_C_pF", "lf_error_pct", "f_3dB_Hz", "f_1dB_Hz", "peaking_dB",
                    "rise_ns", "overshoot_pct", "c_in_pF"])
        w.writerows(rows)
    print(f"\nwrote {path}")


def cmd_sweep(ctx):
    sw = ctx.cfg["sweep"]
    variants = [("uniform", v) for v in sw["top_cap"]]
    variants += [("graded", g) for g in sw.get("graded_base", ["47p"])]
    rows = []
    print(f"{'top caps':<22}{'Rser':>6}{'C_in pF':>9}{'LF err %':>9}{'f-3dB':>11}{'f+/-1dB':>11}"
          f"{'peak dB':>8}{'HFpk dB':>8}{'rise ns':>8}{'ovs %':>7}{'bottom C':>11}{'max seg V':>10}")
    for kind, v in variants:
        for rs in sw["series_r"]:
            if kind == "uniform":
                p = ctx.probe(top_caps=parse_si(v), series_r=rs)
                label = f"{v} x{p.nseg}"
            else:
                p = ctx.probe(series_r=rs)
                caps = Ctx.apply_grading(p, parse_si(v), ctx.cfg["optimize"]["series"])
                label = f"graded {fmt_c(caps[0])}..{fmt_c(caps[-1])}"
            extra, err = optimized(p)
            m = pm.response_metrics(p)
            s = pm.step_metrics(p)
            vhf = p.node_voltages(ctx.cfg["stress"]["hf_frequency"])
            vv = np.array([vhf[n].real for n in p.ladder])
            seg_max = float(np.max(vv[:-1] - vv[1:]) / (vv[0] - vv[-1]) * p.nseg)
            total = p.bottom_fixed_c() + p.p.trim + extra
            print(f"{label:<22}{rs:6.0f}{m['c_in_1M_pF']:9.2f}{m['lf_error_pct']:9.2f}{fmt_f(m['f_3dB']):>11}"
                  f"{fmt_f(m['f_1dB']):>11}{m['peaking_dB']:8.1f}{m['hf_peak_dB']:8.1f}{s['rise_ns']:8.2f}{s['overshoot_pct']:7.0f}"
                  f"{fmt_c(total):>11}{seg_max:9.2f}x")
            rows.append(dict(variant=label, kind=kind, base=v, series_r=rs, c_in_pF=m["c_in_1M_pF"],
                             lf_error_pct=m["lf_error_pct"], f_3dB=m["f_3dB"], f_1dB=m["f_1dB"],
                             peaking_dB=m["peaking_dB"], hf_peak_dB=m["hf_peak_dB"], rise_ns=s["rise_ns"], overshoot_pct=s["overshoot_pct"],
                             bottom_C_pF=total * 1e12, seg_stress_factor=seg_max))
    path = os.path.join(ctx.res, ctx.a.name + ".csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    panels = ((axs[0], "lf_error_pct", "LF error after compensation [%]"),
              (axs[1], "f_1dB", "±1 dB bandwidth [Hz]"),
              (axs[2], "overshoot_pct", "step overshoot [%]"))
    for n, rs in enumerate(sw["series_r"]):
        for kind, marker in (("uniform", "o"), ("graded", "D")):
            sel = [r for r in rows if r["series_r"] == rs and r["kind"] == kind]
            if not sel:
                continue
            x = [r["c_in_pF"] for r in sel]
            for ax, key, lab in panels:
                y = [min(r[key], 3e9) for r in sel]
                ax.plot(x, y, marker, ms=8 if marker == "o" else 7, color=SERIES[n],
                        label=f"series R {rs:g} Ω" + (", graded" if kind == "graded" else ""),
                        markeredgecolor="#fcfcfb", markeredgewidth=1.5)
                ax.set(xlabel="input capacitance [pF]", ylabel=lab)
    axs[1].set_yscale("log")
    axs[0].set_yscale("log")
    # label each design once (LF error does not depend on the series resistor)
    k_graded = 0
    for r in rows:
        if r["series_r"] != sw["series_r"][-1]:
            continue
        if r["kind"] == "graded":
            txt, off = f"graded from {r['base']}", (8, -14 + 12 * k_graded)
            k_graded += 1
        else:
            txt, off = r["base"], (5, 4)
        axs[0].annotate(txt, (r["c_in_pF"], r["lf_error_pct"]), fontsize=7, color=INK2,
                        xytext=off, textcoords="offset points")
    for ax in axs:
        ax.legend(fontsize=8)
    cab = ctx.params()
    fig.suptitle(f"Design trade-off: top capacitor value vs. DUT loading ({cab.cable_type} {cab.length} m, "
                 f"tip R {cab.tip_r:g} Ω, segment R {cab.seg_r:g} Ω)")
    fig.tight_layout()
    path2 = os.path.join(ctx.res, ctx.a.name + ".png")
    fig.savefig(path2, dpi=120)
    plt.close(fig)
    print(f"wrote {path} and {path2}")


# ----------------------------------------------------------------------------
# LTspice

def _spice_val(v):
    return f"{v:.6g}"


def write_spice(p, path, meas_f):
    lines = [f"* HV probe model generated by tools/probe_sim.py from {os.path.basename(p.board.path)}",
             f"* cable {p.p.cable_type} {p.p.length} m, back termination {p.series_r_value():g} ohm, pot {p.p.pot:.1f} ohm",
             "* parasitic capacitances from the 3D field solver (Cp_*)",
             "* dispersion/loss of the dielectrics is NOT included: parasitics at their field-solve value", ""]
    n_r = n_c = 0
    for kind, a, b, prm in p.net.elems:
        ref = re.sub(r"[^A-Za-z0-9_]", "_", prm.get("ref", ""))
        if kind == "R":
            n_r += 1
            lines.append(f"R{ref or n_r}_{n_r} {a} {b} {_spice_val(prm['r'])}")
        elif kind == "CMAT":
            # expand the nodal matrix into pair and ground capacitors
            M, nodes = prm["maxwell"], prm["nodes"]
            for i in range(len(nodes)):
                cg = M[i].sum()
                if cg > 1e-17:
                    lines.append(f"Cp_{nodes[i]}_GND {nodes[i]} 0 {_spice_val(cg)}")
                for j in range(i + 1, len(nodes)):
                    if -M[i, j] > 1e-17:
                        lines.append(f"Cp_{nodes[i]}_{nodes[j]} {nodes[i]} {nodes[j]} {_spice_val(-M[i, j])}")
        elif kind == "C":
            n_c += 1
            extra = f" Lser={_spice_val(prm['esl'])}" if prm.get("esl") else ""
            extra += f" Rser={_spice_val(prm['esr'])}" if prm.get("esr") else ""
            lines.append(f"C{ref or n_c}_{n_c} {a} {b} {_spice_val(prm['c'])}{extra}")
        elif kind == "TL":
            Lp = prm["z0"] / (prm["vf"] * pm.C0)
            Cp = 1 / (prm["z0"] * prm["vf"] * pm.C0)
            lines.append(f"O1 {a} 0 {b} 0 CABLE")
            lines.append(f".model CABLE LTRA(len={prm['length']} R={prm['r_dc']} L={_spice_val(Lp)} C={_spice_val(Cp)})")
            lines.append("* note: LTRA has no skin effect; probe_sim.py models R ~ sqrt(f)")
    src = list(p.fixed)[0]
    lines += ["", f"V1 {src} 0 AC 1 PULSE(0 1 0 1n 1n 1 2)",
              ".ac dec 100 1 1G", "* .tran 0 2u 0 0.1n", ".options plotwinsize=0", ""]
    for i, fm in enumerate(meas_f):
        lines.append(f".meas AC g{i} FIND mag(V(SCOPE)) AT {fm:g}")
    lines += [".end", ""]
    with open(path, "w") as f:
        f.write("\n".join(lines))


def find_ltspice(cfg):
    cands = [cfg.get("ltspice", {}).get("exe", ""),
             os.path.expandvars(r"%LOCALAPPDATA%\Programs\ADI\LTspice\LTspice.exe"),
             r"C:\Program Files\ADI\LTspice\LTspice.exe",
             r"C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe"]
    return next((c for c in cands if c and os.path.exists(c)), None)


def cmd_damping(ctx):
    """Damping study for one divider design: segment resistor x tip resistor."""
    sw = ctx.cfg["sweep"]
    base = ctx.probe()                     # applies --graded / --top-cap once
    caps = base.top_cap_values()
    rows = []
    print(f"divider capacitors (HV end first): {', '.join(fmt_c(c) for c in caps)}")
    if hasattr(base, "grading_kappa"):
        print(f"grading strength kappa = {base.grading_kappa:.3f} (1 = parasitics at the field-solve frequency)")
    print(f"{'seg R':>6}{'tip R':>7}{'C_in pF':>9}{'LF err %':>9}{'f-3dB':>11}{'f+/-1dB':>11}"
          f"{'HFpk dB':>8}{'rise ns':>8}{'ovs %':>7}")
    for tip in sw.get("tip_r", [0.0]):
        for seg in sw.get("seg_r", [0.0]):
            p = ctx.probe(top_caps=caps, seg_r=float(seg), tip_r=float(tip))
            optimized(p)
            m = pm.response_metrics(p)
            s = pm.step_metrics(p)
            print(f"{seg:6g}{tip:7g}{m['c_in_1M_pF']:9.2f}{m['lf_error_pct']:9.2f}{fmt_f(m['f_3dB']):>11}"
                  f"{fmt_f(m['f_1dB']):>11}{m['hf_peak_dB']:8.1f}{s['rise_ns']:8.2f}{s['overshoot_pct']:7.0f}")
            rows.append(dict(seg_r=seg, tip_r=tip, c_in_pF=m["c_in_1M_pF"], lf_error_pct=m["lf_error_pct"],
                             f_3dB=m["f_3dB"], f_1dB=m["f_1dB"], hf_peak_dB=m["hf_peak_dB"],
                             rise_ns=s["rise_ns"], overshoot_pct=s["overshoot_pct"]))
    path = os.path.join(ctx.res, ctx.a.name + ".csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    for n, tip in enumerate(sw.get("tip_r", [0.0])):
        sel = [r for r in rows if r["tip_r"] == tip]
        x = [r["seg_r"] for r in sel]
        for ax, key, lab in ((axs[0], "f_3dB", "-3 dB bandwidth [MHz]"), (axs[1], "rise_ns", "rise time 10-90 % [ns]"),
                             (axs[2], "overshoot_pct", "step overshoot [%]")):
            y = [r[key] / (1e6 if key == "f_3dB" else 1) for r in sel]
            ax.plot(x, y, "-o", ms=6, color=SERIES[n], label=f"tip R {tip:g} Ω")
            ax.set(xlabel="damping resistor per segment [Ω]", ylabel=lab)
    for ax in axs:
        ax.legend(fontsize=8)
    cab = ctx.params()
    fig.suptitle(f"Damping study ({cab.cable_type} {cab.length} m, source edge {ctx.cfg.get('analysis', {}).get('source_rise', 1e-9)*1e9:g} ns)")
    fig.tight_layout()
    path2 = os.path.join(ctx.res, ctx.a.name + ".png")
    fig.savefig(path2, dpi=120)
    plt.close(fig)
    print(f"wrote {path} and {path2}")


def cmd_spice(ctx):
    p = ctx.probe()
    if ctx.a.optimize:
        extra, _ = optimized(p)
        print(f"bottom capacitance changed by {fmt_c(extra)} (Cextra in the netlist)")
    path = os.path.join(ctx.res, ctx.a.name + ".cir")
    meas_f = [1e2, 1e4, 1e6, 1e7, 3e7, 1e8]
    write_spice(p, path, meas_f)
    print(f"wrote {path}")
    if not ctx.a.run:
        return
    exe = find_ltspice(ctx.cfg)
    if not exe:
        print("LTspice not found; set [ltspice] exe = '...' in probe.toml")
        return
    subprocess.run([exe, "-b", os.path.abspath(path)], check=False, timeout=600)
    log = os.path.splitext(path)[0] + ".log"
    raw = open(log, "rb").read()
    txt = raw.decode("utf-16-le", errors="ignore") if raw.count(b"\x00") > len(raw) // 4 else raw.decode("latin-1")
    p.p.dielectric_loss = False        # same network as the netlist
    p.net = p._build()
    h = p.H(meas_f)
    print("(probe_sim evaluated without dielectric loss, like the netlist)")
    print(f"{'f':>10}{'LTspice':>14}{'probe_sim':>14}{'diff':>9}")
    for i, fm in enumerate(meas_f):
        m = re.search(rf"g{i}: mag\(v\(scope\)\)=\(([-0-9.eE+]+)dB", txt, re.I)
        if not m:
            print(f"{fm:10.3g}  (no LTspice result)")
            continue
        lt = float(m.group(1))
        py = 20 * np.log10(abs(h[i]))
        print(f"{fm:10.3g}{lt:12.3f}dB{py:12.3f}dB{py-lt:8.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["info", "response", "optimize", "sweep", "damping", "spice"])
    ap.add_argument("--config", default=os.path.join(os.path.dirname(__file__), "probe.toml"))
    ap.add_argument("--capacitance", help="capacitance json (default from config)")
    ap.add_argument("--length", type=float, help="cable length [m]")
    ap.add_argument("--cable", help="cable type from probe.toml")
    ap.add_argument("--series-r", type=float, help="back-termination resistance [ohm] (replaces R29's value)")
    ap.add_argument("--tip-r", type=float, help="damping resistor in series with the HV tip [ohm]")
    ap.add_argument("--seg-r", type=float, help="damping resistor in series with each divider capacitor [ohm]")
    ap.add_argument("--scope-bw", type=parse_si, help="scope bandwidth (1st order), e.g. 100Meg; 0 = ideal")
    ap.add_argument("--trim", type=parse_si, help="trimmer capacitance, e.g. 30p")
    ap.add_argument("--c-scope", type=parse_si, help="scope input capacitance, e.g. 20p")
    ap.add_argument("--top-cap", help="replace every divider capacitor with this value, e.g. 100p")
    ap.add_argument("--graded", help="grade the divider capacitors for a linear HF distribution, "
                                     "starting from this value at the bottom, e.g. 47p")
    ap.add_argument("--round", help="round graded values to this E-series (E6/E12/E24)")
    ap.add_argument("--kappa", type=float, help="grading strength (default: chosen automatically, see README)")
    ap.add_argument("--no-parasitics", action="store_true")
    ap.add_argument("--no-loss", action="store_true", help="ignore dielectric loss of the parasitics")
    ap.add_argument("--no-esr", action="store_true", help="ignore capacitor ESR")
    ap.add_argument("--optimize", action="store_true", help="optimise the bottom capacitance first")
    ap.add_argument("--rise", type=float, default=1e-9, help="source rise time for step plots [s]")
    ap.add_argument("--t-fast", type=float, default=200e-9, help="fast step plot window [s]")
    ap.add_argument("--name", default=None, help="output file name (without extension)")
    ap.add_argument("--run", action="store_true", help="spice: run LTspice and compare")
    a = ap.parse_args()
    if a.name is None:
        a.name = {"spice": "probe"}.get(a.command, a.command)
    ctx = Ctx(a)
    {"info": cmd_info, "response": cmd_response, "optimize": cmd_optimize,
     "sweep": cmd_sweep, "damping": cmd_damping, "spice": cmd_spice}[a.command](ctx)


if __name__ == "__main__":
    main()
