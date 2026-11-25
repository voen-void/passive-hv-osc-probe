"""Resin permittivity from a sweep of a cast sample on the HM8118 bridge.

    python resin_fit.py                      # fit + plot from the corrected sweep in probe.toml
    python resin_fit.py --part RAW.csv --ref REF_range5.csv REF_range4.csv [--ref-auto REF_auto.csv]
                                             # correct a raw sweep with a C0G reference first

Raw files are the per-frequency summaries of hm8118_sweep.py (f_set_Hz,
Cp_F_mean, Cp_F_sem, D_mean, D_sem, range, flag). The reference is a C0G
capacitor (D ~ 1e-5) read in the same lead layout; what the bridge shows on
it is the bridge's own error. Three errors are taken off the sample's sweep:

  * D: the reference's D on the same range and frequency is subtracted. Up to
    about 30 kHz that error does not depend on the part's impedance. Above, it
    grows with the impedance, and the sample (~40 pF) is 6x the reference
    (225 pF): the correction is applied once and the possible rest (taken as
    3x the correction) goes into the point's uncertainty.
  * C, spot frequencies: at 0.4, 0.45, 0.75, 1.5, 3, 9, 12, 25 and 45 kHz the
    bridge adds a fixed ~0.35 pF on every range (its stored open correction).
    The amount is read off the reference and subtracted.
  * C, range: range 6 reads 0.14 % below range 5 on the reference.

The reference's rise of C above 30 kHz (+0.45 % at 100 kHz) is not a factor
common to all parts: applied as one, it breaks the sample's Kramers-Kronig
consistency. It behaves like a series impedance, so it is scaled by
C_sample / C_reference (0.08 % at 100 kHz) and half the full value is kept as
uncertainty.

Left out of the fit: below --f-min (mains hum on a 45 pF part), above --f-max
(range 4 fails on parts this small: D collapses towards 200 kHz), odd mains
harmonics, and readings the sweep flagged.

Outputs: the corrected sweep (--out, default from probe.toml), a fit report,
and results/resin_fit.png. The models are in dielectric.py.
"""

import argparse
import csv
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import dielectric
from extract_caps import load_config
from probe_sim import SERIES, INK2

JUMP = (400, 450, 750, 1500, 3000, 9000, 12000, 25000, 45000)     # fixed-capacitance artefacts [Hz]
MAINS = 50.0


def hum(f):
    """Odd mains harmonics up to 1 kHz: steady offsets that averaging does not remove."""
    return f <= 1000 and f % MAINS == 0 and (f / MAINS) % 2 == 1


def read_raw(path):
    rows = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            ccol = next(k for k in r if k.endswith("_F_mean"))
            pre = ccol[:-len("_mean")]
            rows[float(r["f_set_Hz"])] = dict(
                f=float(r["f_set_Hz"]), C=float(r[ccol]) * 1e12, Csem=float(r[pre + "_sem"]) * 1e12,
                D=float(r["D_mean"]), Dsem=float(r["D_sem"]), rng=r["range"].strip(), flag=r.get("flag", "").strip())
    return rows


def correct(part, refs, f_min, f_max):
    """Apply the reference corrections; returns the rows with C_c, D_c, sigmas, use."""
    by_range = {}
    for ref in refs:
        for f, r in ref.items():
            by_range.setdefault(r["rng"], {})[f] = r

    def smooth_c(rng, f):
        """reference C on a range at f, without the spot-frequency steps and hum"""
        pts = sorted((k, v["C"]) for k, v in by_range[rng].items()
                     if k not in JUMP and not hum(k) and k >= 80 and k != 150000)
        return float(np.interp(math.log(f), np.log([p[0] for p in pts]), [p[1] for p in pts]))

    # anchor: the reference at 1 kHz on range 5 (the bridge's specified point)
    c_anchor = by_range["5"][1000.0]["C"] if "5" in by_range and 1000.0 in by_range["5"] else None
    gain = {}
    for rng in by_range:
        fs = [f for f in by_range[rng] if f not in JUMP and 100 <= f <= 10000 and not hum(f)]
        gain[rng] = (np.median([by_range[rng][f]["C"] for f in fs]) / c_anchor) if (c_anchor and fs) else 1.0
    c_ref = c_anchor or np.median([r["C"] for ref in refs for r in ref.values()])

    out = []
    for f in sorted(part):
        r = dict(part[f])
        rng = r["rng"]
        ref = by_range.get(rng, {})
        # D error on this range; range 6 needs none (right to 1e-4 away from mains harmonics)
        dD = ref[f]["D"] if (f in ref and rng != "6") else 0.0
        # spot-frequency step: the same fixed capacitance on every range
        jump = 0.0
        if f in JUMP:
            for rr in ([rng] if rng in by_range else []) + [k for k in by_range if k != rng]:
                if f in by_range[rr]:
                    jump = by_range[rr][f]["C"] - smooth_c(rr, f)
                    break
        hf = 0.0
        if rng in by_range and rng != "6":
            base = np.median([v["C"] for k, v in by_range[rng].items() if k not in JUMP and k <= 10000])
            hf = smooth_c(rng, f) / base - 1.0
        r["C_c"] = (r["C"] - jump) / gain.get(rng, 1.0) / (1 + hf * r["C"] / c_ref)
        r["D_c"] = r["D"] - dD
        r["sC"] = math.sqrt(r["Csem"] ** 2 + (0.0015 * r["C_c"]) ** 2 + (0.5 * hf * r["C_c"]) ** 2)
        r["sD"] = math.sqrt(r["Dsem"] ** 2 + 0.0004 ** 2 + ((3.0 * dD) ** 2 if f > 30e3 else 0.0))
        r["use"] = int(f_min <= f <= f_max and not hum(f) and not r["flag"])
        r["why"] = ("" if r["use"] else "below f-min" if f < f_min else "above f-max" if f > f_max
                    else "mains harmonic" if hum(f) else r["flag"])
        out.append(r)
    return out


def write_corrected(rows, path, header):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        for line in header:
            fh.write("# " + line + "\n")
        w = csv.writer(fh)
        w.writerow(["f_Hz", "C_pF", "tan_d", "sigma_C_pF", "sigma_tan_d", "use", "range", "C_raw_pF", "tan_d_raw", "note"])
        for r in rows:
            w.writerow([f"{r['f']:g}", f"{r['C_c']:.4f}", f"{r['D_c']:.6f}", f"{r['sC']:.4f}", f"{r['sD']:.6f}",
                        r["use"], r["rng"], f"{r['C']:.4f}", f"{r['D']:.6f}", r["why"]])


def kk_check(f, c, d):
    """Local Kramers-Kronig check per band: tan d against (pi/2) * -dlnC/dlnf."""
    out = []
    lf, lc = np.log(f), np.log(c)
    edges = np.logspace(math.log10(f.min()), math.log10(f.max()), 5)
    for a, b in zip(edges, edges[1:]):
        m = (f >= a * 0.999) & (f <= b * 1.001)
        if m.sum() >= 4:
            out.append((a, b, -(math.pi / 2) * np.polyfit(lf[m], lc[m], 1)[0], float(d[m].mean())))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(os.path.dirname(__file__), "probe.toml"))
    ap.add_argument("--part", help="raw sweep of the sample")
    ap.add_argument("--ref", nargs="+", help="raw sweeps of the C0G reference with a range held (range 5, range 4)")
    ap.add_argument("--ref-auto", help="raw autorange sweep of the reference (for range 6, below 800 Hz)")
    ap.add_argument("--out", help="corrected sweep to write (default: [materials] resin_data)")
    ap.add_argument("--f-min", type=float, default=200.0)
    ap.add_argument("--f-max", type=float, default=120e3)
    ap.add_argument("--name", default="resin_fit")
    a = ap.parse_args()

    cfg = load_config(a.config)
    mc = cfg["materials"]
    out = a.out or (os.path.normpath(os.path.join(cfg["_base"], mc["resin_data"])) if mc.get("resin_data") else None)
    if out is None:
        raise SystemExit("no [materials] resin_data in the config and no --out")

    if a.part:
        if not a.ref:
            raise SystemExit("--part needs --ref")
        refs = [read_raw(p) for p in a.ref]
        if a.ref_auto:
            refs.append({f: r for f, r in read_raw(a.ref_auto).items() if r["rng"] == "6"})
        rows = correct(read_raw(a.part), refs, a.f_min, a.f_max)
        write_corrected(rows, out, [
            f"sample sweep {os.path.basename(a.part)}, corrected by resin_fit.py",
            "reference (C0G): " + ", ".join(os.path.basename(p) for p in a.ref + ([a.ref_auto] if a.ref_auto else [])),
            f"use = 1: {a.f_min:g} Hz..{a.f_max:g} Hz, no odd mains harmonic, not flagged by the sweep",
            "C_pF, tan_d: corrected; sigma_*: 1-sigma used as fit weights; *_raw: as read"])
        print(f"wrote {out}: {sum(r['use'] for r in rows)} of {len(rows)} points used")

    # --- fit and report ---------------------------------------------------
    f, c, d, sc, sd = dielectric.read_sweep(out)
    print(f"\n{len(f)} points, {f.min():g} Hz .. {f.max():g} Hz")
    print("local Kramers-Kronig check of the corrected data (a lossless series gap does not change it):")
    for lo, hi, dkk, dm in kk_check(f, c, d):
        print(f"  {lo:8.0f}-{hi:<8.0f} Hz: tan d from the slope of C {dkk:.4f}, measured {dm:.4f}, ratio {dkk/dm:.3f}")

    er, fs, td = mc["resin_er"], mc.get("resin_spec_freq", 50.0), mc.get("resin_tan_d", 0.0)
    geom = dict(area=mc.get("resin_sample_area", 0.0), thickness=mc.get("resin_sample_thickness", 0.0))
    variants = [("electrodes in full contact", "none"), ("air gap matching the datasheet tan d", "datasheet")]
    models = {}
    freqs = (10, 50, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
    for label, gap in variants:
        m = models[gap] = dielectric.measured_model(out, er, fs, td, gap, **geom)
        ft = m.fit
        print(f"\n{label}: " + (f"gap {ft['gap_pF']:.0f} pF in series" if ft["gap_pF"] else "no gap"))
        print(f"  eps_inf {m.eps_inf:.3f}, d_eps {m.d_eps:.3f}, tau {m.tau:.4g} s (loss peak {m.f_peak:.3g} Hz), "
              f"beta {m.beta:.4f}")
        print(f"  fit: chi2/dof {ft['chi2']:.2f}, rms error C {ft['rms_c_pct']:.2f} %, tan d {ft['rms_tan_d']:.4f}; "
              f"scale from {ft['scale_from']} ({ft['pF_per_eps']:.2f} pF per unit eps)")
        print("  " + "  ".join(f"{dielectric._ff(x):>6}" for x in freqs) + "  Hz")
        print("  " + "  ".join(f"{m.er(x):6.2f}" for x in freqs) + "  eps'")
        print("  " + "  ".join(f"{m.tan_d(x):6.3f}" for x in freqs) + "  tan d")
        if "t_eff_mm" in ft and "sample geometry" not in ft["scale_from"]:
            line = (f"  sample check: with {ft['area_mm2']:g} mm2 of electrode this scale stands for a flat sample "
                    f"{ft['t_eff_mm']:.2f} mm thick")
            if "t_range_mm" in ft:
                lo, hi = ft["t_range_mm"]
                inside = "inside" if lo <= ft["t_eff_mm"] <= hi else "OUTSIDE"
                line += (f": {inside} the sample's {lo:g}-{hi:g} mm "
                         f"(spread evenly over that range it would act as {ft['t_even_mm']:.2f} mm)")
            if "gap_um" in ft:
                line += f"; the gap is {ft['gap_um']:.0f} um of air on average"
            print(line)
    used = mc.get("resin_electrode_gap", "none")
    print(f"\nprobe.toml uses resin_electrode_gap = {used!r}")

    # --- plot -------------------------------------------------------------
    ff = np.logspace(1, 8, 400)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    old = dielectric.Dielectric.fit(er, td, fs, mc.get("resin_er_hf", 4.0), m1=math.log10(mc.get("resin_loss_f_lo", 0.01)))
    for n, (label, gap) in enumerate(variants):
        m = models[gap]
        cs, ds = dielectric.without_gap(c, d, m.fit["gap_pF"])
        k = 1 / m.fit["pF_per_eps"]
        axs[0].semilogx(f, cs * k, "o", ms=3.5, color=SERIES[n], mec="none")
        axs[0].semilogx(ff, m.er(ff), "-", color=SERIES[n], label=label)
        axs[1].semilogx(f, ds, "o", ms=3.5, color=SERIES[n], mec="none")
        axs[1].semilogx(ff, m.tan_d(ff), "-", color=SERIES[n], label=label)
        axs[2].semilogx(f, 100 * (m.er(f) / (cs * k) - 1), "o-", ms=3, lw=0.8, color=SERIES[n], label=f"C, {label}")
    axs[0].semilogx(ff, old.er(ff), "--", color=INK2, lw=1.2, label="old model (datasheet point, eps_hf 4 assumed)")
    axs[1].semilogx(ff, old.tan_d(ff), "--", color=INK2, lw=1.2, label="old model")
    axs[0].plot([fs], [er], "s", color=INK2, ms=6, label="datasheet point")
    if td:
        axs[1].plot([fs], [td], "s", color=INK2, ms=6, label="datasheet point")
    for ax in axs[:2]:
        ax.axvspan(f.min(), f.max(), color="#000000", alpha=0.04, lw=0)
    axs[0].set(xlabel="frequency [Hz]", ylabel="eps'", title="Permittivity (shaded: measured band)")
    axs[1].set(xlabel="frequency [Hz]", ylabel="tan d", title="Loss tangent")
    axs[2].set(xlabel="frequency [Hz]", ylabel="model / measured - 1 [%]", title="Fit residual of C")
    for ax in axs:
        ax.legend(fontsize=7)
    fig.suptitle(f"Resin model from {os.path.basename(out)}")
    fig.tight_layout()
    path = os.path.join(cfg["files"]["results"], a.name + ".png")
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
