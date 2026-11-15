"""Where should the shield tube end? Sweeps the position of the tube's HV-side
end and reports, for every position:

  * electric-field stress: peak and bulk |E| per material at the stress voltage,
    plus the tangential field along the board surface;
  * shielding: the outer box is kept as a separate conductor "ENV" (the
    surroundings). The coupling capacitance from ENV to every divider node is
    put into the circuit model, which gives the interference at the scope as
    an equivalent tip voltage per volt on the surroundings (50 Hz and 1 MHz).

    python tube_sweep.py                          # default range, grid scale 1.5
    python tube_sweep.py --from 205 --to 245 --step 2.5 --scale 1.2
    python tube_sweep.py --maps                   # also save a stress map per position

Stress values at sharp copper edges depend on the grid, so compare positions
on the same --scale rather than reading absolute numbers. All positions are
forced grid lines, so every position uses the same grid. The tube's LV end
(x0) and everything else stays as drawn on User.1.
"""

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")   # batched small solves: ~80x faster single-threaded
import argparse
import csv
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import kicad_board
import fieldsolver as fs
import field_plots
import probe_model as pm
from extract_caps import load_config
from probe_sim import SERIES, INK2


def part_x_range(board, refs):
    xs = [x for ref, fp, layer, net, poly in board.pad_shapes if ref in refs for x, _ in poly]
    return (min(xs), max(xs)) if xs else (None, None)


_job = None


def _init(board, cfg, xs, scale, node_v, labels, maps, volts):
    global _job
    _job = dict(board=board, cfg=cfg, xs=xs, scale=scale, node_v=node_v, labels=labels, maps=maps, volts=volts)


def _solve_position(x1):
    """One tube-end position: stress report and the coupling of every net to the
    surroundings (ENV). Runs in a worker process (or in-process with 1 worker)."""
    j = _job
    t0 = time.time()
    board, cfg, labels = j["board"], j["cfg"], j["labels"]
    tube0 = board.tube
    board.tube = dict(tube0, x1=float(x1))
    # every position is a grid line of every model: the grid stays the same,
    # so the differences between positions come from the tube alone
    m = fs.build_model(board, cfg, scale=j["scale"], verbose=False, separate_env=True, x_forced=j["xs"])
    board.tube = tube0
    solver = fs.Solver(m, tol=cfg["solver"]["tol"], verbose=False)
    v = np.array([j["node_v"].get(labels.get(n), 0.0) if n not in ("GND", "ENV") else 0.0 for n in m.names])
    phi = solver.solve(v)
    rep, E = field_plots.stress_report(m, phi, cfg)
    if j["maps"]:
        field_plots.plot_stress(m, phi, E, j["maps"].format(x1), f"Tube end at x = {x1:.1f} mm", j["volts"])
    ve = np.zeros(len(m.names))
    ve[m.names.index("ENV")] = 1.0
    q = solver.charges(solver.solve(ve))
    c_env = {labels[n]: -q[k] for k, n in enumerate(m.names) if n not in ("GND", "ENV") and n in labels}
    return float(x1), rep, c_env, time.time() - t0


def pickup(probe, c_env, freqs):
    """Equivalent tip voltage per volt on the surroundings (|V_scope| / |H_signal|)."""
    net = probe._build()
    for lab, c in c_env.items():
        if c > 1e-18:
            net.add("C", "ENV", lab, c=c, ref=f"Cenv_{lab}")
            # the parasitic matrix (box merged with GND) already counts this
            # coupling as capacitance to ground: move it instead of adding it
            net.add("C", lab, "0", c=-c, ref=f"Cenv_{lab}_gnd")
    s = 2j * np.pi * np.asarray(freqs, float)
    src = list(probe.fixed)[0]
    V, _ = net.solve(s, {src: 0.0, "ENV": 1.0})
    v_out = V[:, net.idx["SCOPE"]]
    if probe.p.scope_bw > 0:        # same scope roll-off as probe.H
        v_out = v_out / (1 + s / (2 * np.pi * probe.p.scope_bw))
    return np.abs(v_out) / np.abs(probe.H(freqs))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(os.path.dirname(__file__), "probe.toml"))
    ap.add_argument("--from", dest="x_from", type=float, help="first tube end position x [mm]")
    ap.add_argument("--to", dest="x_to", type=float, help="last tube end position x [mm]")
    ap.add_argument("--step", type=float, default=4.0, help="step [mm]")
    ap.add_argument("--scale", type=float, default=1.5, help="grid scale (1 = extract_caps default)")
    ap.add_argument("--dist", choices=["dc", "hf"], default="dc",
                    help="node voltages from the DC or the fast-edge distribution")
    ap.add_argument("--voltage", type=float, help="HV tip voltage [V] (default [stress] voltage)")
    ap.add_argument("--maps", action="store_true", help="save a stress map for every position")
    ap.add_argument("--name", default="tube_sweep")
    ap.add_argument("--workers", type=int, help="parallel processes, one position each (default [solver] tube_workers)")
    a = ap.parse_args()

    cfg = load_config(a.config)
    res = cfg["files"]["results"]
    board = kicad_board.load_board(cfg["files"]["pcb"])
    capdata = pm.load_capdata(cfg["files"]["capacitance"])
    probe = pm.Probe(cfg, board, capdata, pm.params_from_config(cfg))
    labels = probe.labels
    volts = a.voltage or cfg["stress"]["voltage"]

    # conductor potentials (the distribution hardly depends on the tube end)
    f = 1e-3 if a.dist == "dc" else cfg["stress"]["hf_frequency"]
    V, _ = probe.solve(np.array([2j * np.pi * f]))
    node_v = {lab: volts * V[0, i].real for lab, i in probe.net.idx.items()}

    tube0 = dict(board.tube)
    seg1 = probe.top.get(0, []) + probe.top_r.get(0, [])
    _, x_seg1 = part_x_range(board, seg1)
    x_resin = max(p[0] for p in board.resin)
    xs = np.arange(a.x_from if a.x_from is not None else round(x_seg1 - 12),
                   (a.x_to if a.x_to is not None else round(x_resin - 14)) + 1e-6, a.step)
    print(f"HV-end divider segment ({'+'.join(seg1)}) ends at x = {x_seg1:.1f} mm, resin ends at {x_resin:.1f} mm, "
          f"tube currently ends at {tube0['x1']:.2f} mm")
    print(f"{len(xs)} positions, grid scale {a.scale}, {a.dist.upper()} distribution, {volts/1e3:g} kV")

    freqs = [50.0, 1e6]
    rows = []
    hdr = (f"{'tube end':>9}{'resin pk':>10}{'resin bulk':>11}{'FR4 pk':>8}{'FR4 bulk':>9}{'surf tan':>9}"
           f"{'  @x':>7}{'C div-ENV':>11}{'pickup 50Hz':>13}{'pickup 1MHz':>13}{'time':>6}")
    print(hdr)
    workers = max(1, min(a.workers if a.workers is not None else cfg["solver"].get("tube_workers", 1), len(xs)))
    init = (board, cfg, xs, a.scale, node_v, labels,
            os.path.join(res, a.name + "_{:.1f}.png") if a.maps else None, volts)
    if workers > 1:
        import multiprocessing as mp
        pool = mp.get_context("spawn").Pool(workers, _init, init)
        results = pool.imap(_solve_position, xs)          # in order, computed in parallel
    else:
        pool = None
        _init(*init)
        results = map(_solve_position, xs)
    for x1, rep, c_env, dt in results:
        c_div = sum(c for lab, c in c_env.items() if lab != "HV")
        pk = pickup(probe, c_env, freqs)
        hot = max(("resin_peak", "FR4_peak"), key=lambda k: rep[k]["E_kV_mm"])
        row = dict(tube_end=float(x1),
                   resin_peak=rep["resin_peak"]["E_kV_mm"], resin_bulk=rep["resin_bulk"]["E_kV_mm"],
                   fr4_peak=rep["FR4_peak"]["E_kV_mm"], fr4_bulk=rep["FR4_bulk"]["E_kV_mm"],
                   surface_tangential=rep["board_surface_tangential"]["E_kV_mm"],
                   hotspot_x=rep[hot]["x"], hotspot_y=rep[hot]["y"], hotspot_z=rep[hot]["z"],
                   c_divider_env_fF=c_div * 1e15, pickup_50Hz=pk[0], pickup_1MHz=pk[1])
        rows.append(row)
        print(f"{x1:9.1f}{row['resin_peak']:10.2f}{row['resin_bulk']:11.2f}{row['fr4_peak']:8.2f}"
              f"{row['fr4_bulk']:9.2f}{row['surface_tangential']:9.2f}{row['hotspot_x']:7.1f}"
              f"{c_div*1e15:9.2f}fF{pk[0]:13.3g}{pk[1]:13.3g}{dt:5.0f}s", flush=True)
    if pool:
        pool.close()
        pool.join()

    path = os.path.join(res, a.name + ".csv")
    with open(path, "w", newline="") as fcsv:
        w = csv.DictWriter(fcsv, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # --- plot -------------------------------------------------------------
    x = np.array([r["tube_end"] for r in rows])
    fig, axs = plt.subplots(3, 1, figsize=(11, 11), sharex=True, gridspec_kw={"height_ratios": [1.2, 1, 0.8]})
    ax = axs[0]
    for n, (key, lab) in enumerate((("resin_bulk", "resin, bulk"), ("fr4_bulk", "FR4, bulk"),
                                    ("surface_tangential", "board surface, tangential"))):
        ax.plot(x, [r[key] for r in rows], "-o", ms=4, color=SERIES[n], label=lab)
    for n, (key, lab) in enumerate((("resin_peak", "resin, at copper edges"), ("fr4_peak", "FR4, at copper edges"))):
        ax.plot(x, [r[key] for r in rows], "--", color=SERIES[n], lw=1.2, label=lab)
    ax.set(ylabel="|E| [kV/mm]", title=f"Peak field at {volts/1e3:g} kV ({a.dist.upper()} distribution, grid scale {a.scale})")
    ax.legend(fontsize=8)

    ax = axs[1]
    for n, (key, lab) in enumerate((("pickup_50Hz", "50 Hz"), ("pickup_1MHz", "1 MHz"))):
        ax.semilogy(x, [r[key] for r in rows], "-o", ms=4, color=SERIES[n + 3], label=lab)
    ax.set(ylabel="tip-referred V per V\non the surroundings",
           title="Interference pickup through the uncovered part of the divider (lower = better shielded)")
    ax.legend(fontsize=8)

    ax = axs[2]
    for loop in [board.outline] + board.holes:
        p = np.array(loop + loop[:1])
        ax.plot(p[:, 0], p[:, 1], color=INK2, lw=0.6)
    for ref, fp, layer, net, poly in board.pad_shapes:
        p = np.array(poly + poly[:1])
        ax.fill(p[:, 0], p[:, 1], color=SERIES[0] if net != probe.hv else SERIES[1], lw=0)
    for net_layer, cu in board.copper.items():
        if net_layer[0] == probe.hv:
            for poly in cu.polys:
                p = np.array(poly)
                ax.fill(p[:, 0], p[:, 1], color=SERIES[1], lw=0, alpha=0.7)
            for p1, p2, wdt in cu.lines:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=SERIES[1], lw=1.5)
    if board.resin:
        p = np.array(board.resin + board.resin[:1])
        ax.plot(p[:, 0], p[:, 1], color=SERIES[2], lw=1, ls="--", label="resin outline")
    ax.set(xlabel="tube end position x [mm]", ylabel="y [mm]", title="Board (HV net in orange)")
    ax.set_ylim(tube0["yc"] + tube0["r_in"] + 1, tube0["yc"] - tube0["r_in"] - 1)
    for axx in axs:
        axx.axvline(tube0["x1"], color=INK2, ls=":", lw=1)
        axx.axvline(x_seg1, color=SERIES[0], ls=":", lw=1)
    axs[0].text(tube0["x1"], axs[0].get_ylim()[1], " current tube end", fontsize=7, color=INK2, va="top")
    axs[0].text(x_seg1, axs[0].get_ylim()[1] * 0.93, " end of divider", fontsize=7, color=SERIES[0], va="top")
    axs[2].set_xlim(x.min() - 15, x.max() + 15)
    fig.tight_layout()
    path2 = os.path.join(res, a.name + ".png")
    fig.savefig(path2, dpi=120)
    plt.close(fig)
    print(f"wrote {path} and {path2}")


if __name__ == "__main__":
    main()
