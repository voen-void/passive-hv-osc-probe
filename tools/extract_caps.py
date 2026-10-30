"""Extract the parasitic capacitance matrix of the probe from the KiCad PCB
with a 3D electrostatic field solver, and optionally map the E-field stress.

    python extract_caps.py                 # full run with probe.toml
    python extract_caps.py --scale 2       # quick, coarse grid
    python extract_caps.py --stress-only   # reuse capacitance.json, redo stress maps

Outputs (in the results directory):
    capacitance.json          Maxwell capacitance matrix + metadata
    capacitance_summary.txt   human readable per-node summary
    model_check.png           how the PCB was imported (check this first!)
    grid.png                  grid step sizes
    stress_dc.png / stress_hf.png and stress.json   E-field maps and peak values
"""

import argparse
import json
import os
import time
import tomllib

import numpy as np

import kicad_board
import fieldsolver as fs
import field_plots
import probe_model


def load_config(path):
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    base = os.path.dirname(os.path.abspath(path))
    for k in ("pcb", "capacitance", "results"):
        cfg["files"][k] = os.path.normpath(os.path.join(base, cfg["files"][k]))
    os.makedirs(cfg["files"]["results"], exist_ok=True)
    return cfg


def summary_text(capdata, labels, ladder):
    names = capdata["names"]
    M = np.array(capdata["maxwell"])
    ix = {n: i for i, n in enumerate(names)}
    cp = lambda a, b: -M[ix[a], ix[b]] * 1e12
    lines = ["Parasitic capacitances [pF] (excluding the discrete parts)", ""]
    lines.append(f"{'node':>6} {'to GND':>8} {'to prev':>8} {'to next':>8} {'to HV':>8} {'to rest':>8}")
    for k, net in enumerate(ladder):
        prev = ladder[k - 1] if k > 0 else None
        nxt = ladder[k + 1] if k + 1 < len(ladder) else None
        g = cp(net, "GND")
        p_ = cp(net, prev) if prev else 0
        n_ = cp(net, nxt) if nxt else 0
        hv = cp(net, ladder[0]) if k > 1 else 0
        rest = sum(cp(net, o) for o in names if o not in (net, prev, nxt, "GND", ladder[0]))
        lines.append(f"{labels[net]:>6} {g:8.3f} {p_:8.3f} {n_:8.3f} {hv:8.3f} {rest:8.3f}")
    tot_g = sum(cp(n, "GND") for n in ladder[1:-1])
    lines += ["", f"sum of internal node capacitances to ground: {tot_g:.2f} pF",
              f"HV tip to ground (shield + surroundings):    {cp(ladder[0], 'GND'):.2f} pF"]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(os.path.dirname(__file__), "probe.toml"))
    ap.add_argument("--scale", type=float, default=1.0, help="grid step multiplier (2 = coarse/fast)")
    ap.add_argument("--resin-er", type=float, help="override resin permittivity")
    ap.add_argument("--out", help="capacitance json path (default from config)")
    ap.add_argument("--no-stress", action="store_true", help="skip the E-field stress analysis")
    ap.add_argument("--stress-only", action="store_true", help="reuse the existing capacitance json")
    ap.add_argument("--voltage", type=float, help="HV tip voltage for the stress maps [V]")
    a = ap.parse_args()

    cfg = load_config(a.config)
    if a.resin_er:
        cfg["materials"]["resin_er"] = a.resin_er
    res = cfg["files"]["results"]
    out = a.out or cfg["files"]["capacitance"]

    board = kicad_board.load_board(cfg["files"]["pcb"])
    cc = cfg["circuit"]
    labels = kicad_board.net_labels(board, cc["input_net"], cc["output_net"])
    ladder = kicad_board.ladder_order(board, cc["input_net"], cc["output_net"])

    t0 = time.time()
    m = fs.build_model(board, cfg, scale=a.scale)
    field_plots.plot_model(m, os.path.join(res, "model_check.png"), labels)
    field_plots.plot_grid(m, os.path.join(res, "grid.png"))
    solver = fs.Solver(m, tol=cfg["solver"]["tol"])

    if a.stress_only:
        capdata = probe_model.load_capdata(out)
    else:
        C = solver.capacitance_matrix()
        capdata = {
            "names": m.names,
            "labels": {n: labels.get(n, n) for n in m.names},
            "maxwell": C.tolist(),
            "units": "F",
            "meta": {
                "pcb": os.path.relpath(cfg["files"]["pcb"], os.path.dirname(out)),
                "grid": list(m.shape), "scale": a.scale,
                "resin_er": cfg["materials"]["resin_er"], "fr4_er": cfg["materials"]["fr4_er"],
                "board_offset_z": cfg["geometry"].get("board_offset_z", 0.0),
                "date": time.strftime("%Y-%m-%d %H:%M"), "runtime_s": round(time.time() - t0, 1),
            },
        }
        with open(out, "w") as f:
            json.dump(capdata, f, indent=1)
        txt = summary_text(capdata, labels, ladder)
        with open(os.path.join(res, "capacitance_summary.txt"), "w") as f:
            f.write(txt + "\n")
        print(txt)
        print(f"\nwrote {out}")

    if a.no_stress:
        return

    # --- E-field stress: conductor potentials from the circuit model --------
    volts = a.voltage or cfg["stress"]["voltage"]
    probe = probe_model.Probe(cfg, board, capdata, probe_model.params_from_config(cfg))
    report = {"voltage": volts}
    for tag, f, title in (("dc", 1e-3, "DC / slow voltage distribution"),
                          ("hf", cfg["stress"]["hf_frequency"],
                           f"Fast-edge voltage distribution ({cfg['stress']['hf_frequency']/1e6:g} MHz)")):
        V, _ = probe.solve(np.array([2j * np.pi * f]))
        v = np.zeros(len(m.names))
        for k, name in enumerate(m.names):
            lab = labels.get(name)
            if name != "GND" and lab in probe.net.idx:
                v[k] = volts * V[0, probe.net.idx[lab]].real
        print(f"stress solve ({tag})")
        phi = solver.solve(v)
        rep, E = field_plots.stress_report(m, phi, cfg)
        rep["node_voltages"] = {labels.get(n, n): float(v[k]) for k, n in enumerate(m.names)}
        report[tag] = rep
        field_plots.plot_stress(m, phi, E, os.path.join(res, f"stress_{tag}.png"), title, volts)
        print(f"  {'region':<28}{'|E| kV/mm':>10}{'of strength':>12}   location x,y,z [mm]")
        for k, r in rep.items():
            if k == "node_voltages":
                continue
            u = f"{100*r['utilisation']:.0f} %" if "utilisation" in r else ""
            print(f"  {k:<28}{r['E_kV_mm']:10.2f}{u:>12}   {r['x']:.1f}, {r['y']:.1f}, {r['z']:.1f}")
    with open(os.path.join(res, "stress.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(f"total time {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
