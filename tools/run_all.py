"""Run the whole analysis pipeline after a design change.

    python run_all.py                # everything at the default grid (~10 min)
    python run_all.py --quick        # coarse grids (~2.5 min), for checking a change
    python run_all.py --skip-field   # reuse results/capacitance.json (no new field solve)
    python run_all.py --no-tube      # skip the tube-end sweep

The field solve runs first; then the tube sweep and the circuit analyses run
concurrently (--jobs). Every step writes its console output to results/<step>.txt.
The proposed divider design used for the "graded" analyses is set by PROPOSED below.
"""

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))

# the divider design proposed in the README: capacitors graded from 47 pF at
# the bottom (E12 values), 10 ohm damping resistor in series with each one
PROPOSED = ["--graded", "47p", "--round", "E12", "--seg-r", "10"]


def stages(a):
    """Steps in stages: a stage starts when the previous one has finished; the
    steps within a stage are independent and run concurrently."""
    field_scale = ["--scale", "2"] if a.quick else []
    tube_scale = ["--scale", "3", "--step", "5"] if a.quick else ["--scale", "1.5", "--step", "2.5"]
    out = []
    if not a.skip_field:
        # stress maps for the proposed design (the capacitances don't depend on it)
        out.append([("extract_caps", ["extract_caps.py"] + field_scale + PROPOSED)])
    # everything below only reads results/capacitance.json
    later = []
    if not a.no_tube:
        later.append(("tube_sweep", ["tube_sweep.py", *tube_scale]))     # longest: start it first
    later += [
        ("info", ["probe_sim.py", "info"]),
        ("optimize", ["probe_sim.py", "optimize"]),
        ("response_as_built", ["probe_sim.py", "response", "--optimize", "--name", "response_as_built"]),
        ("response_proposed", ["probe_sim.py", "response", "--optimize", *PROPOSED, "--name", "response_proposed"]),
        ("optimize_proposed", ["probe_sim.py", "optimize", *PROPOSED, "--name", "optimize_proposed"]),
        ("damping", ["probe_sim.py", "damping", *PROPOSED[:4], "--name", "damping"]),
        ("sweep", ["probe_sim.py", "sweep", "--name", "sweep"]),
        ("sweep_damped", ["probe_sim.py", "sweep", "--seg-r", "10", "--name", "sweep_damped"]),
        ("spice_check", ["probe_sim.py", "spice", "--run"]),
    ]
    out.append(later)
    return out


def run_step(name, cmd, res, env):
    t0 = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] {name}: {' '.join(cmd)}", flush=True)
    with open(os.path.join(res, name + ".txt"), "w", encoding="utf-8") as log:
        r = subprocess.run([sys.executable, *cmd], cwd=HERE, stdout=log, stderr=subprocess.STDOUT, env=env)
    print(f"[{time.strftime('%H:%M:%S')}] {name}: "
          f"{'ok' if r.returncode == 0 else 'FAILED (see results/' + name + '.txt)'} in {time.time()-t0:.0f} s",
          flush=True)
    return r.returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--skip-field", action="store_true")
    ap.add_argument("--no-tube", action="store_true")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1,
                    help="steps run at the same time (the tube sweep and the field solve also use "
                         "[solver] workers processes each)")
    a = ap.parse_args()
    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", PYTHONIOENCODING="utf-8")
    t_all = time.time()
    for stage in stages(a):
        with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as ex:
            codes = dict(zip((n for n, _ in stage), ex.map(lambda st: run_step(*st, res, env), stage)))
        if codes.get("extract_caps", 0) != 0:
            sys.exit(1)
    print(f"done in {(time.time()-t_all)/60:.1f} min")


if __name__ == "__main__":
    main()
