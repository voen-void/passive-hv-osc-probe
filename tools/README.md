# Probe modelling tools

Python tools that read the KiCad design and predict how the real probe behaves:
- **Parasitics:** the capacitances of the divider inside its grounded, resin-filled tube, from a 3D field solver.
- **Response:** frequency and step response with a real cable and scope, and the loading on the DUT.
- **Compensation:** the bottom capacitance needed, and a graded divider design.
- **Stress and shielding:** the electric-field stress in resin, FR4 and air, and how far the shield tube should extend.
- **Images:** README images regenerated straight from KiCad.

```
hv-probe-passive.kicad_pcb ─► extract_caps.py ─► results/capacitance.json ─► probe_sim.py ─► response, compensation,
  copper, outline, parts       3D field solver     capacitance matrices        circuit model    grading, damping, sweeps,
  shield (User.1)                │                 (total + per material)                       LTspice netlist
  resin  (User.2)                └─► stress maps                               tube_sweep.py ─► stress vs. shielding
```

## Quick start

```
pip install -r requirements.txt
python run_all.py            # whole pipeline, ~10 min; results/*.png, *.csv, *.txt
python run_all.py --quick    # coarse grids, ~2.5 min, for checking a design change
python render_kicad.py       # schematic + PCB images for the README (../images/kicad/)
```

Python ≥ 3.11 is required (for `tomllib`). The only non-standard dependency is `pyamg`. If your matplotlib needs numpy < 2, pin the versions: `pip install "numpy<2" "scipy<1.12" "pyamg==5.2.1"`.

All assumptions live in **`probe.toml`**: materials, part heights, ESL/ESR, cable, scope, sweep ranges.

## Files

| file | what it does |
|---|---|
| `run_all.py` | Runs the whole pipeline; each step logs to `results/<step>.txt` |
| `render_kicad.py` | Regenerates schematic PDF/SVG, layout SVG, assembly drawings and 3D renders with `kicad-cli` |
| `kicad_board.py` | Reads the `.kicad_pcb` directly. Run it alone to see what was imported |
| `extract_caps.py` | 3D field solve: capacitance matrices, model-check plots, E-field stress maps |
| `tube_sweep.py` | Moves the tube end: stress against shielding |
| `probe_sim.py` | Circuit analyses: `info`, `response`, `optimize`, `sweep`, `damping`, `spice` |
| `fieldsolver.py`, `field_plots.py` | Solver and its plots |
| `probe_model.py` | Circuit model: nodal analysis, transmission line, inverse Laplace, optimisers, grading |
| `dielectric.py` | Causal wideband (Djordjevic–Sarkar) model of resin and FR4 |
| `sexpr.py` | S-expression reader for KiCad files |
| `legacy/cap_calc.py` | The original closed-form coaxial-cylinder estimate, for reference |

## What the PCB has to contain

The tools read everything from the `.kicad_pcb`. They rely on these conventions:

| layer / item | meaning |
|---|---|
| **User.1** | Axial cross-section of all grounded metal around the board. Two long rectangles, one on each side of the axis, are the **tube wall**; the axis is midway. Any other closed shape (the **end cap**, a flared or rolled rim) is extra metal. Everything is revolved around the axis. |
| **User.2** | Closed outline of the **resin body** in the same view, revolved around the axis. |
| `input_net` | `IN` in `probe.toml`. KiCad stores local labels as `/IN`; both forms work. |
| `output_net` | `OUT`, the BNC centre pin. |
| board position | The PCB mid-plane contains the tube axis (`board_offset_z` otherwise). |

**How the divider is found.** The tools start at `IN` and follow the resistor chain while every node has exactly two resistors. The first node where the chain branches is the divider bottom, `BOT`. The nodes are labelled HV, N01…N19, BOT, OUT.

**Other parts:**
- A resistor between BOT and OUT (R29) is the back termination.
- RV2 is set automatically so that the DC ratio, with the scope's 1 MΩ, is exactly 1000:1.
- The trimmer is recognised from its value range ("5.5pF-65pF").
- Other two-terminal parts (the TVS D1) are included through `[circuit.part_capacitance]`.

`python kicad_board.py ../hv-probe-passive/hv-probe-passive.kicad_pcb` prints all of this. `python probe_sim.py info` warns about parts without a model.

## The steps

### Field solve: `extract_caps.py`

```
python extract_caps.py --scale 2   # coarse, ~1 min
python extract_caps.py             # default grid, 6.5 M nodes, ~5 min with 2 workers, ~6 GB RAM peak
```

Check `results/model_check.png` first. It shows how copper, part terminations, board, slots, resin, tube and end cap were discretised.

Outputs:
- `results/capacitance_summary.txt`: each node's capacitance to the shield, its neighbours and the tip.
- `results/capacitance.json`: the full matrices.
- **E-field stress maps** at `[stress] voltage` (20 kV): `stress_dc.png` uses the resistive (slow) voltage distribution, `stress_hf.png` the distribution a fast edge sees. `stress.json` lists peak |E| per material with locations.

`--stress-only` redoes just the stress part; `--voltage` changes the voltage. The divider is compensated first; `--graded 47p --round E12 --seg-r 10` maps the stress for that design instead of the PCB's capacitors (`run_all.py` does this for the proposed design). `--resin-er` changes the resin permittivity the solve uses (it sets `resin_er_hf`).

### Circuit analyses: `probe_sim.py`

```
python probe_sim.py info                         # summary, dielectric model, required compensation C
python probe_sim.py response --optimize          # response.png (see below)
python probe_sim.py optimize                     # compensation capacitor per cable length
python probe_sim.py sweep                        # divider capacitor value / grading vs. DUT loading
python probe_sim.py damping --graded 47p --round E12   # damping resistors for one divider design
python probe_sim.py spice --run                  # LTspice netlist (results/probe.cir) + cross-check
```

`response` plots six panels:
- frequency response;
- compensation flatness up to 10 MHz;
- fast and slow step response;
- voltage share of every segment (DC vs. fast edge);
- input impedance (DUT loading).

`optimize` reports the C34 + C35 + C13 total that gives the flattest response from 10 Hz to 1 MHz, and standard-value combinations that put the trimmer mid-range.

Options override `probe.toml` for one run:

| option | meaning |
|---|---|
| `--length 2 --cable RG174` | cable length [m] and type |
| `--series-r 47` | back-termination resistance; replaces R29's value |
| `--trim 30p` | trimmer setting |
| `--top-cap 100p` | every divider capacitor set to this value |
| `--graded 47p --round E12` | graded divider starting from 47 pF at the bottom, rounded to E12 (see below) |
| `--seg-r 10 --tip-r 100` | damping resistor in series with every divider capacitor / with the tip |
| `--scope-bw 100Meg --c-scope 20p` | scope bandwidth (1st order) and input capacitance |
| `--no-parasitics`, `--no-esr`, `--no-loss` | leave out the parasitics / capacitor ESR / dielectric dispersion and loss |
| `--optimize` | first find the bottom capacitance that makes the response flattest |
| `--name x` | output file name |

### Tube end: `tube_sweep.py`

```
python tube_sweep.py --from 205 --to 250 --step 2.5     # grid scale 1.5 by default
```

This moves the HV-side end of the tube and re-solves at each position. It reports:
- **Stress:** peak and bulk |E| per material, and the tangential field along the board.
- **Shielding:** the surroundings (the outer box) are kept as a separate conductor. Their coupling to every node goes through the circuit model, giving an *equivalent tip voltage per volt on the surroundings* at 50 Hz and 1 MHz.

Compare positions on the same `--scale`, since peaks at copper edges depend on the grid.

### Images: `render_kicad.py`

This uses `kicad-cli` (KiCad ≥ 8; the newest install is picked, or set `[kicad] cli` in `probe.toml`). It writes to `images/kicad/`:
- `schematic.pdf`, `schematic.svg`;
- `layout.svg`: copper, outline, shield and resin;
- `assembly-top.pdf`, `assembly-bottom.pdf`;
- `pcb-top.png`, `pcb-bottom.png`, `pcb-iso.png`: raytraced, a few minutes; `--fast` for quick previews.

## How it works

### Field solver

Electrostatics, ∇·(ε∇φ) = 0, is solved with the finite integration technique on a graded rectilinear grid. The grid is 0.15 mm through the board and 0.4 mm in-plane, growing towards a grounded box 30 mm outside the tube.

**Conductors:**
- copper from zone fills, pads, tracks and vias, as zero-thickness sheets on F.Cu and B.Cu;
- plated holes;
- SMD terminations, as boxes of `[geometry.part_height]`;
- all User.1 metal.

**Dielectrics:** FR4 in the board outline (slots included), resin in the revolved User.2 outline, and air elsewhere (the unpotted LV end, around the tip, outside the tube). On the rectangular grid the tube wall is a staircase. Inside the tube the resin is extended into the cells between the ideal circle and that staircase, so no artificial air layer forms against the shield.

**Capacitance matrix.** Each conductor is set to 1 V in turn; the charges give one column of the Maxwell matrix. One AMG hierarchy (pyamg) is reused for all 23 solves.

**Per-material matrices.** Each element is an energy integral, C_ij = Σ_m ε_m ∫_m ∇φ_i·∇φ_j dV, so it splits exactly into its resin, FR4 and air parts. These are saved as `partial` in `capacitance.json`; their sum equals the total to ~5e-7 (the CG tolerance, `[solver] tol = 1e-6`). The circuit model scales each part by its material's ε(s). That makes **dielectric parameters changeable without a new field solve** (except `solve_freq`). About 83 % of a node's capacitance to ground sits in resin, 17 % in FR4.

**Why not BEM?** A boundary-element solver is the natural choice for conductors in one homogeneous dielectric. Here there are three dielectrics, open dielectric boundaries, and thin conductors on both sides of a board. A volume method handles that directly, and its field solution doubles as the stress map. [FasterCap](https://www.fastfieldsolvers.com/fastercap.htm) (multi-dielectric BEM) and [Elmer FEM](https://www.elmerfem.org/) (`StatElecSolver` outputs a capacitance matrix) are good independent cross-checks.

**Accuracy.** The node-based discretisation slightly shrinks conductors, so values converge from below. On an earlier revision of this board:

| grid | nodes | node to shield | across a segment | tip to ground |
|---|---|---|---|---|
| `--scale 2` | 1.3 M | 1.13–1.19 pF | 0.25–0.34 pF | 2.76 pF |
| default | 6.5 M | 1.21–1.28 pF | 0.34–0.40 pF | 2.75 pF |

The default grid should be within 5–10 % for the modelled geometry. That's below the uncertainty in the inputs: resin permittivity, real part geometry, board centring, voids.

### Dielectrics

Loss and dispersion are linked (Kramers–Kronig): a material with a high loss tangent must have a permittivity that falls with frequency. Each dielectric is therefore a **causal wideband (Djordjevic–Sarkar) model**:

  ε(s) = ε∞ + k·log10((ω2 + s)/(ω1 + s))

It is fitted to one datasheet point (`*_er`, `*_tan_d` at `*_spec_freq`) and, for the resin, an assumed high-frequency permittivity `resin_er_hf` that sets where the loss band ends.

For Arathane CW 5620 (ε′ = 6, tanδ = 0.11 at 50 Hz) with ε′ → 4 at high frequency:

| | 50 Hz | 1 kHz | 10 kHz | 100 kHz | ≥ 1 MHz |
|---|---|---|---|---|---|
| ε′ | 6.0 | 4.7 | 4.1 | 4.0 | 4.0 |
| tanδ | 0.11 | 0.12 | 0.05 | 0.006 | < 0.001 |

**The high-frequency permittivity of 4 is an assumption.** It's typical for cured PU, but no data was found for this resin. It matters: the parasitics above ~10 kHz scale with it. A measurement on a cast sample (LCR meter, 1 kHz–1 MHz) would pin it down.

Below the datasheet point the fit keeps the same loss down to `resin_loss_f_lo` (0.01 Hz), so ε′ rises further: 6.7 at 10 Hz (the bottom of the compensation band) and 9.5 at DC. That extrapolation is an assumption too.

The field solve uses the permittivities at `solve_freq` (1 MHz). The circuit scales every parasitic per material by ε(s)/ε(solve_freq). This scaling is exact if both materials change by the same factor, and first-order otherwise. Checked on the `--scale 2` grid with the resin at εr = 6.7 (10 Hz) instead of 4: the scaled matrix overestimates the summed node-to-shield capacitance by 1.3 % and the tip-to-ground capacitance by 2.1 %; the capacitances across the segments are within 0.1 %. For the self capacitances (diagonal) the error can only be an overestimate, because the old field is not the energy minimum for the new permittivities. Near 10 Hz it is about as large as the remaining ~1 % gain error, so a second field solve at the low-frequency permittivity, interpolated against this one, would be the next refinement.

### Circuit model

Nodal analysis at complex frequency s. It contains:
- all discrete parts from the PCB, with ESL and ESR per footprint;
- the parasitic capacitance matrices, per material, with dispersion and loss;
- the back termination R29 and the TVS capacitance;
- the cable as an exact lossy transmission line (skin effect ∝ √f fitted to the 100 MHz loss);
- the scope as 1 MΩ ∥ C_in, with an optional bandwidth limit;
- optional tip and segment damping resistors.

**Step responses** come from a numerical inverse Laplace transform (damped FFT with a window). This handles nanosecond cable reflections and millisecond RC tails in one framework. Step metrics use a Gaussian 1 ns source edge (`[analysis] source_rise`).

**Valid range.** The model is quasi-static and lumped, good up to roughly 100–200 MHz. The LC resonances it shows above that come from the capacitor ESL against the stray capacitance to the shield. They are real in principle, but their frequency and Q depend on ESL values and losses. Treat anything above ~150 MHz as qualitative.

**LTspice.** `spice --run` writes the same network (ESR as `Rser`, ESL as `Lser`; dispersion is left out) and compares it with a matching probe_sim run. They agree to < 0.01 dB below 1 MHz; the small differences near cable resonances come from LTspice's LTRA line having no skin effect.

**Performance note.** The circuit solves are thousands of tiny (≈27×27) complex systems. Multi-threaded OpenBLAS is ~80× slower on these than single-threaded, so the scripts set `OPENBLAS_NUM_THREADS=1`.

**Where the time goes.** A full `run_all.py` takes ~10 min on an 8-core laptop (Ryzen 7 4700U), down from ~31 min before these optimisations:

- **Field solve (~5 min).** 23 conductor solves, AMG-preconditioned CG at ~1.3 s per iteration, single-threaded. They run in `[solver] workers` processes, each with its own AMG hierarchy. The stress maps are linear combinations of the 23 unit solutions, so they need no further solves; the per-material matrices are formed in banded row chunks. The CG tolerance of 1e-6 changes the matrix by < 0.01 fF against 1e-8.
- **Tube sweep (~5.5 min).** One model and AMG hierarchy per tube position, in `[solver] tube_workers` processes. It runs alongside the circuit analyses.
- **Circuit analyses (~1.2 min together).** They run concurrently. The bottom-capacitance optimiser solves the network only twice per call and gets the response for every trial capacitance by superposition (see `optimize_bottom`), which makes the grading and rounding searches cheap. Most of the remaining time is in the step responses (2¹⁴ frequency points each).

Memory, not cores, limits the parallelism on the default grid; see `[solver]` in `probe.toml`.

### Graded divider and its optimiser

`graded_caps` computes the segment capacitances that make the fast-edge voltage distribution linear. It starts from the bottom segment and works up: each segment carries the current of the one below plus all parasitic currents leaving its node.

With a dispersive resin the parasitics differ across the compensation band. `optimize_grading` therefore scales the grading strength (κ) and the bottom capacitance together for the flattest response, then rounds to the chosen E-series and re-optimises the bottom capacitance.

**Objective.** "Flattest" means the smallest *peak* gain error over `[optimize] f_lo`–`f_hi`, the same number that is reported as the LF error. The bottom capacitance is found by a bounded scalar search; the peak error is quasi-convex in it, so the search finds the minimax optimum.

The peak is often set by a point the bottom capacitance cannot move. In the proposed design that point is the −0.96 % dip at 166 Hz, which comes from the resin dispersion. The peak is then nearly flat: within 1 % of its minimum for any bottom capacitance from 2547 to 2556 pF. The strict minimum, 2548 pF, lies near the end of that range where the rest of the response is worst. So among the capacitances within `[optimize] peak_tolerance` (1 %) of the minimum peak, the one with the smallest RMS error over the same log-spaced points is used:

| bottom C chosen for | C | peak | RMS | step error 0.5–300 µs | error at 10 MHz |
|---|---|---|---|---|---|
| smallest peak only | 2548 pF | 0.96 % | 0.62 % | 0.73 % | 1.71 % |
| **smallest RMS within 1 % of the smallest peak** | 2556 pF | 0.97 % | 0.55 % | 0.53 % | 1.40 % |
| smallest RMS only | 2560 pF | 1.08 % | 0.54 % | 0.62 % | 1.26 % |

The trimmer controls the first ~300 µs of the step response, about three time constants of the bottom RC. The tie-break lands where minimising the error in that window would (2555 pF), roughly where trimming by eye ends up. The largest step error, a −0.9 % sag around 1 ms, is the resin's slow polarisation, and no bottom capacitance changes it. `peak_tolerance = 0` gives the strict minimax.

The grading and rounding searches compare designs by the strict minimax error. RMS alone is a poor design objective because it rewards any small average gain. With RMS as the objective, the rounding put four capacitors one E12 step higher, which lowered the RMS error only from 0.54 to 0.53 %. The peak error rose from 0.96 to 1.07 % and the −3 dB bandwidth fell from 72 to 69 MHz. The fast-edge distribution also became less even: the worst segment took 1.32× its share instead of 1.26×.

**Rounding.** Each capacitor is rounded down or up, not simply to the nearest value. Starting from nearest rounding, one capacitor at a time is switched to the series value on the other side of its unrounded value, as long as that lowers the peak error. Every value stays within one series step of the computed grading, so the fast-edge distribution stays close to linear. For the proposed design this changes two capacitors and lowers the error from 1.00 % to 0.96 %. Most of the remaining ~1 % comes from the resin dispersion, which no fixed set of capacitors can follow.

### Stress maps

The field is computed with every node at its real potential, for both the DC and the fast-edge distribution. The potentials come from the circuit model of the compensated divider: the PCB's capacitors, or with `--graded/--round/--seg-r` the graded design (`run_all.py` passes the proposed design). |E| is reported per material, with and without cells touching a conductor.

Sharp copper edges make the "peak" values grid-dependent. Use them to locate and rank hotspots; the "bulk" values and same-grid comparisons are the meaningful numbers. The air field around the exposed tip depends on the brass terminal (not modelled) and on the distance to grounded surroundings (the 30 mm box).

The maps are electrostatic (permittivity-weighted). Under steady DC the field divides between resin, FR4 and air by conductivity instead, which isn't modelled. The utilisation against bulk breakdown strength is not a lifetime criterion: partial discharge in voids and at interfaces is (see [Reading the stress numbers](../README.md#reading-the-stress-numbers)).

`tube_sweep.py` forces every swept tube-end position as a grid line, so all positions share one grid and the differences between them come from the tube alone, not from the copper snapping differently.

## Limitations

- Part bodies are treated as resin; only terminations are conductive. THT parts at the LV end are simplified to pad columns. Copper thickness is ignored.
- The resin's high-frequency permittivity is assumed (see above).
- ESR is frequency-independent, and radiation loss is ignored, so the resonances above ~150 MHz are probably underdamped in the model.
- The brass tip terminal and the input fuse are not modelled.
