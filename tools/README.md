# Probe modelling tools

Python tools that take the KiCad design and estimate how the real probe behaves:
- the parasitic capacitances of the divider inside the grounded, resin-filled tube;
- the frequency and step response with a real cable and scope;
- the compensation capacitance needed;
- the electric-field stress in the resin and the board.

```
hv-probe-passive.kicad_pcb ──► extract_caps.py ──► results/capacitance.json ──► probe_sim.py ──► response, compensation,
  copper, outline, parts,      3D electrostatic       Maxwell capacitance        circuit model     design sweeps, LTspice
  tube (User.1),               field solver           matrix (23 conductors)                       netlist
  resin (User.2)               └──► E-field stress maps (stress_dc.png, stress_hf.png)
```

## Files

| file | what it does |
|---|---|
| `probe.toml` | **All assumptions in one place**: materials, part heights, mesh, ESL values, cable, scope, sweep ranges |
| `kicad_board.py` | Reads the `.kicad_pcb` file directly (no KiCad API needed). Run it on its own to see what was imported |
| `fieldsolver.py` | 3D finite-volume electrostatic solver (graded grid, algebraic multigrid) |
| `extract_caps.py` | Runs the solver: capacitance matrix, model check plots, E-field stress |
| `field_plots.py` | Plots for the solver |
| `probe_model.py` | Circuit model: nodal analysis, lossy transmission line, step response, optimiser |
| `probe_sim.py` | Command line front end for the circuit analyses and the LTspice export |
| `sexpr.py` | S-expression reader for KiCad files |
| `legacy/cap_calc.py` | The original closed-form estimate (coaxial-cylinder formula), kept for reference |
| `results/` | Outputs of the last run |

## Install

Python ≥ 3.11 is needed (for `tomllib`).

```
pip install -r requirements.txt
```

`pyamg` is the only extra dependency. If you have an older matplotlib that needs numpy < 2, pin the versions, e.g. `pip install "numpy<2" "scipy<1.12" "pyamg==5.2.1"`.

## Workflow

All commands are run from this directory. Results are written to `results/`.

### 1. Check the import

```
python kicad_board.py ../hv-probe-passive/hv-probe-passive.kicad_pcb
```

This lists the nets, the divider chain (HV, N01…N19, OUT), the parts, the tube (axis, ID, wall, extent) and the resin outline. It's worth running after every PCB change.

The importer relies on these drawing conventions in the PCB:
- **User.1** holds two rectangles: the upper and lower tube wall in the axial cross-section. The tube axis is midway between them.
- **User.2** holds the closed outline of the resin body in the same view. It is revolved around the tube axis.
- The **board is centred in the tube**, i.e. the PCB mid-plane contains the axis. Change `board_offset_z` if it isn't.

### 2. Field solve: capacitances and E-field stress

```
python extract_caps.py --scale 2      # coarse, ~2.5 min, 1.3 M nodes
python extract_caps.py                # default grid, ~13 min, 6.5 M nodes, ~3 GB RAM
```

Open `results/model_check.png` first. It shows how the copper, part terminations, board, slots, resin and tube were discretised. Then read `results/capacitance_summary.txt`, which lists each node's parasitic capacitance to the shield, to its neighbours and to the HV tip. The full matrix goes to `results/capacitance.json`.

The same run produces the **E-field stress maps** at the `[stress] voltage` set in the config (20 kV by default):
- `stress_dc.png`: node voltages from the resistive (DC / slow) distribution.
- `stress_hf.png`: node voltages from the capacitive distribution a fast edge sees, including the parasitics.
- `stress.json`: peak |E| per material with locations, the tangential field along the board surface, and the node voltages used.

Use `--stress-only` to redo just the stress part from an existing `capacitance.json`, and `--voltage` to change the applied voltage.

### 3. Circuit simulation

```
python probe_sim.py info                       # summary, DC ratio, required compensation capacitance
python probe_sim.py response --optimize        # response.png: frequency/step response, voltage per segment, Z_in
python probe_sim.py optimize                   # compensation capacitor per cable length (sweep.cable_lengths)
python probe_sim.py sweep                      # design study: top capacitor value / grading vs. DUT loading
python probe_sim.py spice --run                # LTspice netlist (results/probe.cir) + cross-check
```

Options override the config for a single run:

| option | meaning |
|---|---|
| `--length 2 --cable RG174` | cable length [m] and type |
| `--series-r 47` | back-termination resistor between the divider output and the cable |
| `--trim 30p` | trimmer setting |
| `--top-cap 100p` | replace every divider capacitor by this value |
| `--graded 47p --round E12` | grade the divider capacitors (see below), starting at 47p at the bottom |
| `--tip-r 100 --seg-r 10` | damping resistor at the tip / in series with every divider capacitor |
| `--scope-bw 100Meg --c-scope 20p` | scope bandwidth (1st order) and input capacitance |
| `--no-parasitics` | ideal divider, for comparison |
| `--optimize` | first find the bottom capacitance that makes the response flattest |
| `--name x` | output file name |

What `optimize` reports:
- **The required fixed + trimmer capacitance** for each cable length. This is the total of C34 + C35 + C13 that gives the flattest response between `f_lo` and `f_hi` (10 Hz to 1 MHz by default).
- **Combinations of standard values** that put the trimmer near mid-range.

## How it works

### Field solver

Electrostatics, ∇·(ε∇φ) = 0, solved with the finite integration technique. The grid is a non-uniform rectilinear grid that is fine near the board (0.15 mm through the thickness, 0.4 mm in-plane) and grows geometrically towards a grounded box 30 mm outside the tube.

**Conductors**, each held at a fixed potential:
- copper polygons from the zone fills, pads, tracks and vias, as zero-thickness sheets on the F.Cu and B.Cu planes;
- plated holes;
- SMD terminations, as conductive boxes of the height given in `[geometry.part_height]`;
- the tube wall.

**Dielectrics:** FR4 for the board (slots included), resin inside the revolved User.2 outline, and air everywhere else. That covers the unpotted LV end, the region around the exposed tip, and outside the tube.

**Capacitance matrix:** each conductor is set to 1 V in turn with all the others at 0 V. The charges on every conductor give one column of the Maxwell capacitance matrix. A single AMG hierarchy is reused for all 22 solves.

**Why not BEM?** A boundary-element solver is the natural choice for conductors in one homogeneous dielectric. This problem has three dielectrics (resin, FR4, air), a dielectric boundary at the open end and at the unpotted LV end, and thin conductors on both sides of a board. A volume method handles all of that without extra effort, and its field solution doubles as the stress map. [FasterCap](https://www.fastfieldsolvers.com/fastercap.htm) (free, multi-dielectric BEM) and [Elmer FEM](https://www.elmerfem.org/) (free; its `StatElecSolver` outputs a capacitance matrix) are good independent cross-checks if you export the geometry.

**Accuracy.** The node-based discretisation slightly shrinks conductors, so values converge from below:

| grid | nodes | node to shield | across a segment | HV tip to ground |
|---|---|---|---|---|
| `--scale 2` | 1.3 M | 1.13–1.19 pF | 0.25–0.34 pF | 2.76 pF |
| default | 6.5 M | 1.21–1.28 pF | 0.34–0.40 pF | 2.75 pF |

Expect the default grid to be within about 5–10 % of the true values for the modelled geometry. That error is smaller than the uncertainty in the inputs: resin εr against temperature and frequency, real part geometry, how well the board is centred, and bubbles.

### Circuit model

The circuit is solved by nodal analysis at complex frequency s. It contains:
- **Discrete parts from the PCB:** 20 × (2 MΩ ∥ 47 pF), the 56k/56k/RV2 bottom resistors, C34, C35 and trimmer C13. RV2 is set automatically so that the DC ratio, including the scope's 1 MΩ, is exactly 1000:1.
- **ESL on all capacitors** (`[circuit.esl]`).
- **Every parasitic capacitance from the matrix:** node to shield, node to node and node to tip.
- **The cable as an exact lossy transmission line**, with skin-effect resistance ∝ √f fitted to the datasheet loss at 100 MHz.
- **The scope** as 1 MΩ ∥ C_in, with an optional bandwidth limit.
- Optional source, tip and segment damping resistors, and a back-termination resistor.

Step responses use a numerical inverse Laplace transform: damped FFT with a Hanning window. That gives clean results for both the nanosecond cable reflections and the millisecond RC tails. The step metrics use a Gaussian 1 ns source edge (`[analysis] source_rise`).

The LTspice export contains exactly the same network. `spice --run` runs LTspice in batch mode and compares |H| at a few frequencies. It agrees to < 0.01 dB. Near cable resonances there are small differences, because LTspice's LTRA line has no skin effect.

**Valid range.** The model is quasi-static and lumped. It is good up to roughly 100–200 MHz. Above that, the 15 cm divider behaves as a distributed structure. The LC resonances the model shows there are real in principle (ESL of the capacitor string against the shield capacitance), but their frequency and Q depend on ESL values that are assumptions, and on losses the model doesn't include. Treat anything above ~150 MHz as qualitative.

### Graded divider (`--graded`)

`probe_model.graded_caps` computes the segment capacitances that make the fast-edge voltage distribution exactly linear, given the parasitic matrix. It starts from the bottom segment and works upward. Each segment must carry the current of the segment below it plus all parasitic currents leaving its node. `--round E12` snaps the result to standard values.

### Stress maps

The field is computed with every divider node at its real potential, using both the DC and the fast-edge distribution from the circuit model. |E| is evaluated per cell and reported per material, with and without the cells touching a conductor.

Copper edges are mathematically sharp, so the "peak" values there grow with grid refinement: 9 → 13 kV/mm in resin between the two grids. Use them to rank and locate hotspots. The "bulk" values, and the comparison between design variants on the same grid, are the more meaningful numbers.

The air field around the exposed tip depends on the brass terminal (not modelled) and on the distance to grounded surroundings: the 30 mm box.

## Limitations / not modelled

- Part bodies are treated as resin, and only the terminations are conductive. Leads of the THT parts at the LV end are simplified to pad columns.
- Copper thickness is ignored.
- Real resin εr: check the datasheet value at your frequency and temperature, and set `materials.resin_er`. All parasitics scale roughly with it.
- No dielectric loss, so resonances are underdamped.
- The HV resistors' own body capacitance is included only through their terminations.
- The brass tip terminal and the fuse are not in the model.
