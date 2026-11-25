# Probe modelling tools

Python tools that read the KiCad design and predict how the real probe behaves:
- **Parasitics:** the capacitances of the divider inside its grounded, resin-filled tube, from a 3D field solver.
- **Response:** frequency and step response with a real cable and scope, and the loading on the DUT.
- **Compensation:** the bottom capacitance needed, and a graded divider design.
- **Stress and shielding:** the electric-field stress in resin, FR4 and air, and how far the shield tube should extend.
- **Copper pours:** the zones that connect each divider step, drawn into the board from its pads.
- **Resin:** a dielectric model fitted to a bridge measurement of a cast sample.
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
python make_pours.py         # after moving divider parts: redraw the pours and refill the zones (KiCad closed)
python run_all.py            # whole pipeline; results/*.png, *.csv, *.txt
python run_all.py --quick    # coarse grids, for checking a design change
python render_kicad.py       # schematic + PCB images for the README (../images/kicad/)
python draw_schematics.py    # circuit diagrams for the README (../images/schematics/)
```

Python ≥ 3.11 is required (for `tomllib`). The only non-standard dependencies are `pyamg` and, for the circuit diagrams alone, `schemdraw`. If your matplotlib needs numpy < 2, pin the versions: `pip install "numpy<2" "scipy<1.12" "pyamg==5.2.1"`. On a machine with several Pythons, use the one that has pyamg for every command (`py -3.12 run_all.py`); the field solver says so if it is missing.

All assumptions live in **`probe.toml`**: materials, part heights, ESL/ESR, cable, scope, sweep ranges.

## Files

| file | what it does |
|---|---|
| `run_all.py` | Runs the whole pipeline; each step logs to `results/<step>.txt` |
| `render_kicad.py` | Regenerates schematic PDF/SVG, layout SVG, assembly drawings and 3D renders with `kicad-cli` |
| `draw_schematics.py` | Draws the circuit diagrams of the main README with `schemdraw`, with the values on the board |
| `kicad_board.py` | Reads the `.kicad_pcb` directly. Run it alone to see what was imported |
| `make_pours.py` | Draws the copper pours that connect each divider step and refills all zones (runs in KiCad's Python) |
| `resin_fit.py` | Corrects a bridge sweep of a resin sample with a C0G reference and fits the resin model |
| `extract_caps.py` | 3D field solve: capacitance matrices, model-check plots, E-field stress maps |
| `tube_sweep.py` | Moves the tube end: stress against shielding |
| `probe_sim.py` | Circuit analyses: `info`, `response`, `optimize`, `sweep`, `damping`, `spice` |
| `fieldsolver.py`, `field_plots.py` | Solver and its plots |
| `probe_model.py` | Circuit model: nodal analysis, transmission line, inverse Laplace, optimisers, grading |
| `dielectric.py` | Causal wideband models: Cole–Cole fitted to the measured resin, Djordjevic–Sarkar for FR4 |
| `measurements/` | The resin sample's sweep, corrected (`huntsman_cw5620_2026-11-17.csv`) and as read with its C0G reference (`hm8118/`) |
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

**Damping resistors.** A resistor in series with a capacitor, joined at a net that nothing else connects to, is a damping resistor (R30–R49). It is left out of the chain above. The net between it and its capacitor is treated as part of the ladder node on the resistor's other side: one conductor in the field solve, and the resistance in series with the segment capacitor in the circuit. 33 Ω into that net's ~0.3 pF of stray capacitance is 10 ps, so nothing is lost below a few GHz, and the field solve needs 23 conductors instead of 43.

**Zones must be filled.** The tools read the filled polygons KiCad saves. A zone that was drawn but not filled has none; its outline is then used and a warning printed. `make_pours.py` refills everything.

**The capacitance file belongs to one board file.** `capacitance.json` records a hash of the `.kicad_pcb` it was solved for. The circuit tools say so when the board has changed since, and warn if divider nets are missing from it.

**Other parts:**
- A resistor between BOT and OUT (R29) is the back termination.
- RV2 is set automatically so that the DC ratio, with the scope's 1 MΩ, is exactly 1000:1.
- The trimmer is recognised from its value range ("5.5pF-65pF").
- Other two-terminal parts (the TVS D1) are included through `[circuit.part_capacitance]`.

`python kicad_board.py ../hv-probe-passive/hv-probe-passive.kicad_pcb` prints all of this. `python probe_sim.py info` warns about parts without a model.

## The steps

### Copper pours: `make_pours.py`

```
python make_pours.py                       # into the project's .kicad_pcb (close KiCad first)
python make_pours.py --out copy.kicad_pcb --bar 1.0 --stem 1.0     # try a variant on a copy
python make_pours.py --remove              # take the generated zones out again
```

It reads the pads of every divider net and draws one zone per net on their copper layer, then refills all zones and saves. It needs KiCad's own Python for the `pcbnew` module; started with another Python it restarts itself with the newest KiCad it finds (or `[kicad] python` in `probe.toml`).

| zone | net | shape |
|---|---|---|
| `s0` | HV (`IN`) | convex hull of the first resistor pad and the first damping-resistor pad, grown 1 mm towards the tip and the tube wall, running into a 2.5 mm strip to the tip zone |
| `s1` … `s19` | divider nodes | a T: a bar over the node's two resistor pads and a stem down to its capacitor and damping-resistor pads |
| `s20` | divider bottom | the last resistor pad and capacitor pad joined by a narrow slanted band, plus a track to the low-voltage circuit's dangling track end |
| `s0c` … `s19c` | between a damping resistor and its capacitor | a bar from pad centre to pad centre. Electrically part of step `s<n>` |

All zones are solid-connected to their pads (no thermal spokes), have a 1 mm corner fillet, priority 2, and keep the settings of a hand-drawn zone, so they can be edited in KiCad afterwards. Zones with these names are replaced on every run; other zones are left as drawn.

**Rules behind the shapes:**
- **Never outside the pads.** A node's pour stays inside the outline of the pads it joins. The copper therefore comes no closer to the shield or to the next step than the pads already do: the smallest copper-to-copper distance between adjacent steps is 2.42 mm with and without the pours.
- **No added points.** Pads are joined solid, every corner the pour adds is filleted, and where a pour's corner coincides with a pad's corner the pad's own rounded corner is the outline.
- **The HV pour is the exception.** Its capacitance to the shield only loads the source and does not affect the compensation, so it is made larger and rounder than its pads, on the sides facing the tip and the tube wall. That keeps the terminations of the first parts, which sit in the tube mouth, in its shadow.

**Why this shape.** Each variant was built on a copy of the board and put through the field solver: 20 kV, linear voltage distribution, the default grid, resin 4.25 and FR4 4.5 for all. Fields are the largest values in cells that do not touch a conductor.

| node pours (HV pour flush, 1 mm strip) | node → shield, total | resin field, divider body |
|---|---|---|
| none (pads only; the lower bound) | 24.59 pF | 7.18 kV/mm |
| slim: 1 mm bar, 1 mm stem | 24.96 pF (+1.5 %) | 7.06 kV/mm (−2 %) |
| **flush with the pads (used)** | **25.90 pF (+5.3 %)** | **6.50 kV/mm (−9.5 %)** |
| flush, grown 0.55 mm towards the tube wall | 26.50 pF (+7.8 %) | 6.78 kV/mm (−6 %) |
| flush, grown 1.1 mm towards the tube wall | 27.36 pF (+11.2 %) | 7.05 kV/mm (−2 %) |

| HV pour and strip (node pours flush) | hottest resin, anywhere | HV → ground |
|---|---|---|
| hull flush with the pads, 1 mm strip | 10.7 kV/mm, at the tube mouth | 2.9 pF |
| flush, 2.5 mm strip | 10.2 kV/mm | 3.1 pF |
| **grown 1 mm towards tip and wall, 2.5 mm strip (used)** | **8.4 kV/mm**, in the divider body | **3.2 pF** |
| grown 0.8 mm on every side, 4 mm strip | 8.0 kV/mm | 3.5 pF |
| grown 1.5 mm on every side, 4 mm strip | 8.0 kV/mm | 3.6 pF |

- **Copper between the pads costs capacitance and buys field:** filling the gaps removes edges. Flush pours were chosen because the graded capacitors absorb a few percent of stray capacitance, while the field sets the life of the potting. `--bar 1.0 --stem 1.0` gives the slim version: 3.6 % less capacitance than flush, and nearly all of the field reduction given back.
- **Copper beyond the pads towards the tube wall costs both:** more capacitance and a higher field than flush.
- **At the HV end the gain is large and nearly free.** Beyond about 1 mm of growth the hottest spot is no longer at the HV pour, and more copper there does not help.
- A margin smaller than one grid step (0.55 mm) cannot be compared this way: the grid lines sit on the pad edges, and the pour's edge would fall between two of them.

### Field solve: `extract_caps.py`

```
python extract_caps.py --scale 1.5   # coarse, ~2 min
python extract_caps.py               # default grid, 4.7 M nodes, 2 workers, ~4.5 GB RAM peak
```

Check `results/model_check.png` first. It shows how copper, part terminations, board, slots, resin, tube and end cap were discretised.

Outputs:
- `results/capacitance_summary.txt`: each node's capacitance to the shield, its neighbours and the tip.
- `results/capacitance.json`: the full matrices.
- **E-field stress maps** at `[stress] voltage` (20 kV): `stress_dc.png` uses the resistive (slow) voltage distribution, `stress_hf.png` the distribution a fast edge sees. `stress.json` lists peak |E| per material with locations.

`--stress-only` redoes just the stress part; `--voltage` changes the voltage. The divider is compensated first; `--graded 47p --round E12` maps the stress for the grading the tools compute instead of the PCB's capacitors (`run_all.py` passes it; for a board that carries that grading the two are the same). `--resin-er` sets the resin permittivity the solve uses, instead of the model's value at `solve_freq`.

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
| `--seg-r 22 --tip-r 100` | damping resistor in series with every divider capacitor (replaces the PCB's) / with the tip |
| `--scope-bw 100Meg --c-scope 20p` | scope bandwidth (1st order) and input capacitance |
| `--no-parasitics`, `--no-esr`, `--no-loss` | leave out the parasitics / capacitor ESR / dielectric dispersion and loss |
| `--optimize` | first find the bottom capacitance that makes the response flattest |
| `--name x` | output file name |

### Tube end: `tube_sweep.py`

```
python tube_sweep.py --from 265 --to 310 --step 2.5     # default grid
```

This moves the HV-side end of the tube and re-solves at each position. It reports:
- **Stress:** peak and bulk |E| per material, and the tangential field along the board.
- **Shielding:** the surroundings (the outer box) are kept as a separate conductor. Their coupling to every node goes through the circuit model, giving an *equivalent tip voltage per volt on the surroundings* at 50 Hz and 1 MHz.

Compare positions on the same `--scale`, since peaks at copper edges depend on the grid.

### Resin: `resin_fit.py`

```
python resin_fit.py                    # fit and plot from the corrected sweep named in probe.toml
python resin_fit.py --part measurements/hm8118/hm8118_huntsman_sample_2026-11-17.csv \
       --ref measurements/hm8118/hm8118_c0g220_2026-11-17b_range5.csv measurements/hm8118/hm8118_c0g220_2026-11-17b_range4.csv \
       --ref-auto measurements/hm8118/hm8118_c0g220_2026-11-17.csv      # redo the correction from the raw sweeps
```

The second form turns a raw sweep of a sample capacitor into the corrected file the model is fitted to; see [Resin measurement](#resin-measurement). Both print the fit and write `results/resin_fit.png`. After a new measurement, only `run_all.py --skip-field` is needed.

### Images: `render_kicad.py`

This uses `kicad-cli` (KiCad ≥ 8; the newest install is picked, or set `[kicad] cli` in `probe.toml`). It writes to `images/kicad/`:
- `schematic.pdf`, `schematic.svg`;
- `layout.svg`: copper, outline, shield and resin;
- `assembly-top.pdf`, `assembly-bottom.pdf`;
- `pcb-top.png`, `pcb-bottom.png`, `pcb-iso.png`: raytraced and cropped to the board with a 3 % margin; `--fast` for quick previews.

### Circuit diagrams: `draw_schematics.py`

This draws the three diagrams of the main README with [`schemdraw`](https://schemdraw.readthedocs.io/) and writes them to `images/schematics/`: the plain divider, the compensated divider, and the divider as built with the stray capacitance to the shield. The segment resistor, the first and last divider capacitor, the damping resistor and the number of segments are read from the board file.

Each diagram comes as `*-light.svg` and `*-dark.svg` with a transparent background, and the README picks one by the reader's theme with a `<picture>` element. Text is stored as outlines, so the diagrams look the same on every machine. The README shows them at their natural size, where the lettering is 11 pt. The script prints that width; if a layout change alters it, the `width` in the README needs the new number. `--png DIR` also writes previews on the page colours.

## How it works

### Field solver

Electrostatics, ∇·(ε∇φ) = 0, is solved with the finite integration technique on a graded rectilinear grid. The grid is 0.2 mm through the board and 0.55 mm in-plane, with grid lines on the pad edges of the divider nets, growing towards a grounded box 30 mm outside the tube.

**Conductors:**
- copper from zone fills, pads, tracks and vias, as zero-thickness sheets on F.Cu and B.Cu;
- plated holes;
- SMD terminations, as boxes of `[geometry.part_height]`;
- all User.1 metal.

**Dielectrics:** FR4 in the board outline (slots included), resin in the revolved User.2 outline, and air elsewhere (the unpotted LV end, around the tip, outside the tube). On the rectangular grid the tube wall is a staircase. Inside the tube the resin is extended into the cells between the ideal circle and that staircase, so no artificial air layer forms against the shield.

**Capacitance matrix.** Each conductor is set to 1 V in turn; the charges give one column of the Maxwell matrix. One AMG hierarchy (pyamg) is reused for all 23 solves. The nets between a damping resistor and its capacitor are merged with their ladder node first (see [Damping resistors](#what-the-pcb-has-to-contain)).

**Per-material matrices.** Each element is an energy integral, C_ij = Σ_m ε_m ∫_m ∇φ_i·∇φ_j dV, so it splits exactly into its resin, FR4 and air parts. These are saved as `partial` in `capacitance.json`; their sum equals the total to ~5e-7 (the CG tolerance, `[solver] tol = 1e-6`). The circuit model scales each part by its material's ε(s). That makes **dielectric parameters changeable without a new field solve** (except `solve_freq`). About 83 % of a node's capacitance to ground sits in resin, 17 % in FR4.

**Why not BEM?** A boundary-element solver is the natural choice for conductors in one homogeneous dielectric. Here there are three dielectrics, open dielectric boundaries, and thin conductors on both sides of a board. A volume method handles that directly, and its field solution doubles as the stress map. [FasterCap](https://www.fastfieldsolvers.com/fastercap.htm) (multi-dielectric BEM) and [Elmer FEM](https://www.elmerfem.org/) (`StatElecSolver` outputs a capacitance matrix) are good independent cross-checks.

**Accuracy.** A conductor is the set of grid nodes inside its outline. An edge that falls between two grid lines is therefore moved inwards by up to one step: the pads come out too small and the gaps between steps too wide, by a different amount for every step depending on where the grid falls. On this board that makes nominally identical steps differ by ±6 % and the total read about 8 % low. The grid therefore has lines on the pad edges of all divider nets (`[mesh] snap_pads`), so pads and terminations keep their true size on any grid. Lines closer than 0.4 of a step are merged, so no tiny cells appear. Node-to-shield capacitance, flush pours, resin 4.25, FR4 4.5:

| grid | in-plane step | without snapping: total, per node | with snapping: total, per node |
|---|---|---|---|
| `--scale 2` | 1.1 mm | 22.95 pF, 1.13–1.27 pF | 25.43 pF, 1.335–1.339 pF |
| `--scale 1.5` | 0.8 mm | 23.60 pF, 1.22–1.26 pF | 25.41 pF, 1.333–1.338 pF |
| default | 0.55 mm | not run | 25.90 pF, 1.36 pF |

On the default grid the inner nodes agree within 1.5 %, most within 0.2 %: the two pad rows have slightly different pitches (10.929 and 10.922 mm), so at some steps two pad edges fall within 0.4 of a step and share one grid line. What remains beyond that is the field singularity at the copper edges, which no grid resolves: the values still rise slowly with refinement (+1.9 % from 0.8 to 0.55 mm steps). The default grid is therefore probably a few percent low; a finer one needs more than 16 GB of memory. That is of the order of the other unknowns: resin permittivity, real part geometry, board centring, voids.

### Dielectrics

Loss and dispersion are linked (Kramers–Kronig): a material with a high loss tangent must have a permittivity that falls with frequency. Each dielectric is therefore a causal model that is analytic in s, so it goes straight into the circuit and the inverse Laplace transform:

- **FR4: Djordjevic–Sarkar,** ε(s) = ε∞ + k·log10((ω2 + s)/(ω1 + s)), fitted to one datasheet point (`fr4_er`, `fr4_tan_d` at `fr4_spec_freq`): a constant loss tangent over many decades.
- **Resin: Cole–Cole,** ε(s) = ε∞ + Δε / (1 + (sτ)^β), fitted to a bridge sweep of a cast sample (`[materials] resin_data`). One broad relaxation: β = 1 would be a single Debye pole, the resin has β = 0.24.

| resin (Cole–Cole fit) | 10 Hz | 50 Hz | 1 kHz | 10 kHz | 100 kHz | 1 MHz | 100 MHz |
|---|---|---|---|---|---|---|---|
| ε′ | 6.51 | 6.00 | 5.19 | 4.74 | 4.44 | 4.25 | 4.06 |
| tanδ | 0.076 | 0.078 | 0.069 | 0.054 | 0.038 | 0.025 | 0.009 |

The dispersion is spread over the whole band and beyond: at 100 kHz the permittivity is still 4.5 % above its 1 MHz value and the loss tangent is 0.04. The 200 Hz–120 kHz part is measured; below and above it the Cole–Cole curve is an extrapolation (`results/resin_fit.png`).

Without `resin_data` the resin falls back to the Djordjevic–Sarkar form with the datasheet point and an assumed high-frequency permittivity (`resin_er_hf`). A single datasheet point does not give the shape: with 4.0 assumed, that model puts all the dispersion below 10 kHz.

#### Resin measurement

**The sample.** A cast layer of cured Arathane CW 5620 between two aluminium-foil electrodes, on an R&S HM8118 bridge (1 V, parallel C–D), 69 frequencies from 20 Hz to 200 kHz, 17 November 2026. The layer was made to try the resin, not to be measured: it has the shape of a half ellipse, 65 × 43 mm, and its thickness varies from 1 to about 4 mm. The electrodes end 2 mm inside its edge (61 × 39 mm, about 1870 mm²). The files as read are in `measurements/hm8118/`, together with the sweeps of a 220 pF C0G capacitor taken on the same leads that day.

**The bridge's errors.** A C0G capacitor has a loss factor of about 10⁻⁵, so whatever the bridge shows on it is the bridge. `resin_fit.py` takes three errors off the sample's sweep:

| error | size on the 45 pF sample | correction |
|---|---|---|
| D reads low above 1 kHz, by a range-dependent amount | −0.0002 to −0.0006 on range 5, up to −0.0012 on range 4 at 100 kHz | the reference's D on the same range and frequency is subtracted |
| a fixed +0.35 pF at 0.4, 0.45, 0.75, 1.5, 3, 9, 12, 25 and 45 kHz, on every range | +0.8 % | read off the reference and subtracted |
| range 6 reads C 0.14 % below range 5 | −0.14 % below 4.5 kHz | divided out |

- The reference's rise of C above 30 kHz (+0.45 % at 100 kHz) is **not** applied as a factor. On the sample it would break the Kramers–Kronig consistency, which the readings have without it. It behaves like a series impedance, so it is scaled by C_sample/C_reference (0.08 % at 100 kHz) and half the full value goes into the point's uncertainty.
- Above about 30 kHz the D error grows with the part's impedance, and the sample is at six times the reference's. The correction is applied once, and three times its size is carried as that point's uncertainty.
- **Not used:** below 200 Hz (mains hum on a 45 pF part), above 120 kHz (range 4 fails on parts this small: D collapses towards 200 kHz), odd mains harmonics, and readings the sweep flagged. 46 of 69 points remain.

**Checks.**
- Kramers–Kronig, locally: tanδ ≈ −(π/2)·d ln C / d ln f. Over four bands from 200 Hz to 120 kHz the loss factor implied by the slope of C agrees with the measured one within +4.3, +2.2, −0.4 and +0.3 %.
- The Cole–Cole fit follows the 46 points to 0.15 % rms in C and 0.0015 rms in tanδ.

**What the sweep does not give.** Its shape over frequency is the material's. Two things are not in it, and `probe.toml` has to settle them:

1. **The scale.** Capacitance becomes permittivity through the sample's thickness and electrode area. A layer that is 1 to 4 mm thick cannot set it, so the datasheet point does: ε′ = 6.0 at 50 Hz. The uneven thickness costs nothing else: every part of the sample has the same ε(f), so the shape of C over frequency and tanδ are those of a flat sample. `resin_sample_thickness = [min, max]` records the range and is only used for the check below; a single number would make the geometry set the scale instead.
2. **Air under the electrodes.** Air between foil and sample is a lossless capacitance in series. It lowers C, tanδ and the dispersion all by the same factor, and nothing in the sweep itself can show it. Three things suggest there is some:
   - an earlier sweep of the same sample read 29.1 pF and tanδ 0.046 at 1 kHz, against 44.6 pF and 0.070 here: both lower by the same factor 0.65, as a larger gap would do;
   - a tube sample measured with the same foil electrodes that evening read tanδ 30–40 % below a later measurement with adhesive aluminium tape;
   - the fitted tanδ at 50 Hz is 0.078, against the datasheet's 0.11.

   `resin_electrode_gap = "none"` takes the sweep as it is. `"datasheet"` solves for the gap that brings the fitted tanδ at 50 Hz to the datasheet's 0.11: 177 pF, in series with a sample of 59 pF at 1 kHz.

| resin, by electrode gap | 10 Hz | 50 Hz | 1 kHz | 10 kHz | 100 kHz | 1 MHz | 100 MHz |
|---|---|---|---|---|---|---|---|
| ε′, `"none"` (used) | 6.51 | 6.00 | 5.19 | 4.74 | 4.44 | 4.25 | 4.06 |
| ε′, `"datasheet"` | 6.73 | 6.00 | 4.92 | 4.38 | 4.02 | 3.81 | 3.60 |
| tanδ, `"none"` (used) | 0.076 | 0.078 | 0.069 | 0.054 | 0.038 | 0.025 | 0.009 |
| tanδ, `"datasheet"` | 0.109 | 0.110 | 0.092 | 0.070 | 0.048 | 0.031 | 0.011 |

The two differ by 10 % in the permittivity at 1 MHz and by 30 % in loss. Through the circuit model the choice matters less: the graded capacitor values come out the same but for one (C31: 56 or 47 pF), and the remaining gain error is 2.8 % with `"none"` and 3.3 % with `"datasheet"`. The list on the board, kept under the `"datasheet"` resin, gives 3.4 %.

**What the sample's dimensions say.** An uneven sample has the capacitance of a flat one with the thickness t_eff = A / ∫ dA/t, in which the thin places count most. `resin_fit.py` prints the t_eff that each reading of the sweep stands for, with ε′ = 6 at 50 Hz and 1870 mm² of electrode:

| reading | t_eff | |
|---|---|---|
| full contact (`"none"`, used) | 1.93 mm | |
| 177 pF gap (`"datasheet"`) | 1.37 mm | the gap is 94 µm of air on average |
| the sample, if its thickness is spread evenly from 1 to 4 mm | 2.16 mm | 3 mm / ln 4 |

- Fringing at the edge of the unguarded electrodes is left out. A 2D solve of two equal foils on a slab with a 2 mm rim puts it at 3 % of the capacitance if the sample is 1 mm thick at the electrode's edge, 6 % at 2 mm and 11 % at 4 mm. That moves the first two rows up by as much: about 2.0–2.1 mm and 1.4–1.5 mm.
- **Full contact fits:** it agrees within about 10 % with an evenly spread sample. Read the other way, the geometry alone gives ε′ ≈ 6.0–6.7 at 50 Hz, against the datasheet's 6.0.
- **The gap reading is not excluded.** It needs most of the electrode area to lie over resin thinner than about 1.5 mm. Both values are inside 1–4 mm, so the dimensions lean towards `"none"` without deciding it.

**To settle it:**
- measure the sample's thickness at a dozen points under the electrodes. That gives t_eff to perhaps 10 %, enough to tell 1.4 mm from 2.0 mm;
- or cast a flat sample (between two plates, on spacers) and use electrodes that leave no air: conductive paint, or foil on a thin film of grease;
- above 120 kHz there is no data. An S11 sweep of the same sample on a VNA from 100 kHz to 100 MHz would replace the extrapolation that the bandwidth figures rest on.

#### Scaling the parasitics

The field solve uses the permittivities at `solve_freq`. The circuit scales every parasitic per material by ε(s)/ε(solve_freq). This is the tangent to the true dependence: exact at `solve_freq` and exact if all materials change by the same factor. Elsewhere it overestimates the self capacitances, by an amount that grows with the square of the distance in permittivity (the old field is not the energy minimum for the new permittivities).

`solve_freq` is 1 kHz, where the resin (5.19) lies midway between its values at 10 Hz (6.51) and 100 MHz (4.06). Against direct field solves at those permittivities, the scaled node-to-shield capacitance is 0.29 % high at 10 Hz, 0.25 % at 1 MHz and 0.38 % at 100 MHz, and the tip-to-ground capacitance 0.3–0.4 %. Both ends err the same way, so the difference across the compensation band, which is what the flatness depends on, is 0.04 %. A solve at 1 MHz instead would be 1.3–2.1 % off at 10 Hz, about as large as the gain error being computed.

### Circuit model

Nodal analysis at complex frequency s. It contains:
- all discrete parts from the PCB, with ESL and ESR per footprint;
- the damping resistors on the PCB, in series with their divider capacitors (`--seg-r` or `[circuit] segment_r` replaces them);
- the series inductance of each divider segment and of the tip strip (below);
- the parasitic capacitance matrices, per material, with dispersion and loss;
- the back termination R29 and the TVS capacitance;
- the cable as an exact lossy transmission line (skin effect ∝ √f fitted to the 100 MHz loss);
- the scope as 1 MΩ ∥ C_in, with an optional bandwidth limit;
- optional tip and segment damping resistors.

**Step responses** come from a numerical inverse Laplace transform (damped FFT with a window). This handles nanosecond cable reflections and millisecond RC tails in one framework. Step metrics use a Gaussian 1 ns source edge (`[analysis] source_rise`).

**Inductance.** The datasheet ESL of an 1808 capacitor, 1 nH, is the part's own 4.5 mm over a return a fraction of a millimetre away. In the probe the fast current runs the whole 10.9 mm pitch of a segment and returns on the tube wall, 14 mm away. `[circuit] segment_l = "auto"` therefore gives each segment the inductance of a coaxial line of its length, (µ0/2π)·pitch·ln(r_tube/r), with r = `segment_radius` (0.8 mm) the equivalent radius of the parts and pours that carry the current: 6.1 nH for the 10.9 mm pitch, 122 nH along the divider. `tip_l = "auto"` adds the partial inductance of the board's strip from the tip to the divider, and `source_l` whatever connects the probe to the circuit under test.

With these the divider is a lossy LC ladder with a characteristic impedance of about 70 Ω and a delay of about 1.8 ns end to end, not twenty independent RC sections. For the board's divider with its 33 Ω resistors, a 1 ns edge then overshoots by 11 %, where the datasheet ESL alone gives 5 %; with 22 Ω per capacitor it would be 29 % against 12 % (`results/damping.txt`). The inductance of the bottom capacitors counts as well: with half the 0.74 nH assumed for C34 and C35 the overshoot is 4 %, with twice 23 %. Both inductances are estimates from the geometry, good to perhaps ±30 %. `segment_l = "parts"` and `tip_l = 0` leave both out.

**Valid range.** The model is quasi-static and lumped: the inductances above are single lumped values, there is no magnetic coupling between segments and no radiation. It should be right in shape up to roughly 100 MHz; treat the frequency and height of anything above that as qualitative. In use, the lead and ground return to the circuit under test (about 1 nH per mm) will usually dominate what happens above some tens of MHz; set `source_l` to see it.

**LTspice.** `spice --run` writes the same network (ESR as `Rser`, ESL as `Lser`; dispersion is left out) and compares it with a matching probe_sim run. They agree to < 0.01 dB below 1 MHz; the small differences near cable resonances come from LTspice's LTRA line having no skin effect.

**Performance note.** The circuit solves are thousands of tiny (≈27×27) complex systems. Multi-threaded OpenBLAS is ~80× slower on these than single-threaded, so the scripts set `OPENBLAS_NUM_THREADS=1`.

**Where the time goes.** A full `run_all.py` takes about 13 min on an 8-core laptop (Ryzen 7 4700U) with 3–5 GB of memory free: 4.5 min for the field solve, then 8 min for the tube sweep on the two processes that fit, with the circuit analyses (1.5 min together) alongside. With memory for four tube-sweep processes the second stage halves.

- **Field solve.** 23 conductor solves, AMG-preconditioned CG, single-threaded. They run in `[solver] workers` processes, each with its own AMG hierarchy. The board with its damping resistors has 43 nets; merging each damping resistor's net into its ladder node keeps the count at 23. The stress maps are linear combinations of the 23 unit solutions, so they need no further solves; the per-material matrices are formed in banded row chunks. The CG tolerance of 1e-6 changes the matrix by < 0.01 fF against 1e-8.
- **Tube sweep.** One model and AMG hierarchy per tube position, in `[solver] tube_workers` processes. It runs alongside the circuit analyses.
- **Circuit analyses.** They run concurrently. The bottom-capacitance optimiser solves the network only twice per call and gets the response for every trial capacitance by superposition (see `optimize_bottom`), which makes the grading and rounding searches cheap. Most of the remaining time is in the step responses (2¹⁴ frequency points each).
- **Comparing layouts.** For a change of copper, two solves are enough: all divider nets at 1 V give every node's capacitance to the shield, and one solve at the working voltages gives the stress. That is how the pour variants above were compared, at 40–90 s each.

Memory, not cores, limits the parallelism on the default grid; see `[solver]` in `probe.toml`.

### Graded divider and its optimiser

`graded_caps` computes the segment capacitances that make the fast-edge voltage distribution linear. It starts from the bottom segment and works up: each segment carries the current of the one below plus all parasitic currents leaving its node.

With a dispersive resin the parasitics differ across the compensation band. `optimize_grading` therefore scales the grading strength (κ) and the bottom capacitance together for the flattest response, then rounds to the chosen E-series and re-optimises the bottom capacitance.

For the board as it is, the remaining error is 2.8 % at a grading strength of 1.55, with no capacitor above 270 pF. `run_all.py` runs the grading itself as well (`response_proposed`, `optimize_proposed`); on a board that carries the computed list these repeat the results for the board.

**Largest value.** `[optimize] cap_max` is the largest divider capacitor that can be bought. Values above 270 pF are hard to source in 3 kV C0G 1808, so it is set to 270 pF. `graded_caps` cuts the grading off there, and the rounding uses no larger value. The segments at the HV end then get less than the grading asks for, and the strength scan finds a stronger grading that makes up for it further down:

| largest value | κ | list, C14 → C33 [pF] | gain error | bottom C (1.5 m) | fast-edge share per segment |
|---|---|---|---|---|---|
| none | 1.15 | 470 390 390 330 330 270 220 180 180 150 120 120 100 100 82 68 56 56 47 47 | 3.0 % | 2.76 nF | 0.74–1.32 |
| 330 pF | 1.25 | 330 330 330 330 330 270 270 220 220 150 150 120 120 82 82 68 56 56 47 47 | 3.0 % | 2.79 nF | |
| **270 pF (on the board)** | **1.55** | **270 270 270 270 270 270 270 270 220 220 180 150 120 120 100 68 68 56 47 47** | **2.9 %** | **3.01 nF** | **0.70–1.43** |
| 220 pF | 2.15 | 220 220 220 220 220 220 220 220 220 220 220 180 150 150 82 82 68 56 47 47 | 3.6 % | 2.87 nF | |

- **The limit costs no flatness down to 270 pF.** The flatness is set by the resin's dispersion, and many capacitor lists reach it.
- **It costs evenness on fast edges.** The worst segment (C32, near the bottom) takes 1.43 times its DC share instead of 1.32: 1.43 kV at 20 kV, on a 3 kV part.
- **Cutting a finished list off is not the same.** The list of the first row with everything above 270 pF replaced by 270 pF gives 3.6 %, and 4.8–6.0 % if the parasitics are 5–10 % larger than modelled. The grading has to be redone with the limit.

**Objective.** "Flattest" means the smallest *peak* gain error over `[optimize] f_lo`–`f_hi`, the same number that is reported as the LF error. The bottom capacitance is found by a bounded scalar search; the peak error is quasi-convex in it, so the search finds the minimax optimum.

The peak can be set by a point the bottom capacitance cannot move, such as a low-frequency dip from the resin dispersion. It is then nearly flat over a range of capacitances, and the strict minimum may lie at the end of that range where the rest of the response is worst. So among the capacitances within `[optimize] peak_tolerance` (1 %) of the minimum peak, the one with the smallest RMS error over the same log-spaced points is used. `peak_tolerance = 0` gives the strict minimax.

On the board as it is the optimum is sharp instead. The peak is shared by −2.8 % at 6 kHz and +2.8 % at 1 MHz, which the bottom capacitance moves in opposite directions: the peak stays within 1 % of its minimum only from 3009 to 3011 pF, and the tie-break moves the result by less than 1 pF (3010 pF, RMS error 1.6 %).

The trimmer controls the first ~300 µs of the step response, about three time constants of the bottom RC, which is roughly where trimming by eye ends up. The largest slow error of the step, a −1.5 % sag around 1 ms, is the resin's slow polarisation, and no bottom capacitance changes it.

The grading and rounding searches compare designs by the strict minimax error. RMS alone is a poor design objective because it rewards any small average gain: it buys a slightly lower average with a higher peak, less bandwidth and a less even fast-edge distribution.

**Rounding.** Each capacitor is rounded down or up, not simply to the nearest value. Starting from nearest rounding, one capacitor at a time is switched to the series value on the other side of its unrounded value, if that lowers the peak error by more than `[optimize] peak_tolerance` (1 %). Every value stays within one series step of the computed grading, so the fast-edge distribution stays close to linear. For the current design this switches four capacitors (C23, C24, C27 and C28, each one step up) and lowers the error from 3.05 % to 2.88 %. The rest comes from the resin dispersion, which no fixed set of capacitors can follow.

**The list is not unique.** Lists that differ by one E12 step in a few places give the same error within a few hundredths of a percent, and which one the search returns depends on the parasitics at the percent level. If the search took every gain above 0.1 %, a change of 0.3 % in the parasitics would change two values and make the list non-monotonic for the same error, which is why it asks for 1 %. What counts is how a list, once bought, holds if the real parasitics differ from the model's (peak error with all parasitics scaled and the bottom capacitance re-trimmed):

| all parasitics × | 0.95 | 1.00 | 1.05 | 1.10 |
|---|---|---|---|---|
| **list on the board (largest 270 pF)** | 3.5 % | **2.9 %** | 3.1 % | 3.3 % |
| the same with C28 = 82 pF instead of 100 pF | 3.4 % | 2.9 % | 3.1 % | 3.3 % |
| nearest E12 of a weaker grading (κ = 1.5) | 2.9 % | 3.1 % | 3.3 % | 3.6 % |
| graded without the limit (470…47 pF) | 3.9 % | 3.0 % | 3.2 % | 3.4 % |

Any of the first four is a sound choice. The field solver's capacitances are more likely low than high, which is the side where the refined lists are better than nearest rounding.

### Stress maps

The field is computed with every node at its real potential, for both the DC and the fast-edge distribution. The potentials come from the circuit model of the compensated divider: the PCB's capacitors, or with `--graded/--round/--seg-r` the grading the tools compute. |E| is reported per material, with and without cells touching a conductor.

Sharp copper edges make the "peak" values grid-dependent. Use them to locate and rank hotspots; the "bulk" values and same-grid comparisons are the meaningful numbers. The air field around the exposed tip depends on the brass terminal (not modelled) and on the distance to grounded surroundings (the 30 mm box).

The maps are electrostatic (permittivity-weighted). Under steady DC the field divides between resin, FR4 and air by conductivity instead, which isn't modelled. The utilisation against bulk breakdown strength is not a lifetime criterion: partial discharge in voids and at interfaces is (see [Reading the stress numbers](../README.md#reading-the-stress-numbers)).

`tube_sweep.py` forces every swept tube-end position as a grid line, so all positions share one grid and the differences between them come from the tube alone.

## Limitations

- Part bodies are treated as resin; only terminations are conductive. The ceramic of the 1808 capacitors and the alumina of the resistors have a higher permittivity than the resin, so the capacitance across each segment is underestimated. THT parts at the LV end are simplified to pad columns. Copper thickness is ignored.
- The resin's permittivity is measured from 200 Hz to 120 kHz only. Its scale is the datasheet's, and the electrode air gap is not decided (see [Resin measurement](#resin-measurement)).
- The inductances are geometric estimates, without coupling between segments. ESR is frequency-independent, and radiation loss is ignored, so resonances above ~100 MHz are probably underdamped in the model.
- The field solve's capacitances are probably a few percent low on the default grid (see Accuracy).
- The brass tip terminal, the input fuse and the connection to the circuit under test are not modelled.
