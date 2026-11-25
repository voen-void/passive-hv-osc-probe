"""Draw the explanatory circuit diagrams of the README with schemdraw.

    python draw_schematics.py              # SVGs into ../images/schematics/
    python draw_schematics.py --png DIR    # PNG previews on the page colours as well

Three diagrams, each in a light and a dark variant (transparent background, so
the README can pick the one that fits the reader's theme):

    divider              resistive divider
    divider-compensated  the same with the compensation capacitors
    divider-ladder       the divider as built: segments, damping resistors and
                         the stray capacitance to the shield

The values in the diagrams are read from the board file: segment resistor,
first and last divider capacitor, damping resistor, number of segments.
"""

import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
import schemdraw
import schemdraw.elements as elm
from schemdraw.segments import SegmentText

import kicad_board
from extract_caps import load_config

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# line and text colour, colour of the stray capacitances, page colour (previews only)
THEMES = {
    "light": {"fg": "#1f2328", "stray": "#cf222e", "page": "#ffffff"},
    "dark": {"fg": "#e6edf3", "stray": "#ff7b72", "page": "#0d1117"},
}
FONTS = ["Segoe UI", "Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]
FONT_SIZE = 11        # pt; 1 pt = 1.33 px at the size the README shows the diagrams
UNIT = 2.4            # length of a two-terminal element, in drawing units of 0.5 inch
RISE = 1.25           # height of the capacitor branch above its resistor
DOT_R = 0.075         # radius of schemdraw's dots


def setup_fonts():
    have = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in FONTS if f in have), "DejaVu Sans")
    matplotlib.rcParams.update({
        "font.family": font,
        "mathtext.fontset": "custom", "mathtext.rm": font, "mathtext.it": f"{font}:italic",
        "mathtext.bf": f"{font}:bold",
        "svg.fonttype": "path",          # text as outlines: the same on every machine
        "svg.hashsalt": "hv-probe",      # stable ids, so a re-run gives the same file
    })
    return font


def eng(x, unit):
    """2e6, 'Ω' -> '2 MΩ'"""
    for scale, prefix in ((1e9, "G"), (1e6, "M"), (1e3, "k"), (1, ""), (1e-3, "m"), (1e-6, "µ"), (1e-9, "n"), (1e-12, "p")):
        if abs(x) >= scale * 0.9995:
            return f"{x / scale:.3g} {prefix}{unit}"
    return f"{x:g} {unit}"


def board_values(cfg):
    board = kicad_board.load_board(cfg["files"]["pcb"])
    ladder = kicad_board.ladder_order(board, kicad_board.resolve_net(board, cfg["circuit"]["input_net"]))
    damping, _ = kicad_board.segment_damping(board, ladder)
    value = lambda ref: kicad_board.parse_si(board.parts[ref].value)
    seg_r = []
    for a, b in zip(ladder, ladder[1:]):
        seg_r.append(next(value(p.ref) for p in board.parts.values()
                          if p.ref.startswith("R") and not p.ref.startswith("RV") and set(p.nets) == {a, b}))
    n = len(seg_r)
    first, last = damping[0][0], damping[n - 1][0]        # (damping resistor, capacitor, net between them)
    r_top = sum(seg_r)
    return {
        "n": n,
        "r_seg": eng(seg_r[0], "Ω"),
        "r_top": eng(r_top, "Ω"),
        "r_bottom": eng(r_top / (cfg["circuit"]["target_ratio"] - 1), "Ω").replace(".0 ", " "),
        "ratio": f"{cfg['circuit']['target_ratio']:.0f}",
        "c_first": eng(value(first[1]), "F"), "c_last": eng(value(last[1]), "F"),
        "rd_first": eng(value(first[0]), "Ω"), "rd_last": eng(value(last[0]), "Ω"),
        "c_refs": (first[1], last[1]),
    }


class Text(elm.Element):
    """A line of text at a point, aligned as given (horizontal, vertical)."""

    def __init__(self, text, align=("center", "center")):
        super().__init__()
        self.segments.append(SegmentText((0, 0), text, align=align, fontsize=FONT_SIZE))


def drawing(theme):
    d = schemdraw.Drawing(show=False)
    d.config(unit=UNIT, fontsize=FONT_SIZE, font=matplotlib.rcParams["font.family"][0], color=theme["fg"],
             lw=1.3, margin=0.08)
    return d


def terminal(d, at, text, side, lead):
    """An open circle with its name beside it and a lead of the given length
    towards the circuit. Returns the circuit end of the lead."""
    sign = 1 if side == "left" else -1
    d.add(elm.Dot(open=True).at(at).fill(False))
    d.add(Text(text, ("right" if side == "left" else "left", "center")).at((at[0] - sign * 0.3, at[1])))
    end = (at[0] + sign * lead, at[1])
    d.add(elm.Line().at((at[0] + sign * DOT_R, at[1])).to(end))       # from the circle's edge: it has no fill
    return end


def branch(d, start, r_label, c_label=None, rd_label=None, dots=(True, True), shared=False):
    """A resistor to the right of `start`, with a capacitor (and a damping
    resistor in series with it) bridging it. Returns the right-hand end.
    `shared`: the segment before ends on the same node, so its capacitor
    branch meets this one at the top."""
    r = elm.ResistorIEC().right().at(start)
    d.add(r.label(r_label, loc="bottom"))
    if c_label is not None:
        up = elm.Line().up().at(r.start).length(RISE)
        d.add(up)
        if shared:
            d.add(elm.Dot().at(up.end))
        if rd_label is None:
            d.add(elm.Capacitor().right().label(c_label, loc="top"))
        else:
            d.add(elm.Capacitor().right().length(UNIT * 0.42).label(c_label, loc="top"))
            d.add(elm.ResistorIEC().right().length(UNIT * 0.58).label(rd_label, loc="top"))
        d.add(elm.Line().down().length(RISE))
        for on, at in zip(dots, (r.start, r.end)):
            if on:
                d.add(elm.Dot().at(at))
    return r.end


def shunt(d, at, kind, label, color=None, rail=False):
    """An element from `at` down to the ground rail UNIT below."""
    e = (elm.ResistorIEC() if kind == "R" else elm.Capacitor()).down().at(at).label(label, loc="bottom")
    if color:
        e = e.color(color)
    d.add(e)
    d.add(elm.Dot().at(at))
    if rail:
        d.add(elm.Dot().at(e.end))
    return e


def draw_divider(theme, v, compensated):
    d = drawing(theme)
    start = terminal(d, (0, 0), "input", "left", 0.6)
    end = branch(d, start, f"$R_\\mathrm{{top}}$\n{v['r_top']}", "$C_\\mathrm{top}$" if compensated else None)
    node = (end[0] + 0.9, 0)
    out = (node[0] + 2.3, 0) if compensated else node
    d.add(elm.Line().at(end).to(out))
    d.add(elm.Ground().at(shunt(d, node, "R", f"$R_\\mathrm{{bottom}}$\n{v['r_bottom']}").end))
    if compensated:
        d.add(elm.Ground().at(shunt(d, out, "C", "$C_\\mathrm{bottom}$").end))
    terminal(d, (out[0] + 1.0 + DOT_R, 0), "output", "right", 1.0 + DOT_R)
    return d


def draw_ladder(theme, v):
    """The first two segments and the last one, with the stray capacitance of
    the nodes between segments to the shield, which is the ground rail."""
    d = drawing(theme)
    rail = -UNIT
    stray = lambda at, label="": shunt(d, at, "C", label, color=theme["stray"], rail=True)

    start = terminal(d, (0, 0), "input", "left", 0.6)
    n1 = branch(d, start, v["r_seg"], v["c_first"], v["rd_first"])
    stray(n1, "$C_\\mathrm{s}$")
    n2 = branch(d, n1, v["r_seg"], "", "", dots=(False, True), shared=True)
    stray(n2)
    # the segments left out
    lead, skip = 0.35, 1.9
    a, b = (n2[0] + lead, 0), (n2[0] + lead + skip, 0)
    n3 = (b[0] + lead, 0)
    d.add(elm.Line().at(n2).to(a))
    d.add(elm.Line().at(a).to(b).linestyle(":"))
    d.add(elm.Line().at(b).to(n3))
    d.add(Text(f"{v['n']} in total", ("center", "bottom")).at(((a[0] + b[0]) / 2, 0.22)))
    stray(n3)
    n4 = branch(d, n3, v["r_seg"], v["c_last"], v["rd_last"], dots=(False, True))

    # bottom of the divider
    xr, xc = n4[0] + 0.9, n4[0] + 3.2
    out = terminal(d, (xc + 1.0 + DOT_R, 0), "output", "right", 1.0 + DOT_R)
    d.add(elm.Line().at(n4).to(out))
    shunt(d, (xr, 0), "R", "$R_\\mathrm{bottom}$", rail=True)
    shunt(d, (xc, 0), "C", "$C_\\mathrm{bottom}$", rail=True)

    x0 = n1[0] - 0.9
    d.add(elm.Line().at((x0, rail)).to((out[0], rail)))
    d.add(elm.Ground().at((out[0], rail)))
    d.add(Text("shield", ("right", "center")).at((x0 - 0.2, rail)))
    return d


def save_svg(d, path):
    d.save(path, transparent=True)
    with open(path, encoding="utf-8") as f:
        s = f.read()
    s = re.sub(r"\s*<metadata>.*?</metadata>", "", s, flags=re.S)     # holds the date of the run
    s = "\n".join(line.rstrip() for line in s.splitlines()) + "\n"
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(s)
    m = re.search(r'<svg[^>]*?width="([\d.]+)pt"', s)
    return round(float(m.group(1)) * 96 / 72) if m else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(HERE, "probe.toml"))
    ap.add_argument("--out", default=os.path.join(ROOT, "images", "schematics"))
    ap.add_argument("--png", metavar="DIR", help="also write PNG previews on the page colour of each theme")
    a = ap.parse_args()

    schemdraw.use("matplotlib")
    font = setup_fonts()
    v = board_values(load_config(a.config))
    os.makedirs(a.out, exist_ok=True)
    line = (f"font: {font}; from the board: {v['n']} x {v['r_seg']}, {v['c_refs'][0]} {v['c_first']} ... "
            f"{v['c_refs'][1]} {v['c_last']}, damping {v['rd_first']}")
    print(line.replace("Ω", "ohm").replace("µ", "u"))      # consoles with a legacy code page
    figures = {
        "divider": lambda t: draw_divider(t, v, False),
        "divider-compensated": lambda t: draw_divider(t, v, True),
        "divider-ladder": lambda t: draw_ladder(t, v),
    }
    for name, make in figures.items():
        for tname, theme in THEMES.items():
            width = save_svg(make(theme), os.path.join(a.out, f"{name}-{tname}.svg"))
            if a.png:
                d = make(theme)
                d.config(bgcolor=theme["page"])
                d.save(os.path.join(a.png, f"{name}-{tname}.png"), transparent=False, dpi=192)
        print(f"  {name}: {width} px wide at its natural size")


if __name__ == "__main__":
    main()
