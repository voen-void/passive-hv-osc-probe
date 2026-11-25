"""Draw the copper pours that connect each divider step, and refill all zones.

    python make_pours.py                 # writes the zones into the project's .kicad_pcb
    python make_pours.py --out copy.kicad_pcb --margin 0.5   # try a variant on a copy
    python make_pours.py --remove        # take the generated zones out again

Close the board in KiCad first: the file is rewritten. It needs KiCad's own
Python (the `pcbnew` module); started with any other Python it restarts itself
with the newest KiCad it finds, or with [kicad] python = "..." in probe.toml.

What it draws, all on the layer the divider parts sit on:

    s0            the HV net: first resistor + first damping resistor, and the
                  strip to the tip zone
    s1 .. s19     one pour per divider node: the two resistor pads, the
                  capacitor pad and the damping-resistor pad of that node
    s20           the divider bottom, plus a track to the low-voltage circuit
    s0c .. s19c   the short link between a damping resistor and its capacitor
                  (electrically part of step s<n>: 10 ohm away from it)

Shape rules (see tools/README.md, "Copper pours"):
  * a pour never reaches outside the pads it joins, so the copper comes no
    closer to the shield or to the neighbouring step than the pads already do;
  * pads are joined solid (no thermal spokes) and every corner the pour adds
    is filleted, so the pours add no points for the field to concentrate on;
  * --margin grows the pours towards the tube wall only, with corners of
    --radius (rounder outline, slightly more capacitance to the shield).

Zones with the names above are replaced on every run, so the pours follow the
parts when they are moved. Other zones are left as drawn; all are refilled.
"""

import argparse
import glob
import os
import re
import subprocess
import sys
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
OWN_ZONE = re.compile(r"^s\d+c?$")


# ----------------------------------------------------------------------------
# KiCad's Python

def find_kicad_python():
    try:
        with open(os.path.join(HERE, "probe.toml"), "rb") as f:
            exe = tomllib.load(f).get("kicad", {}).get("python", "")
        if exe and os.path.exists(exe):
            return exe
    except OSError:
        pass

    def version(p):   # newest install first: "10.0" > "9.0" numerically
        v = os.path.basename(os.path.dirname(os.path.dirname(p)))
        return tuple(int(x) if x.isdigit() else 0 for x in v.split("."))
    cands = sorted(glob.glob(r"C:\Program Files\KiCad\*\bin\python.exe"), key=version, reverse=True)
    cands += ["/Applications/KiCad/KiCad.app/Contents/Frameworks/Python.framework/Versions/Current/bin/python3",
              "/usr/bin/python3"]
    return next((c for c in cands if os.path.exists(c)), None)


try:
    import pcbnew
except ImportError:
    pcbnew = None


# ----------------------------------------------------------------------------
# geometry [mm], KiCad board coordinates (y down)

def bbox(polys):
    xs = [p[0] for poly in polys for p in poly]
    ys = [p[1] for poly in polys for p in poly]
    return min(xs), max(xs), min(ys), max(ys)


def hull(points):
    """Convex hull (monotone chain), in drawing order."""
    pts = sorted(set((round(x, 6), round(y, 6)) for x, y in points))
    if len(pts) < 3:
        return pts

    def half(seq):
        out = []
        for p in seq:
            while len(out) >= 2 and ((out[-1][0] - out[-2][0]) * (p[1] - out[-2][1])
                                     - (out[-1][1] - out[-2][1]) * (p[0] - out[-2][0])) <= 0:
                out.pop()
            out.append(p)
        return out
    lower, upper = half(pts), half(reversed(pts))
    return lower[:-1] + upper[:-1]


def box_corners(b):
    x0, x1, y0, y1 = b
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def tee(upper, lower, stem=0.0):
    """Two pad groups, one above the other, joined by a stem. The stem is as
    wide as the overlap of the two groups (or `stem`, centred in it), so
    neither side of it comes closer to a neighbouring step than the pads do.
    Returns None without overlap."""
    ux0, ux1, uy0, uy1 = upper
    lx0, lx1, ly0, ly1 = lower
    sx0, sx1 = max(ux0, lx0), min(ux1, lx1)
    if sx1 - sx0 < 0.5 or uy1 >= ly0:
        return None
    if 0 < stem < sx1 - sx0:
        mid = 0.5 * (sx0 + sx1)
        sx0, sx1 = mid - stem / 2, mid + stem / 2
    pts = [(ux0, uy0), (ux1, uy0), (ux1, uy1)]
    if ux1 - sx1 > 1e-6:
        pts.append((sx1, uy1))
    if lx1 - sx1 > 1e-6:
        pts += [(sx1, ly0), (lx1, ly0)]
    pts += [(lx1, ly1), (lx0, ly1)]
    if sx0 - lx0 > 1e-6:
        pts += [(lx0, ly0), (sx0, ly0)]
    if sx0 - ux0 > 1e-6:
        pts += [(sx0, uy1), (ux0, uy1)]
    return pts


def zee(upper, lower, width=0.6):
    """Two pad groups that do not overlap in x, joined by a slanted band
    between their facing corners. The band is narrow enough that its outer
    edge stays behind the line between the two pad corners, so it comes no
    closer to the next step than those corners do."""
    ux0, ux1, uy0, uy1 = upper
    lx0, lx1, ly0, ly1 = lower
    if uy1 >= ly0:
        return None
    w = min(width, 0.5 * (ux1 - ux0), 0.5 * (lx1 - lx0))
    if lx1 <= ux0 + 0.5:        # lower group to the left of the upper one
        return [(ux0, uy0), (ux1, uy0), (ux1, uy1), (ux0 + w, uy1), (lx1, ly0), (lx1, ly1),
                (lx0, ly1), (lx0, ly0), (lx1 - w, ly0), (ux0, uy1)]
    if lx0 >= ux1 - 0.5:        # lower group to the right
        return [(ux0, uy0), (ux1, uy0), (ux1, uy1), (lx0 + w, ly0), (lx1, ly0), (lx1, ly1),
                (lx0, ly1), (lx0, ly0), (ux1 - w, uy1), (ux0, uy1)]
    return None


def grow_outward(b, yc, margin):
    """Move the edge of a pad group that faces away from the axis outwards."""
    x0, x1, y0, y1 = b
    return (x0, x1, y0 - margin, y1) if 0.5 * (y0 + y1) < yc else (x0, x1, y0, y1 + margin)


def node_outline(r_pads, c_pads, yc, margin, bar=0.0, stem=0.0):
    """Pour of one divider node: its resistor pads and its capacitor-row pads.
    bar > 0 limits the height of the bar between the resistor pads (it then
    sits on their side facing the capacitor row); stem > 0 the stem's width."""
    groups = [grow_outward(bbox(g), yc, margin) for g in (r_pads, c_pads) if g]
    if len(groups) == 1:
        return box_corners(groups[0])
    if r_pads and c_pads and 0 < bar < groups[0][3] - groups[0][2]:
        x0, x1, y0, y1 = groups[0]
        groups[0] = (x0, x1, y1 - bar, y1) if y0 < groups[1][2] else (x0, x1, y0, y0 + bar)
    upper, lower = sorted(groups, key=lambda b: b[2])
    return tee(upper, lower, stem) or zee(upper, lower) or hull(box_corners(upper) + box_corners(lower))


def link_outline(pads):
    """Bridge between two pads in one row: from pad centre to pad centre, as
    high as the lower of the two pads, so it stays inside both."""
    (ax0, ax1, ay0, ay1), (bx0, bx1, by0, by1) = sorted((bbox([p]) for p in pads), key=lambda b: b[0])
    return box_corners((0.5 * (ax0 + ax1), 0.5 * (bx0 + bx1), max(ay0, by0), min(ay1, by1)))


def hv_outline(r_pads, c_pads, tip_x, yc, width, margin, taper=1.5):
    """HV net: the convex hull of the first resistor pad, the first
    damping-resistor pad and the root of the strip, with the strip running
    along the axis to the tip zone. The hull is convex, so with its corners
    filleted the HV end of the divider has no corner pointing at the shield.
    `margin` grows the pads before the hull is taken, towards the tip and
    towards the tube wall but not towards the next step: the HV net's
    capacitance to the shield only loads the source, so here a larger, rounder
    electrode that keeps the part terminations in its shadow costs nothing."""
    boxes = [bbox(g) for g in (r_pads, c_pads) if g]
    sgn = 1.0 if tip_x > max(b[1] for b in boxes) else -1.0
    groups = [(x0 - (margin if sgn < 0 else 0.0), x1 + (margin if sgn > 0 else 0.0), y0 - margin, y1 + margin)
              for x0, x1, y0, y1 in boxes]
    corners = [p for g in groups for p in box_corners(g)]
    xe = max(sgn * p[0] for p in corners) * sgn + sgn * taper          # root of the strip
    root = [(round(xe, 6), round(yc - width / 2, 6)), (round(xe, 6), round(yc + width / 2, 6))]
    body = hull(corners + root)
    i, j = body.index(root[0]), body.index(root[1])
    far = {root[0]: (tip_x, yc - width / 2), root[1]: (tip_x, yc + width / 2)}
    if (i + 1) % len(body) == j:        # ... root0, root1 ...
        first, second = i, j
    elif (j + 1) % len(body) == i:      # ... root1, root0 ...
        first, second = j, i
    else:
        raise ValueError("the strip root is not an edge of the HV pour's hull")
    out = []
    for k, p in enumerate(body):
        out.append(p)
        if k == first:
            out += [far[body[first]], far[body[second]]]
    return out


# ----------------------------------------------------------------------------
# what to draw

def plan(board, cfg, a):
    """{zone name: (net, outline)} and the bottom track, from the board."""
    import kicad_board as kb
    hv = kb.resolve_net(board, cfg["circuit"]["input_net"])
    ladder = kb.ladder_order(board, hv)
    damping, merged = kb.segment_damping(board, ladder)
    yc = board.tube["yc"]

    seg_of = {frozenset(pair) for pair in zip(ladder, ladder[1:])}
    ladder_r = {p.ref for p in board.parts.values() if kb._is_resistor(p) and frozenset(p.nets) in seg_of}
    pads = {}                                     # net -> [(ref, polygon)]
    layers = set()
    for ref, fp, layer, net, poly in board.pad_shapes:
        pads.setdefault(net, []).append((ref, poly))
        if net in ladder or net in merged:
            layers.add(layer)
    if len(layers) != 1:
        raise SystemExit(f"the divider parts sit on {sorted(layers) or 'no layer'}: expected one copper layer")
    layer = layers.pop()

    # parts of the capacitor row: divider capacitors and their damping resistors
    row = {r for brs in damping.values() for rr, cc, _ in brs for r in (rr, cc)}
    row |= {p.ref for p in board.parts.values() if kb._is_capacitor(p) and frozenset(p.nets) in seg_of}

    zones = {}
    for k, net in enumerate(ladder):
        r_pads = [poly for ref, poly in pads.get(net, []) if ref in ladder_r]
        c_pads = [poly for ref, poly in pads.get(net, []) if ref in row]
        if not r_pads and not c_pads:
            continue
        if k == 0:
            tip = a.tip_x if a.tip_x is not None else tip_zone_edge(board, hv, yc, bbox(r_pads + c_pads))
            outline = hv_outline(r_pads, c_pads, tip, yc, a.strip, a.hv_margin) if tip is not None \
                else node_outline(r_pads, c_pads, yc, a.margin)
        else:
            last = k == len(ladder) - 1
            outline = node_outline(r_pads, c_pads, yc, 0.0 if last else a.margin,
                                   0.0 if last else a.bar, 0.0 if last else a.stem)
        zones[f"s{k}"] = (net, outline)
    for mid, to in merged.items():
        zones[f"s{ladder.index(to)}c"] = (mid, link_outline([poly for _, poly in pads[mid]]))

    # divider bottom -> the dangling track end of the low-voltage circuit
    track = None
    bot = ladder[-1]
    r_bot = [poly for ref, poly in pads.get(bot, []) if ref in ladder_r]
    if r_bot:
        x0, x1, y0, y1 = bbox(r_bot)
        ends = [e for (net, lay), cu in board.copper.items() if net == bot and lay == layer
                for p1, p2, w in cu.lines for e in ((p1, w), (p2, w))]
        bot_pads = [poly for _, poly in pads.get(bot, [])] + [poly for net, poly in board.barrels if net == bot]
        inside = lambda p: any(_in_box(p, bbox([poly])) for poly in bot_pads)
        dangling = [(p, w) for p, w in ends if not inside(p) and sum(1 for q, _ in ends if _close(p, q)) == 1]
        if dangling:
            (px, py), w = min(dangling, key=lambda e: abs(e[0][0] - 0.5 * (x0 + x1)))
            if y0 < py < y1:
                track = (bot, layer, (px, py), (0.5 * (x0 + x1), py), w)
    return layer, zones, track


def _in_box(p, b):
    return b[0] - 1e-3 <= p[0] <= b[1] + 1e-3 and b[2] - 1e-3 <= p[1] <= b[3] + 1e-3


def _close(p, q):
    return abs(p[0] - q[0]) < 1e-3 and abs(p[1] - q[1]) < 1e-3


def tip_zone_edge(board, hv, yc, body):
    """x where the strip from the divider must end: just inside the HV net's
    other copper (the tip zone), on the side away from the divider."""
    x0, x1, _, _ = body
    best = None
    for (net, _), cu in board.copper.items():
        if net != hv:
            continue
        for poly in cu.polys:
            px0, px1, py0, py1 = bbox([poly])
            if not (py0 <= yc <= py1) or (px1 - px0) < 3.0:
                continue
            if px0 > x1:                    # tip beyond the divider's high-x end
                cand = px0 + 0.5
                best = cand if best is None else min(best, cand)
            elif px1 < x0:
                cand = px1 - 0.5
                best = cand if best is None else max(best, cand)
    return best


# ----------------------------------------------------------------------------
# KiCad side

def find_net(brd, name):
    net = brd.FindNet(name)
    if net is None:
        raise SystemExit(f"net {name} not found on the board")
    return net


def apply(a, cfg):
    sys.path.insert(0, HERE)
    import kicad_board as kb

    mm = pcbnew.FromMM
    brd = pcbnew.LoadBoard(a.pcb)
    for z in list(brd.Zones()):
        if OWN_ZONE.match(z.GetZoneName() or ""):
            brd.Remove(z)
    if a.remove:
        pcbnew.ZONE_FILLER(brd).Fill(brd.Zones())
        pcbnew.SaveBoard(a.out, brd, True)        # True: leave the project settings files alone
        print(f"removed the generated zones, wrote {a.out}")
        return

    # geometry from the same reader the field solver uses (pads and nets only,
    # so zones from an earlier run that are still in the file do not matter)
    layer, zones, track = plan(kb.load_board(a.pcb), cfg, a)
    lid = brd.GetLayerID(layer)

    made = []          # own references: zones added from Python come back from brd.Zones() untyped
    for name in sorted(zones, key=lambda n: (int(re.sub(r"\D", "", n)), n)):
        net, outline = zones[name]
        z = pcbnew.ZONE(brd)
        made.append((name, z))
        z.SetLayer(lid)
        z.SetNet(find_net(brd, net))
        z.SetZoneName(name)
        z.SetAssignedPriority(a.priority)
        z.SetPadConnection(pcbnew.ZONE_CONNECTION_FULL)
        z.SetLocalClearance(mm(a.clearance))
        z.SetMinThickness(mm(0.25))
        z.SetCornerSmoothingType(pcbnew.ZONE_SETTINGS.CORNER_SMOOTHING_FILLET)
        z.SetCornerRadius(mm(a.radius))
        z.SetIslandRemovalMode(pcbnew.ISLAND_REMOVAL_MODE_ALWAYS)
        z.SetBorderDisplayStyle(pcbnew.ZONE_BORDER_DISPLAY_STYLE_DIAGONAL_EDGE, mm(0.5), True)
        for x, y in outline:
            z.AppendCorner(pcbnew.VECTOR2I(mm(x), mm(y)), -1)     # -1: the outer outline
        brd.Add(z)

    if track:
        net, lay, p1, p2, w = track
        t = pcbnew.PCB_TRACK(brd)
        t.SetStart(pcbnew.VECTOR2I(mm(p1[0]), mm(p1[1])))
        t.SetEnd(pcbnew.VECTOR2I(mm(p2[0]), mm(p2[1])))
        t.SetWidth(mm(w))
        t.SetLayer(brd.GetLayerID(lay))
        t.SetNet(find_net(brd, net))
        brd.Add(t)

    pcbnew.ZONE_FILLER(brd).Fill(brd.Zones())
    unfilled = [name for name, z in made if not z.HasFilledPolysForLayer(lid)]
    pcbnew.SaveBoard(a.out, brd, True)            # True: leave the project settings files alone
    print(f"{len(zones)} zones on {layer} (fillet {a.radius:g} mm, margin {a.margin:g} mm, "
          f"HV strip {a.strip:g} mm, HV margin {a.hv_margin:g} mm)" + (", 1 track at the divider bottom" if track else ""))
    if unfilled:
        print("WARNING: no copper after the fill in:", ", ".join(unfilled))
    print(f"wrote {a.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(HERE, "probe.toml"))
    ap.add_argument("--pcb", help="board to read (default: [files] pcb in probe.toml)")
    ap.add_argument("--out", help="board to write (default: the board that was read)")
    ap.add_argument("--radius", type=float, default=1.0, help="fillet radius of the pour corners [mm]")
    ap.add_argument("--margin", type=float, default=0.0,
                    help="grow the pours beyond the pads towards the tube wall [mm]")
    ap.add_argument("--bar", type=float, default=0.0,
                    help="height of the bar between a node's two resistor pads [mm]; 0 = as high as the pads")
    ap.add_argument("--stem", type=float, default=0.0,
                    help="width of the stem between a node's resistor and capacitor rows [mm]; 0 = as wide as "
                         "the capacitor-row pads")
    ap.add_argument("--strip", type=float, default=2.5, help="width of the HV strip to the tip zone [mm]")
    ap.add_argument("--hv-margin", type=float, default=1.0,
                    help="grow the HV pour beyond its pads towards the tip and the tube wall [mm]")
    ap.add_argument("--tip-x", type=float, help="x where the HV strip ends (default: inside the tip zone)")
    ap.add_argument("--clearance", type=float, default=0.5, help="zone clearance [mm]")
    ap.add_argument("--priority", type=int, default=2, help="zone priority (above the ground and tip zones)")
    ap.add_argument("--remove", action="store_true", help="only remove the generated zones")
    a = ap.parse_args()

    if pcbnew is None:
        exe = find_kicad_python()
        if not exe or os.path.samefile(exe, sys.executable):
            sys.exit("KiCad's Python (module pcbnew) not found; set [kicad] python = '...' in probe.toml")
        sys.exit(subprocess.run([exe, os.path.abspath(__file__), *sys.argv[1:]]).returncode)

    with open(a.config, "rb") as f:
        cfg = tomllib.load(f)
    base = os.path.dirname(os.path.abspath(a.config))
    a.pcb = os.path.abspath(a.pcb or os.path.join(base, cfg["files"]["pcb"]))
    a.out = os.path.abspath(a.out or a.pcb)
    lock = os.path.join(os.path.dirname(a.out), "~" + os.path.basename(a.out) + ".lck")
    if os.path.exists(lock):
        sys.exit(f"{os.path.basename(a.out)} is open in KiCad (lock file found): close it first")
    apply(a, cfg)


if __name__ == "__main__":
    main()
