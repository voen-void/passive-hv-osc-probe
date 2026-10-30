"""Read the probe geometry and netlist straight from a .kicad_pcb file.

Coordinates are KiCad board coordinates in mm (x right, y down). The probe
axis runs along x. Extracted:

  * copper polygons per (net, layer): zone fills, pads, tracks, vias
  * plated holes (pads/vias spanning both copper layers)
  * board outline (outer loop + internal slots) from Edge.Cuts
  * the shield tube from User.1 (two rectangles = tube wall in cross section)
  * the resin outline from User.2 (revolved around the tube axis)
  * two-terminal parts (R, C, trimmer, pot) with their nets and values
"""

import math
import re
from dataclasses import dataclass, field

import numpy as np

import sexpr
from sexpr import child, children, num, value, xy


# ----------------------------------------------------------------------------
# small geometry helpers

def rotate(px, py, deg):
    """KiCad rotation (y axis points down, positive angle = CCW on screen)."""
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return px * c + py * s, -px * s + py * c


def arc_points(start, mid, end, seg_len=0.25):
    """Polyline approximating the circular arc start -> mid -> end."""
    (x1, y1), (x2, y2), (x3, y3) = start, mid, end
    d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < 1e-12:
        return [start, end]
    ux = ((x1**2 + y1**2) * (y2 - y3) + (x2**2 + y2**2) * (y3 - y1) + (x3**2 + y3**2) * (y1 - y2)) / d
    uy = ((x1**2 + y1**2) * (x3 - x2) + (x2**2 + y2**2) * (x1 - x3) + (x3**2 + y3**2) * (x2 - x1)) / d
    r = math.hypot(x1 - ux, y1 - uy)
    a1 = math.atan2(y1 - uy, x1 - ux)
    a2 = math.atan2(y2 - uy, x2 - ux)
    a3 = math.atan2(y3 - uy, x3 - ux)
    # sweep from a1 to a3 passing through a2
    def norm(a):
        return a % (2 * math.pi)
    ccw = norm(a2 - a1) < norm(a3 - a1)
    sweep = norm(a3 - a1) if ccw else -norm(a1 - a3)
    n = max(4, int(abs(sweep) * r / seg_len))
    return [(ux + r * math.cos(a1 + sweep * t / n), uy + r * math.sin(a1 + sweep * t / n)) for t in range(n + 1)]


def circle_poly(cx, cy, r, n=24):
    return [(cx + r * math.cos(2 * math.pi * i / n), cy + r * math.sin(2 * math.pi * i / n)) for i in range(n)]


def stadium_poly(cx, cy, length, width, angle_deg, n=8):
    """Oval/stadium with total length `length` (along local x) and width."""
    r = width / 2
    half = max(length / 2 - r, 0)
    pts = []
    for i in range(n + 1):
        a = -math.pi / 2 + math.pi * i / n
        pts.append((half + r * math.cos(a), r * math.sin(a)))
    for i in range(n + 1):
        a = math.pi / 2 + math.pi * i / n
        pts.append((-half + r * math.cos(a), r * math.sin(a)))
    out = []
    for px, py in pts:
        rx, ry = rotate(px, py, angle_deg)
        out.append((cx + rx, cy + ry))
    return out


def rect_poly(cx, cy, w, h, angle_deg):
    out = []
    for px, py in ((-w / 2, -h / 2), (w / 2, -h / 2), (w / 2, h / 2), (-w / 2, h / 2)):
        rx, ry = rotate(px, py, angle_deg)
        out.append((cx + rx, cy + ry))
    return out


def capsule_poly(p1, p2, width, n=8):
    (x1, y1), (x2, y2) = p1, p2
    length = math.hypot(x2 - x1, y2 - y1)
    ang = -math.degrees(math.atan2(y2 - y1, x2 - x1))
    return stadium_poly((x1 + x2) / 2, (y1 + y2) / 2, length + width, width, ang, n)


def chain_loops(segments, tol=1e-3):
    """Join polyline pieces [(p0..pn), ...] into closed loops."""
    pieces = [list(s) for s in segments]
    loops = []
    while pieces:
        loop = pieces.pop(0)
        changed = True
        while changed and pieces:
            changed = False
            end = loop[-1]
            for i, p in enumerate(pieces):
                if math.dist(end, p[0]) < tol:
                    loop += p[1:]
                elif math.dist(end, p[-1]) < tol:
                    loop += p[::-1][1:]
                else:
                    continue
                pieces.pop(i)
                changed = True
                break
        if math.dist(loop[0], loop[-1]) < tol:
            loop = loop[:-1]
        loops.append(loop)
    return loops


def poly_area(pts):
    a = 0.0
    for (x1, y1), (x2, y2) in zip(pts, pts[1:] + pts[:1]):
        a += x1 * y2 - x2 * y1
    return abs(a) / 2


# ----------------------------------------------------------------------------

SI = {"f": 1e-15, "p": 1e-12, "n": 1e-9, "u": 1e-6, "µ": 1e-6, "m": 1e-3,
      "k": 1e3, "K": 1e3, "M": 1e6, "G": 1e9, "R": 1.0, "": 1.0}


def parse_si(s):
    """'47p' -> 47e-12, '2M' -> 2e6, '4k7' -> 4700, '5.5pF' -> 5.5e-12."""
    s = s.strip().replace("Ω", "").replace("ohm", "")
    m = re.fullmatch(r"([0-9]*\.?[0-9]+)\s*([fpnuµmkKMGR]?)([0-9]*)\s*[FfHh]?", s)
    if not m:
        raise ValueError(f"cannot parse value '{s}'")
    mant, pre, tail = m.groups()
    if tail:  # 4k7 notation
        mant = f"{mant}.{tail}"
    return float(mant) * SI[pre]


@dataclass
class Part:
    ref: str
    value: str
    footprint: str
    layer: str
    pads: dict            # pad number -> net name
    mpn: str = ""

    @property
    def nets(self):
        return sorted(set(self.pads.values()))


@dataclass
class Copper:
    net: str
    layer: str            # "F.Cu" or "B.Cu"
    polys: list           # filled regions (lists of (x, y))
    lines: list = field(default_factory=list)  # tracks: (p1, p2, width)


@dataclass
class Board:
    path: str
    thickness: float
    outline: list          # outer loop
    holes: list            # internal slots/cutouts
    copper: dict           # (net, layer) -> Copper
    barrels: list          # (net, polygon) spanning both copper layers
    pad_shapes: list       # (ref, footprint, layer, net, polygon) for SMD pads
    parts: dict            # ref -> Part
    tube: dict             # x0, x1, yc, r_in, wall
    resin: list            # closed outline polygon in (x, y)
    texts: list

    @property
    def nets(self):
        s = set()
        for net, _ in self.copper:
            s.add(net)
        for net, _ in self.barrels:
            s.add(net)
        return sorted(s)


def _graphic_polyline(g):
    kind = g[0]
    if kind == "gr_line":
        return [xy(g, "start"), xy(g, "end")]
    if kind == "gr_arc":
        return arc_points(xy(g, "start"), xy(g, "mid"), xy(g, "end"))
    if kind == "gr_rect":
        (x1, y1), (x2, y2) = xy(g, "start"), xy(g, "end")
        return [(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)]
    if kind == "gr_circle":
        (cx, cy), (ex, ey) = xy(g, "center"), xy(g, "end")
        pts = circle_poly(cx, cy, math.hypot(ex - cx, ey - cy), 64)
        return pts + pts[:1]
    if kind == "gr_poly":
        pts = [(num(p[1]), num(p[2])) for p in child(g, "pts")[1:] if p[0] == "xy"]
        return pts + pts[:1]
    return None


def _pad_polygon(pad, fx, fy, fangle):
    lx, ly = num(child(pad, "at")[1]), num(child(pad, "at")[2])
    at = child(pad, "at")
    pangle = num(at[3]) if len(at) > 3 else 0.0  # absolute angle in file
    ox, oy = rotate(lx, ly, fangle)
    cx, cy = fx + ox, fy + oy
    w, h = num(child(pad, "size")[1]), num(child(pad, "size")[2])
    shape = pad[3]
    if shape == "circle":
        return circle_poly(cx, cy, w / 2)
    if shape == "oval":
        if w >= h:
            return stadium_poly(cx, cy, w, h, pangle)
        return stadium_poly(cx, cy, h, w, pangle + 90)
    # rect, roundrect, trapezoid, custom -> bounding rectangle
    return rect_poly(cx, cy, w, h, pangle)


def load_board(path):
    pcb = sexpr.load(path)

    thickness = num(value(child(pcb, "general"), "thickness", "1.6"))

    copper = {}

    def add_poly(net, layer, poly):
        key = (net, layer)
        copper.setdefault(key, Copper(net, layer, []))
        copper[key].polys.append(poly)

    def add_line(net, layer, p1, p2, w):
        key = (net, layer)
        copper.setdefault(key, Copper(net, layer, []))
        copper[key].lines.append((p1, p2, w))

    # --- zones (use the filled polygons KiCad saved) ------------------------
    for z in children(pcb, "zone"):
        net_node = child(z, "net")
        if net_node is None or child(z, "keepout") is not None:
            continue
        net = net_node[1]
        for fp in children(z, "filled_polygon"):
            layer = value(fp, "layer")
            if layer not in ("F.Cu", "B.Cu"):
                continue
            pts = [(num(p[1]), num(p[2])) for p in child(fp, "pts")[1:] if p[0] == "xy"]
            add_poly(net, layer, pts)

    # --- tracks & vias ------------------------------------------------------
    barrels = []
    for s in children(pcb, "segment"):
        add_line(child(s, "net")[1], value(s, "layer"), xy(s, "start"), xy(s, "end"), num(value(s, "width")))
    for s in children(pcb, "arc"):
        pts = arc_points(xy(s, "start"), xy(s, "mid"), xy(s, "end"))
        for p1, p2 in zip(pts, pts[1:]):
            add_line(child(s, "net")[1], value(s, "layer"), p1, p2, num(value(s, "width")))
    for v in children(pcb, "via"):
        net = child(v, "net")[1]
        (x, y) = xy(v, "at")
        poly = circle_poly(x, y, num(value(v, "size")) / 2, 12)
        for layer in ("F.Cu", "B.Cu"):
            add_poly(net, layer, poly)
        barrels.append((net, poly))

    # --- footprints ---------------------------------------------------------
    parts = {}
    pad_shapes = []
    for fp in children(pcb, "footprint"):
        props = {p[1]: p[2] for p in children(fp, "property")}
        ref = props.get("Reference", "?")
        at = child(fp, "at")
        fx, fy = num(at[1]), num(at[2])
        fangle = num(at[3]) if len(at) > 3 else 0.0
        fp_layer = value(fp, "layer")
        pads = {}
        for pad in children(fp, "pad"):
            net_node = child(pad, "net")
            net = net_node[1] if net_node else None
            layers = child(pad, "layers")[1:]
            poly = _pad_polygon(pad, fx, fy, fangle)
            if net is None:
                continue
            pads.setdefault(pad[1], net)
            through = pad[2] in ("thru_hole",) or "*.Cu" in layers
            cu_layers = ["F.Cu", "B.Cu"] if through else [l for l in layers if l in ("F.Cu", "B.Cu")]
            for l in cu_layers:
                add_poly(net, l, poly)
            if through:
                barrels.append((net, poly))
            else:
                pad_shapes.append((ref, fp[1], cu_layers[0], net, poly))
        if ref[0] in "RC" or ref.startswith("RV"):
            parts[ref] = Part(ref, props.get("Value", ""), fp[1], fp_layer, pads, props.get("Description", ""))

    # --- graphics: edge cuts, shield, resin ---------------------------------
    edge, shield, resin, texts = [], [], [], []
    for g in pcb[1:]:
        if not (isinstance(g, list) and g and g[0].startswith("gr_")):
            continue
        layer = value(g, "layer")
        if g[0] == "gr_text":
            texts.append((layer, g[1]))
            continue
        pl = _graphic_polyline(g)
        if pl is None:
            continue
        if layer == "Edge.Cuts":
            edge.append(pl)
        elif layer == "User.1":
            shield.append(pl)
        elif layer == "User.2":
            resin.append(pl)

    loops = chain_loops(edge)
    loops.sort(key=poly_area, reverse=True)
    outline, holes = loops[0], loops[1:]

    tube = _tube_from_shield(shield)
    resin_loops = chain_loops(resin)
    resin_loops.sort(key=poly_area, reverse=True)
    resin_outline = resin_loops[0] if resin_loops else []

    return Board(path, thickness, outline, holes, copper, barrels, pad_shapes,
                 parts, tube, resin_outline, texts)


def _tube_from_shield(shield):
    """The shield is drawn as two rectangles (upper and lower tube wall)."""
    if len(shield) < 2:
        raise ValueError("expected two rectangles on User.1 describing the tube wall")
    boxes = []
    for pl in shield:
        xs = [p[0] for p in pl]
        ys = [p[1] for p in pl]
        boxes.append((min(xs), max(xs), min(ys), max(ys)))
    boxes.sort(key=lambda b: b[2])
    top, bot = boxes[0], boxes[-1]
    r_in = (bot[2] - top[3]) / 2
    return {
        "x0": min(top[0], bot[0]),
        "x1": max(top[1], bot[1]),
        "yc": (bot[2] + top[3]) / 2,
        "r_in": r_in,
        "wall": ((top[3] - top[2]) + (bot[3] - bot[2])) / 2,
    }


# ----------------------------------------------------------------------------
# divider topology helpers

def ladder_order(board, hv_net, out_net):
    """Walk the divider from the HV net to the output net through R parts.

    Returns the list of nets [hv, n1, ..., out] in order.
    """
    adj = {}
    for p in board.parts.values():
        if p.ref.startswith("R") and not p.ref.startswith("RV") and len(p.nets) == 2:
            a, b = p.nets
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
    order = [hv_net]
    prev = None
    cur = hv_net
    while cur != out_net:
        nxt = [n for n in adj.get(cur, ()) if n != prev and n not in order]
        if len(nxt) != 1:
            raise ValueError(f"divider chain is ambiguous or broken at {cur}: {nxt}")
        prev, cur = cur, nxt[0]
        order.append(cur)
    return order


def net_labels(board, hv_net, out_net):
    """Friendly short names: HV, N01..Nxx, OUT, GND, others sanitised."""
    labels = {"GND": "GND"}
    try:
        order = ladder_order(board, hv_net, out_net)
    except ValueError:
        order = [hv_net, out_net]
    labels[hv_net] = "HV"
    labels[out_net] = "OUT"
    for i, n in enumerate(order[1:-1], 1):
        labels[n] = f"N{i:02d}"
    for n in board.nets + [n for p in board.parts.values() for n in p.nets]:
        if n not in labels:
            labels[n] = re.sub(r"[^A-Za-z0-9_]", "_", n).strip("_")
    return labels


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Print what the tools read from a .kicad_pcb file")
    ap.add_argument("pcb")
    ap.add_argument("--hv-net", default="Net-(C14-Pad1)")
    ap.add_argument("--out-net", default="Net-(J2-In)")
    a = ap.parse_args()
    b = load_board(a.pcb)
    t = b.tube
    print(f"board thickness {b.thickness} mm, outline {len(b.outline)} pts, {len(b.holes)} slots")
    print(f"tube: x {t['x0']:.2f}..{t['x1']:.2f} mm, axis y={t['yc']:.3f}, ID {2*t['r_in']:.2f} mm, wall {t['wall']:.2f} mm")
    xs = [p[0] for p in b.resin]
    print(f"resin outline: x {min(xs):.2f}..{max(xs):.2f} mm")
    labels = net_labels(b, a.hv_net, a.out_net)
    print("nets with copper:")
    for n in b.nets:
        print(f"  {labels[n]:>8}  {n}")
    print("parts:")
    for p in sorted(b.parts.values(), key=lambda p: (p.ref[0], len(p.ref), p.ref)):
        print(f"  {p.ref:>5} {p.value:>12} {p.footprint:<45} {' - '.join(labels[n] for n in p.nets)}")
