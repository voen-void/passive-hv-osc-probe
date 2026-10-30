"""3D electrostatic finite-volume solver for the probe.

The domain is a box discretised by a non-uniform rectilinear grid. Potentials
live on grid nodes, permittivity on cells (node-based finite integration
technique). Conductors are sets of nodes with a fixed potential:

  * copper (zero-thickness sheets on the F.Cu / B.Cu planes)
  * plated holes / THT pins (columns of nodes through the board)
  * SMD part terminations (boxes above the pads)
  * the shield tube and the outer box (ground)

Dielectrics: FR4 board (with its slots), resin body (User.2 outline revolved
around the tube axis), air everywhere else.

Coordinates: x, y as in KiCad [mm]; z is perpendicular to the board, z = 0
on the tube axis, F.Cu on the +z side.
"""

import time
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from matplotlib.path import Path

from kicad_board import capsule_poly

EPS0 = 8.8541878128e-12

AIR, RESIN, FR4 = 0, 1, 2
MATERIAL_NAMES = {AIR: "air", RESIN: "resin", FR4: "FR4"}


# ----------------------------------------------------------------------------
# grid

def make_axis(lo, hi, regions, forced, h_max, growth):
    """Graded 1D grid.

    regions: [(a, b, h)] intervals with target step h. Away from a region
    the step grows linearly with distance (≈ geometric growth per cell),
    capped at h_max. Every value in `forced` becomes a grid line.
    """
    s = np.linspace(lo, hi, 200001)
    h = np.full_like(s, h_max)
    for a, b, hr in regions:
        d = np.maximum(a - s, 0) + np.maximum(s - b, 0)
        h = np.minimum(h, hr + (growth - 1) * d)
    pts = sorted(set([lo, hi] + [f for f in forced if lo < f < hi]))
    out = [pts[0]]
    for a, b in zip(pts, pts[1:]):
        m = (s >= a) & (s <= b)
        ss = np.concatenate(([a], s[m], [b]))
        hh = np.interp(ss, s, h)
        cum = np.concatenate(([0], np.cumsum(np.diff(ss) * 0.5 * (1 / hh[1:] + 1 / hh[:-1]))))
        n = max(1, int(np.ceil(cum[-1] - 1e-9)))
        inner = np.interp(np.arange(1, n) * cum[-1] / n, cum, ss)
        out.extend(inner.tolist())
        out.append(b)
    return np.array(out)


@dataclass
class Model:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    mat: np.ndarray        # (nx-1, ny-1, nz-1) uint8 material id
    eps: np.ndarray        # (nx-1, ny-1, nz-1) relative permittivity
    label: np.ndarray      # (nx, ny, nz) int16: -1 free, 0 GND, k conductor
    names: list            # names[k] = net name for label k (names[0] = "GND")
    z_f: float
    z_b: float
    tube: dict
    board: object

    @property
    def shape(self):
        return (len(self.x), len(self.y), len(self.z))


def _nodes_in_poly(poly, x, y, pad=0.0):
    """Indices (i, j) of grid nodes inside a polygon."""
    poly = np.asarray(poly)
    xmin, ymin = poly.min(0) - pad - 1e-9
    xmax, ymax = poly.max(0) + pad + 1e-9
    i0, i1 = np.searchsorted(x, xmin), np.searchsorted(x, xmax, "right")
    j0, j1 = np.searchsorted(y, ymin), np.searchsorted(y, ymax, "right")
    if i1 <= i0 or j1 <= j0:
        return np.empty(0, int), np.empty(0, int)
    X, Y = np.meshgrid(x[i0:i1], y[j0:j1], indexing="ij")
    inside = Path(poly).contains_points(np.c_[X.ravel(), Y.ravel()], radius=1e-6 + 2 * pad).reshape(X.shape)
    ii, jj = np.nonzero(inside)
    return ii + i0, jj + j0


def _nearest(arr, v):
    i = np.searchsorted(arr, v)
    i = np.clip(i, 1, len(arr) - 1)
    return np.where(np.abs(arr[i - 1] - v) <= np.abs(arr[i] - v), i - 1, i)


def _in_any_poly(polys, px, py):
    m = np.zeros(px.shape, bool)
    pts = np.c_[px.ravel(), py.ravel()]
    for p in polys:
        m |= Path(np.asarray(p)).contains_points(pts).reshape(px.shape)
    return m


def build_model(board, cfg, scale=1.0, verbose=True):
    gcfg, mcfg, mat_cfg = cfg["geometry"], cfg["mesh"], cfg["materials"]
    t = board.thickness
    off = gcfg.get("board_offset_z", 0.0)
    z_f, z_b = off + t / 2, off - t / 2
    tube = board.tube
    yc, r_in, wall = tube["yc"], tube["r_in"], tube["wall"]
    r_out = r_in + wall

    part_h = {}
    for ref, fp, layer, net, poly in board.pad_shapes:
        for key, h in gcfg.get("part_height", {}).items():
            if key in fp:
                part_h[ref] = h
    h_parts = max(part_h.values(), default=0.0)

    hx, hy = mcfg["h_x"] * scale, mcfg["h_y"] * scale
    hzb, hzp = mcfg["h_z_board"] * scale, mcfg["h_z_parts"] * scale
    ht, hmax, gr, margin = mcfg["h_tube"] * scale, mcfg["h_max"], mcfg["growth"], mcfg["margin"]

    bx = [p[0] for p in board.outline]
    by = [p[1] for p in board.outline]
    rx = [p[0] for p in board.resin] or bx
    x_lo = min(min(bx), tube["x0"], min(rx)) - margin
    x_hi = max(max(bx), tube["x1"], max(rx)) + margin
    x = make_axis(x_lo, x_hi, [(min(bx) - 1, max(bx) + 1, hx), (min(rx), max(rx), ht)],
                  [tube["x0"], tube["x1"]], hmax, gr)
    y = make_axis(yc - r_out - margin, yc + r_out + margin,
                  [(min(by) - 0.5, max(by) + 0.5, hy), (yc - r_out, yc + r_out, ht)],
                  [yc, yc - r_in, yc + r_in, yc - r_out, yc + r_out], hmax, gr)
    zp = h_parts + 0.6
    z_forced = [z_f, z_b, -r_in, r_in, -r_out, r_out]
    for h in set(part_h.values()):
        z_forced += [z_f + h, z_b - h]
    z = make_axis(-r_out - margin, r_out + margin,
                  [(z_b - 0.4, z_f + 0.4, hzb), (z_b - zp, z_f + zp, hzp), (-r_out, r_out, ht)],
                  z_forced, hmax, gr)
    nx, ny, nz = len(x), len(y), len(z)
    if verbose:
        print(f"grid {nx} x {ny} x {nz} = {nx*ny*nz/1e6:.2f} M nodes")

    # --- materials on cells ---------------------------------------------------
    xc, ycc, zc = (x[1:] + x[:-1]) / 2, (y[1:] + y[:-1]) / 2, (z[1:] + z[:-1]) / 2
    mat = np.full((nx - 1, ny - 1, nz - 1), AIR, np.uint8)

    if board.resin:
        # resin: revolve the outline around the axis -> test (x, yc - r)
        R = np.sqrt((ycc[:, None] - yc) ** 2 + zc[None, :] ** 2)            # (ny-1, nz-1)
        rq = np.unique(np.round(R, 4))
        path = Path(np.asarray(board.resin))
        Xg, Rg = np.meshgrid(xc, rq, indexing="ij")
        inside = path.contains_points(np.c_[Xg.ravel(), (yc - Rg).ravel()]).reshape(Xg.shape)
        inside |= path.contains_points(np.c_[Xg.ravel(), (yc + Rg).ravel()]).reshape(Xg.shape)
        idx = np.searchsorted(rq, np.round(R, 4))
        mat[inside[:, idx]] = RESIN

    Xb, Yb = np.meshgrid(xc, ycc, indexing="ij")
    on_board = _in_any_poly([board.outline], Xb, Yb) & ~_in_any_poly(board.holes, Xb, Yb)
    kz = (zc > z_b) & (zc < z_f)
    mat[on_board[:, :, None] & kz[None, None, :]] = FR4

    eps_of = np.array([mat_cfg["air_er"], mat_cfg["resin_er"], mat_cfg["fr4_er"]])
    eps = eps_of[mat]

    # --- conductors on nodes -------------------------------------------------
    names = ["GND"] + [n for n in board.nets if n != "GND"]
    lab_of = {n: i for i, n in enumerate(names)}
    label = np.full((nx, ny, nz), -1, np.int16)
    conflicts = 0

    def mark(ii, jj, kk, lab):
        nonlocal conflicts
        cur = label[ii, jj, kk]
        conflicts += int(np.count_nonzero((cur >= 0) & (cur != lab)))
        label[ii, jj, kk] = lab

    # shield tube (ground)
    ix = np.nonzero((x >= tube["x0"] - 1e-9) & (x <= tube["x1"] + 1e-9))[0]
    Rn = np.sqrt((y[:, None] - yc) ** 2 + z[None, :] ** 2)
    jw, kw = np.nonzero((Rn >= r_in - 1e-6) & (Rn <= r_out + 1e-6))
    label[ix[:, None], jw[None, :], kw[None, :]] = 0

    kF, kB = int(_nearest(z, z_f)), int(_nearest(z, z_b))
    k_of_layer = {"F.Cu": kF, "B.Cu": kB}

    for (net, layer), cu in board.copper.items():
        lab = lab_of[net]
        k = k_of_layer[layer]
        for poly in cu.polys:
            ii, jj = _nodes_in_poly(poly, x, y)
            mark(ii, jj, np.full_like(ii, k), lab)
        for p1, p2, w in cu.lines:
            ii, jj = _nodes_in_poly(capsule_poly(p1, p2, w), x, y)
            mark(ii, jj, np.full_like(ii, k), lab)
            # always keep the centre line connected even if narrower than the grid
            n = int(np.hypot(p2[0] - p1[0], p2[1] - p1[1]) / (min(hx, hy) / 3)) + 2
            sx = np.linspace(p1[0], p2[0], n)
            sy = np.linspace(p1[1], p2[1], n)
            mark(_nearest(x, sx), _nearest(y, sy), np.full(n, k), lab)

    for net, poly in board.barrels:
        ii, jj = _nodes_in_poly(poly, x, y)
        if len(ii) == 0:
            c = np.mean(poly, axis=0)
            ii, jj = np.array([_nearest(x, c[0])]), np.array([_nearest(y, c[1])])
        for k in range(kB, kF + 1):
            mark(ii, jj, np.full_like(ii, k), lab_of[net])

    for ref, fp, layer, net, poly in board.pad_shapes:
        h = part_h.get(ref, 0.0)
        if h <= 0:
            continue
        ii, jj = _nodes_in_poly(poly, x, y)
        if layer == "F.Cu":
            ks = np.nonzero((z >= z_f - 1e-9) & (z <= z_f + h + 1e-9))[0]
        else:
            ks = np.nonzero((z <= z_b + 1e-9) & (z >= z_b - h - 1e-9))[0]
        for k in ks:
            mark(ii, jj, np.full_like(ii, k), lab_of[net])

    # outer box = ground (far environment)
    label[[0, -1], :, :] = 0
    label[:, [0, -1], :] = 0
    label[:, :, [0, -1]] = 0

    if verbose:
        counts = np.bincount(label[label >= 0].ravel(), minlength=len(names))
        empty = [names[i] for i in range(1, len(names)) if counts[i] == 0]
        print(f"conductors: {len(names)-1} nets + GND, {conflicts} node conflicts")
        if empty:
            print("WARNING: nets with no grid nodes (grid too coarse?):", empty)

    return Model(x, y, z, mat, eps, label, names, z_f, z_b, tube, board)


# ----------------------------------------------------------------------------
# matrix

def assemble(m):
    """Symmetric positive (semi)definite FIT Laplacian. Units: eps0 * mm."""
    x, y, z, eps = m.x, m.y, m.z, m.eps
    nx, ny, nz = m.shape
    dx, dy, dz = np.diff(x), np.diff(y), np.diff(z)

    def pad2(a, axes):
        pw = [(0, 0)] * 3
        for ax in axes:
            pw[ax] = (1, 1)
        return np.pad(a, pw)

    # x edges: (nx-1, ny, nz)
    W = pad2(eps * dy[None, :, None] * dz[None, None, :] / 4, (1, 2))
    gx = (W[:, :-1, :-1] + W[:, 1:, :-1] + W[:, :-1, 1:] + W[:, 1:, 1:]) / dx[:, None, None]
    # y edges: (nx, ny-1, nz)
    W = pad2(eps * dx[:, None, None] * dz[None, None, :] / 4, (0, 2))
    gy = (W[:-1, :, :-1] + W[1:, :, :-1] + W[:-1, :, 1:] + W[1:, :, 1:]) / dy[None, :, None]
    # z edges: (nx, ny, nz-1)
    W = pad2(eps * dx[:, None, None] * dy[None, :, None] / 4, (0, 1))
    gz = (W[:-1, :-1, :] + W[1:, :-1, :] + W[:-1, 1:, :] + W[1:, 1:, :]) / dz[None, None, :]

    N = nx * ny * nz
    sx, sy, sz = ny * nz, nz, 1
    ox = gx.ravel()                                            # length N - sx
    oy = np.zeros((nx, ny, nz)); oy[:, :-1, :] = gy; oy = oy.ravel()[: N - sy]
    oz = np.zeros((nx, ny, nz)); oz[:, :, :-1] = gz; oz = oz.ravel()[: N - sz]
    diag = np.zeros(N)
    diag[:-sx] += ox; diag[sx:] += ox
    diag[:-sy] += oy; diag[sy:] += oy
    diag[:-sz] += oz; diag[sz:] += oz
    L = sp.diags([diag, -ox, -ox, -oy, -oy, -oz, -oz], [0, sx, -sx, sy, -sy, sz, -sz],
                 shape=(N, N), format="csr")
    return L


class Solver:
    """Holds the reduced system and its AMG hierarchy for repeated solves."""

    def __init__(self, model, tol=1e-8, verbose=True):
        import pyamg
        self.m = model
        self.tol = tol
        self.verbose = verbose
        t0 = time.time()
        L = assemble(model)
        lab = model.label.ravel()
        self.free = np.nonzero(lab < 0)[0]
        self.fixed = np.nonzero(lab >= 0)[0]
        self.fixed_lab = lab[self.fixed]
        self.L = L
        self.A = L[self.free][:, self.free].tocsr()
        self.Afc = L[self.free][:, self.fixed].tocsr()
        t1 = time.time()
        self.ml = pyamg.smoothed_aggregation_solver(self.A, symmetry="hermitian", max_coarse=1000)
        if verbose:
            print(f"assembly {t1-t0:.1f} s, AMG setup {time.time()-t1:.1f} s, "
                  f"{len(self.free)/1e6:.2f} M unknowns, {len(self.ml.levels)} levels")

    def solve(self, conductor_volts):
        """conductor_volts: array over labels (len(names)). Returns node potentials."""
        phi_c = np.asarray(conductor_volts, float)[self.fixed_lab]
        b = -self.Afc @ phi_c
        res = []
        x = self.ml.solve(b, tol=self.tol, accel="cg", maxiter=300, residuals=res)
        phi = np.empty(self.L.shape[0])
        phi[self.free] = x
        phi[self.fixed] = phi_c
        if self.verbose:
            print(f"    {len(res)-1} CG iterations, rel. residual {res[-1]/max(res[0],1e-300):.1e}")
        return phi

    def charges(self, phi):
        """Charge on each conductor [C] for potential vector phi [V]."""
        r = self.L @ phi
        q = np.bincount(self.fixed_lab, weights=r[self.fixed], minlength=len(self.m.names))
        return q * EPS0 * 1e-3

    def capacitance_matrix(self):
        """Maxwell capacitance matrix [F] including the GND row/column."""
        n = len(self.m.names)
        C = np.zeros((n, n))
        for k in range(1, n):
            t0 = time.time()
            v = np.zeros(n); v[k] = 1.0
            if self.verbose:
                print(f"  conductor {k}/{n-1}: {self.m.names[k]}", end="", flush=True)
            phi = self.solve(v)
            C[:, k] = self.charges(phi)
            if self.verbose:
                print(f"    ({time.time()-t0:.1f} s)")
        # GND column from charge neutrality of each row
        C[:, 0] = -C[:, 1:].sum(axis=1)
        return 0.5 * (C + C.T)


# ----------------------------------------------------------------------------
# electric field

def efield(m, phi):
    """|E| and components on cell centres [V/mm]."""
    P = phi.reshape(m.shape)
    dx, dy, dz = np.diff(m.x), np.diff(m.y), np.diff(m.z)
    ex = np.diff(P, axis=0)
    ex = (ex[:, :-1, :-1] + ex[:, 1:, :-1] + ex[:, :-1, 1:] + ex[:, 1:, 1:]) / 4 / dx[:, None, None]
    ey = np.diff(P, axis=1)
    ey = (ey[:-1, :, :-1] + ey[1:, :, :-1] + ey[:-1, :, 1:] + ey[1:, :, 1:]) / 4 / dy[None, :, None]
    ez = np.diff(P, axis=2)
    ez = (ez[:-1, :-1, :] + ez[1:, :-1, :] + ez[:-1, 1:, :] + ez[1:, 1:, :]) / 4 / dz[None, None, :]
    return -ex, -ey, -ez


def conductor_adjacent_cells(m):
    """Cells that touch at least one conductor node (field singularities live here)."""
    c = (m.label >= 0)
    a = c[:-1, :-1, :-1] | c[1:, :-1, :-1] | c[:-1, 1:, :-1] | c[:-1, :-1, 1:]
    a |= c[1:, 1:, :-1] | c[1:, :-1, 1:] | c[:-1, 1:, 1:] | c[1:, 1:, 1:]
    return a


def cell_is_inside_conductor(m):
    c = (m.label >= 0)
    return c[:-1, :-1, :-1] & c[1:, :-1, :-1] & c[:-1, 1:, :-1] & c[:-1, :-1, 1:] & \
        c[1:, 1:, :-1] & c[1:, :-1, 1:] & c[:-1, 1:, 1:] & c[1:, 1:, 1:]
