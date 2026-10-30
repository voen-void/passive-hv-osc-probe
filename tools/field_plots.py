"""Plots for the field solver: model check and E-field stress maps."""

import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap

import fieldsolver as fs

MAT_COLORS = ListedColormap(["#ffffff", "#9ec5e8", "#c8d98a"])   # air, resin, FR4


def _outlines_xy(ax, m):
    b, t = m.board, m.tube
    for loop in [b.outline] + b.holes:
        p = np.array(loop + loop[:1])
        ax.plot(p[:, 0], p[:, 1], color="k", lw=0.6)
    if b.resin:
        p = np.array(b.resin + b.resin[:1])
        ax.plot(p[:, 0], p[:, 1], color="tab:blue", lw=0.8, ls="--")
    for yy in (t["yc"] - t["r_in"], t["yc"] + t["r_in"]):
        ax.plot([t["x0"], t["x1"]], [yy, yy], color="0.3", lw=2)


def _outlines_xz(ax, m):
    b, t = m.board, m.tube
    bx = [p[0] for p in b.outline]
    ax.plot([min(bx), max(bx), max(bx), min(bx), min(bx)], [m.z_b, m.z_b, m.z_f, m.z_f, m.z_b], color="k", lw=0.6)
    for zz in (-t["r_in"], t["r_in"]):
        ax.plot([t["x0"], t["x1"]], [zz, zz], color="0.3", lw=2)
    if b.resin:
        p = np.array(b.resin + b.resin[:1])
        ax.plot(p[:, 0], p[:, 1] - t["yc"], color="tab:blue", lw=0.8, ls="--")


def _view_limits(m):
    b, t = m.board, m.tube
    bx = [p[0] for p in b.outline]
    r = t["r_in"] + t["wall"] + 2
    return (min(bx) - 5, max(bx) + 5), (t["yc"] - r, t["yc"] + r), (-r, r)


def plot_model(m, path, labels):
    """Three views of the discretised model, to check the import."""
    xl, yl, zl = _view_limits(m)
    kB = int(np.argmin(np.abs(m.z - m.z_b)))
    fig, axs = plt.subplots(3, 1, figsize=(16, 13), gridspec_kw={"height_ratios": [1, 1, 1.3]})

    lab = m.label[:, :, kB].astype(float)
    lab[lab < 0] = np.nan
    ax = axs[0]
    ax.pcolormesh(m.x, m.y, lab.T, cmap="tab20", shading="nearest", rasterized=True)
    _outlines_xy(ax, m)
    ax.set(xlim=xl, ylim=yl[::-1], aspect="equal", title="Conductors on the B.Cu plane (colour = net)",
           xlabel="x [mm]", ylabel="y [mm]")
    for k, name in enumerate(m.names):
        ii, jj = np.nonzero(m.label[:, :, kB] == k)
        if len(ii) and k > 0:
            ax.text(m.x[int(np.median(ii))], m.y[int(np.median(jj))], labels.get(name, name),
                    fontsize=6, ha="center", va="center")

    jc = int(np.argmin(np.abs(m.y - m.tube["yc"])))
    j_cap = int(np.argmin(np.abs(m.y - 36.09)))   # through the 1808 capacitor row
    ax = axs[1]
    yc_cell = min(jc, len(m.y) - 2)
    ax.pcolormesh(m.x, m.z, m.mat[:, yc_cell, :].T, cmap=MAT_COLORS, vmin=-0.5, vmax=2.5,
                  shading="flat", rasterized=True)
    cond = (m.label[:, j_cap, :] > 0).astype(float)
    cond[cond == 0] = np.nan
    ax.pcolormesh(m.x, m.z, cond.T, cmap="autumn", shading="nearest", rasterized=True)
    gnd = (m.label[:, jc, :] == 0).astype(float)
    gnd[gnd == 0] = np.nan
    ax.pcolormesh(m.x, m.z, gnd.T, cmap="Greys", vmin=0, vmax=1.2, shading="nearest", rasterized=True)
    ax.set(xlim=xl, ylim=zl, aspect="equal", xlabel="x [mm]", ylabel="z [mm]",
           title="Side view: materials at the axis (white air, blue resin, green FR4), "
                 "conductors at the capacitor row (orange), ground (grey)")

    ax = axs[2]
    seg_x = np.median([p[0] for p in m.board.outline])
    i = int(np.argmin(np.abs(m.x - seg_x)))
    ax.pcolormesh(m.y, m.z, m.mat[min(i, len(m.x) - 2)].T, cmap=MAT_COLORS, vmin=-0.5, vmax=2.5,
                  shading="flat", rasterized=True)
    lab = m.label[i].astype(float)
    lab[lab < 0] = np.nan
    ax.pcolormesh(m.y, m.z, lab.T, cmap="tab20", shading="nearest", rasterized=True)
    Y, Z = np.meshgrid(m.y, m.z, indexing="ij")
    ax.set(xlim=yl, ylim=zl, aspect="equal", xlabel="y [mm]", ylabel="z [mm]",
           title=f"Cross section at x = {m.x[i]:.1f} mm")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_grid(m, path):
    fig, axs = plt.subplots(1, 3, figsize=(15, 3.5))
    for ax, a, n in zip(axs, (m.x, m.y, m.z), "xyz"):
        ax.semilogy((a[1:] + a[:-1]) / 2, np.diff(a), ".", ms=2)
        ax.set(xlabel=f"{n} [mm]", ylabel="step [mm]", title=f"{n}: {len(a)} nodes")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------
# stress

def stress_report(m, phi, cfg):
    """Peak fields per material [kV/mm], with and without conductor-edge cells."""
    ex, ey, ez = fs.efield(m, phi)
    E = np.sqrt(ex**2 + ey**2 + ez**2) / 1000.0
    inside = fs.cell_is_inside_conductor(m)
    adj = fs.conductor_adjacent_cells(m)
    xc, yc, zc = [(a[1:] + a[:-1]) / 2 for a in (m.x, m.y, m.z)]
    mats = cfg["materials"]
    strength = {fs.AIR: mats["air_strength"], fs.RESIN: mats["resin_strength"], fs.FR4: mats["fr4_strength"]}
    rep = {}
    for mid, name in fs.MATERIAL_NAMES.items():
        sel = (m.mat == mid) & ~inside
        if not sel.any():
            continue
        for tag, mask in (("peak", sel), ("bulk", sel & ~adj)):
            if not mask.any():
                continue
            Em = np.where(mask, E, 0)
            i = np.unravel_index(np.argmax(Em), Em.shape)
            rep[f"{name}_{tag}"] = {
                "E_kV_mm": float(Em[i]), "x": float(xc[i[0]]), "y": float(yc[i[1]]), "z": float(zc[i[2]]),
                "utilisation": float(Em[i] / strength[mid]),
            }
    # tangential field along the board surfaces (tracking / creepage)
    kz = np.nonzero((zc > m.z_f) & (zc < m.z_f + 0.5) | (zc < m.z_b) & (zc > m.z_b - 0.5))[0]
    Et = np.sqrt(ex**2 + ey**2)[:, :, kz] / 1000.0
    on_board = (m.mat[:, :, (m.mat.shape[2] // 2)] == fs.FR4)
    sel = on_board[:, :, None] & ~adj[:, :, kz] & ~inside[:, :, kz]
    Etm = np.where(sel, Et, 0)
    i = np.unravel_index(np.argmax(Etm), Etm.shape)
    rep["board_surface_tangential"] = {
        "E_kV_mm": float(Etm[i]), "x": float(xc[i[0]]), "y": float(yc[i[1]]), "z": float(zc[kz[i[2]]]),
    }
    return rep, E


def plot_stress(m, phi, E, path, title, volts):
    warnings.filterwarnings("ignore", "All-NaN slice", RuntimeWarning)
    E = E.copy()
    E[fs.cell_is_inside_conductor(m)] = np.nan
    xl, yl, zl = _view_limits(m)
    vmax = np.nanpercentile(E, 99.99)
    norm = LogNorm(vmin=max(vmax / 300, 1e-3), vmax=vmax)

    fig = plt.figure(figsize=(16, 16))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 0.6, 1, 1.1], width_ratios=[1, 1])

    ax = fig.add_subplot(gs[0, :])
    top = np.nanmax(E, axis=2)
    pc = ax.pcolormesh(m.x, m.y, top.T, norm=norm, cmap="inferno", shading="flat", rasterized=True)
    _outlines_xy(ax, m)
    ax.set(xlim=xl, ylim=yl[::-1], aspect="equal", xlabel="x [mm]", ylabel="y [mm]",
           title=f"{title}: max |E| through the thickness (top view)")
    fig.colorbar(pc, ax=ax, label="|E| [kV/mm]", fraction=0.02, pad=0.01)

    ax = fig.add_subplot(gs[1, :])
    xc = (m.x[1:] + m.x[:-1]) / 2
    for mid, col in ((fs.RESIN, "tab:blue"), (fs.FR4, "tab:green"), (fs.AIR, "0.4")):
        Em = np.where(m.mat == mid, E, np.nan)
        prof = np.nanmax(np.nanmax(Em, axis=2), axis=1) if np.isfinite(Em).any() else None
        if prof is not None:
            ax.semilogy(xc, prof, color=col, lw=1, label=fs.MATERIAL_NAMES[mid])
    ax.axvline(m.tube["x0"], color="0.5", ls=":")
    ax.axvline(m.tube["x1"], color="0.5", ls=":", label="tube ends")
    ax.set(xlim=xl, xlabel="x [mm]", ylabel="max |E| [kV/mm]", title="Peak field along the probe, per material")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)

    ax = fig.add_subplot(gs[2, :])
    side = np.nanmax(E, axis=1)
    pc = ax.pcolormesh(m.x, m.z, side.T, norm=norm, cmap="inferno", shading="flat", rasterized=True)
    _outlines_xz(ax, m)
    ax.set(xlim=xl, ylim=zl, aspect="equal", xlabel="x [mm]", ylabel="z [mm]",
           title="max |E| across the board width (side view)")
    fig.colorbar(pc, ax=ax, label="|E| [kV/mm]", fraction=0.02, pad=0.01)

    # cross sections: at the tube mouth and at the hottest resin cell
    P = phi.reshape(m.shape)
    Er = np.where(m.mat == fs.RESIN, E, np.nan)
    Er[fs.conductor_adjacent_cells(m)] = np.nan
    ihot = np.unravel_index(np.nanargmax(Er), Er.shape)[0]
    Er[xc > m.tube["x1"] - 10] = np.nan
    iin = np.unravel_index(np.nanargmax(Er), Er.shape)[0]
    for col, (i, what) in enumerate(((ihot, "hottest resin"), (iin, "hottest resin >10 mm inside tube"))):
        ax = fig.add_subplot(gs[3, col])
        pc = ax.pcolormesh(m.y, m.z, E[i].T, norm=norm, cmap="inferno", shading="flat", rasterized=True)
        cs = ax.contour(m.y, m.z, P[i].T / 1000, levels=12, colors="w", linewidths=0.5)
        ax.clabel(cs, fontsize=6, fmt="%.1f kV")
        th = np.linspace(0, 2 * np.pi, 200)
        r = m.tube["r_in"]
        ax.plot(m.tube["yc"] + r * np.cos(th), r * np.sin(th), color="0.6", lw=1)
        ax.set(xlim=yl, ylim=zl, aspect="equal", xlabel="y [mm]", ylabel="z [mm]",
               title=f"Cross section at x = {xc[i]:.1f} mm ({what})")
        fig.colorbar(pc, ax=ax, label="|E| [kV/mm]", fraction=0.04)

    fig.suptitle(f"{title} — HV tip at {volts/1000:.1f} kV", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
