#!/usr/bin/env python3
"""
Contour and profile figures comparing VSR cases (reconstructed, un-rotated, ascii, with C written;
for the coarse cases use <case>/post from Allpost; profiles CSVs from profiles.py):

    contours.py <outDir> label=<caseDir>:<time>:<profiles.csv> [label=...]

Figures (one row per case, identical colour scales, so the cases can be compared directly):
  midplane_T.png, midplane_G.png, midplane_U.png, midplane_p.png   mid-plane (z = 25 mm) fields
  porous_faces_T.png                                                porous front (x = 75 mm) and back (x = 125 mm) faces
  profiles.png                                                      bulk/centreline T, wall T, wall radiative flux along x
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.tri import Triangulation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, patch_geometry

ALL_REGIONS = ['inletRegion', 'interface', 'porousRegion', 'exitRegion']
FACE_PATCHES = {'front': 'porousRegion_to_inletRegion', 'back': 'porousRegion_to_exitRegion'}   # corrected model
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
ORANGE = LinearSegmentedColormap.from_list('orange_seq', ['#fdf0e9', '#f9c9b0', '#f29a72', '#eb6834', '#c24e1f', '#8f3612', '#5c2008'])
BLUE = LinearSegmentedColormap.from_list('blue_seq', ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'])
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']   # categorical slots 1-4, fixed order, never cycled
STYLES = ['-', (0, (5, 3)), (0, (1, 2)), (0, (6, 2, 1, 2))]   # secondary encoding, one per case

plt.rcParams.update({'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
                     'text.color': INK, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.edgecolor': '#b8b7b2', 'font.size': 10, 'axes.titlesize': 11})


def frame(case):
    """Channel frame taken from the mesh itself (handles any residual rigid rotation, roll included):
    axis e from the inlet-plane normal, e3 from the normals of the z = 50 mm wall, e2 = e3 x e;
    origin such that the inlet plane is x = 0 and the channel spans y, z = 0..50 mm."""
    geo = patch_geometry(f'{case}/constant/inletRegion/polyMesh')
    Cf = np.vstack([geo[p][0] for p in ('inlet', 'opticalWindow')])
    Sf = np.vstack([geo[p][1] for p in ('inlet', 'opticalWindow')])
    A = np.linalg.norm(Sf, axis=1)
    e = -Sf.sum(0); e /= np.linalg.norm(e)
    P0 = (Cf * A[:, None]).sum(0) / A.sum()                      # centre of the inlet plane
    Sw = geo['inletWall'][1]
    nw = Sw / np.linalg.norm(Sw, axis=1)[:, None]
    top = nw @ np.array([0.0, 0.0, 1.0]) > 0.9                     # faces of the z = 50 mm wall
    e3 = nw[top].mean(0); e3 -= (e3 @ e) * e; e3 /= np.linalg.norm(e3)
    e2 = np.cross(e3, e)
    return lambda X: np.stack([(X - P0) @ e, (X - P0) @ e2 + 0.025, (X - P0) @ e3 + 0.025], axis=-1) * 1e3


def load(case, t):
    to_local = frame(case)
    cells = {k: [] for k in ('x', 'y', 'z', 'T', 'G', 'U', 'p')}
    faces = {}
    for r in [r for r in ALL_REGIONS if os.path.isdir(f'{case}/constant/{r}/polyMesh')]:
        L = to_local(read_field(f'{case}/{t}/{r}/C')[0])
        T, Tb = read_field(f'{case}/{t}/{r}/T')
        G = read_field(f'{case}/{t}/{r}/G')[0]
        U = read_field(f'{case}/{t}/{r}/U')[0]
        p = read_field(f'{case}/{t}/{r}/p_rgh')[0]
        for k, v in zip(('x', 'y', 'z'), L.T):
            cells[k].append(v)
        cells['T'].append(T); cells['G'].append(np.broadcast_to(G, T.shape))
        cells['U'].append(np.linalg.norm(U, axis=1)); cells['p'].append(np.broadcast_to(p, T.shape))
        if r == 'porousRegion':
            geo = patch_geometry(f'{case}/constant/{r}/polyMesh')
            for name, xs in (('front', L[:, 0].min()), ('back', L[:, 0].max())):
                if FACE_PATCHES[name] in geo:          # opaque face: its own temperature
                    F = to_local(geo[FACE_PATCHES[name]][0])
                    Tf = np.broadcast_to(Tb[FACE_PATCHES[name]], len(F))
                    faces[name] = (F[:, 1], F[:, 2], Tf, F[:, 0].mean(), 'face')
                else:                                   # reference layout: first/last porous cell layer
                    m = np.abs(L[:, 0] - xs) < 1e-3
                    faces[name] = (L[m, 1], L[m, 2], T[m], xs, 'cell layer')
    return {k: np.concatenate(v) for k, v in cells.items()}, faces


def midplane(c):
    """Cells of the two layers bracketing z = 25 mm, averaged onto the mid-plane."""
    d = np.abs(c['z'] - 25.0)
    band = 1.5 * d.min()
    m = d <= band
    ix, _ = cluster(c['x'][m])
    iy, _ = cluster(c['y'][m])
    _, inv = np.unique(ix * (iy.max() + 1) + iy, return_inverse=True)
    out = {}
    for f in ('x', 'y', 'T', 'G', 'U', 'p'):
        out[f] = np.bincount(inv, weights=c[f][m]) / np.bincount(inv)
    return out


def cluster(v, tol=0.02):
    """Index of each value's cluster (values within tol mm of their neighbour belong together) and the cluster centres."""
    order = np.argsort(v)
    vs = v[order]
    cid_sorted = np.cumsum(np.concatenate([[True], np.diff(vs) > tol])) - 1
    cid = np.empty_like(cid_sorted)
    cid[order] = cid_sorted
    return cid, np.bincount(cid, weights=v) / np.bincount(cid)


def grid2d(a, b, v):
    """Values at points (a, b) of a tensor grid -> 2-D arrays for pcolormesh (NaN where missing)."""
    ia, ua = cluster(a)
    ib, ub = cluster(b)
    V = np.full((len(ub), len(ua)), np.nan)
    V[ib, ia] = v
    A, B = np.meshgrid(ua, ub)
    return A, B, np.ma.masked_invalid(V)


def outline(ax):
    kw = dict(color=INK, lw=0.8, ls=(0, (4, 3)))
    for x in (75.0, 125.0):
        ax.axvline(x, **kw)
    ax.plot([0, 0], [12.5, 37.5], color=INK, lw=3, solid_capstyle='butt')
    ax.set_xlim(0, 275); ax.set_ylim(0, 50); ax.set_aspect('equal')
    ax.set_ylabel('y [mm]')


def field_figure(cases, field, title, unit, cmap, norm, path):
    fig, axes = plt.subplots(len(cases), 1, figsize=(12, 1.0 + 2.25 * len(cases)), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, (label, mp) in zip(axes, cases):
        X, Y, V = grid2d(mp['x'], mp['y'], mp[field])
        if field == 'G':
            neg = np.ma.masked_where(~(V < 0), np.ones_like(V))
            V = np.ma.masked_where(V <= 0, V)
        im = ax.pcolormesh(X, Y, V, cmap=cmap, norm=norm, shading='gouraud', rasterized=True)
        if field == 'G' and neg.count():
            ax.pcolormesh(X, Y, neg, cmap=LinearSegmentedColormap.from_list('neg', ['#8f8e89', '#8f8e89']),
                          shading='nearest', rasterized=True)
        outline(ax)
        v = mp[field]
        ax.set_title(f'{label}    min {v.min():.4g}   max {v.max():.4g} {unit}', loc='left')
    axes[-1].set_xlabel('x from the inlet [mm]')
    cb = fig.colorbar(im, ax=axes, shrink=0.9, pad=0.01)
    cb.set_label(f'{title} [{unit}]')
    note = '   gray: G < 0 (non-physical)' if field == 'G' else ''
    fig.suptitle(f'{title}, mid-plane z = 25 mm\nThick line at x = 0: window.  Dashed: porous zone (75-125 mm).{note}',
                 x=0.01, ha='left', fontsize=10.5)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def faces_figure(cases, path):
    vals = np.concatenate([np.concatenate([np.asarray(f[k][2]) for k in ('front', 'back')]) for _, f in cases])
    norm = Normalize(300, max(vals.max(), 301))
    fig, axes = plt.subplots(len(cases), 2, figsize=(8.5, 0.6 + 3.6 * len(cases)), constrained_layout=True)
    axes = np.atleast_2d(axes)
    for row, (label, f) in zip(axes, cases):
        for ax, k in zip(row, ('front', 'back')):
            y, z, T, xs, kind = f[k]
            name = f'{k} face' if kind == 'face' else ('first' if k == 'front' else 'last') + ' porous cell layer'
            A, B, V = grid2d(y, z, T)
            im = ax.pcolormesh(A, B, V, cmap=ORANGE, norm=norm, shading='gouraud', rasterized=True)
            ax.plot([12.5, 37.5, 37.5, 12.5, 12.5], [12.5, 12.5, 37.5, 37.5, 12.5], color=INK, lw=0.8, ls=(0, (4, 3)))
            ax.set_aspect('equal'); ax.set_xlim(0, 50); ax.set_ylim(0, 50)
            ax.set_title(f'{label}\n{name} (x = {xs:.1f} mm): {T.min():.0f}-{T.max():.0f} K', loc='left', fontsize=10)
            ax.set_xlabel('y [mm]'); ax.set_ylabel('z [mm]')
    cb = fig.colorbar(im, ax=axes, shrink=0.8); cb.set_label('Temperature [K]')
    fig.suptitle('Porous zone, front and back faces (dashed: beam footprint)', x=0.01, ha='left', fontsize=12)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def profiles_figure(profs, path):
    fig, axes = plt.subplots(3, 1, figsize=(11, 10.5), sharex=True, constrained_layout=True)
    if len(profs) > len(STYLES):
        sys.exit(f'profiles.png: at most {len(STYLES)} cases')
    panels = [
        ('Fluid temperature [K]', [('bulkT', 'bulk (mass-flow weighted)'), ('Tcentre', 'centreline')]),
        ('Wall temperature [K]', [('wallT_mean', 'wall, perimeter mean'), ('wallT_max', 'wall, maximum')]),
        ('Radiation absorbed by the walls [kW/m2]', [('wall_qr', 'wall, perimeter mean')]),
    ]
    for ax, (ylabel, cols) in zip(axes, panels):
        for ci, (col, cname) in enumerate(cols):
            for si, (label, d) in enumerate(profs):
                v = d[col] / (1e3 if col == 'wall_qr' else 1)
                ax.plot(d['x_mm'], v, color=SERIES[ci], ls=STYLES[si], lw=2, label=f'{cname} - {label}')
        for x in (75, 125):
            ax.axvline(x, color=INK2, lw=0.8, ls=(0, (4, 3)))
        ax.grid(True, color=GRID, lw=0.6); ax.set_ylabel(ylabel)
        ax.legend(frameon=False, fontsize=8, loc='upper right')
    axes[0].set_title('Axial profiles (dashed verticals: porous zone 75-125 mm)', loc='left')
    axes[-1].set_xlabel('x from the inlet [mm]'); axes[-1].set_xlim(0, 275)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    mids, faces, profs = [], [], []
    for arg in sys.argv[2:]:
        label, rest = arg.split('=', 1)
        case, t, csv = rest.split(':')
        c, f = load(case, t)
        mids.append((label, midplane(c)))
        faces.append((label, f))
        d = np.genfromtxt(csv, delimiter=',', names=True)
        profs.append((label, d))
    allv = lambda k: np.concatenate([m[k] for _, m in mids])
    field_figure(mids, 'T', 'Temperature', 'K', ORANGE, Normalize(300, allv('T').max()), f'{out}/midplane_T.png')
    G = allv('G')
    field_figure(mids, 'G', 'Incident radiation G', 'W/m2', BLUE, LogNorm(1e3, max(G.max(), 1e4)), f'{out}/midplane_G.png')
    field_figure(mids, 'U', 'Velocity magnitude', 'm/s', BLUE, Normalize(0, allv('U').max()), f'{out}/midplane_U.png')
    field_figure(mids, 'p', 'Pressure p_rgh (outlet = 0)', 'Pa', BLUE, Normalize(min(allv('p').min(), 0), allv('p').max()), f'{out}/midplane_p.png')
    faces_figure(faces, f'{out}/porous_faces_T.png')
    profiles_figure(profs, f'{out}/profiles.png')
    print('figures written to', out)


if __name__ == '__main__':
    main()
