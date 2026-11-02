#!/usr/bin/env python3
"""
Angular (ray-set) error of a corrected-model run, estimated from the run alone:

    symmetryCheck.py [--figure out.png] label=<case>/post:<time> [label=...]

The case is symmetric under the 8 symmetries of the square cross-section (reflections and
quarter turns of y, z about the channel centre): uniform inflow over x = 0, centred square window,
beam along x, symmetric mesh. The fvDOM ray set is not (it is rotated so that one ray carries the
beam), so the part of a solution that is not symmetric is angular-discretisation error.
Prints, per region, the cell T and per wall/join patch the face T and the absorbed radiation qr,
minus their symmetric average (area/volume-weighted rms, max and where), and the rms of the
half-differences to the y-mirror, z-mirror and y<->z (diagonal) images separately. With an odd
nTheta (beam ray on the equator at an azimuthal band centre, channel rotated about z only) both
mirrors are symmetries of the ray set too, so only the diagonal part can remain.
With --figure: the inlet walls unwrapped (x vs perimeter), T minus its symmetric average, one row
per run, and the wall temperature around the perimeter on the face row next to the front face.
"""
import os, sys
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, patch_geometry
from contours import frame, grid2d, plt, SERIES, STYLES, INK, INK2, GRID
from matplotlib.colors import LinearSegmentedColormap, Normalize

REGIONS = ['inletRegion', 'porousRegion', 'exitRegion']
WALLS = {'inletRegion': 'inletWall', 'porousRegion': 'porousWall', 'exitRegion': 'exitWall'}
D4 = [lambda y, z: (y, z), lambda y, z: (50 - y, z), lambda y, z: (y, 50 - z), lambda y, z: (50 - y, 50 - z),
      lambda y, z: (z, y), lambda y, z: (50 - z, y), lambda y, z: (z, 50 - y), lambda y, z: (50 - z, 50 - y)]
DIVERGING = LinearSegmentedColormap.from_list('blue_red', ['#184f95', '#3987e5', '#9ec5f4', '#f0efec',
                                                           '#f4b6b5', '#e34948', '#9f2626'])


def images(X):
    """Index of the symmetric image of every point under each of the 8 symmetries."""
    tree = cKDTree(X)
    idx = []
    for g in D4:
        y, z = g(X[:, 1], X[:, 2])
        dist, i = tree.query(np.stack([X[:, 0], y, z], 1))
        if dist.max() > 2e-3:
            sys.exit(f'mesh is not symmetric: image {dist.max():.3g} mm away')
        idx.append(i)
    return idx


def asym(v, idx):
    return v - np.mean([v[i] for i in idx], axis=0)


def load(post, t):
    to_local = frame(post)
    out = {}
    for r in REGIONS:
        geo = patch_geometry(f'{post}/constant/{r}/polyMesh')
        T, Tb = read_field(f'{post}/{t}/{r}/T')
        qr = read_field(f'{post}/{t}/{r}/qr')[1]
        sets = {'cells': dict(X=to_local(read_field(f'{post}/{t}/{r}/C')[0]), w=read_field(f'{post}/{t}/{r}/V')[0], T=T)}
        for p in [WALLS[r]] + [p for p in geo if p.startswith(r + '_to_')]:
            Cf, Sf = geo[p]
            sets[p] = dict(X=to_local(Cf), w=np.linalg.norm(Sf, axis=1),
                           T=np.broadcast_to(Tb[p], len(Cf)), qr=np.broadcast_to(qr[p], len(Cf)))
        out[r] = sets
    return out


def perimeter(X):
    """Unwrapped position around the cross-section [mm]: y=0 wall, z=50, y=50, z=0."""
    y, z = X[:, 1], X[:, 2]
    return np.select([y < 0.01, z > 49.99, y > 49.99], [z, 50 + y, 150 - z], 200 - y)


def figure(runs, path):
    if len(runs) > len(SERIES):
        sys.exit(f'figure: at most {len(SERIES)} runs (categorical colours are never cycled)')
    fig, axes = plt.subplots(len(runs) + 1, 1, figsize=(11, 2.6 * len(runs) + 3.4), constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[1] * len(runs) + [1.1]))
    norm = Normalize(-30, 30)
    for ax, (label, s) in zip(axes, runs):
        X, d = s['X'], s['dT']
        A, B, V = grid2d(X[:, 0], perimeter(X), d)
        im = ax.pcolormesh(A, B, V, cmap=DIVERGING, norm=norm, shading='nearest', rasterized=True)
        for b in (50, 100, 150):
            ax.axhline(b, color=INK2, lw=0.6)
        ax.set_yticks([25, 75, 125, 175], ['y = 0', 'z = 50', 'y = 50', 'z = 0'])
        ax.set_xlim(0, 75); ax.set_ylim(0, 200)
        ax.set_title(f'{label}:  rms {s["rms"]:.1f} K,  max |dT| {s["max"]:.0f} K', loc='left')
    axes[len(runs) - 1].set_xlabel('x from the inlet [mm]   (front face of the porous zone at 75 mm)')
    cb = fig.colorbar(im, ax=axes[:len(runs)], shrink=0.9, pad=0.01, extend='both')
    cb.set_label('T minus its symmetric average [K]')
    ax = axes[-1]
    for k, (label, s) in enumerate(runs):
        row = s['X'][:, 0] > s['X'][:, 0].max() - 1e-3
        p = perimeter(s['X'][row]); o = np.argsort(p)
        ax.plot(p[o], s['T'][row][o], color=SERIES[k], ls=STYLES[k], lw=2, label=label)
    for b in (50, 100, 150):
        ax.axvline(b, color=INK2, lw=0.8, ls=(0, (4, 3)))
    ax.set_xticks([25, 75, 125, 175], ['y = 0 wall', 'z = 50 wall', 'y = 50 wall', 'z = 0 wall'])
    ax.set_xlim(0, 200); ax.set_ylabel('Wall T [K]')
    ax.set_title(f'Inlet-wall temperature on the face row next to the front face (x = {s["X"][:, 0].max():.2f} mm)', loc='left')
    ax.grid(True, color=GRID, lw=0.6); ax.legend(frameon=False, fontsize=9)
    fig.suptitle('Inlet walls, unwrapped: departure from the square symmetry of the case = angular-discretisation error',
                 x=0.01, ha='left', fontsize=11)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    args = sys.argv[1:]
    fig = None
    if args and args[0] == '--figure':
        fig, args = args[1], args[2:]
    runs = []
    idx = {}
    for arg in args:
        label, rest = arg.split('=', 1)
        post, t = rest.rsplit(':', 1)
        data = load(post, t)
        print(f'== {label}  ({post} @ {t})')
        for r in REGIONS:
            for name, s in data[r].items():
                key = (r, name)
                if key not in idx:
                    idx[key] = images(s['X'])
                for f in ['T', 'qr'] if 'qr' in s else ['T']:
                    sc, unit = (1e3, 'kW/m2') if f == 'qr' else (1, 'K')
                    v = np.asarray(s[f], float) / sc
                    d = asym(v, idx[key])
                    i = np.argmax(np.abs(d))
                    rms = np.sqrt((s['w'] * d**2).sum() / s['w'].sum())
                    part = [np.sqrt((s['w'] * (0.5 * (v - v[idx[key][g]]))**2).sum() / s['w'].sum()) for g in (1, 2, 4)]
                    print(f'  {r:13s} {name:28s} {f:2s}  rms {rms:8.3g}  max {abs(d[i]):8.3g} {unit:5s} at x,y,z ({s["X"][i, 0]:6.2f} {s["X"][i, 1]:5.2f} {s["X"][i, 2]:5.2f})'
                          f'   rms by image: y-mirror {part[0]:7.3g}  z-mirror {part[1]:7.3g}  diagonal {part[2]:7.3g}')
                    if (r, name, f) == ('inletRegion', 'inletWall', 'T'):
                        runs.append((label, dict(X=s['X'], T=np.asarray(s['T'], float), dT=d, rms=rms, max=abs(d[i]))))
    if fig:
        figure(runs, fig)
        print('figure written to', fig)


if __name__ == '__main__':
    main()
