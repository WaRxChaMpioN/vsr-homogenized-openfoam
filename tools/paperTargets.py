#!/usr/bin/env python3
"""
Compare a run (runs/<case>) with the values reported in the papers for G 5e6, Re 300:

    paperTargets.py <case> <time> <profiles.csv> [figure.png]

<case>/post/<time> is the un-rotated copy written by Allpost (cell centres C written), the CSV comes
from profiles.py. Paper values are read from paper 2 (Table 4, Figs 12-14, text) with the read-off
uncertainty noted; paper 1 is the pore-scale result. With figure.png, the bulk and top-wall
temperature profiles are also plotted against the paper points.
"""
import os, sys, glob
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, patch_geometry
from contours import frame, plt, SERIES, INK2, GRID

BULK = [(50, 398), (69, 473), (75, 1355), (100, 1277), (125, 1126), (150, 1055), (200, 966), (275, 914)]
BULK_PORE = [(50, 427), (75, 721), (90, 1469), (125, 1214), (275, 1018.9)]      # paper 1 Figs 13/14b, read +-8-10 K
WALL_MAX = (86, 1400)                                                           # paper 2 Fig 14b
WALL_MEANS = [(0, 75, 660), (75, 125, 1040), (125, 275, 490)]                   # paper 2 Fig 14b section means


def figure(prof, xw, Tw, label, path):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True, constrained_layout=True)
    ax = axes[0]
    ax.fill_between(prof['x_mm'], prof['Tmin'], prof['Tmax'], color=SERIES[0], alpha=0.12, lw=0,
                    label='model: cross-section range (min-max)')
    ax.plot(prof['x_mm'], prof['bulkT'], color=SERIES[0], lw=2, label=f'model: bulk ({label})')
    ax.plot(*zip(*BULK), 'o', color=SERIES[1], ms=8, label='paper 2 (homogenised), Figs 12/13')
    ax.plot(*zip(*BULK_PORE), 's', color=SERIES[2], ms=8, label='paper 1 (pore scale)')
    ax.set_ylabel('Temperature [K]')
    ax.set_title('Bulk (mass-flow weighted) temperature', loc='left')
    ax = axes[1]
    ax.plot(xw, Tw, color=SERIES[0], lw=2, label='model: top wall, centreline (z = 50 mm, y = 25 mm)')
    ax.plot(*WALL_MAX, 'o', color=SERIES[1], ms=8, label='paper 2 Fig 14b: maximum')
    for i, (a, b, T) in enumerate(WALL_MEANS):
        ax.plot([a, b], [T, T], color=SERIES[1], lw=2, ls=(0, (5, 3)),
                label='paper 2 Fig 14b: section means' if i == 0 else None)
    ax.set_ylabel('Temperature [K]')
    ax.set_title('Wall temperature', loc='left')
    ax.set_xlabel('x from the inlet [mm]')
    ax.set_xlim(0, 280)                     # outlet markers at 275 mm fully visible
    for ax in axes:
        for x in (75, 125):
            ax.axvline(x, color=INK2, lw=0.8, ls=(0, (4, 3)))
        ax.grid(True, color=GRID, lw=0.6)
        ax.legend(frameon=False, fontsize=9, loc='upper right')
    fig.savefig(path, dpi=150)
    plt.close(fig)


def last_value(case, region, fo):
    d = sorted(glob.glob(f'{case}/postProcessing/{region}/{fo}/*'), key=lambda p: float(os.path.basename(p)))[-1]
    rows = [l.split() for l in open(glob.glob(f'{d}/*.dat')[0]) if l.strip() and not l.startswith('#')]
    return float(rows[-1][1])


def main():
    case, t, csv = sys.argv[1], sys.argv[2], sys.argv[3]
    post = f'{case}/post'
    prof = np.genfromtxt(csv, delimiter=',', names=True)
    to_local = frame(post)
    rows = []

    Tout = last_value(case, 'exitRegion', 'outletBulkT')
    rows.append(('Outlet bulk T [K]', 914, Tout, 'Table 4 (+-2 %); pore scale 1018.9'))

    # front / back face temperatures (patch values of the opaque faces)
    geo = patch_geometry(f'{post}/constant/porousRegion/polyMesh')
    Tb = read_field(f'{post}/{t}/porousRegion/T')[1]
    for name, patch, target, note in [('Front face max T [K]', 'porousRegion_to_inletRegion', 2250, 'text, Fig 8a ~2300'),
                                      ('Back face max T [K]', 'porousRegion_to_exitRegion', 1040, 'Fig 11b core')]:
        rows.append((name, target, float(np.max(Tb[patch])), note))

    # bulk profile
    for x, target in BULK:
        v = np.interp(x, prof['x_mm'], prof['bulkT'])
        rows.append((f'Bulk T at x = {x} mm [K]', target, v, 'Figs 12/13, read +-5-15 K'))

    # walls: z = 50 mm ("top") wall, maximum and section means
    xs, ys, Ts = [], [], []
    for r, w in [('inletRegion', 'inletWall'), ('porousRegion', 'porousWall'), ('exitRegion', 'exitWall')]:
        Cf, Sf = patch_geometry(f'{post}/constant/{r}/polyMesh')[w]
        T = np.broadcast_to(read_field(f'{post}/{t}/{r}/T')[1][w], len(Cf))
        L = to_local(Cf)
        top = L[:, 2] > 49.9
        xs.append(L[top, 0]); ys.append(L[top, 1]); Ts.append(T[top])
    xs, Ts = np.concatenate(xs), np.concatenate(Ts)
    if len(sys.argv) > 4:
        ys = np.concatenate(ys)
        c = np.abs(ys - 25.0) < np.abs(ys - 25.0).min() + 1e-3        # face column nearest y = 25 mm
        o = np.argsort(xs[c])
        figure(prof, xs[c][o], Ts[c][o], os.path.basename(os.path.abspath(case)), sys.argv[4])
    i = np.argmax(Ts)
    rows.append(('Top wall max T [K]', 1400, Ts[i], f'Fig 14b, at x ~ 85-87 mm (model: {xs[i]:.1f} mm)'))
    for name, a, b, target in [('inlet', 0, 75, 660), ('porous', 75, 125, 1040), ('exit', 125, 275, 490)]:
        m = (xs > a) & (xs < b)
        rows.append((f'Top wall mean T, {name} section [K]', target, Ts[m].mean(), 'Fig 14b pixel mean, loose'))

    # upstream wall heating: radiation absorbed by the inlet-section walls
    rows.append(('Inlet-wall radiation absorbed [W]', None, last_value(case, 'inletRegion', 'Qr_inletWall'),
                 ''))

    w = max(len(r[0]) for r in rows)
    print(f'{"quantity":{w}s}  {"paper 2":>8s}  {"model":>8s}  {"diff":>7s}   note')
    for name, target, model, note in rows:
        tgt = f'{target:8.0f}' if target is not None else f'{"-":>8s}'
        diff = f'{100 * (model - target) / target:+6.1f}%' if target else f'{"":>7s}'
        print(f'{name:{w}s}  {tgt}  {model:8.1f}  {diff}   {note}')


if __name__ == '__main__':
    main()
