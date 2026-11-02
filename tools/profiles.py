#!/usr/bin/env python3
"""
Axial profiles along the receiver channel for one case (reconstructed, un-rotated, ascii,
with cell centres C and volumes V written; for the coarse cases use <case>/post from Allpost):

    profiles.py <caseDir> <time> <out.csv>

The channel axis is taken from the inlet-plane normal (so a residual tilt of the mesh does not
matter) and every bin is one axial cell layer of the mesh, bounded by its own point planes.

Per cell layer (x = distance from the inlet along the axis, mm):
  bulkT      mass-flow-weighted fluid temperature  sum(Ux T V)/sum(Ux V)   [K]
             (cell-centre velocity; the exact face-plane bulk temperature of Figs 12-13,
             from the face mass fluxes, is bulkFaces.py)
  meanT      volume-mean temperature                                     [K]
  Tmin, Tmax temperature range in the cross-section                      [K]
  Tcentre    temperature of the cells nearest the channel axis (y = z = 25 mm) [K]
  wallT_mean area-mean wall temperature (all four walls)                 [K]
  wallT_max  maximum wall temperature                                    [K]
  wall_qr    mean radiative flux absorbed by the walls (qr > 0: into the wall) [W/m2]
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, patch_geometry, read_points

ALL_REGIONS = ['inletRegion', 'interface', 'porousRegion', 'exitRegion']
WALLS = {'inletRegion': 'inletWall', 'interface': 'interfaceWall', 'porousRegion': 'porousWall', 'exitRegion': 'exitWall'}
regions, walls = ALL_REGIONS, WALLS


def case_regions(case):
    """Regions present in the case: four with the 0.5 mm slab (reference layout), three in the corrected model."""
    return [r for r in ALL_REGIONS if os.path.isdir(f'{case}/constant/{r}/polyMesh')]


def load_cells(case, t):
    X, T, Ux, V, R = [], [], [], [], []
    for r in regions:
        C = read_field(f'{case}/{t}/{r}/C')[0]
        X.append(C)
        T.append(read_field(f'{case}/{t}/{r}/T')[0])
        U = read_field(f'{case}/{t}/{r}/U')[0]
        Ux.append(U[:, 0] if U is not None and not np.isscalar(U) else np.zeros(len(C)))
        Vf = f'{case}/{t}/{r}/V'
        V.append(read_field(Vf)[0] if os.path.exists(Vf) else np.ones(len(C)))
        R.append(np.full(len(C), r))
    return np.vstack(X), np.concatenate(T), np.concatenate(Ux), np.concatenate(V), np.concatenate(R)


def load_walls(case, t):
    xs, Ts, qs, As = [], [], [], []
    for r, w in walls.items():
        Cf, Sf = patch_geometry(f'{case}/constant/{r}/polyMesh')[w]
        Tb = read_field(f'{case}/{t}/{r}/T')[1][w]
        qb = read_field(f'{case}/{t}/{r}/qr')[1].get(w, 0.0)
        n = len(Cf)
        xs.append(Cf)
        Ts.append(np.full(n, Tb) if np.isscalar(Tb) else Tb)
        qs.append(np.full(n, qb) if np.isscalar(qb) else qb)
        As.append(np.linalg.norm(Sf, axis=1))
    return np.vstack(xs), np.concatenate(Ts), np.concatenate(qs), np.concatenate(As)


def axis_frame(case):
    """Channel axis e (unit, from the inlet towards the outlet) and the inlet-plane centroid P0."""
    geo = patch_geometry(f'{case}/constant/inletRegion/polyMesh')
    Cf = np.vstack([geo[p][0] for p in ('inlet', 'opticalWindow')])
    Sf = np.vstack([geo[p][1] for p in ('inlet', 'opticalWindow')])
    A = np.linalg.norm(Sf, axis=1)
    e = -Sf.sum(0)
    e /= np.linalg.norm(e)
    P0 = (Cf * A[:, None]).sum(0) / A.sum()
    return e, P0


def layer_edges(case, e, P0):
    """Axial positions (m) of the mesh's own point planes, all regions."""
    s = []
    for r in regions:
        P = read_points(f'{case}/constant/{r}/polyMesh')
        s.append((P - P0) @ e)
    s = np.sort(np.concatenate(s))
    planes = [s[0]]
    for v in s[1:]:
        if v - planes[-1] > 2e-5:       # > 0.02 mm apart: a new plane
            planes.append(v)
    return np.array(planes)


def main():
    global regions, walls
    case, t, out = sys.argv[1], sys.argv[2], sys.argv[3]
    regions = case_regions(case)
    walls = {r: WALLS[r] for r in regions}
    e, P0 = axis_frame(case)
    C, T, Ux, V, _ = load_cells(case, t)
    D = C - P0
    s = D @ e
    rc = np.linalg.norm(D - np.outer(s, e), axis=1)
    U = Ux                                   # x-component; replaced by the axial component below
    Ufull = []
    for r in regions:
        Ur = read_field(f'{case}/{t}/{r}/U')[0]
        Ufull.append(Ur)
    Ua = np.vstack(Ufull) @ e
    Cw, Tw, qw, Aw = load_walls(case, t)
    sw = (Cw - P0) @ e
    edges = layer_edges(case, e, P0)
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (s > a) & (s < b)
        if not m.any():
            continue
        w = Ua[m] * V[m]
        bulk = np.sum(w * T[m]) / np.sum(w) if np.sum(w) > 0 else np.nan
        mean = np.average(T[m], weights=V[m])
        rmin = rc[m].min()
        centre = T[m][rc[m] <= rmin + 1e-6].mean()
        mw = (sw > a) & (sw < b)
        wT = np.average(Tw[mw], weights=Aw[mw]) if mw.any() else np.nan
        wTmax = Tw[mw].max() if mw.any() else np.nan
        wq = np.average(qw[mw], weights=Aw[mw]) if mw.any() else np.nan
        rows.append((0.5 * (a + b) * 1e3, bulk, mean, T[m].min(), T[m].max(), centre, wT, wTmax, wq))
    rows = np.array(rows)
    hdr = 'x_mm,bulkT,meanT,Tmin,Tmax,Tcentre,wallT_mean,wallT_max,wall_qr'
    np.savetxt(out, rows, delimiter=',', header=hdr, comments='', fmt='%.4f')
    print(f'{len(rows)} cell layers written to {out}  (axis tilt from x: {np.degrees(np.arccos(e[0])):.3f} deg)')
    print('  x[mm]   bulkT   meanT    Tmin    Tmax  Tcentre  wallT_mean wallT_max  wall_qr[W/m2]')
    for xs in [2.5, 12.5, 25, 37.5, 50, 62.5, 72, 75.5, 80, 90, 100, 110, 124, 130, 150, 175, 200, 225, 250, 272]:
        r = rows[np.argmin(abs(rows[:, 0] - xs))]
        print('  ' + '  '.join(f'{v:7.1f}' for v in r))


if __name__ == '__main__':
    main()
