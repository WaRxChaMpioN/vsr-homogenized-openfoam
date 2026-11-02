#!/usr/bin/env python3
"""
Bulk (mass-flow-weighted) temperature on every plane of cell faces normal to the channel axis, from
the face mass fluxes, for one case (un-rotated copy written by Allpost):

    bulkFaces.py <case>/post <time> <out.csv>

bulk(x) = sum(mdot_f T_f) / sum(mdot_f) over the faces of the plane at x, with mdot_f the solver's
face mass flux (phi; exact, also where the reconstructed cell-centre velocity wiggles next to the
porous faces) and T_f the face temperature: linear interpolation between the two cells inside a
region, and the boundary value on the inlet plane, on the two faces of the porous zone (both sides,
same value) and at the outlet (zero gradient: the cell value). At the outlet this is the bulk
temperature of the energy books. profiles.py instead weights the cells with their cell-centre
velocity and puts each value at the cell centre; with the 9.2 mm last cell and the zero-gradient
outlet its last value (= the outlet value) sits 4.6 mm upstream of the outlet.
Output: x_mm (from the inlet), bulkT, mdot_rel (mass flow of the plane / median; 2 on the porous
faces, which both regions book).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, read_points, read_faces, read_boundary, _list_after
from contours import frame, cluster

REGIONS = ['inletRegion', 'porousRegion', 'exitRegion']


def labels(mesh, name):
    buf = open(os.path.join(mesh, name), 'rb').read()
    return _list_after(buf, buf.index(b'}') + 1)[0].astype(int)


def region_planes(post, t, r, to_local):
    mesh = f'{post}/constant/{r}/polyMesh'
    P = read_points(mesh)
    F = np.array(read_faces(mesh))                                   # block mesh: quads
    Cf = to_local(P[F].mean(1))
    Sf = 0.5 * np.cross(P[F[:, 2]] - P[F[:, 0]], P[F[:, 3]] - P[F[:, 1]])
    own, nei = labels(mesh, 'owner'), labels(mesh, 'neighbour')
    nint = len(nei)
    C = to_local(read_field(f'{post}/{t}/{r}/C')[0])
    T, Tb = read_field(f'{post}/{t}/{r}/T')
    phi_i, phi_b = read_field(f'{post}/{t}/{r}/phi')
    nx = np.abs(Sf[:, 0]) / np.linalg.norm(Sf, axis=1)               # un-rotated copy: channel axis = x
    xs, mT, m = [], [], []
    # internal faces normal to the axis
    sel = np.where(nx[:nint] > 0.99)[0]
    xo, xn, xf = C[own[sel], 0], C[nei[sel], 0], Cf[sel, 0]
    wO = np.abs(xn - xf) / np.abs(xn - xo)
    Tf = wO * T[own[sel]] + (1 - wO) * T[nei[sel]]
    mdot = np.abs(phi_i[sel])
    xs.append(xf); mT.append(mdot * Tf); m.append(mdot)
    # boundary planes: inlet, window, porous faces, outlet
    for p, (typ, n, start) in read_boundary(mesh).items():
        if p not in phi_b:
            continue
        f = np.arange(start, start + n)
        if not np.all(nx[f] > 0.99):
            continue                                                  # walls
        ph = np.abs(np.broadcast_to(phi_b[p], n))
        Tp = np.broadcast_to(Tb[p], n) if p in Tb else T[own[f]]      # zeroGradient outlet: cell value
        xs.append(Cf[f, 0]); mT.append(ph * Tp); m.append(ph)
    return np.concatenate(xs), np.concatenate(mT), np.concatenate(m)


def main():
    post, t, out = sys.argv[1], sys.argv[2], sys.argv[3]
    to_local = frame(post)
    X, MT, M = [], [], []
    for r in REGIONS:
        x, mt, m = region_planes(post, t, r, to_local)
        X.append(x); MT.append(mt); M.append(m)
    X, MT, M = np.concatenate(X), np.concatenate(MT), np.concatenate(M)
    ix, ux = cluster(X, 1e-3)
    bulk = np.bincount(ix, weights=MT) / np.bincount(ix, weights=M)
    flow = np.bincount(ix, weights=M)
    o = np.argsort(ux)
    with open(out, 'w') as f:
        f.write('x_mm,bulkT,mdot_rel\n')
        for k in o:
            f.write(f'{ux[k]:.4f},{bulk[k]:.4f},{flow[k] / np.median(flow):.6f}\n')
    print(f'{len(o)} face planes written to {out}; mass flow per plane / median: {flow.min() / np.median(flow):.6f} .. {flow.max() / np.median(flow):.6f}')


if __name__ == '__main__':
    main()
