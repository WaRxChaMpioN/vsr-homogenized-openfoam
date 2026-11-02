#!/usr/bin/env python3
"""
Bulk (mass-flow-weighted) temperature on every plane of cell faces normal to the channel axis, from
the face mass fluxes, for one case (un-rotated copy written by Allpost):

    bulkFaces.py <case>/post <time> <out.csv>

bulk(x) = sum(phi_f T_f) / sum(phi_f) over the faces of the plane at x, with phi_f the solver's face
mass flux and T_f the face value of the solver's convection term (div(phi,h): bounded Gauss
linearUpwind with the gradient cellLimited Gauss linear 1; with constant c_p the same face value
holds for T): the upwind cell value extrapolated to the face with its limited gradient inside a
region, and the boundary value on the inlet plane, on the porous faces (both sides, same value) and
at the outlet (zero gradient: the cell value). Every plane is therefore the convective enthalpy flux
of the solver divided by mdot c_p, and the outlet value is the outlet bulk temperature of the
energy books. profiles.py instead weights the cells with their cell-centre velocity and puts each
value at the cell centre; with the 9.2 mm last cell and the zero-gradient outlet its last value
(= the outlet value) sits 4.6 mm upstream of the outlet.
Needs the cell centres C and volumes V of <time> (Allpost writes them).
Output: x_mm (from the inlet), bulkT, mdot_rel (mass flow of the plane / median; 2 on the porous
faces, which both regions book).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, read_points, read_faces, read_boundary, _list_after, frame, cluster

REGIONS = ['inletRegion', 'porousRegion', 'exitRegion']


def labels(mesh, name):
    buf = open(os.path.join(mesh, name), 'rb').read()
    return _list_after(buf, buf.index(b'}') + 1)[0].astype(int)


def face_values(T, Tb, C, V, Cf, Sf, own, nei, phi):
    """Face values of T on the internal faces as the solver's 'bounded Gauss linearUpwind limitedGrad'
    builds them, limitedGrad = cellLimited Gauss linear 1. Tb: values on the boundary faces."""
    n = len(nei)
    o = own[:n]
    nf = Sf[:n] / np.linalg.norm(Sf[:n], axis=1)[:, None]
    w = np.einsum('ij,ij->i', nf, C[nei] - Cf[:n]) / np.einsum('ij,ij->i', nf, C[nei] - C[o])
    g = np.zeros((len(T), 3))                                         # Gauss linear gradient
    Tl = w * T[o] + (1 - w) * T[nei]
    np.add.at(g, o, Sf[:n] * Tl[:, None])
    np.add.at(g, nei, -Sf[:n] * Tl[:, None])
    np.add.at(g, own[n:], Sf[n:] * Tb[:, None])
    g /= V[:, None]
    hi, lo = T.copy(), T.copy()                                       # extremes over neighbours and boundary values
    for cells, vals in ((o, T[nei]), (nei, T[o]), (own[n:], Tb)):
        np.maximum.at(hi, cells, vals)
        np.minimum.at(lo, cells, vals)
    hi -= T
    lo -= T
    lim = np.ones(len(T))                                             # cellLimited, k = 1
    for cells, faces in ((o, slice(0, n)), (nei, slice(0, n)), (own[n:], slice(n, None))):
        d = np.einsum('ij,ij->i', Cf[faces] - C[cells], g[cells])
        r = np.ones_like(d)
        up, dn = d > hi[cells], d < lo[cells]
        r[up] = hi[cells][up] / d[up]
        r[dn] = lo[cells][dn] / d[dn]
        np.minimum.at(lim, cells, r)
    g *= lim[:, None]
    u = np.where(phi >= 0, o, nei)                                    # upwind cell
    return T[u] + np.einsum('ij,ij->i', Cf[:n] - C[u], g[u])


def region_planes(post, t, r, to_local):
    mesh = f'{post}/constant/{r}/polyMesh'
    P = read_points(mesh)
    F = np.array(read_faces(mesh))                                   # block mesh: planar quads
    Cf = P[F].mean(1)
    Sf = 0.5 * np.cross(P[F[:, 2]] - P[F[:, 0]], P[F[:, 3]] - P[F[:, 1]])
    own, nei = labels(mesh, 'owner'), labels(mesh, 'neighbour')
    nint = len(nei)
    C = read_field(f'{post}/{t}/{r}/C')[0]
    V = read_field(f'{post}/{t}/{r}/V')[0]
    T, Tb = read_field(f'{post}/{t}/{r}/T')
    phi_i, phi_b = read_field(f'{post}/{t}/{r}/phi')
    Tbf, phib = np.empty(len(F) - nint), np.zeros(len(F) - nint)
    for p, (typ, n, start) in read_boundary(mesh).items():
        f = slice(start - nint, start - nint + n)
        Tbf[f] = np.broadcast_to(Tb[p], n) if p in Tb else T[own[start:start + n]]   # zero gradient: the cell value
        phib[f] = np.broadcast_to(phi_b.get(p, 0.0), n)
    Tf = np.concatenate([face_values(T, Tbf, C, V, Cf, Sf, own, nei, phi_i), Tbf])
    phi = np.abs(np.concatenate([phi_i, phib]))
    axial = np.abs(Sf[:, 0]) > 0.99 * np.linalg.norm(Sf, axis=1)    # un-rotated copy: channel axis = x
    return to_local(Cf[axial])[:, 0], phi[axial] * Tf[axial], phi[axial]


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
