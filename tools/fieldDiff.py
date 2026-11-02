#!/usr/bin/env python3
"""
Cell-by-cell and face-by-face difference between two corrected-model runs on the same flow split
(two angular sets, two G values, a warm start against a cold start, ...). Both must be un-rotated
copies made by Allpost (cell centres C and volumes V written):

    fieldDiff.py <caseA>/post:<tA> <caseB>/post:<tB>

The two runs may be rotated onto different beam rays: cells and faces are matched by index (same
split of the same flow mesh), and the match is checked in the channel frame.
Prints, per region, the cell T and G (max |B-A| and where, volume-weighted RMS of B-A), and per
wall and join patch the face T and the absorbed radiation qr (max |B-A| and where, area-weighted
RMS, patch totals). Positions are x, y, z in mm, x from the inlet.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, patch_geometry
from contours import frame

REGIONS = ['inletRegion', 'porousRegion', 'exitRegion']
WALLS = {'inletRegion': 'inletWall', 'porousRegion': 'porousWall', 'exitRegion': 'exitWall'}


def load(post, t):
    to_local = frame(post)
    out = {}
    for r in REGIONS:
        geo = patch_geometry(f'{post}/constant/{r}/polyMesh')
        T, Tb = read_field(f'{post}/{t}/{r}/T')
        G = read_field(f'{post}/{t}/{r}/G')[0]
        qr = read_field(f'{post}/{t}/{r}/qr')[1]
        patches = {}
        for p in [WALLS[r]] + [p for p in geo if p.startswith(r + '_to_')]:
            Cf, Sf = geo[p]
            n = len(Cf)
            patches[p] = dict(X=to_local(Cf), A=np.linalg.norm(Sf, axis=1),
                              T=np.broadcast_to(Tb[p], n), qr=np.broadcast_to(qr[p], n))
        out[r] = dict(X=to_local(read_field(f'{post}/{t}/{r}/C')[0]), V=read_field(f'{post}/{t}/{r}/V')[0],
                      T=T, G=np.broadcast_to(G, T.shape), patches=patches)
    return out


def where(X, i):
    return '({:5.1f} {:4.1f} {:4.1f})'.format(*X[i])


def row(name, a, b, w, X, unit, scale=1.0):
    d = (b - a) / scale
    i = np.argmax(np.abs(d))
    rms = np.sqrt((w * d**2).sum() / w.sum())
    print(f'  {name:34s} max|B-A| {abs(d[i]):10.4g} {unit:5s} at {where(X, i)}  (A {a[i] / scale:.5g}, B {b[i] / scale:.5g})'
          f'   rms {rms:9.3g}   range A {a.min() / scale:.5g}..{a.max() / scale:.5g}')


def main():
    (pa, ta), (pb, tb) = [s.rsplit(':', 1) for s in sys.argv[1:3]]
    A, B = load(pa, ta), load(pb, tb)
    print(f'A = {pa} @ {ta}\nB = {pb} @ {tb}')
    for r in REGIONS:
        a, b = A[r], B[r]
        dX = np.abs(a['X'] - b['X']).max()
        print(f'{r}: {len(a["T"])} cells, centres match to {dX:.2e} mm')
        if dX > 1e-3:
            sys.exit('cells do not match: not the same flow split')
        row('cell T', a['T'], b['T'], a['V'], a['X'], 'K')
        row('cell G', a['G'], b['G'], a['V'], a['X'], 'kW/m2', 1e3)
        for p, fa in a['patches'].items():
            fb = b['patches'][p]
            QA, QB = (fa['qr'] * fa['A']).sum(), (fb['qr'] * fb['A']).sum()
            print(f'  {p} ({len(fa["A"])} faces): absorbed qr A {QA:.2f} W, B {QB:.2f} W (B-A {QB - QA:+.2f} W)')
            row('face T', fa['T'], fb['T'], fa['A'], fa['X'], 'K')
            row('face qr', fa['qr'], fb['qr'], fa['A'], fa['X'], 'kW/m2', 1e3)


if __name__ == '__main__':
    main()
