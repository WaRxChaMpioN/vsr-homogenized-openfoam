#!/usr/bin/env python3
"""Axial/transverse cell sizes of the single-region blockMesh (replicates
blockMesh lineDivide.C multi-grading: secnDivs = label(frac*nDiv+0.5), expFact = r^(1/(n-1)))."""
import math
def section_sizes(L, n, r):
    if abs(r-1) < 1e-12: return [L/n]*n
    q = r**(1.0/(n-1))
    first = L*(1-q)/(1-q**n)
    return [first*q**i for i in range(n)]
def multi(L, nDiv, secs):
    out=[]
    for lf, nf, r in secs:
        n=int(nf*nDiv+0.5); out += section_sizes(L*lf, n, r)
    return out
if __name__ == '__main__':
    inlet = multi(75.0, 48, [(70/75, 35/48, 1/3.5), (5/75, 13/48, 0.08)])
    porous = section_sizes(50.0, 24, 8.0)
    exitz = section_sizes(150.0, 30, 4.0)
    for name, s in [('inlet', inlet), ('porous', porous), ('exit', exitz)]:
        print(f"{name:7s} n={len(s):3d} sum={sum(s):8.4f} mm first={s[0]:.4f} last={s[-1]:.4f} max ratio={max(max(a/b,b/a) for a,b in zip(s,s[1:])):.3f}")
    print('inlet last 6 cells (mm):', ' '.join(f'{x:.4f}' for x in inlet[-6:]))
    print('size jump at x=75: %.4f -> %.4f (x%.1f); at x=125: %.3f -> %.3f'%(inlet[-1],porous[0],porous[0]/inlet[-1],porous[-1],exitz[0]))
    wall = section_sizes(12.5, 12, 2.5)
    print('wall band: first %.4f last %.4f; window 1.25 mm'%(wall[0], wall[-1]))
    nax = len(inlet)+len(porous)+len(exitz); ntr = (12+20+12)**2
    print('cells: axial %d x %d = %d (inlet %d, porous %d, exit %d)'%(nax, ntr, nax*ntr, len(inlet)*ntr, len(porous)*ntr, len(exitz)*ntr))
    # Peclet of the last upstream cell, alpha = kappa/(rho cp)
    alpha = 0.026/(1.225*1006.0)
    for Re, U in [(300,0.1068),(600,0.2136),(900,0.3204)]:
        for f in (1.0, 1.3):
            print(f'Re {Re}: U_face {f:.1f} U_in -> P_cell = {f*U*inlet[-1]*1e-3/alpha:.2f}  (bounded if <= 2); old 1.4 mm cell: {f*U*1.4e-3/alpha:.1f}')
