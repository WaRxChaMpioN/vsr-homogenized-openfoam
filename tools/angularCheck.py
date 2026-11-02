#!/usr/bin/env python3
"""
Angular-discretisation check for the beam-aligned (rotated) channel.

For each (nPhi, nTheta): rotate the channel exactly as Allrun does
(transformPoints -rotate '((1 0 0) e)', e = dAve^ of the beam ray), then for the six
wall/plane orientations of the channel compute the hemispherical weight that fvDOM uses
for diffuse emission and reflection,
    S(n) = sum over rays with d.n > 0 of dAve.n          (exact value: pi)
S/pi - 1 is the relative error ("overhang") of the emitted/reflected/absorbed wall flux.

Usage: angularCheck.py nPhi:nTheta [nPhi:nTheta ...]
"""
import math, sys
import numpy as np
from beamAlign import rays, select, unit


def rotation(a, b):
    """Minimal rotation taking unit vector a to unit vector b (as transformPoints -rotate)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    v, c = np.cross(a, b), float(a @ b)
    if np.linalg.norm(v) < 1e-15:
        return np.eye(3)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K / (1 + c)


print(f"{'nPhi x nTheta':>13} {'rays':>5} {'beam ray':>9} {'tilt':>7}   overhang S/pi-1 for faces with outward normal (rotated)")
print(f"{'':>13} {'':>5} {'':>9} {'[deg]':>7}   {'+x(outlet)':>11} {'-x(inlet)':>10} {'+-y walls':>10} {'+-z walls':>10}")
for arg in sys.argv[1:]:
    nPhi, nTheta = map(int, arg.split(':'))
    rs = rays(nPhi, nTheta)
    beam = rs[select(rs, (1.0, 0.0, 0.0))]
    e = unit(beam['dAve'])
    R = rotation((1, 0, 0), e)
    D = np.array([r['d'] for r in rs]); A = np.array([r['dAve'] for r in rs])
    err = {}
    for name, n0 in [('+x', (1, 0, 0)), ('-x', (-1, 0, 0)), ('+y', (0, 1, 0)), ('-y', (0, -1, 0)), ('+z', (0, 0, 1)), ('-z', (0, 0, -1))]:
        n = R @ np.array(n0, float)
        S = float(np.sum((A @ n)[D @ n > 0]))
        err[name] = S / math.pi - 1
    ey = max(abs(err['+y']), abs(err['-y'])); ez = max(abs(err['+z']), abs(err['-z']))
    print(f"{nPhi:>6} x {nTheta:<4} {len(rs):>5} {beam['id']:>9} {math.degrees(math.acos(e[0])):7.3f}   "
          f"{100*err['+x']:+10.3f}% {100*err['-x']:+9.3f}% {100*ey:9.3f}% {100*ez:9.3f}%")
