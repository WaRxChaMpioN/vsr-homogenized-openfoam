#!/usr/bin/env python3
"""Probe points in the ROTATED frame (or unrotated with axis (1 0 0), flow stage).
transformPoints -rotate '(n1 n2)' applies rotationTensor(n1, n2) (OpenFOAM transform.H):
    R = s I + (1 - s) n3 n3/|n3|^2 + (n2 n1 - n1 n2),  s = n1.n2,  n3 = n1 x n2
Usage: rotPoints.py '(ex ey ez)' > system/probePoints     (axis from beamAlign.py --shell)
Writes a 'probeLocations' list: the centreline (y = z = 25 mm) at the x stations below."""
import sys, re
def rot(n1, n2):
    s = sum(a*b for a, b in zip(n1, n2))
    n3 = (n1[1]*n2[2]-n1[2]*n2[1], n1[2]*n2[0]-n1[0]*n2[2], n1[0]*n2[1]-n1[1]*n2[0])
    m3 = sum(a*a for a in n3)
    I = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
    if m3 < 1e-30:                       # n1 parallel to n2: identity (unrotated flow stage)
        return [[float(I[i][j]) for j in range(3)] for i in range(3)]
    return [[s*I[i][j] + (1-s)*n3[i]*n3[j]/m3 + (n2[i]*n1[j] - n1[i]*n2[j]) for j in range(3)] for i in range(3)]
def apply(R, p):
    return tuple(sum(R[i][j]*p[j] for j in range(3)) for i in range(3))
if __name__ == '__main__':
    e = tuple(float(x) for x in re.findall(r'[-+0-9.eE]+', sys.argv[1]))
    R = rot((1.0, 0.0, 0.0), e)
    xs_mm = [1, 10, 20, 30, 40, 50, 60, 70, 74, 74.9, 75.5, 80, 90, 100, 110, 120, 124.5,
             125.5, 130, 150, 175, 200, 225, 250, 274]
    print('// x [mm] along the channel axis at y = z = 25 mm, rotated with rotationTensor((1 0 0), e)')
    print('// R = ' + ' '.join(f'{R[i][j]:.10g}' for i in range(3) for j in range(3)))
    print('probeLocations\n(')
    for x in xs_mm:
        p = apply(R, (x*1e-3, 0.025, 0.025))
        print(f'    ({p[0]:.9g} {p[1]:.9g} {p[2]:.9g})    // x = {x} mm')
    print(');')
