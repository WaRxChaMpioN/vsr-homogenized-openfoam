#!/usr/bin/env python3
"""1-D check of the upstream side of the front face (air flows OUT of inletRegion through a
patch whose T is ~Dirichlet at the hot face temperature Tf).
Discretisation as in OpenFOAM: gaussConvectionScheme boundary flux = phi_b * h_b with h_b the
mixed-BC face value (gaussConvectionScheme.C:107-108, mixedFvPatchField.C:194-213), upwind
interior faces, laplacian with half-cell distance at the boundary.
Prints the last-cell temperatures and the enthalpy a 300 K limitTemperature clip would create."""
import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from meshcalc import multi
rho, cp, kap = 1.225, 1006.0, 0.026
def solve(sizes_mm, U, Tin=300.0, Tf=2250.0, vf=0.9925):
    dx = np.array(sizes_mm)*1e-3; n = len(dx)
    xc = np.cumsum(dx) - dx/2
    F = rho*cp*U
    A = np.zeros((n, n)); b = np.zeros(n)
    for i in range(n):
        # west face
        if i == 0:
            Dw = kap/(dx[0]/2); A[i,i] += Dw; b[i] += F*Tin + Dw*Tin
        else:
            Dw = kap/(xc[i]-xc[i-1]); A[i,i] += Dw; A[i,i-1] -= Dw + F   # inflow from i-1 (upwind)
        # east face
        if i == n-1:
            De = kap/(dx[-1]/2)
            # mixed: face value Tb = vf*Tf + (1-vf)*Tc ; outflow F*Tb ; diffusion De*(Tb - Tc)
            A[i,i] += F*(1-vf) + De*vf; b[i] += -F*vf*Tf + De*vf*Tf
        else:
            De = kap/(xc[i+1]-xc[i]); A[i,i] += De + F; A[i,i+1] -= De
    T = np.linalg.solve(A, b)
    return T, F
r = 0.4; n = 32; q = r**(1/(n-1)); first = 74.5*(1-q)/(1-q**n); oldmesh = [first*q**i for i in range(n)]
newmesh = multi(75.0, 48, [(70/75, 35/48, 1/3.5), (5/75, 13/48, 0.08)])
for Re, U0 in [(300, 0.1068), (900, 0.3204)]:
    U = 1.3*U0
    for name, m in [('old 32 cells, last 1.4 mm', oldmesh), ('adopted 48 cells, last 0.081 mm', newmesh)]:
        T, F = solve(m, U)
        clipE = F*np.clip(300.0 - T, 0, None).max()
        print(f'Re {Re} U_face={U:.3f} {name:34s}: last cells T = ' + ' '.join(f'{t:8.1f}' for t in T[-4:]) +
              f' | min T {T.min():8.1f} K')
