#!/usr/bin/env python3
"""Check the arithmetic of radiativeInterfaceT / radiativeWallTc in Python:
c weights on the rotated 6x10 ray set (as fvDOM builds it), Newton convergence, tangent = root."""
import math, numpy as np
nPhi, nTheta = 6, 10
dPhi = math.pi/(2*nPhi); dTh = math.pi/nTheta
d, dA = [], []
for n in range(1, nTheta+1):
    for m in range(1, 4*nPhi+1):
        th = (2*n-1)*dTh/2; ph = (2*m-1)*dPhi/2
        f = math.sin(0.5*dPhi)*(dTh - math.cos(2*th)*math.sin(dTh))
        d.append([math.sin(th)*math.sin(ph), math.sin(th)*math.cos(ph), math.cos(th)])
        dA.append([math.sin(ph)*f, math.cos(ph)*f, 0.5*dPhi*math.sin(2*th)*math.sin(dTh)])
d = np.array(d); dA = np.array(dA)
R = np.array([0.9793728467, 0.1289368813, -0.1555766943, -0.1289368813, 0.9916010168,
              0.01013430785, 0.1555766943, 0.01013430785, 0.9877718299]).reshape(3, 3)
def weights(n):
    nd = d @ n; ndA = dA @ n
    return -ndA[nd < 0].sum()/math.pi, ndA[nd > 0].sum()/math.pi
for name, n in [('front face, upstream side (n=+e)', R @ [1, 0, 0]), ('front face, porous side (n=-e)', -(R @ [1, 0, 0])),
                ('wall y=0 (n=-y)', -(R @ [0, 1, 0])), ('wall y=50', R @ [0, 1, 0]), ('wall z=0', -(R @ [0, 0, 1])), ('wall z=50', R @ [0, 0, 1])]:
    c1, c2 = weights(n)
    print(f'{name:34s} c_own = {c1:.5f}  c_other_side = {c2:.5f}')
sigma, eps = 5.670374419e-8, 0.92
def interface(k1, k2, Tc, TcN, qin, qinN, c1, c2):
    K = k1+k2; Q = (1-c1*(1-eps))*max(qin, 0) + (1-c2*(1-eps))*max(qinN, 0); E = (c1+c2)*eps*sigma
    C = k1*Tc + k2*TcN + Q
    T = min(C/K, (C/E)**0.25); hist = [T]
    for _ in range(50):
        dT = (C - K*T - E*T**4)/(-K - 4*E*T**3); T -= dT; hist.append(T)
        if abs(dT) < 1e-9*T: break
    a = k2 + 4*E*T**3; b = k2*TcN + Q + 3*E*T**4
    Tb = (b + k1*Tc)/(a + k1)             # mixed face value with the current cell value
    return T, Tb, a/(a+k1), C - K*T - E*T**4, hist
c = weights(R @ [1, 0, 0])[0]
for case in [('cold start (porous 300 K)', 650, 81700, 300, 300, 5e6, 0),
             ('G5e6 hot, upstream side', 650, 81700, 400, 2200, 5.2e6, 1.4e6),
             ('G5e6 hot, porous side', 81700, 650, 2200, 400, 1.4e6, 5.2e6),
             ('G1e7, porous side', 81700, 650, 2900, 450, 3.5e6, 1.05e7),
             ('back face, exit side', 22.6, 10000, 900, 1050, 3e4, 9e4)]:
    T, Tb, vf, res, hist = interface(*case[1:], c, c)
    mono = all(x >= y for x, y in zip(hist, hist[1:]))
    print(f'{case[0]:28s}: Tf = {T:8.2f} K, tangent face value = {Tb:8.2f} K, valueFraction = {vf:.4f}, '
          f'balance residual = {res:.2e} W/m2, Newton its = {len(hist)-1}, monotone = {mono}')
