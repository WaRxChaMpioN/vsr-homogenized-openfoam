#!/usr/bin/env python3
"""
Exact geometry rotation for a collimated beam in OpenFOAM fvDOM (3D).

The beam is given to ONE discrete ordinate by greyDiffusiveRadiation (qRadExt, qRadExtDir).
fvDOM transports every ray with the face flux Ji = dAve & Sf (radiativeIntensityRay.C),
so the beam only travels without lateral (false) scattering if the channel axis is
parallel to that ray's dAve -- not to its nominal direction d.

This script rebuilds the ray set exactly as fvDOM does, picks the ray closest to the
channel axis (+x before rotation), and prints:
  - the transformPoints rotation that maps +x onto dAve/|dAve| (and its exact inverse),
  - qRadExtDir that selects that ray unambiguously (the 4 nearest rays tie for (1 0 0)),
  - the inlet velocity vector in the rotated frame,
  - the boundary sink rays of the rotated channel (see below).

Boundary sink rays: greyDiffusiveRadiation decides from d whether a ray leaves a boundary
face (-n.d > 0: fixed emitted value), fvDOM moves it with dAve. A ray with -n.d > 0 but
n.dAve > 0 carries that emitted value OUT through the face, a sink in every adjacent cell;
along a boundary the deficit accumulates into large negative intensities (thin cells make it
worse). With nTheta EVEN the beam ray is off the equator, the channel is tilted by dTheta/2 and
the faces normal to the axis get such rays. With nTheta ODD the beam ray is on the equator,
the rotation is about z only, every boundary normal is horizontal or vertical, and
sign(n.d) == sign(n.dAve) for every ray: no sinks. Use an odd nTheta.

Usage:  beamAlign.py <nPhi> <nTheta>          or   beamAlign.py <file with nPhi/nTheta, e.g. system/caseParameters>
        add --shell to print only shell variable assignments (used by Allrun)
"""
import math, re, sys


def rays(nPhi, nTheta):
    """Ray set of fvDOM::initialise (3D) and radiativeIntensityRay (d, dAve)."""
    dPhi = math.pi / (2 * nPhi)
    dTheta = math.pi / nTheta
    out = []
    for n in range(1, nTheta + 1):
        for m in range(1, 4 * nPhi + 1):
            th = (2 * n - 1) * dTheta / 2.0
            ph = (2 * m - 1) * dPhi / 2.0
            st, ct, sp, cp = math.sin(th), math.cos(th), math.sin(ph), math.cos(ph)
            d = (st * sp, st * cp, ct)
            f = math.sin(0.5 * dPhi) * (dTheta - math.cos(2.0 * th) * math.sin(dTheta))
            dAve = (sp * f, cp * f, 0.5 * dPhi * math.sin(2.0 * th) * math.sin(dTheta))
            out.append(dict(id=len(out), theta=th, phi=ph, d=d, dAve=dAve))
    return out


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def unit(a):
    m = math.sqrt(dot(a, a))
    return tuple(x / m for x in a)


def select(rs, q):
    """Same loop as greyDiffusiveRadiationMixedFvPatchScalarField (strict '>')."""
    best, maxRay = -1, -1e300
    for r in rs:
        v = dot(q, r['d'])
        if v > maxRay:
            maxRay, best = v, r['id']
    return best


def rotation(a, b):
    """Rows of rotationTensor(a, b) (unit a, b), the rotation of transformPoints -rotate '(a b)'."""
    s, n3 = dot(a, b), (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    m = dot(n3, n3)
    return [[s * (i == j) + (1 - s) * n3[i] * n3[j] / m + b[i] * a[j] - a[i] * b[j] for j in range(3)]
            for i in range(3)]


def sinks(rs, e, tol=1e-12):
    """Sink rays (-n.d > 0, n.dAve > tol) on the six boundary normals of the rotated channel."""
    R = rotation((1.0, 0.0, 0.0), e)
    out = {}
    for name, v in [('+x', (1, 0, 0)), ('-x', (-1, 0, 0)), ('+y', (0, 1, 0)), ('-y', (0, -1, 0)),
                    ('+z', (0, 0, 1)), ('-z', (0, 0, -1))]:
        n = tuple(dot(row, v) for row in R)
        out[name] = [r['id'] for r in rs if -dot(n, r['d']) > 0 and dot(n, r['dAve']) > tol]
    return out


def read_nphi_ntheta(path):
    txt = re.sub(r'//.*', '', open(path).read())
    return int(re.search(r'\bnPhi\s+(\d+)', txt).group(1)), int(re.search(r'\bnTheta\s+(\d+)', txt).group(1))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    shell = '--shell' in sys.argv
    if len(args) == 1:
        nPhi, nTheta = read_nphi_ntheta(args[0])
    else:
        nPhi, nTheta = int(args[0]), int(args[1])

    rs = rays(nPhi, nTheta)
    beam = rs[select(rs, (1.0, 0.0, 0.0))]
    e = unit(beam['dAve'])                       # rotated channel axis
    assert select(rs, e) == beam['id']           # qRadExtDir = e picks the same ray, no tie
    second = sorted((dot(e, r['d']) for r in rs), reverse=True)[1]
    misalign_d = math.degrees(math.acos(min(1.0, dot(e, beam['d']))))
    fmt = lambda v: '(%.15g %.15g %.15g)' % v
    sink = sinks(rs, e)
    nsink = sum(len(v) for v in sink.values())

    if shell:
        print(f"BEAM_RAY={beam['id']}")
        print(f"BEAM_AXIS='{fmt(e)}'")
        print(f"ROTATE_FWD='((1 0 0) {fmt(e)})'")
        print(f"ROTATE_INV='({fmt(e)} (1 0 0))'")
        print(f"SINK_RAYS={nsink}")
        return

    print(f"fvDOM nPhi={nPhi} nTheta={nTheta}: {len(rs)} rays, dPhi={math.degrees(math.pi/(2*nPhi)):.4f} deg, dTheta={180/nTheta:.4f} deg")
    print(f"beam ray            : {beam['id']}  (ILambda_{beam['id']}_0)  theta={math.degrees(beam['theta']):.4f}  phi={math.degrees(beam['phi']):.4f}")
    print(f"ray direction d     : {fmt(beam['d'])}")
    print(f"transport dir dAve^ : {fmt(e)}   (d vs dAve: {misalign_d:.4f} deg)")
    print(f"axis tilt from +x   : {math.degrees(math.acos(e[0])):.4f} deg"
          f"  (azimuth {math.degrees(math.atan2(e[1], e[0])):.4f}, elevation {math.degrees(math.asin(e[2])):.4f})")
    print(f"selection margin    : next ray has e.d = {second:.6f} (beam ray {dot(e, beam['d']):.9f})")
    print(f"boundary sink rays  : {nsink}" + ''.join(f"  {k}: {v}" for k, v in sink.items() if v)
          + ('' if nsink else '  (none)'))
    if nsink:
        print(f"  WARNING: nTheta {nTheta} is even - use an odd nTheta (e.g. {nTheta + 1}) so that the beam ray"
              " is on the equator and no ray is a boundary sink")
    print()
    print("transformPoints (forward, before the run):")
    print(f"  transformPoints -allRegions -rotate '((1 0 0) {fmt(e)})'")
    print("transformPoints (exact inverse, post-processing copy only):")
    print(f"  transformPoints -allRegions -rotateFields -rotate '({fmt(e)} (1 0 0))'")
    print()
    print(f"opticalWindow ray BC: qRadExtDir {fmt(e)};")
    print(f"inlet velocity      : surfaceNormalFixedValue refValue uniform -U  (or U*{fmt(e)})")


if __name__ == '__main__':
    main()
