#!/usr/bin/env python3
"""
Negative ray intensities in a reconstructed fvDOM time (all ILambda_<i>_0 fields, every region):

    negativeRays.py <case> <time>

A correct run has none. Negative intensities next to the porous faces come from boundary sink
rays (beamAlign.py): rays that greyDiffusiveRadiation emits from a face by their direction d
while fvDOM transports them out through it by dAve. They appear when nTheta is even.
Exit status 1 if any ray is negative anywhere.
"""
import glob, re, sys
import numpy as np

REGIONS = ['inletRegion', 'porousRegion', 'exitRegion']
HEADER = re.compile(rb'internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(')


def internal(path):
    buf = open(path, 'rb').read()
    m = HEADER.search(buf)
    if m is None:                                   # uniform field
        return np.array([float(re.search(rb'internalField\s+uniform\s+([-+0-9.eE]+)', buf).group(1))])
    v = np.array(buf[m.end():buf.index(b')', m.end())].split(), dtype=float)
    assert len(v) == int(m.group(1)), path
    return v


def main():
    case, t = sys.argv[1], sys.argv[2]
    bad = 0
    for r in REGIONS:
        files = sorted(glob.glob(f'{case}/{t}/{r}/ILambda_*_0'), key=lambda p: int(p.split('_')[-2]))
        if not files:
            sys.exit(f'{case}/{t}/{r}: no ILambda_*_0 fields (give the case itself and a fully reconstructed time, not post/)')
        neg, lo = [], (np.inf, None)
        for f in files:
            v = internal(f)
            ray = int(f.split('_')[-2])
            if v.min() < 0:
                neg.append((v.min(), ray, int((v < 0).sum())))
            if v.min() < lo[0]:
                lo = (v.min(), ray)
        bad += len(neg)
        print(f'{r:13s} {len(files)} rays, {len(neg)} with negative values; minimum {lo[0]:.4g} W/m2/sr (ray {lo[1]})')
        for mn, ray, n in sorted(neg)[:12]:
            print(f'    ray {ray:3d}: min {mn:12.4g} W/m2/sr in {n} cells')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
