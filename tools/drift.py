#!/usr/bin/env python3
"""
Change of the temperature between two reconstructed times of a heat-transfer case (run frame;
T is a scalar, so the rotation does not matter):

    drift.py <case> <t0> <t1>

Prints max |T(t1) - T(t0)| over the cells and the wall faces of each region. A converged run
changes by well under 1 K over its last write interval.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field

WALLS = {'inletRegion': 'inletWall', 'porousRegion': 'porousWall', 'exitRegion': 'exitWall'}


def main():
    case, t0, t1 = sys.argv[1:4]
    worst = 0.0
    for r, w in WALLS.items():
        a, ab = read_field(f'{case}/{t0}/{r}/T'); b, bb = read_field(f'{case}/{t1}/{r}/T')
        dc = float(np.abs(b - a).max()); dw = float(np.abs(np.asarray(bb[w]) - np.asarray(ab[w])).max())
        worst = max(worst, dc, dw)
        print(f'{r:13s} max |dT| cells {dc:8.3f} K, wall faces {dw:8.3f} K')
    print(f'drift {t0} -> {t1}: max {worst:.3f} K')


if __name__ == '__main__':
    main()
