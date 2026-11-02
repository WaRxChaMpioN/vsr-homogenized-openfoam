#!/usr/bin/env python3
"""Rebuild phi on the inter-region patches after 'splitMeshRegions -cellZonesOnly -overwrite'.

splitMeshRegions (v2206) sets every field, phi included, to 0 on the patches it adds
(splitMeshRegions.C:244-253 'Hack: set value to 0 for introduced patches'). The signed base
flux is recovered from constant/<R>/polyMesh/faceRegionAddressing (splitMeshRegions.C:925-961):
a = +(oldFace+1) if the region face keeps the base orientation, -(oldFace+1) if reversed, so
    phi_R[k] = sign(a) * phi_base[|a| - 1]      (interface faces are base INTERNAL faces).

Usage:  restorePhi.py <case> <time> [--dry-run] [--out <dir>] [--tol 1e-8] [--peclet <dx_m>]
        --peclet 8.12e-5 also checks max(u_n)*dx/alpha_air <= 2 on inletRegion_to_porousRegion
  <case>/<time>/phi                     base-mesh phi (the single-region flow solution)
  <case>/constant/<R>/polyMesh/{boundary,faceRegionAddressing}
  <case>/<time>/<R>/phi                 rewritten in place (or into --out/<R>/phi)
Regions are found from constant/*/polyMesh/faceRegionAddressing; interface patches are the
patches named <R>_to_<S>. Then it checks, and exits non-zero if a check fails:
  1. |a|-1 < nInternalFaces(base) for every interface face (they must be base internal faces)
  2. <R>_to_<S> and <S>_to_<R> list the same base faces in the same order, phi equal and opposite
  3. net boundary mass flow of every region ~ 0 (relative to the inlet mass flow)
"""
import os, re, sys, glob

def strip_comments(t):
    t = re.sub(r'/\*.*?\*/', '', t, flags=re.S)
    return re.sub(r'//[^\n]*', '', t)

def parse_list(txt, pos):
    """Parse an OpenFOAM ascii list starting at txt[pos:] ('N ( ... )' or 'N{v}'); return (list of str, end)."""
    m = re.compile(r'\s*(\d+)\s*([({])').match(txt, pos)
    if not m:
        raise ValueError('no list at %d: %r' % (pos, txt[pos:pos+40]))
    n = int(m.group(1))
    if m.group(2) == '{':
        e = txt.index('}', m.end())
        return [txt[m.end():e].strip()]*n, e+1
    e = txt.index(')', m.end())
    items = txt[m.end():e].split()
    if len(items) != n:
        raise ValueError('list size %d != %d' % (len(items), n))
    return items, e+1

def read_labelIOList(path):
    txt = strip_comments(open(path).read())
    body = txt[txt.index('}', txt.index('FoamFile'))+1:]
    items, _ = parse_list(body, 0)
    return [int(x) for x in items]

def read_boundary(path):
    """Return {name: (startFace, nFaces, type)} from a polyMesh/boundary file."""
    txt = strip_comments(open(path).read())
    body = txt[txt.index('}', txt.index('FoamFile'))+1:]
    out = {}
    for m in re.finditer(r'(\w+)\s*\{([^{}]*)\}', body):
        d = m.group(2)
        out[m.group(1)] = (int(re.search(r'startFace\s+(\d+)', d).group(1)),
                           int(re.search(r'nFaces\s+(\d+)', d).group(1)),
                           re.search(r'\btype\s+(\w+)', d).group(1))
    return out

def block_span(txt, name, start):
    """Find 'name { ... }' at brace depth 1 inside txt[start:]; return (open, close) indices of the braces."""
    for m in re.finditer(r'(^|\s)(' + re.escape(name) + r')\s*\{', txt[start:]):
        o = start + m.end() - 1
        depth = 0
        for i in range(o, len(txt)):
            if txt[i] == '{':
                depth += 1
            elif txt[i] == '}':
                depth -= 1
                if depth == 0:
                    return o, i
    raise KeyError(name)

def read_internal(path):
    txt = open(path).read()
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*', txt)
    items, _ = parse_list(txt, m.end())
    return [float(x) for x in items]

def patch_values(txt, name, size):
    o, c = block_span(txt, name, txt.index('boundaryField'))
    blk = txt[o+1:c]
    m = re.search(r'\bvalue\s+uniform\s+([^;]+);', blk)
    if m:
        return [float(m.group(1))]*size
    m = re.search(r'\bvalue\s+nonuniform\s+List<scalar>\s*', blk)
    if not m:
        return [0.0]*size                     # e.g. empty patch
    items, _ = parse_list(blk, m.end())
    return [float(x) for x in items]

def set_patch_value(txt, name, vals):
    o, c = block_span(txt, name, txt.index('boundaryField'))
    lst = 'nonuniform List<scalar> \n%d\n(\n%s\n)\n' % (len(vals), '\n'.join('%.14g' % v for v in vals))
    blk = '\n    {\n        type            calculated;\n        value           ' + lst + ';\n    }'
    return txt[:o] + blk.lstrip() + txt[c+1:]

def read_points(path):
    txt = strip_comments(open(path).read())
    body = txt[txt.index('}', txt.index('FoamFile'))+1:]
    m = re.match(r'\s*(\d+)\s*\(', body)
    pts = re.findall(r'\(([^()]+)\)', body[m.end():])
    return [tuple(map(float, q.split())) for q in pts[:int(m.group(1))]]

def read_faces(path, start, n):
    txt = strip_comments(open(path).read())
    body = txt[txt.index('}', txt.index('FoamFile'))+1:]
    m = re.match(r'\s*(\d+)\s*\(', body)
    fs = re.findall(r'\d+\(([^()]+)\)', body[m.end():])
    return [list(map(int, f.split())) for f in fs[start:start+n]]

def face_area(pts, f):
    p0 = pts[f[0]]; ax = ay = az = 0.0
    for i in range(1, len(f)-1):
        a = [pts[f[i]][k]-p0[k] for k in range(3)]; b = [pts[f[i+1]][k]-p0[k] for k in range(3)]
        ax += a[1]*b[2]-a[2]*b[1]; ay += a[2]*b[0]-a[0]*b[2]; az += a[0]*b[1]-a[1]*b[0]
    return 0.5*(ax*ax+ay*ay+az*az)**0.5

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    case, time = args[0], args[1]
    dry = '--dry-run' in sys.argv
    tol = float(sys.argv[sys.argv.index('--tol')+1]) if '--tol' in sys.argv else 1e-8
    outdir = sys.argv[sys.argv.index('--out')+1] if '--out' in sys.argv else None
    phi_base = read_internal(os.path.join(case, time, 'phi'))
    nInt = len(phi_base)
    regions = sorted(os.path.basename(os.path.dirname(os.path.dirname(p)))
                     for p in glob.glob(os.path.join(case, 'constant', '*', 'polyMesh', 'faceRegionAddressing')))
    ok = True
    restored = {}
    nets = {}
    inflow = None
    for R in regions:
        bnd = read_boundary(os.path.join(case, 'constant', R, 'polyMesh', 'boundary'))
        addr = read_labelIOList(os.path.join(case, 'constant', R, 'polyMesh', 'faceRegionAddressing'))
        phifile = os.path.join(case, time, R, 'phi')
        txt = open(phifile).read()
        for name, (start, n, typ) in bnd.items():
            if not re.fullmatch(r'%s_to_\w+' % re.escape(R), name):
                continue
            a = addr[start:start+n]
            if any(abs(x)-1 >= nInt for x in a):
                print('FAIL %s/%s: a face is not a base internal face' % (R, name)); ok = False
            vals = [(1.0 if x > 0 else -1.0)*phi_base[abs(x)-1] for x in a]
            restored[(R, name)] = (a, vals)
            txt = set_patch_value(txt, name, vals)
            print('%-14s %-28s %6d faces  sum(phi) = %+.10e kg/s' % (R, name, n, sum(vals)))
            if '--peclet' in sys.argv and R == 'inletRegion':
                # cell Peclet of the upstream face cells: must be <= 2 for a bounded outflow
                # through the near-Dirichlet face (see outflow_face_1d.py)
                dx = float(sys.argv[sys.argv.index('--peclet')+1])      # face-cell thickness [m]
                pts = read_points(os.path.join(case, 'constant', R, 'polyMesh', 'points'))
                fcs = read_faces(os.path.join(case, 'constant', R, 'polyMesh', 'faces'), start, n)
                un = max(v/(1.225*face_area(pts, f)) for v, f in zip(vals, fcs))
                Pe = un*dx/(0.026/(1.225*1006.0))
                print('    max face velocity %.4f m/s -> cell Peclet %.2f %s' % (un, Pe, 'OK' if Pe <= 2 else 'FAIL (> 2)'))
                ok = ok and Pe <= 2
        net = 0.0
        for name, (start, n, typ) in bnd.items():
            if typ in ('empty', 'processor'):
                continue
            vals = restored[(R, name)][1] if (R, name) in restored else patch_values(txt, name, n)
            net += sum(vals)
            if name in ('inlet', 'opticalWindow'):
                inflow = (inflow or 0.0) - sum(vals)
        print('%-14s net boundary mass flow = %+.3e kg/s' % (R, net))
        nets[R] = net
        if not dry:
            dst = os.path.join(outdir, R, 'phi') if outdir else phifile
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            open(dst, 'w').write(txt)
    print('inlet mass flow (inlet + opticalWindow) = %.10e kg/s' % inflow)
    for R, net in nets.items():
        if abs(net) > tol*inflow:
            print('FAIL %s: net mass flow %.3e > %g x inflow' % (R, net, tol)); ok = False
    for (R, name), (a, vals) in restored.items():
        S = name[len(R)+4:]
        other = restored.get((S, '%s_to_%s' % (S, R)))
        if other is None:
            print('FAIL: no partner patch for %s/%s' % (R, name)); ok = False; continue
        if [abs(x) for x in a] != [abs(x) for x in other[0]]:
            print('FAIL %s/%s: face order differs from the partner patch' % (R, name)); ok = False
        if max(abs(u + w) for u, w in zip(vals, other[1])) > 1e-14*max(1e-30, max(map(abs, vals))):
            print('FAIL %s/%s: phi not equal and opposite face by face' % (R, name)); ok = False
        if abs(abs(sum(vals)) - inflow) > tol*inflow:
            print('FAIL %s/%s: |sum phi| %.10e != inflow %.10e' % (R, name, abs(sum(vals)), inflow)); ok = False
    print('ALL CHECKS PASSED' if ok else 'CHECKS FAILED')
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
