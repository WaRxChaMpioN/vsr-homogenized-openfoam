"""Minimal OpenFOAM ascii reader: polyMesh (points/faces/boundary) and vol fields (internal + patch values)."""
import os, re
import numpy as np

_num = re.compile(rb'[-+0-9.eE]+')


def _list_after(buf, start):
    """Parse 'N ( ... )' starting at index start; return (np array of floats, end index, N)."""
    m = re.compile(rb'(\d+)\s*\(').search(buf, start)
    n = int(m.group(1))
    i = m.end()
    depth, j = 1, i
    while depth:
        c = buf[j:j + 1]
        if c == b'(':
            depth += 1
        elif c == b')':
            depth -= 1
        j += 1
    return np.array(_num.findall(buf[i:j - 1]), dtype=float), j, n


def read_points(mesh):
    buf = open(os.path.join(mesh, 'points'), 'rb').read()
    v, _, n = _list_after(buf, buf.index(b'}') + 1)
    return v.reshape(n, 3)


def read_faces(mesh):
    buf = open(os.path.join(mesh, 'faces'), 'rb').read()
    start = buf.index(b'}') + 1
    m = re.compile(rb'(\d+)\s*\(').search(buf, start)
    n = int(m.group(1))
    body = buf[m.end():]
    faces = [list(map(int, f.split())) for f in re.findall(rb'\d+\(([^)]*)\)', body)[:n]]
    return faces


def read_boundary(mesh):
    txt = open(os.path.join(mesh, 'boundary')).read()
    txt = txt[txt.index('}') + 1:]
    out = {}
    for name, body in re.findall(r'(\w+)\s*\{([^}]*)\}', txt):
        d = dict(re.findall(r'(\w+)\s+([^;]+);', body))
        out[name] = (d.get('type'), int(d['nFaces']), int(d['startFace']))
    return out


def patch_geometry(mesh):
    """Face centres, area vectors for every boundary patch: {patch: (Cf, Sf)}."""
    P = read_points(mesh)
    F = read_faces(mesh)
    geo = {}
    for name, (typ, n, s) in read_boundary(mesh).items():
        Cf = np.zeros((n, 3)); Sf = np.zeros((n, 3))
        for k in range(n):
            q = P[F[s + k]]
            c = q.mean(0)
            a = 0.5 * sum(np.cross(q[i] - c, q[(i + 1) % len(q)] - c) for i in range(len(q)))
            Cf[k], Sf[k] = c, a
        geo[name] = (Cf, Sf)
    return geo


def read_field(path):
    """Return (internal array, {patch: array or scalar})."""
    buf = open(path, 'rb').read()
    ncomp = 3 if b'volVectorField' in buf[:2000] else 1
    i = buf.index(b'internalField')
    if buf[i:i + 60].split()[1] == b'uniform':
        internal = float(_num.findall(buf[i + 22:i + 80])[0]) if ncomp == 1 else None
        j = i
    else:
        internal, j, n = _list_after(buf, i)
        internal = internal.reshape(n, ncomp) if ncomp == 3 else internal
    patches = {}
    bf = buf.index(b'boundaryField', j)
    pos = buf.index(b'{', bf) + 1
    pat = re.compile(rb'\n\s{4}(\w+)\s*\n\s{4}\{')
    heads = list(pat.finditer(buf, pos))
    for k, h in enumerate(heads):
        end = heads[k + 1].start() if k + 1 < len(heads) else len(buf)
        seg = buf[h.end():end]
        vm = re.search(rb'\n\s+value\s+(uniform|nonuniform)', seg)
        if not vm:
            continue
        if vm.group(1) == b'uniform':
            if ncomp == 1:
                patches[h.group(1).decode()] = float(_num.findall(seg[vm.end():vm.end() + 60])[0])
            else:
                patches[h.group(1).decode()] = np.array(_num.findall(seg[vm.end():vm.end() + 120])[:3], dtype=float)
        else:
            v, _, n = _list_after(seg, vm.end())
            patches[h.group(1).decode()] = v.reshape(n, ncomp) if ncomp == 3 else v
    return internal, patches


def frame(case):
    """Channel frame taken from the mesh itself (handles any residual rigid rotation, roll included):
    axis e from the inlet-plane normal, e3 from the normals of the z = 50 mm wall, e2 = e3 x e;
    origin such that the inlet plane is x = 0 and the channel spans y, z = 0..50 mm."""
    geo = patch_geometry(f'{case}/constant/inletRegion/polyMesh')
    Cf = np.vstack([geo[p][0] for p in ('inlet', 'opticalWindow')])
    Sf = np.vstack([geo[p][1] for p in ('inlet', 'opticalWindow')])
    A = np.linalg.norm(Sf, axis=1)
    e = -Sf.sum(0); e /= np.linalg.norm(e)
    P0 = (Cf * A[:, None]).sum(0) / A.sum()                      # centre of the inlet plane
    Sw = geo['inletWall'][1]
    nw = Sw / np.linalg.norm(Sw, axis=1)[:, None]
    top = nw @ np.array([0.0, 0.0, 1.0]) > 0.9                     # faces of the z = 50 mm wall
    e3 = nw[top].mean(0); e3 -= (e3 @ e) * e; e3 /= np.linalg.norm(e3)
    e2 = np.cross(e3, e)
    return lambda X: np.stack([(X - P0) @ e, (X - P0) @ e2 + 0.025, (X - P0) @ e3 + 0.025], axis=-1) * 1e3


def cluster(v, tol=0.02):
    """Index of each value's cluster (values within tol mm of their neighbour belong together) and the cluster centres."""
    order = np.argsort(v)
    vs = v[order]
    cid_sorted = np.cumsum(np.concatenate([[True], np.diff(vs) > tol])) - 1
    cid = np.empty_like(cid_sorted)
    cid[order] = cid_sorted
    return cid, np.bincount(cid, weights=v) / np.bincount(cid)
