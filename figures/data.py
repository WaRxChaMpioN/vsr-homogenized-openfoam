"""
Data for the paper-2 figures, read from the run directory of Allrun.matrix ($VSR_RUNS, default
<repository>/runs):
  G<G>_Re<Re>_6x11   heat-transfer runs; un-rotated copy in post/<time> (Allpost), profiles CSV in
                     figures/, energy balance in log.energyBalance
  flow_Re<Re>        single-region flow runs (un-rotated by construction, metres)
Coordinates are channel-frame mm: x from the inlet (0-275), y and z across the channel (0-50).
The mid-plane is y = 25 mm (the x-z plane of the manuscript's Figs 5, 7, 9 and 10), the top wall
is z = 50 mm (Fig 14). Slices are structured arrays that include the boundary values (walls,
opaque porous faces, inlet, outlet), so contours reach the walls.
"""
import os, re, sys, glob
from functools import lru_cache
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = os.path.abspath(os.environ.get('VSR_RUNS', os.path.join(REPO, 'runs')))   # run directory of Allrun.matrix
sys.path.insert(0, os.path.join(REPO, 'tools'))
from foamio import read_field, patch_geometry, read_points      # noqa: E402
from contours import frame, cluster                              # noqa: E402

RAYS = '6x11'
G_LEVELS = ['1e6', '5e6', '1e7']
RE_LEVELS = [300, 600, 900]
U_IN = {300: 0.1068, 600: 0.2136, 900: 0.3204}
REGIONS = ['inletRegion', 'porousRegion', 'exitRegion']
X_RANGE = {'inletRegion': (0.0, 75.0), 'porousRegion': (75.0, 125.0), 'exitRegion': (125.0, 275.0)}
WALL = {'inletRegion': 'inletWall', 'porousRegion': 'porousWall', 'exitRegion': 'exitWall'}
FACE = {'front': 'porousRegion_to_inletRegion', 'back': 'porousRegion_to_exitRegion'}
G_LABEL = {'1e6': r'1\times10^{6}', '5e6': r'5\times10^{6}', '1e7': r'1\times10^{7}'}
BEAM_W = 3125.0 / 5e6          # beam power per unit G: (25 mm)^2


def case_dir(G, Re):
    return os.path.join(CASES, f'G{G}_Re{Re}_{RAYS}')


def available(G, Re):
    return os.path.isfile(os.path.join(case_dir(G, Re), 'post', 'post.foam'))


def post_time(G, Re):
    post = os.path.join(case_dir(G, Re), 'post')
    return max((d for d in os.listdir(post) if re.fullmatch(r'\d+', d)), key=int)


def _scalars(internal, bnd):
    """Internal and patch values as scalars: vector fields become magnitudes (a uniform patch
    value of a vector field is a single 3-vector)."""
    internal = np.asarray(internal, float)
    if not (internal.ndim == 2 and internal.shape[1] == 3):
        return internal, {p: np.asarray(v, float) for p, v in bnd.items()}
    return (np.linalg.norm(internal, axis=1),
            {p: np.linalg.norm(np.asarray(v, float), axis=-1) for p, v in bnd.items()})


# --------------------------------------------------------------------------- heat-transfer cases
@lru_cache(maxsize=None)
def _geometry(G, Re, region):
    post = os.path.join(case_dir(G, Re), 'post')
    to_local = frame(post)
    t = post_time(G, Re)
    X = to_local(read_field(f'{post}/{t}/{region}/C')[0])
    patches = {p: (to_local(Cf), np.linalg.norm(Sf, axis=1))
               for p, (Cf, Sf) in patch_geometry(f'{post}/constant/{region}/polyMesh').items()}
    return X, patches


@lru_cache(maxsize=None)
def _field(G, Re, region, name):
    post = os.path.join(case_dir(G, Re), 'post')
    return _scalars(*read_field(f'{post}/{post_time(G, Re)}/{region}/{name}'))


def _grid(a, b, v, tol=1e-3):
    """Points (a, b, v) of a tensor grid -> 2-D arrays; duplicates averaged, gaps (corners) filled."""
    ia, ua = cluster(np.asarray(a, float), tol)
    ib, ub = cluster(np.asarray(b, float), tol)
    S = np.zeros((len(ub), len(ua))); N = np.zeros_like(S)
    np.add.at(S, (ib, ia), v); np.add.at(N, (ib, ia), 1)
    V = np.where(N > 0, S / np.maximum(N, 1), np.nan)
    for _ in range(4):                                   # fill isolated gaps from their neighbours
        if not np.isnan(V).any():
            break
        P = np.pad(V, 1, constant_values=np.nan)
        nb = np.stack([P[:-2, 1:-1], P[2:, 1:-1], P[1:-1, :-2], P[1:-1, 2:]])
        V = np.where(np.isnan(V), np.nanmean(np.where(np.isnan(nb), np.nan, nb), axis=0), V)
    A, B = np.meshgrid(ua, ub)
    return A, B, V


def _extend(A, B, V, a_range, b_range, fill=None):
    """Add the boundary columns/rows a grid lacks (patches that write no value: zeroGradient,
    noSlip, surfaceNormalFixedValue). fill: {'a0', 'a1', 'b0', 'b1'} -> constant, or edge copy."""
    fill = fill or {}
    a, b = A[0], B[:, 0]
    if a[0] > a_range[0] + 1e-3:
        a = np.r_[a_range[0], a]; V = np.c_[V[:, :1] * 0 + fill['a0'] if 'a0' in fill else V[:, :1], V]
    if a[-1] < a_range[1] - 1e-3:
        a = np.r_[a, a_range[1]]; V = np.c_[V, V[:, -1:] * 0 + fill['a1'] if 'a1' in fill else V[:, -1:]]
    if b[0] > b_range[0] + 1e-3:
        b = np.r_[b_range[0], b]; V = np.r_[V[:1] * 0 + fill['b0'] if 'b0' in fill else V[:1], V]
    if b[-1] < b_range[1] - 1e-3:
        b = np.r_[b, b_range[1]]; V = np.r_[V, V[-1:] * 0 + fill['b1'] if 'b1' in fill else V[-1:]]
    A, B = np.meshgrid(a, b)
    return A, B, V


def _plane_points(G, Re, region, name, keep_cells, keep_faces):
    X, patches = _geometry(G, Re, region)
    vi, vb = _field(G, Re, region, name)
    vi = np.broadcast_to(vi, len(X))
    m = keep_cells(X)
    pts = [(X[m], vi[m])]
    for p, (Cf, _) in patches.items():
        if p not in vb:
            continue
        f = keep_faces(Cf)
        if f.any():
            pts.append((Cf[f], np.broadcast_to(vb[p], len(Cf))[f]))
    P = np.vstack([q for q, _ in pts]); V = np.concatenate([v for _, v in pts])
    return P, V


def midplane(G, Re, name, y0=25.0, half=0.7):
    """Field on the mid-plane y = y0, one block per region: [(region, X, Z, V)]; vectors as magnitude.
    The two cell layers either side of the plane are averaged; boundary faces on the plane add the
    wall, face, inlet and outlet values."""
    blocks = []
    wall = {'b0': 0.0, 'b1': 0.0} if name == 'U' else {}           # noSlip writes no value
    for r in REGIONS:
        P, V = _plane_points(G, Re, r, name, lambda C: np.abs(C[:, 1] - y0) < half,
                             lambda C: np.abs(C[:, 1] - y0) < half)
        blocks.append((r,) + _extend(*_grid(P[:, 0], P[:, 2], V), X_RANGE[r], (0.0, 50.0), wall))
    return blocks


def face(G, Re, which, name='T'):
    """Field on the front (x = 75 mm) or back (x = 125 mm) face of the porous zone: (Y, Z, V),
    padded to the walls with the edge values."""
    X, patches = _geometry(G, Re, 'porousRegion')
    Cf, _ = patches[FACE[which]]
    v = np.broadcast_to(_field(G, Re, 'porousRegion', name)[1][FACE[which]], len(Cf))
    Y, Z, V = _grid(Cf[:, 1], Cf[:, 2], v)
    y = np.concatenate([[0.0], Y[0], [50.0]]); z = np.concatenate([[0.0], Z[:, 0], [50.0]])
    Y, Z = np.meshgrid(y, z)
    return Y, Z, np.pad(V, 1, mode='edge')


def top_wall(G, Re, name, z0=50.0):
    """Field on the top wall z = z0, one block per region: [(region, X, Y, V)], padded to the region
    ends and the side walls with the edge values."""
    blocks = []
    for r in REGIONS:
        _, patches = _geometry(G, Re, r)
        Cf, _ = patches[WALL[r]]
        v = np.broadcast_to(_field(G, Re, r, name)[1][WALL[r]], len(Cf))
        m = np.abs(Cf[:, 2] - z0) < 0.01
        X, Y, V = _grid(Cf[m, 0], Cf[m, 1], v[m])
        x = np.concatenate([[X_RANGE[r][0]], X[0], [X_RANGE[r][1]]])
        y = np.concatenate([[0.0], Y[:, 0], [50.0]])
        X, Y = np.meshgrid(x, y)
        blocks.append((r, X, Y, np.pad(V, 1, mode='edge')))
    return blocks


def profile(G, Re):
    """Axial profiles of profiles.py (x_mm, bulkT, meanT, Tmin, Tmax, Tcentre, wallT_mean, ...)."""
    t = post_time(G, Re)
    return np.genfromtxt(os.path.join(case_dir(G, Re), 'figures', f'profiles_{t}.csv'),
                         delimiter=',', names=True)


def energy(G, Re):
    """Numbers of log.energyBalance: outlet bulk T, air enthalpy gain, wall losses, closure."""
    txt = open(os.path.join(case_dir(G, Re), 'log.energyBalance')).read()
    num = lambda pat: float(re.search(pat, txt).group(1))
    return {'T_out': num(r'bulk T out ([\d.]+) K'),
            'Q_air': num(r'air enthalpy gain\s+([-\d.]+) W'),
            'Q_wall': num(r'wall losses \(external\)\s+([-\d.]+) W'),
            'closure_pct': num(r'CLOSURE\s+[-\d.]+ W\s+=\s+([-+\d.]+) %'),
            'beam': num(r'beam in\s+([\d.]+) W')}


# --------------------------------------------------------------------------------- flow cases
FLOW = lambda Re: os.path.join(CASES, f'flow_Re{Re}')


def flow_available(Re):
    return os.path.isfile(os.path.join(FLOW(Re), 'FLOW_TIME'))


@lru_cache(maxsize=None)
def _flow_geometry():
    """Cell centres and patches of the flow mesh (the same block mesh for every Re), in mm."""
    d = FLOW(300)
    t = open(os.path.join(d, 'FLOW_TIME')).read().strip()
    C = 1e3 * read_field(f'{d}/{t}/C')[0]
    patches = {p: (1e3 * Cf, np.linalg.norm(Sf, axis=1) * 1e6)
               for p, (Cf, Sf) in patch_geometry(f'{d}/constant/polyMesh').items()}
    return C, patches


@lru_cache(maxsize=None)
def _flow_field(Re, name):
    d = FLOW(Re)
    t = open(os.path.join(d, 'FLOW_TIME')).read().strip()
    return _scalars(*read_field(f'{d}/{t}/{name}'))


def flow_midplane(Re, name, y0=25.0, half=0.7):
    """Flow field on the mid-plane y = y0 over the whole channel: (X, Z, V)."""
    C, patches = _flow_geometry()
    vi, vb = _flow_field(Re, name)
    m = np.abs(C[:, 1] - y0) < half
    pts, vals = [C[m]], [np.broadcast_to(vi, len(C))[m]]
    for p, (Cf, _) in patches.items():
        f = np.abs(Cf[:, 1] - y0) < half
        if p in vb and f.any():
            pts.append(Cf[f]); vals.append(np.broadcast_to(vb[p], len(Cf))[f])
    P = np.vstack(pts)
    fill = {'a0': U_IN[Re], 'b0': 0.0, 'b1': 0.0} if name == 'U' else {}   # inlet U_in, noSlip walls
    return _extend(*_grid(P[:, 0], P[:, 2], np.concatenate(vals)), (0.0, 275.0), (0.0, 50.0), fill)


@lru_cache(maxsize=None)
def _flow_faces():
    """Face centres (mm) and area vectors (m2) of the flow mesh, and the number of internal faces."""
    from foamio import read_faces, read_boundary
    mesh = os.path.join(FLOW(300), 'constant', 'polyMesh')
    P = read_points(mesh)
    F = np.array(read_faces(mesh))                       # block mesh: every face is a quad
    Sf = 0.5 * np.cross(P[F[:, 2]] - P[F[:, 0]], P[F[:, 3]] - P[F[:, 1]])
    n_internal = min(s for _, _, s in read_boundary(mesh).values())
    return 1e3 * P[F].mean(1), Sf, n_internal


def flow_midplane_ux(Re, y0=25.0, half=0.7):
    """Axial velocity on the mid-plane from the face mass fluxes, u_x = phi / (rho S_x), on the faces
    normal to x: the exact transported velocity, free of the cell-centre wiggle that the
    reconstructed U shows next to the porous faces. (X, Z, u_x)."""
    Cf, Sf, n = _flow_faces()
    from foamio import read_boundary
    d = FLOW(Re)
    t = open(os.path.join(d, 'FLOW_TIME')).read().strip()
    phi_i, phi_b = read_field(f'{d}/{t}/phi')
    rho = float(np.mean(read_field(f'{d}/{t}/rho')[0]))
    phi = np.zeros(len(Cf)); phi[:n] = phi_i
    for p, (_, nf, start) in read_boundary(os.path.join(FLOW(300), 'constant', 'polyMesh')).items():
        phi[start:start + nf] = np.broadcast_to(phi_b.get(p, 0.0), nf)
    m = (np.abs(Sf[:, 0]) > 0.99 * np.linalg.norm(Sf, axis=1)) & (np.abs(Cf[:, 1] - y0) < half)
    ux = phi[m] / (rho * Sf[m, 0])
    return _extend(*_grid(Cf[m, 0], Cf[m, 2], ux), (0.0, 275.0), (0.0, 50.0), {'b0': 0.0, 'b1': 0.0})


def centreline(Re, name='p_rgh', half=0.7):
    """Field along the centreline y = z = 25 mm (mean of the four cells around it), with the inlet
    and outlet patch values at x = 0 and 275 mm: (x, v)."""
    C, patches = _flow_geometry()
    vi, vb = _flow_field(Re, name)
    m = (np.abs(C[:, 1] - 25) < half) & (np.abs(C[:, 2] - 25) < half)
    ix, ux = cluster(C[m, 0], 1e-3)
    v = np.bincount(ix, weights=np.broadcast_to(vi, len(C))[m]) / np.bincount(ix)
    ends = []
    for p, x in (('opticalWindow', 0.0), ('outlet', 275.0)):
        Cf, _ = patches[p]
        f = (np.abs(Cf[:, 1] - 25) < half) & (np.abs(Cf[:, 2] - 25) < half)
        ends.append((x, float(np.mean(np.broadcast_to(vb[p], len(Cf))[f]))))
    x = np.concatenate([[ends[0][0]], ux, [ends[1][0]]])
    return x, np.concatenate([[ends[0][1]], v, [ends[1][1]]])


@lru_cache(maxsize=None)
def mesh_lines():
    """Grid-line positions of the block mesh (mm): unique x, y and z of the mesh points."""
    P = 1e3 * read_points(os.path.join(FLOW(300), 'constant', 'polyMesh'))
    return tuple(cluster(P[:, k], 1e-4)[1] for k in range(3))
