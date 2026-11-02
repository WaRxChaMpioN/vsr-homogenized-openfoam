#!/usr/bin/env python3
"""
Tables of paper 2 from the corrected model, as LaTeX (booktabs, no other packages) plus a compiled
preview (out/tables_preview.pdf / .png):

  out/table4_outlet_temperatures.tex   Table 4: outlet bulk temperature, pore scale vs macroscopic,
                                        deviation, and receiver efficiency of both models
  out/table3_angular_sensitivity.tex   Table 3: sensitivity of the G 5e6 / Re 300 results to
                                        (a) the mesh, with the adopted ray set, once the mesh-study
                                        runs of studies/Allrun.studies exist, and (b) the angular
                                        discretisation, on the adopted mesh
  out/mesh_study.txt                    the mesh-study numbers and their grid convergence index

    tables.py
"""
import os, re, glob, subprocess
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import data as D
import style as S

T_PORE = {('1e6', 300): 711.9, ('1e6', 600): 679.18, ('1e6', 900): 639.2,        # paper 1,
          ('5e6', 300): 1018.9, ('5e6', 600): 1037.65, ('5e6', 900): 1040.93,    # tab:porous_outlet_temperatures
          ('1e7', 300): 1170.88, ('1e7', 600): 1196.27, ('1e7', 900): 1204.5}
RHO, CP, A_IN, A_WIN, T_IN = 1.225, 1006.0, 2.5e-3, 6.25e-4, 300.0


def r0(v):
    """Round half up to an integer (format() rounds 1204.5 to the even 1204; the paper prints 1205)."""
    return f"{Decimal(str(v)).quantize(Decimal('1'), ROUND_HALF_UP)}"


def r1(v):
    """Round half up to one decimal (as printed in Table 4)."""
    return float(Decimal(str(v)).quantize(Decimal('0.1'), ROUND_HALF_UP))


def num(v, fmt):
    """A signed table number in math mode, so that negatives get a minus sign, not a hyphen."""
    return '$' + format(v, fmt) + '$'


def efficiency(G, Re, T_out):
    """Share of the beam power carried away by the air: mdot cp (T_out - T_in) / (G A_window)."""
    return 100 * RHO * D.U_IN[Re] * A_IN * CP * (T_out - T_IN) / (float(G) * A_WIN)


def table4():
    rows = []
    for G in D.G_LEVELS:
        for k, Re in enumerate(D.RE_LEVELS):
            # temperatures to 0.1 K, and the deviation and efficiencies from the printed values, so
            # that they can be recomputed from the table itself
            tp = r1(T_PORE[(G, Re)])
            if D.available(G, Re):
                tm = r1(D.energy(G, Re)['T_out'])
                cells = [f'{tm:.1f}', num(100 * (tp - tm) / tp, '.1f'), f'{efficiency(G, Re, tp):.1f}', f'{efficiency(G, Re, tm):.1f}']
            else:
                cells = ['--'] * 2 + [f'{efficiency(G, Re, tp):.1f}', '--']
            g = rf'${D.G_LABEL[G]}$' if k == 0 else ''
            rows.append(f'{g} & {Re} & {tp:.1f} & ' + ' & '.join(cells) + r' \\')
        rows.append(r'\addlinespace')
    body = '\n'.join(rows[:-1])
    return rf"""\begin{{table}}[H]
\centering
\small
\setlength{{\tabcolsep}}{{3.2pt}}
\caption{{Bulk mean temperature at the channel outlet ($x = 275$~mm) from the pore-scale \cite{{effective_our}} and macroscopic models for each irradiation $G$ and Reynolds number. Deviation (Dev.) $= (T_{{pore}} - T_{{macro}})/T_{{pore}}$. The efficiency $\eta = \dot m c_p (T_{{out}} - T_{{in}})/(G A_w)$ is the share of the beam power, through the window of area $A_w$, carried away by the air.}}
\label{{tab:Outlet_temperatures_comp}}
\begin{{tabular}}{{@{{}}llrrrrr@{{}}}}
\toprule
 & & \multicolumn{{2}}{{c}}{{$T_{{out}}$ (K)}} & & \multicolumn{{2}}{{c}}{{$\eta$ (\%)}} \\
\cmidrule(lr){{3-4}} \cmidrule(l){{6-7}}
$G$ (W\,m$^{{-2}}$) & $Re$ & Pore scale & Macroscopic & Dev. (\%) & Pore scale & Macroscopic \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{table}}
"""


def _last(c, region, fo):
    """Last value written by a function object (full precision; latest restart folder)."""
    import glob
    d = sorted(glob.glob(os.path.join(c, 'postProcessing', region, fo, '*')), key=lambda p: float(os.path.basename(p)))[-1]
    rows = [l.split() for l in open(glob.glob(os.path.join(d, '*.dat'))[0]) if l.strip() and not l.startswith('#')]
    return float(rows[-1][1])


def _top_wall_max(c):
    """Maximum temperature on the top wall (z = 50 mm), from the un-rotated copy post/<t>."""
    post = os.path.join(c, 'post')
    t = max((d for d in os.listdir(post) if d.isdigit()), key=int)
    to_local, best = D.frame(post), 0.0
    for r, w in D.WALL.items():
        Cf, _ = D.patch_geometry(f'{post}/constant/{r}/polyMesh')[w]
        T = np.broadcast_to(D.read_field(f'{post}/{t}/{r}/T')[1][w], len(Cf))
        best = max(best, float(np.max(T[to_local(Cf)[:, 2] > 49.9])))
    return best


def _bulk_peak(c):
    """Peak of the face-plane bulk temperature (bulkFaces.py, figures/bulk_faces_<t>.csv)."""
    f = max(glob.glob(os.path.join(c, 'figures', 'bulk_faces_*.csv')), key=lambda p: int(re.findall(r'(\d+)\.csv$', p)[0]))
    return float(np.genfromtxt(f, delimiter=',', names=True)['bulkT'].max())


def _case(c):
    """Outlet bulk T and front-face maximum (function objects, full precision), peak of the bulk
    temperature, top-wall maximum (post/) and cost per iteration (log) of one G 5e6 / Re 300 run."""
    # cost from the longest solver log of the run: a continued run keeps its first log as
    # log.chtMultiRegionSimpleFoam.<time>, timed on an otherwise idle machine
    logs = [open(f).read() for f in glob.glob(os.path.join(c, 'log.chtMultiRegionSimpleFoam*'))]
    log = max(logs, key=lambda l: len(re.findall(r'^Time = ', l, re.M)))
    n = len(re.findall(r'^Time = ', log, re.M))
    secs = float(re.findall(r'^ExecutionTime = ([\d.]+) s', log, re.M)[-1])
    return dict(T_out=_last(c, 'exitRegion', 'outletBulkT'), front=_last(c, 'porousRegion', 'Tmax_porousRegion_to_inletRegion'),
                peak=_bulk_peak(c), wall=_top_wall_max(c), cost=secs / n)


STUDIES = os.path.join(D.CASES, 'studies')
MESHES = [('coarse', 33, 80, 87120), ('base', 44, 108, 209088), ('fine', 55, 136, 411400)]   # cross-section, axial cells


def _mesh_dir(m):
    return os.path.join(D.CASES, 'G5e6_Re300_6x11') if m == 'base' else os.path.join(STUDIES, f'mesh_{m}', 'G5e6_Re300_6x11')


def _front_cell(c):
    """Axial size (mm) of the upstream cell next to the front face, from the un-rotated copy."""
    x = np.unique(np.round(D.read_points(os.path.join(c, 'post', 'constant', 'inletRegion', 'polyMesh'))[:, 0] * 1e3, 6))
    return x[-1] - x[-2]


def gci(f1, f2, f3, r21, r32):
    """Grid convergence index of the finest of three solutions (Celik et al. 2008, J. Fluids Eng.
    130, 078001): apparent order p, extrapolated value, GCI_fine (%); p is None for oscillatory or
    no convergence, and GCI_fine is then the fine-medium difference with a safety factor 3."""
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 / e21 <= 0:
        return None, None, 300 * abs(e21 / f1)
    p = 2.0
    for _ in range(100):
        q = np.log((r21**p - 1) / (r32**p - 1))
        p = abs(np.log(abs(e32 / e21)) + q) / np.log(r21)
    ext = (r21**p * f1 - f2) / (r21**p - 1)
    return p, ext, 100 * 1.25 * abs(e21 / f1) / (r21**p - 1)


def mesh_study():
    """Mesh-study values of the three meshes, or None while a run is missing."""
    have = [m for m, *_ in MESHES if os.path.isfile(os.path.join(_mesh_dir(m), 'log.energyBalance'))
            and glob.glob(os.path.join(_mesh_dir(m), 'figures', 'bulk_faces_*.csv'))]
    if len(have) < 3:
        print('mesh study incomplete; runs ready:', have)
        return None
    vals = {m: dict(_case(_mesh_dir(m)), dx=_front_cell(_mesh_dir(m))) for m, *_ in MESHES}
    r21, r32 = (MESHES[2][3] / MESHES[1][3]) ** (1 / 3), (MESHES[1][3] / MESHES[0][3]) ** (1 / 3)
    lines = [f'refinement ratios r21 (fine/base) {r21:.4f}, r32 (base/coarse) {r32:.4f}']
    for m, n, nx, cells in MESHES:
        v = vals[m]
        lines.append(f"{m:6s} {n}x{n}x{nx} {cells:7d} cells, front cell {v['dx']:.4f} mm: T_out {v['T_out']:.2f} K, "
                     f"front max {v['front']:.1f} K, bulk peak {v['peak']:.1f} K, top wall max {v['wall']:.1f} K")
    for k, name in (('T_out', 'outlet bulk T'), ('front', 'front-face maximum'), ('peak', 'bulk-T peak'), ('wall', 'top-wall maximum')):
        p, ext, g = gci(vals['fine'][k], vals['base'][k], vals['coarse'][k], r21, r32)
        pg = 'no monotonic convergence' if p is None else f'apparent order {p:.2f}, extrapolated {ext:.2f} K'
        lines.append(f'{name:20s}: {pg}, GCI_fine {g:.3f} %, base vs fine {100 * (vals["base"][k] - vals["fine"][k]) / vals["fine"][k]:+.3f} %')
    os.makedirs(S.OUT, exist_ok=True)
    open(os.path.join(S.OUT, 'mesh_study.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    return vals


def table3():
    sets = [('6x10', 6, 10, 81), ('12x10', 12, 10, 81), ('6x11', 6, 11, 90), ('8x15', 8, 15, 90)]
    have = [t for t in sets if os.path.isfile(os.path.join(D.CASES, f'G5e6_Re300_{t[0]}', 'post', 'post.foam'))
            and glob.glob(os.path.join(D.CASES, f'G5e6_Re300_{t[0]}', 'figures', 'bulk_faces_*.csv'))]
    if len(have) < 2 or have[-1][0] != '8x15':
        missing = [f'G5e6_Re300_{t[0]}' for t in sets if t not in have]
        print('Table 3 skipped: needs the 8x15 run and at least one other ray set; missing (or no post/):', ', '.join(missing))
        return None
    sets = have
    vals = {s: _case(os.path.join(D.CASES, f'G5e6_Re300_{s}')) for s, *_ in sets}
    meshes = mesh_study()
    if meshes is not None:
        return _table3_combined(sets, vals, meshes)
    ref = vals['8x15']['T_out']
    rows = []
    for s, nphi, nth, th in sets:
        v = vals[s]
        mark = r'$^{\dagger}$' if s == '6x11' else (r'$^{\ast}$' if nth % 2 == 0 else '')
        dev = '--' if s == '8x15' else num(100 * (v['T_out'] - ref) / ref, '+.2f')
        rows.append(f"$({nphi}, {nth})${mark} & {4 * nphi * nth} & {th} & {v['T_out']:.1f} & {dev} & {r0(v['front'])} & {r0(v['wall'])} & {v['cost']:.2f}" + r' \\')
    body = '\n'.join(rows)
    return rf"""\begin{{table}}[H]
\centering
\small
\setlength{{\tabcolsep}}{{3.5pt}}
\caption{{Sensitivity of the results to the angular discretisation for $G = 5 \times 10^{{6}}$~W\,m$^{{-2}}$ and $Re = 300$, on the adopted mesh (209\,088 cells). $\theta_b$ is the polar angle of the ordinate that carries the beam; the deviation of the outlet bulk temperature is relative to $(8, 15)$.}}
\label{{tab:grid_independence_of}}
\begin{{tabular}}{{@{{}}lrrrrrrr@{{}}}}
\toprule
$(n_\phi, n_\theta)$ & Ordinates & $\theta_b$ ($^\circ$) & $T_{{out}}$ (K) & Dev. (\%) & $T_{{i,max}}$ (K) & $T_{{w,max}}$ (K) & s/it. \\
\midrule
{body}
\bottomrule
\multicolumn{{8}}{{@{{}}l}}{{\footnotesize $^{{\dagger}}$ Adopted. $^{{\ast}}$ Even $n_\theta$: beam ordinate off the equator, negative intensities next to the}}\\
\multicolumn{{8}}{{@{{}}l}}{{\footnotesize \phantom{{$^{{\ast}}$}} porous faces. $T_{{i,max}}$: front face; $T_{{w,max}}$: top wall. Cost on 14 cores.}}\\
\end{{tabular}}
\end{{table}}
"""


def _audit(c):
    """Last 'energyAudit global' imbalance of a run, % of the beam."""
    log = open(os.path.join(c, 'log.chtMultiRegionSimpleFoam')).read()
    return abs(float(re.findall(r'^energyAudit global.*\(([-\d.e+]+) % of beam\)', log, re.M)[-1]))


def _mesh_costs(base_cost):
    """Cost per iteration of the three meshes on 14 ranks: the back-to-back timing runs of
    studies/Allrun.timing (log.timing), scaled so that the adopted mesh has its cost on an otherwise
    idle machine. None while the timing has not been run."""
    f = os.path.join(STUDIES, 'log.timing')
    if not os.path.isfile(f):
        return None
    t = dict(re.findall(r'^(coarse|base|fine)\s+([\d.]+) s/it', open(f).read(), re.M))
    if len(t) < 3:
        return None
    return {m: base_cost * float(t[m]) / float(t['base']) for m in t}


def _table3_combined(sets, vals, meshes):
    """Table 3 with (a) the mesh study, ray set (6, 11), and (b) the ray sets on the adopted mesh;
    deviations of the outlet bulk temperature relative to the finest mesh and the finest ray set."""
    even = [_audit(os.path.join(D.CASES, f'G5e6_Re300_{s}')) for s, nphi, nth, th in sets if nth % 2 == 0]
    even_res = ' and '.join(f'{v:.2f}\\%' for v in even)
    rows = [r'\multicolumn{7}{@{}l}{\textit{(a) Mesh, $(n_\phi, n_\theta) = (6, 11)$}} \\']
    ref = meshes['fine']['T_out']
    costs = _mesh_costs(meshes['base']['cost'])
    for m, n, nx, cells in MESHES:
        v = meshes[m]
        label = f'{cells:,}'.replace(',', r'\,') + rf' ($ {n} \times {n} \times {nx} $)' + (r'$^{\dagger}$' if m == 'base' else '')
        dev = '--' if m == 'fine' else num(100 * (v['T_out'] - ref) / ref, '+.2f')
        cost = f"{v['cost']:.2f}" if m == 'base' else (f"{costs[m]:.2f}" if costs else '--')
        rows.append(f"{label} & {v['T_out']:.1f} & {dev} & {r0(v['front'])} & {r0(v['peak'])} & {r0(v['wall'])} & {cost}" + r' \\')
    rows += [r'\addlinespace', r'\multicolumn{7}{@{}l}{\textit{(b) Ordinates $(n_\phi, n_\theta)$, mesh of 209\,088 cells}} \\']
    ref = vals['8x15']['T_out']
    for s, nphi, nth, th in sets:
        v = vals[s]
        mark = r'$^{\dagger}$' if s == '6x11' else (r'$^{\ast}$' if nth % 2 == 0 else '')
        dev = '--' if s == '8x15' else num(100 * (v['T_out'] - ref) / ref, '+.2f')
        rows.append(f"$({nphi}, {nth})$, {4 * nphi * nth}{mark} & {v['T_out']:.1f} & {dev} & {r0(v['front'])} & {r0(v['peak'])} & {r0(v['wall'])} & {v['cost']:.2f}" + r' \\')
    body = '\n'.join(rows)
    return rf"""\begin{{table}}[H]
\centering
\small
\setlength{{\tabcolsep}}{{4pt}}
\caption{{Sensitivity of the results to (a) the mesh and (b) the angular discretisation for $G = 5 \times 10^{{6}}$~W\,m$^{{-2}}$ and $Re = 300$. Dev.\ $= (T_{{out}} - T_{{out,ref}})/T_{{out,ref}}$, with the finest mesh in (a) and $(8, 15)$ in (b) as reference. $T_{{i,max}}$: maximum temperature of the interface; $T_{{b,max}}$: peak bulk temperature; $T_{{w,max}}$: maximum temperature of the top wall; s/it.: cost per iteration on 14 cores.}}
\label{{tab:grid_independence_of}}
\begin{{tabular}}{{@{{}}lrrrrrr@{{}}}}
\toprule
 & $T_{{out}}$ (K) & Dev. (\%) & $T_{{i,max}}$ (K) & $T_{{b,max}}$ (K) & $T_{{w,max}}$ (K) & s/it. \\
\midrule
{body}
\bottomrule
\multicolumn{{7}}{{@{{}}l}}{{\footnotesize $^{{\dagger}}$ Adopted. $^{{\ast}}$ Even $n_\theta$: beam ordinate off the equator, negative intensities next to}}\\
\multicolumn{{7}}{{@{{}}l}}{{\footnotesize \phantom{{$^{{\ast}}$}} the porous faces; these runs end at an energy imbalance of {even_res} of the beam power.}}\\
\end{{tabular}}
\end{{table}}
"""


def preview(snippets):
    tex = os.path.join(S.OUT, 'tables_preview.tex')
    with open(tex, 'w') as f:
        f.write('\\documentclass[12pt]{article}\n\\usepackage{booktabs,float,amsmath}\n'
                '\\renewcommand{\\cite}[1]{[1]}\n\\begin{document}\n' + '\n'.join(snippets) + '\\end{document}\n')
    subprocess.run(['pdflatex', '-interaction=nonstopmode', '-output-directory', S.OUT, tex],
                   check=True, capture_output=True)
    subprocess.run(['pdftoppm', '-r', '150', '-png', '-singlefile', tex[:-4] + '.pdf', tex[:-4]], check=True)
    for ext in ('aux', 'log'):
        os.remove(tex[:-3] + ext)


if __name__ == '__main__':
    os.makedirs(S.OUT, exist_ok=True)
    t4, t3 = table4(), table3()
    for name, t in (('table4_outlet_temperatures', t4), ('table3_angular_sensitivity', t3)):
        if t is not None:
            open(os.path.join(S.OUT, name + '.tex'), 'w').write(t)
    preview([t for t in (t3, t4) if t is not None])
    print('tables written to', S.OUT)
