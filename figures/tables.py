#!/usr/bin/env python3
"""
Tables of paper 2 from the corrected model, as LaTeX (booktabs, no other packages) plus a compiled
preview (out/tables_preview.pdf / .png):

  out/table4_outlet_temperatures.tex   Table 4: outlet bulk temperature, pore scale vs macroscopic,
                                        deviation, and receiver efficiency of both models
  out/table3_angular_sensitivity.tex   replacement for Table 3: sensitivity of the G 5e6 / Re 300
                                        results to the angular discretisation (mesh fixed)

    tables.py
"""
import os, re, subprocess
import numpy as np
import data as D
import style as S

T_PORE = {('1e6', 300): 711.9, ('1e6', 600): 679.18, ('1e6', 900): 639.2,        # paper 1,
          ('5e6', 300): 1018.9, ('5e6', 600): 1037.65, ('5e6', 900): 1040.93,    # tab:porous_outlet_temperatures
          ('1e7', 300): 1170.88, ('1e7', 600): 1196.27, ('1e7', 900): 1204.5}
RHO, CP, A_IN, A_WIN, T_IN = 1.225, 1006.0, 2.5e-3, 6.25e-4, 300.0


def efficiency(G, Re, T_out):
    """Share of the beam power carried away by the air: mdot cp (T_out - T_in) / (G A_window)."""
    return 100 * RHO * D.U_IN[Re] * A_IN * CP * (T_out - T_IN) / (float(G) * A_WIN)


def table4():
    rows = []
    for G in D.G_LEVELS:
        for k, Re in enumerate(D.RE_LEVELS):
            tp = T_PORE[(G, Re)]
            if D.available(G, Re):
                tm = D.energy(G, Re)['T_out']
                cells = [f'{tm:.0f}', f'{100 * (tp - tm) / tp:.1f}', f'{efficiency(G, Re, tp):.1f}', f'{efficiency(G, Re, tm):.1f}']
            else:
                cells = ['--'] * 2 + [f'{efficiency(G, Re, tp):.1f}', '--']
            g = rf'${D.G_LABEL[G]}$' if k == 0 else ''
            rows.append(f'{g} & {Re} & {tp:.0f} & ' + ' & '.join(cells) + r' \\')
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


def _ray_case(rays):
    c = os.path.join(D.CASES, f'G5e6_Re300_{rays}')
    eb = open(os.path.join(c, 'log.energyBalance')).read()
    T_out = float(re.search(r'bulk T out ([\d.]+) K', eb).group(1))
    pc = open(os.path.join(c, 'figures', 'paper_comparison.txt')).read()
    front = float(re.search(r'Front face max T \[K\]\s+\d+\s+([\d.]+)', pc).group(1))
    top = float(re.search(r'Top wall max T \[K\]\s+\d+\s+([\d.]+)', pc).group(1))
    log = open(os.path.join(c, 'log.chtMultiRegionSimpleFoam')).read()
    n = len(re.findall(r'^Time = ', log, re.M))
    secs = float(re.findall(r'^ExecutionTime = ([\d.]+) s', log, re.M)[-1])
    return T_out, front, top, secs / n


def table3():
    sets = [('6x10', 6, 10, 81), ('12x10', 12, 10, 81), ('6x11', 6, 11, 90), ('8x15', 8, 15, 90)]
    sets = [t for t in sets if os.path.isfile(os.path.join(D.CASES, f'G5e6_Re300_{t[0]}', 'figures', 'paper_comparison.txt'))]
    if len(sets) < 2 or sets[-1][0] != '8x15':
        return None                      # needs the 8x15 reference run and at least one other set
    vals = {s: _ray_case(s) for s, *_ in sets}
    ref = vals['8x15'][0]
    rows = []
    for s, nphi, nth, th in sets:
        T_out, front, top, cost = vals[s]
        mark = r'$^{\dagger}$' if s == '6x11' else (r'$^{\ast}$' if nth % 2 == 0 else '')
        dev = '--' if s == '8x15' else f'{100 * (T_out - ref) / ref:+.2f}'
        rows.append(f'$({nphi}, {nth})${mark} & {4 * nphi * nth} & {th} & {T_out:.1f} & {dev} & {front:.0f} & {top:.0f} & {cost:.2f}' + r' \\')
    body = '\n'.join(rows)
    return rf"""\begin{{table}}[H]
\centering
\small
\setlength{{\tabcolsep}}{{3.5pt}}
\caption{{Sensitivity of the results to the angular discretisation for $G = 5 \times 10^{{6}}$~W\,m$^{{-2}}$ and $Re = 300$, on the adopted mesh (209\,088 cells). $\theta_b$ is the polar angle of the ordinate that carries the beam; the deviation of the outlet bulk temperature is relative to $(8, 15)$.}}
\label{{tab:grid_independence_of}}
\begin{{tabular}}{{@{{}}lrrrrrrr@{{}}}}
\toprule
$(n_\phi, n_\theta)$ & Ordinates & $\theta_b$ ($^\circ$) & $T_{{out}}$ (K) & Dev. (\%) & $T_{{f,max}}$ (K) & $T_{{w,max}}$ (K) & s/it. \\
\midrule
{body}
\bottomrule
\multicolumn{{8}}{{@{{}}l}}{{\footnotesize $^{{\dagger}}$ Adopted. $^{{\ast}}$ Even $n_\theta$: beam ordinate off the equator, negative intensities next to the}}\\
\multicolumn{{8}}{{@{{}}l}}{{\footnotesize \phantom{{$^{{\ast}}$}} porous faces. $T_{{f,max}}$: front face; $T_{{w,max}}$: top wall. Cost on 14 cores.}}\\
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
