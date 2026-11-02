#!/usr/bin/env python3
"""
Proof of all paper-2 figures and tables at the manuscript's sizes (article, 12 pt, default
margins), with suggested captions: out/figures_proof.pdf. Uses only graphicx, float, booktabs and
amsmath (the manuscript's subcaption is not installed here, so panels are labelled by hand).

    make_proof.py
"""
import os, subprocess
import style as S

PANEL = r'\par\smallskip{\small (%s)}\par\smallskip'


def fig(body, caption):
    return '\\begin{figure}[H]\n\\centering\n' + body + '\n\\caption{' + caption + '}\n\\end{figure}\n'


def g(name, width=r'\linewidth'):
    return r'\includegraphics[width=%s]{%s.pdf}' % (width, name)


def stacked(names, width=r'\linewidth'):
    return '\n'.join(g(n, width) + PANEL % chr(97 + i) for i, n in enumerate(names))


def side(names, width):
    return r'\hfill'.join(r'\begin{minipage}[b]{%s}\centering %s\par{\small (%s)}\end{minipage}' % (width, g(n, r'\linewidth'), chr(97 + i))
                          for i, n in enumerate(names))


FIGS = [
    (g('fig04_openfoam_mesh'),
     r'Computational mesh (209\,088 cells) on the mid-plane $y = 25$~mm: the whole channel, the front face of the pseudo-porous zone, and the $0.081$~mm cells next to it. Orange lines: faces of the pseudo-porous zone.'),
    (stacked(['fig05a_velocity_Re300', 'fig05b_velocity_Re600', 'fig05c_velocity_Re900']),
     r'Axial velocity $u_x/U_{in}$ on the mid-plane of the channel ($y = 25$~mm), from the face mass fluxes, for (a) $Re = 300$, (b) $Re = 600$ and (c) $Re = 900$. Dashed lines: faces of the pseudo-porous zone ($x = 75$ and $125$~mm); thick line at $x = 0$: window.'),
    (g('fig06_centreline_pressure', r'0.85\linewidth'),
     r'Centreline pressure along the channel, relative to the outlet pressure, for $Re = 300$, 600 and 900, with the pressure drop across the pseudo-porous zone.'),
    (r'\includegraphics[width=0.62\linewidth]{fig07a_irradiation_isometric.pdf}' + PANEL % 'a' + g('fig07b_irradiation_midplane') + PANEL % 'b',
     r'Incident radiation $\mathcal{G}$ for $G = 5 \times 10^{6}$~W\,m$^{-2}$ and $Re = 300$: (a) on the mid-planes $y = 25$~mm and $z = 25$~mm (colour scale of (b); dashed: hidden edges); (b) on the mid-plane $y = 25$~mm. The collimated beam runs from the window to the front face of the pseudo-porous zone.'),
    (side(['fig08a_front_face_temperature', 'fig08b_back_face_temperature'], r'0.47\linewidth'),
     r'Temperature on (a) the front face and (b) the back face of the pseudo-porous zone for $Re = 300$ and $G = 5 \times 10^{6}$~W\,m$^{-2}$. Dashed: beam footprint.'),
    (stacked(['fig09a_midplane_T_G1e6_Re300', 'fig09b_midplane_T_G5e6_Re300', 'fig09c_midplane_T_G1e7_Re300']) + g('fig09_legend'),
     r'Temperature on the mid-plane of the channel for $Re = 300$: (a) $G = 1 \times 10^{6}$, (b) $5 \times 10^{6}$ and (c) $1 \times 10^{7}$~W\,m$^{-2}$ (one colour scale; isolines every 100~K).'),
    (stacked(['fig10a_midplane_T_G5e6_Re300', 'fig10b_midplane_T_G5e6_Re600', 'fig10c_midplane_T_G5e6_Re900']) + g('fig10_legend'),
     r'Temperature on the mid-plane of the channel for $G = 5 \times 10^{6}$~W\,m$^{-2}$: (a) $Re = 300$, (b) $Re = 600$ and (c) $Re = 900$ (one colour scale; isolines every 100~K).'),
    (g('fig11_legend') + r'\par\medskip' + side(['fig11a_outlet_T_G1e6', 'fig11b_outlet_T_G5e6', 'fig11c_outlet_T_G1e7'], r'0.32\linewidth'),
     r'Temperature on the exit face of the pseudo-porous zone for $Re = 300$: (a) $G = 1 \times 10^{6}$, (b) $5 \times 10^{6}$ and (c) $1 \times 10^{7}$~W\,m$^{-2}$ (one colour scale; isolines every 25~K).'),
    (g('fig12_fluid_temperature_profile', r'0.7\linewidth'),
     r'Mass-weighted average (bulk) fluid temperature along the channel for all $G$ and $Re$. The distance is measured from the inlet; shaded: pseudo-porous zone.'),
    (g('fig13_fluid_temperature_comparison', r'0.7\linewidth'),
     r'Bulk fluid temperature from the pore-scale (digitised from \cite{effective_our}) and macroscopic models for $Re = 300$. The distance is measured from the inlet; shaded: pseudo-porous zone.'),
    (stacked(['fig14a_wall_total_heat_flux', 'fig14b_wall_radiative_flux', 'fig14c_wall_temperature']),
     r'(a) Total heat flux into the wall, by conduction from the adjacent medium and by radiation (equal to the heat lost to the surroundings), (b) net radiative heat flux absorbed by the wall and (c) temperature on the top wall of the channel ($z = 50$~mm) for $G = 5 \times 10^{6}$~W\,m$^{-2}$ and $Re = 300$. The distance is measured from the inlet.'),
]


def main():
    parts = [fig(b, c) for b, c in FIGS]
    for t in ('table3_angular_sensitivity', 'table4_outlet_temperatures'):
        p = os.path.join(S.OUT, t + '.tex')
        if os.path.isfile(p):
            parts.append(open(p).read())
    tex = os.path.join(S.OUT, 'figures_proof.tex')
    with open(tex, 'w') as f:
        f.write('\\documentclass[12pt]{article}\n'
                '\\def\\DoNotLoadEpstopdf{}  % epstopdf-base needs grfext, not installed here; no EPS figures\n'
                '\\usepackage{graphicx,float,booktabs,amsmath}\n'
                '\\renewcommand{\\cite}[1]{[paper 1]}\n\\setlength{\\parindent}{0pt}\n\\begin{document}\n'
                '\\section*{Paper 2: figures and tables (proof)}\n'
                + '\n'.join(parts) + '\n\\end{document}\n')
    run = subprocess.run(['pdflatex', '-interaction=nonstopmode', 'figures_proof.tex'], cwd=S.OUT, capture_output=True, text=True)
    log = open(os.path.join(S.OUT, 'figures_proof.log')).read()
    over = log.count('Overfull \\hbox')
    for ext in ('aux', 'log'):
        os.remove(os.path.join(S.OUT, 'figures_proof.' + ext))
    print('proof written:', os.path.join(S.OUT, 'figures_proof.pdf'), f'(overfull boxes: {over}, exit {run.returncode})')


if __name__ == '__main__':
    main()
