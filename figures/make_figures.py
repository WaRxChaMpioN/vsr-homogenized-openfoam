#!/usr/bin/env python3
"""
Figures of paper 2 from the runs (6x11 rays), named as the files in the manuscript's figures/
folder (PDF + PNG in out/):

    make_figures.py [fig04 fig05 ... | all]

A panel whose case has not finished is drawn as 'pending'. Shared colour ranges: every panel of a
figure uses one colour scale, written once to the figure's legend file where the manuscript has
one (Figs 9-11), else into each panel.
"""
import os, sys
import numpy as np
from matplotlib import ticker
from matplotlib.colors import Normalize, LogNorm
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

import style as S
import data as D

S.apply()
W = S.LINEWIDTH
SUB = {300: 'a', 600: 'b', 900: 'c', '1e6': 'a', '5e6': 'b', '1e7': 'c'}


def strip_axes(fig_w=W, cbar=False):
    """Figure with one channel side view (275 x 50 mm), optional vertical colour bar on the right."""
    ax_w = fig_w - (0.62 if cbar else 0.45)
    fig = plt.figure(figsize=(fig_w, ax_w * 50 / 275 + 0.45))
    ax = fig.add_axes([0.40 / fig_w, 0.36 / fig.get_figheight(), ax_w / fig_w, ax_w * 50 / 275 / fig.get_figheight()])
    cax = fig.add_axes([(0.40 + ax_w + 0.08) / fig_w, ax.get_position().y0, 0.07 / fig_w, ax.get_position().height]) if cbar else None
    return fig, ax, cax


def legend_bar(name, norm, label, ticks, fmt=None, minor=None, width=W, log=False):
    """A stand-alone horizontal colour bar (the manuscript's figNN_legend files)."""
    fig = plt.figure(figsize=(width, 0.42))
    cax = fig.add_axes([0.06, 0.62, 0.88, 0.22])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=S.RAINBOW)
    cb = S.colorbar(fig, sm, cax, label, ticks=ticks, fmt=fmt, minor=minor)
    cb.ax.xaxis.set_label_position('bottom')
    S.save(fig, name)


def T_range(cases, getter):
    vals = [getter(*c) for c in cases if D.available(*c)]
    return (min(v[0] for v in vals), max(v[1] for v in vals)) if vals else (300.0, 2300.0)


def mid_T_limits(G, Re):
    blocks = D.midplane(G, Re, 'T')
    return (min(np.nanmin(b[3]) for b in blocks), max(np.nanmax(b[3]) for b in blocks))


# ------------------------------------------------------------------------------------- Fig 4
def fig04():
    """Computational mesh: side view on the mid-plane and two close-ups at the front face."""
    xs, ys, zs = D.mesh_lines()
    fig = plt.figure(figsize=(W, 3.05))
    top = fig.add_axes([0.07, 0.60, 0.91, 0.91 * W * 50 / 275 / 3.05])
    zooms = [(fig.add_axes([0.07, 0.10, 0.40, 0.36]), (62.5, 87.5), (0.0, 15.0)),
             (fig.add_axes([0.575, 0.10, 0.40, 0.36]), (74.0, 76.5), (0.0, 1.5))]
    for ax, (x0, x1), (z0, z1) in [(top, (0, 275), (0, 50))] + [(a, xr, zr) for a, xr, zr in zooms]:
        lw = 0.12 if ax is top else 0.35
        ax.vlines(xs[(xs >= x0) & (xs <= x1)], z0, z1, color='k', lw=lw)
        ax.hlines(zs[(zs >= z0) & (zs <= z1)], x0, x1, color='k', lw=lw)
        for x in (75, 125):
            if x0 <= x <= x1:
                ax.axvline(x, color=S.SERIES[1], lw=0.9 if ax is not top else 0.6)
        ax.set_xlim(x0, x1); ax.set_ylim(z0, z1)
        ax.set_xlabel(r'$x$ (mm)'); ax.set_ylabel(r'$y$ (mm)')
    top.set_aspect('equal')
    top.xaxis.set_major_locator(ticker.MultipleLocator(25)); top.yaxis.set_major_locator(ticker.MultipleLocator(25))
    top.set_xlabel(r'$x$ (mm)', labelpad=1)
    for (ax, (x0, x1), (z0, z1)), host, colour in zip(zooms, (top, zooms[0][0]), (S.SERIES[0], S.SERIES[2])):
        host.add_patch(Rectangle((x0, z0), x1 - x0, z1 - z0, fill=False, ec=colour, lw=1.0, zorder=6))
        for s in ax.spines.values():
            s.set_edgecolor(colour); s.set_linewidth(1.0)
    zooms[0][0].set_title(r'front face of the porous zone ($x = 75$ mm)', fontsize=8)
    zooms[1][0].set_title(r'$0.081$ mm cells at the front face', fontsize=8)
    zooms[1][0].xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    top.text(100, 52, r'porous zone', ha='center', va='bottom', fontsize=8)
    S.save(fig, 'fig04_openfoam_mesh')


# ------------------------------------------------------------------------------------- Fig 5
def fig05():
    """Axial velocity on the mid-plane, u_x/U_in from the face mass fluxes, Re = 300, 600, 900."""
    norm = Normalize(0.0, 2.0)
    for Re in D.RE_LEVELS:
        fig, ax, cax = strip_axes(cbar=True)
        if D.flow_available(Re):
            X, Z, U = D.flow_midplane_ux(Re)
            m = S.field(ax, X, Z, U / D.U_IN[Re], norm, levels=np.arange(0.1, 2.0, 0.1))
            S.colorbar(fig, m, cax, r'$u_x/U_{in}$', ticks=[0, 0.5, 1.0, 1.5, 2.0],
                       orientation='vertical', minor=ticker.MultipleLocator(0.1))
        else:
            S.pending(ax, f'Re = {Re}: run pending')
        S.channel(ax, xlabel=(Re == 900))
        ax.set_title(rf'$Re = {Re}$, $U_{{in}} = {D.U_IN[Re]:.4f}$ m\,s$^{{-1}}$', loc='right', fontsize=8, pad=2)
        S.save(fig, f'fig05{SUB[Re]}_velocity_Re{Re}')


# ------------------------------------------------------------------------------------- Fig 6
def fig06():
    """Centreline pressure relative to the outlet, Re = 300, 600, 900."""
    fig, ax = plt.subplots(figsize=(0.85 * W, 2.7))
    ax.axvspan(75, 125, color='#ececec', lw=0)
    ax.text(106, 0.86, 'porous\nzone', transform=ax.get_xaxis_transform(), ha='center', va='center', fontsize=8,
            linespacing=1.0)
    for k, Re in enumerate(D.RE_LEVELS):
        if not D.flow_available(Re):
            continue
        x, p = D.centreline(Re)
        dp = np.interp(75, x, p) - np.interp(125, x, p)
        ax.plot(x, p, color=S.SERIES[k], ls=S.LINESTYLES[k], marker=S.MARKERS[k], markevery=0.06, ms=3.5,
                mfc='white', mew=0.9, label=rf'$Re = {Re}$ ($\Delta p = {dp:.1f}$ Pa)')
    ax.set_xlim(0, 275); ax.set_ylim(bottom=-2)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(25)); ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.set_xlabel(r'$x$ (mm)'); ax.set_ylabel(r'$p - p_{out}$ (Pa)')
    ax.legend(loc='upper right')
    S.save(fig, 'fig06_centreline_pressure')


# ------------------------------------------------------------------------------------- Fig 7
def fig07():
    """Radiation field G for G = 5e6, Re = 300: (a) two mid-planes in 3-D, (b) the y = 25 mm mid-plane."""
    G, Re = '5e6', 300
    norm = LogNorm(1e4, 1e7)
    if not D.available(G, Re):
        return
    # (b) mid-plane strip
    fig, ax, cax = strip_axes(cbar=True)
    for r, X, Z, V in D.midplane(G, Re, 'G'):
        m = S.field(ax, X, Z, np.clip(V, 1e4, None), norm, levels=10.0 ** np.arange(4.5, 7.01, 0.5))
    S.channel(ax)
    S.colorbar(fig, m, cax, r'$\mathcal{G}$ (W\,m$^{-2}$)', ticks=[1e4, 1e5, 1e6, 1e7], orientation='vertical',
               fmt=S.sci_formatter())
    S.save(fig, 'fig07b_irradiation_midplane')
    # (a) orthographic view from upstream, in front and above: the mid-planes y = 25 mm (vertical)
    # and z = 25 mm (horizontal). Painter's order: back half of the horizontal plane, the vertical
    # plane, front half of the horizontal plane (they only cross along y = z = 25 mm).
    proj, camera = _projector(az=-128.0, el=26.0)
    fig, ax = plt.subplots(figsize=(0.5 * W, 0.5 * W * 0.66))
    line = lambda p, q, **kw: ax.plot(*proj(*np.array([p, q], float).T), **kw)
    hidden = [e for e in _box_edges(0, 275) if not _edge_visible(e, camera)]
    for p, q in hidden:                                  # edges behind the planes: dashed, drawn first
        line(p, q, color='k', lw=0.45, ls=(0, (3, 2)))
    draw = lambda X, Y, Z, V: ax.pcolormesh(*proj(X, Y, Z), np.clip(V, 1e4, 1e7), cmap=S.RAINBOW, norm=norm,
                                            shading='gouraud', rasterized=True)
    horiz = _plane_z(G, Re, 25.0)
    for part in ('back', 'vertical', 'front'):
        if part == 'vertical':
            for r, X, Z, V in D.midplane(G, Re, 'G'):
                draw(X, np.full_like(X, 25.0), Z, V)
            continue
        for r, X, Y, V in horiz:
            rows = Y[:, 0] >= 25.0 - 1e-9 if part == 'back' else Y[:, 0] <= 25.0 + 1e-9
            draw(X[rows], Y[rows], np.full_like(X[rows], 25.0), V[rows])
    for p, q in _box_edges(75, 125):                     # porous zone, thin, over the planes
        line(p, q, color='k', lw=0.4)
    for p, q in _box_edges(0, 275):
        if (p, q) not in hidden:
            line(p, q, color='k', lw=0.7)
    line((0, 12.5, 12.5), (0, 37.5, 12.5), color='k', lw=1.5); line((0, 37.5, 12.5), (0, 37.5, 37.5), color='k', lw=1.5)
    line((0, 37.5, 37.5), (0, 12.5, 37.5), color='k', lw=1.5); line((0, 12.5, 37.5), (0, 12.5, 12.5), color='k', lw=1.5)
    ax.set_aspect('equal'); ax.axis('off')
    at = lambda *p: tuple(float(v) for v in proj(*np.array(p, float)))
    kw = dict(fontsize=7.5, arrowprops=dict(arrowstyle='-', lw=0.5, color='k', shrinkA=1, shrinkB=0))
    ax.annotate('window', xy=at(0, 37.5, 25), xytext=at(0, 75, 70), ha='right', va='bottom', **kw)
    ax.annotate('porous zone', xy=at(100, 0, 0), xytext=at(100, -25, -38), ha='center', va='top', **kw)
    ax.annotate('outlet', xy=at(275, 25, 50), xytext=at(275, 10, 85), ha='center', va='bottom', **kw)
    S.save(fig, 'fig07a_irradiation_isometric')


def _projector(az, el):
    """Orthographic projection; returns proj(x, y, z) -> (u, v) and the unit vector towards the
    camera, (cos el cos az, cos el sin az, sin el)."""
    a, e = np.radians(az), np.radians(el)
    right = np.array([-np.sin(a), np.cos(a), 0.0])
    up = np.array([-np.sin(e) * np.cos(a), -np.sin(e) * np.sin(a), np.cos(e)])
    def proj(x, y, z):
        P = np.stack([np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)], -1)
        return P @ right, P @ up
    return proj, np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])


def _box_edges(x0, x1):
    e = []
    for y in (0, 50):
        for z in (0, 50):
            e.append(((x0, y, z), (x1, y, z)))
    for x in (x0, x1):
        e += [((x, 0, 0), (x, 50, 0)), ((x, 50, 0), (x, 50, 50)), ((x, 50, 50), (x, 0, 50)), ((x, 0, 50), (x, 0, 0))]
    return e


def _edge_visible(edge, camera, box=((0, 275), (0, 50), (0, 50))):
    """An edge of the channel box is hidden when both faces meeting at it face away from the camera."""
    (p, q) = np.array(edge, float)
    faces = []
    for k, (lo, hi) in enumerate(box):
        if p[k] == q[k] and p[k] in (lo, hi):
            faces.append(-camera[k] if p[k] == lo else camera[k])     # outward normal . camera
    return any(f > 0 for f in faces)


def _plane_z(G, Re, z0, half=0.7):
    """Field on the plane z = z0 (like midplane, in x-y)."""
    blocks = []
    for r in D.REGIONS:
        P, V = D._plane_points(G, Re, r, 'G', lambda C: np.abs(C[:, 2] - z0) < half, lambda C: np.abs(C[:, 2] - z0) < half)
        blocks.append((r,) + D._extend(*D._grid(P[:, 0], P[:, 1], V), D.X_RANGE[r], (0.0, 50.0)))
    return blocks


# ------------------------------------------------------------------------------------- Fig 8
def fig08():
    """Temperature on the front and back faces of the porous zone, G = 5e6, Re = 300."""
    G, Re = '5e6', 300
    if not D.available(G, Re):
        return
    for which, name, step in (('front', 'fig08a_front_face_temperature', 100), ('back', 'fig08b_back_face_temperature', 10)):
        Y, Z, V = D.face(G, Re, which)
        lo, hi = np.floor(V.min() / step) * step, np.ceil(V.max() / step) * step
        _face_panel(Y, Z, V, Normalize(lo, hi), np.arange(lo + step, hi, step), name, cbar=True)


def _face_panel(Y, Z, V, norm, levels, name, cbar, w=0.45 * W):
    fig = plt.figure(figsize=(w, w * (0.86 if cbar else 1.0)))
    side = (w - (0.95 if cbar else 0.45))
    ax = fig.add_axes([0.40 / w, 0.34 / fig.get_figheight(), side / w, side / fig.get_figheight()])
    m = S.field(ax, Y, Z, V, norm, levels=levels, lw=0.35)
    ax.plot([12.5, 37.5, 37.5, 12.5, 12.5], [12.5, 12.5, 37.5, 37.5, 12.5], color='k', lw=0.6, ls=(0, (3, 2)))
    ax.set_xlim(0, 50); ax.set_ylim(0, 50); ax.set_aspect('equal')
    for a in (ax.xaxis, ax.yaxis):
        a.set_major_locator(ticker.MultipleLocator(10)); a.set_minor_locator(ticker.MultipleLocator(5))
    ax.set_xlabel(r'$z$ (mm)'); ax.set_ylabel(r'$y$ (mm)')
    if cbar:
        cax = fig.add_axes([(0.40 + side + 0.08) / w, ax.get_position().y0, 0.07 / w, ax.get_position().height])
        span = norm.vmax - norm.vmin
        step = next(s for s in (10, 20, 25, 50, 100, 200, 250, 500) if span / s <= 6)
        S.colorbar(fig, m, cax, r'$T$ (K)', orientation='vertical', fmt=S.int_formatter(),
                   ticks=ticker.MultipleLocator(step), minor=ticker.MultipleLocator(step / 2))
    S.save(fig, name)


# -------------------------------------------------------------------------------- Figs 9, 10
def _midplane_T_figure(cases, names, legend, step, xlabel_last=True):
    lo, hi = 300.0, T_range(cases, mid_T_limits)[1]
    hi = np.ceil(hi / 100) * 100
    norm = Normalize(lo, hi)
    for k, ((G, Re), name) in enumerate(zip(cases, names)):
        fig, ax, _ = strip_axes()
        if D.available(G, Re):
            for r, X, Z, V in D.midplane(G, Re, 'T'):
                S.field(ax, X, Z, V, norm, levels=S.nice_levels(lo + step, hi, step))
        else:
            S.pending(ax, rf'$G = {D.G_LABEL[G]}$ W\,m$^{{-2}}$, $Re = {Re}$: run pending')
        S.channel(ax, xlabel=(k == len(cases) - 1) or not xlabel_last)
        S.save(fig, name)
    ticks = S.nice_levels(lo, hi, 200 if hi - lo < 2500 else 400)
    legend_bar(legend, norm, r'$T$ (K)', ticks=ticks, fmt=S.int_formatter(), minor=ticker.MultipleLocator(step))


def fig09():
    """Mid-plane temperature, Re = 300, G = 1e6, 5e6, 1e7 (one colour scale)."""
    cases = [(G, 300) for G in D.G_LEVELS]
    _midplane_T_figure(cases, [f'fig09{SUB[G]}_midplane_T_G{G}_Re300' for G, _ in cases], 'fig09_legend', 100)


def fig10():
    """Mid-plane temperature, G = 5e6, Re = 300, 600, 900 (one colour scale)."""
    cases = [('5e6', Re) for Re in D.RE_LEVELS]
    _midplane_T_figure(cases, [f'fig10{SUB[Re]}_midplane_T_G5e6_Re{Re}' for _, Re in cases], 'fig10_legend', 100)


# ------------------------------------------------------------------------------------ Fig 11
def fig11():
    """Temperature on the exit (back) face of the porous zone, Re = 300, G = 1e6, 5e6, 1e7."""
    cases = [(G, 300) for G in D.G_LEVELS]
    faces = {c: D.face(*c, 'back') for c in cases if D.available(*c)}
    lo = np.floor(min(v[2].min() for v in faces.values()) / 100) * 100 if faces else 300
    hi = np.ceil(max(v[2].max() for v in faces.values()) / 100) * 100 if faces else 1400
    norm = Normalize(lo, hi)
    for G, Re in cases:
        name = f'fig11{SUB[G]}_outlet_T_G{G}'
        if (G, Re) in faces:
            _face_panel(*faces[(G, Re)], norm, np.arange(lo + 25, hi, 25), name, cbar=False, w=0.32 * W)
        else:
            fig, ax = plt.subplots(figsize=(0.32 * W, 0.32 * W)); S.pending(ax, 'run pending')
            ax.set_xticks([]); ax.set_yticks([]); S.save(fig, name)
    legend_bar('fig11_legend', norm, r'$T$ (K)', ticks=S.nice_levels(lo, hi, 100), fmt=S.int_formatter(),
               minor=ticker.MultipleLocator(25))


# ------------------------------------------------------------------------------------ Fig 12
def bulk_profile(G, Re):
    """Bulk temperature on the face planes, from the face mass fluxes: 300 K on the inlet plane, the
    face temperature on the two porous faces and the outlet bulk temperature of the energy balance at
    x = 275 mm (the Table 4 value)."""
    b = D.bulk_faces(G, Re)
    assert abs(b['bulkT'][-1] - D.energy(G, Re)['T_out']) < 0.01
    return b['x_mm'], b['bulkT']


def fig12():
    """Bulk (mass-weighted) temperature along the channel, all G and Re."""
    fig, ax = plt.subplots(figsize=(0.7 * W, 2.95))
    ax.axvspan(75, 125, color='#ececec', lw=0)
    top = 0.0
    for k, Re in enumerate(D.RE_LEVELS):
        for G in D.G_LEVELS:
            if not D.available(G, Re):
                continue
            x, T = bulk_profile(G, Re)
            top = max(top, T.max())
            ax.plot(x, T, color=S.SERIES[k], ls=S.LINESTYLES[k], lw=1.1, label=rf'$Re = {Re}$')
    _profile_axes(ax, top)
    _group_labels(ax, [(G, [bulk_profile(G, Re) for Re in D.RE_LEVELS if D.available(G, Re)]) for G in D.G_LEVELS])
    handles, labels = ax.get_legend_handles_labels()
    keep = {l: h for h, l in zip(handles, labels)}
    ax.legend(list(keep.values()), list(keep.keys()), loc='upper right', bbox_to_anchor=(1.0, 0.97))
    S.save(fig, 'fig12_fluid_temperature_profile')


def _group_labels(ax, groups, min_gap=0.075):
    """Direct labels right of the axes, level with where each G group ends (mean outlet value),
    pushed apart to at least min_gap (axes fraction) so they never overlap."""
    ends = [(G, np.mean([np.interp(275.0, x, T) for x, T in curves])) for G, curves in groups if curves]
    if not ends:
        return
    lo, hi = ax.get_ylim()
    ys = sorted([(G, (y - lo) / (hi - lo)) for G, y in ends], key=lambda e: e[1])
    for k in range(1, len(ys)):
        if ys[k][1] - ys[k - 1][1] < min_gap:
            ys[k] = (ys[k][0], ys[k - 1][1] + min_gap)
    for G, f in ys:
        ax.text(1.015, f, rf'${D.G_LABEL[G]}$', transform=ax.transAxes, ha='left', va='center', fontsize=7.5, clip_on=False)
    ax.text(1.015, 1.0, r'$G$ (W\,m$^{-2}$)', transform=ax.transAxes, ha='left', va='top', fontsize=7.5, clip_on=False)


def _profile_axes(ax, top):
    ax.set_xlim(0, 275)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(50)); ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(200 if top < 2200 else 400))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(50 if top < 2200 else 100))
    ax.set_ylim(250, np.ceil(1.08 * top / 100) * 100)
    ax.set_xlabel(r'$x$ (mm)'); ax.set_ylabel(r'Bulk temperature, $T_b$ (K)')
    ax.text(100, 0.03, 'porous\nzone', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=7.5,
            linespacing=1.0)


# ------------------------------------------------------------------------------------ Fig 13
def fig13():
    """Bulk temperature, macroscopic (this model) vs pore scale (paper 1, digitised), Re = 300."""
    import os
    csv = os.path.join(S.HERE, 'data', 'pore_scale_bulkT_digitised.csv')
    pore = np.genfromtxt(csv, delimiter=',', names=True, dtype=None, encoding=None) if os.path.isfile(csv) else None
    fig, ax = plt.subplots(figsize=(0.7 * W, 2.95))
    ax.axvspan(75, 125, color='#ececec', lw=0)
    top = 0.0
    groups = []
    for G in D.G_LEVELS:
        group = []
        if pore is not None:
            m = np.isclose(pore['G'].astype(float), float(G)) & (pore['Re'] == 300)
            if m.any():
                o = np.argsort(pore['x_mm_from_inlet'][m])
                x, T = pore['x_mm_from_inlet'][m][o], pore['T_K'][m][o]
                group.append((x, T)); top = max(top, T.max())
                pts = np.where(np.isin(np.round(x, 1), [25, 50] + list(range(75, 130, 5)) + [150, 175, 200, 225, 250, 275]))[0]
                # markers on paper 1's stations, the same on every curve (its vertices; one is hidden under another curve)
                ax.plot(x, T, color=S.SERIES[1], ls=S.LINESTYLES[1], lw=1.0, marker='o', markevery=list(pts), ms=2.6,
                        mfc='white', mec=S.SERIES[1], mew=0.8, label='pore scale')
        if D.available(G, 300):
            x, T = bulk_profile(G, 300)
            group.append((x, T)); top = max(top, T.max())
            ax.plot(x, T, color=S.SERIES[0], ls=S.LINESTYLES[0], lw=1.1, label='macroscopic (this work)')
        groups.append((G, group))
    _profile_axes(ax, top)
    _group_labels(ax, groups)
    handles, labels = ax.get_legend_handles_labels()
    keep = {l: h for h, l in zip(handles, labels)}
    ax.legend(list(keep.values()), list(keep.keys()), loc='upper right', bbox_to_anchor=(1.0, 0.97))
    S.save(fig, 'fig13_fluid_temperature_comparison')


# ------------------------------------------------------------------------------------ Fig 14
def fig14():
    """Top wall (z = 50 mm), G = 5e6, Re = 300: (a) total heat flux into the wall (conduction +
    radiation = loss to the surroundings), (b) net radiative heat flux absorbed by the wall,
    (c) temperature. wallHeatFlux = alpha snGrad(h) - qr is the flux into the channel, so the flux
    into the wall is its negative."""
    G, Re = '5e6', 300
    if not D.available(G, Re):
        return
    panels = [('fig14a_wall_total_heat_flux', [(r, X, Y, -1e-3 * V) for r, X, Y, V in D.top_wall(G, Re, 'wallHeatFlux')],
               r'$q_w$ (kW\,m$^{-2}$)', 20),
              ('fig14b_wall_radiative_flux', [(r, X, Y, 1e-3 * V) for r, X, Y, V in D.top_wall(G, Re, 'qr')],
               r'$q_{rad,w}$ (kW\,m$^{-2}$)', 20),
              ('fig14c_wall_temperature', D.top_wall(G, Re, 'T'), r'$T$ (K)', 50)]
    flux = [b for _, blocks, _, _ in panels[:2] for b in blocks]          # (a) and (b): one colour scale
    rng = {name: (min(np.nanmin(b[3]) for b in (flux if k < 2 else blocks)), max(np.nanmax(b[3]) for b in (flux if k < 2 else blocks)))
           for k, (name, blocks, _, _) in enumerate(panels)}
    figs = []
    for name, blocks, label, step in panels:
        lo = np.floor(rng[name][0] / step) * step
        hi = np.ceil(rng[name][1] / step) * step
        norm = Normalize(lo, hi)
        fig, ax, cax = strip_axes(cbar=True)
        for r, X, Y, V in blocks:
            m = S.field(ax, X, Y, V, norm, levels=[v for v in S.nice_levels(lo + step, hi, step) if abs(v) > 1e-9])
            if lo < 0 < hi:                                  # where the wall turns from gaining to losing
                ax.contour(X, Y, V, levels=[0.0], colors='k', linewidths=0.8, linestyles=(0, (4, 2)))
        S.channel(ax, window=False, ylabel=r'$z$ (mm)', xlabel=name.endswith('temperature'))
        span = hi - lo
        tick = next(s for s in (10, 20, 25, 50, 100, 200, 250, 500) if span / s <= 6)
        S.colorbar(fig, m, cax, label, orientation='vertical', fmt=S.int_formatter(), ticks=ticker.MultipleLocator(tick),
                   minor=ticker.MultipleLocator(tick / 2))
        figs.append((name, fig))
    from matplotlib.transforms import Bbox                     # same canvas for the three stacked panels
    box = Bbox.union([f.get_tightbbox(f.canvas.get_renderer()) for _, f in figs]).padded(0.02)
    for name, fig in figs:
        S.save(fig, name, bbox=box)
    for old in ('fig14a_wall_radiative_flux', 'fig14b_wall_temperature'):      # previous two-panel layout
        for ext in ('.pdf', '.png'):
            p = os.path.join(S.OUT, old + ext)
            if os.path.isfile(p):
                os.remove(p)


FIGS = {f'fig{n:02d}': globals()[f'fig{n:02d}'] for n in (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)}

if __name__ == '__main__':
    todo = sys.argv[1:] or ['all']
    failed = []
    for k in (FIGS if todo == ['all'] else todo):
        print('making', k, flush=True)
        try:
            FIGS[k]()
        except Exception as e:                          # one missing input must not stop the others
            failed.append(k)
            print(f'  {k} FAILED: {type(e).__name__}: {e}', flush=True)
    if failed:
        sys.exit('failed: ' + ' '.join(failed))
