"""
Shared look of the paper-2 figures.

- Text is set by LaTeX in Computer Modern, as in the manuscript; the PDFs embed Type 1 fonts.
- Every field plot uses one colour map: ParaView's 'Rainbow Uniform', the map of the manuscript's
  Figs 9-11 (sampled from figures/fig11_legend.png), with thin isolines in a darker shade of the
  same colours.
- Line plots use the categorical palette in fixed order, with line styles and markers as a second
  cue, and boxed axes with inward ticks.
- Widths are for the manuscript (article class, 12 pt): \\linewidth = 390 pt = 5.40 in.
- save(fig, name) writes out/<name>.pdf (vector; field fills rasterised at 600 dpi) and
  out/<name>.png (300 dpi, rasterised from the PDF with pdftoppm, since matplotlib's usetex PNG
  path needs dvipng).
"""
import os, subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LinearSegmentedColormap, Normalize, LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out')
if os.environ.get('VSR_TEX_BIN'):                       # LaTeX not on PATH: e.g. VSR_TEX_BIN=~/.local/bin (TinyTeX)
    os.environ['PATH'] = os.path.expanduser(os.environ['VSR_TEX_BIN']) + os.pathsep + os.environ.get('PATH', '')

LINEWIDTH = 390 / 72.27                                  # in
RAINBOW = LinearSegmentedColormap.from_list(
    'rainbow_uniform', ['#0562fe', '#058ddd', '#05b2af', '#05d876', '#73f51f', '#fcff04',
                        '#ffcd37', '#ff9a37', '#fe5c36', '#e81549', '#b00d69'], N=512)
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']    # categorical slots 1-4, never cycled
LINESTYLES = ['-', (0, (5, 2)), (0, (1, 1.3)), (0, (6, 1.6, 1, 1.6))]
MARKERS = ['o', 's', '^', 'D']
GREY = '#6b6b6b'

RC = {
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'text.color': 'black', 'axes.labelcolor': 'black', 'axes.edgecolor': 'black',
    'xtick.color': 'black', 'ytick.color': 'black',
    'text.usetex': True, 'font.family': 'serif',
    'text.latex.preamble': r'\usepackage{amsmath}',
    'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9, 'legend.fontsize': 8,
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.linewidth': 0.6, 'lines.linewidth': 1.2, 'lines.markersize': 4,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 3, 'ytick.major.size': 3, 'xtick.minor.size': 1.6, 'ytick.minor.size': 1.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
    'legend.frameon': False, 'legend.handlelength': 2.6, 'axes.unicode_minus': False,
    'figure.dpi': 150, 'savefig.dpi': 600, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
    'pdf.fonttype': 42, 'image.cmap': 'viridis',
}


def apply():
    """(Re)apply the paper style: tools/contours.py, imported by data.py, sets its own rcParams."""
    plt.rcParams.update(RC)


apply()


def iso_colours(levels, norm, cmap=RAINBOW, shade=0.62):
    """Isoline colours: the colour map at each level, darkened."""
    return [tuple(shade * np.array(cmap(norm(v))[:3])) for v in levels]


def field(ax, X, Y, V, norm, levels=None, cmap=RAINBOW, lw=0.3):
    """Smooth (Gouraud) fill of a structured field plus thin isolines at `levels`."""
    V = np.ma.masked_invalid(V)
    m = ax.pcolormesh(X, Y, V, cmap=cmap, norm=norm, shading='gouraud', rasterized=True)
    if levels is not None:
        lv = [v for v in levels if V.min() < v < V.max()]
        if lv:
            ax.contour(X, Y, V, levels=lv, colors=iso_colours(lv, norm, cmap), linewidths=lw)
    return m


def colorbar(fig, mappable, cax, label, ticks=None, orientation='horizontal', fmt=None, minor=None):
    cb = fig.colorbar(mappable, cax=cax, orientation=orientation, ticks=ticks)
    cb.set_label(label)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(direction='in', length=2.6, width=0.5, which='major')
    cb.ax.tick_params(direction='in', length=1.4, width=0.4, which='minor')
    if fmt is not None:
        cb.ax.xaxis.set_major_formatter(fmt) if orientation == 'horizontal' else cb.ax.yaxis.set_major_formatter(fmt)
    if minor is not None:
        (cb.ax.xaxis if orientation == 'horizontal' else cb.ax.yaxis).set_minor_locator(minor)
    return cb


def channel(ax, h=50.0, faces=True, window=True, xlabel=True, ylabel=r'$y$ (mm)'):
    """Axes of a side view of the channel: x 0-275 mm, porous faces dashed, window marked at x = 0."""
    ax.set_xlim(0, 275); ax.set_ylim(0, h); ax.set_aspect('equal')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(25))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(5))
    if faces:
        for x in (75, 125):
            ax.axvline(x, color='k', lw=0.45, ls=(0, (3, 2)))
    if window:
        ax.plot([0, 0], [12.5, 37.5], color='k', lw=2.2, solid_capstyle='butt', clip_on=False, zorder=5)
    if xlabel:
        ax.set_xlabel(r'$x$ (mm)')
    ax.set_ylabel(ylabel)


def int_formatter():
    return ticker.FuncFormatter(lambda v, _: f'{v:.0f}')


def sci_formatter():
    """10^n labels for log colour bars."""
    return ticker.FuncFormatter(lambda v, _: r'$10^{%d}$' % round(np.log10(v)) if v > 0 else '')


def nice_levels(lo, hi, step):
    return np.arange(np.ceil(lo / step) * step, hi + 0.5 * step, step)


def save(fig, name, bbox=None):
    """PDF (tight bounding box, or the given one in inches) and a 300 dpi PNG preview."""
    os.makedirs(OUT, exist_ok=True)
    pdf = os.path.join(OUT, name + '.pdf')
    fig.savefig(pdf, **({'bbox_inches': bbox} if bbox is not None else {}))
    plt.close(fig)
    subprocess.run(['pdftoppm', '-r', '300', '-png', '-singlefile', pdf, os.path.join(OUT, name)], check=True)
    return pdf


def pending(ax, text):
    """Placeholder for a panel whose case has not finished yet."""
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha='center', va='center', color=GREY)
    ax.set_facecolor('#f2f2f2')
