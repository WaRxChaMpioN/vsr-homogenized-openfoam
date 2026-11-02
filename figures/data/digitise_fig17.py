#!/usr/bin/env python3
"""
Wall temperature along the channel of the pore-scale paper (its Fig 17: line-averaged wall
temperature against the distance from the absorber centre, cm; three G groups, one colour per Re),
read from the figure image, for the wall-temperature comparison with the macroscopic model:

    digitise_fig17.py <fig17_wall_temperature_profile.png> <out.csv> [overlay.png]

Axes are calibrated on the tick marks of the frame. Curve pixels are classified by hue (blue
Re 300, red Re 600, orange Re 900). For each colour, the pixel columns that show its three curves
(one per G) as separate runs give the curves first (outliers against the local median dropped);
then, in every column, the remaining runs are matched to the nearest of these curves, which adds
the columns where one curve is hidden under another colour. A running median removes the outlines
of the markers, and the curves are interpolated linearly in between (the figure draws straight
segments between its stations). Output: G, Re, x_mm (from the inlet = 100 mm
+ 10 x cm), T_K, every 0.5 mm from the inlet to the end of the range read. The overlay redraws the read curves on the image as a check.
"""
import sys
import numpy as np
from PIL import Image

COLOURS = {300: (48, 48, 255), 600: (211, 3, 3), 900: (255, 171, 87)}
G_ORDER = ('1e6', '5e6', '1e7')                     # low to high temperature
LEGEND = (slice(225, 370), slice(790, 1000))         # rows, columns of the three legend lines


def calibrate(im):
    """Pixel -> (x cm, T K) from the frame and its tick marks."""
    dark = im.sum(2) < 200
    colsum, rowsum = dark.sum(0), dark.sum(1)
    left = 150 + int(np.argmax(colsum[150:300])); right = 1100 + int(np.argmax(colsum[1100:1250]))
    top = 120 + int(np.argmax(rowsum[120:300])); bottom = 1000 + int(np.argmax(rowsum[1000:1150]))
    yt = [y for y in range(top + 2, bottom - 2) if dark[y, left + 2:left + 12].mean() > 0.8]
    yt = [g.mean() for g in np.split(np.array(yt), np.where(np.diff(yt) > 2)[0] + 1)]
    xt = [x for x in range(left + 2, right - 2) if dark[bottom - 12:bottom - 2, x].mean() > 0.8]
    xt = [g.mean() for g in np.split(np.array(xt), np.where(np.diff(xt) > 2)[0] + 1)]
    # y ticks: 1600, 1400, ..., 400 K from the top; x ticks: -5, 0, 5, 10, 15 cm (-10 and 17.5 on the frame)
    ky, cy = np.polyfit(yt, np.linspace(1600, 400, len(yt)), 1)
    kx, cx = np.polyfit(xt, np.arange(-5, -5 + 5 * len(xt), 5), 1)
    return (left, right, top, bottom), (lambda px: kx * px + cx), (lambda py: ky * py + cy)


def runs_of(mask, c, top, bottom):
    ys = np.where(mask[top + 2:bottom - 2, c])[0] + top + 2
    return [r.mean() for r in np.split(ys, np.where(np.diff(ys) > 3)[0] + 1) if len(r) >= 2] if len(ys) else []


def masks(im):
    R, G, B = im[..., 0], im[..., 1], im[..., 2]
    sat = im.max(2) - im.min(2) > 60
    out = {300: sat & (B - np.maximum(R, G) > 60),
           600: sat & (R - G > 80) & (np.abs(G - B) < 35),
           900: sat & (R - B > 80) & (G - B > 35)}
    for m in out.values():
        m[LEGEND] = False
    return out


def curves(im, frame, xof, Tof):
    left, right, top, bottom = frame
    out = {}
    cols = np.arange(left + 2, right - 1)
    for Re, mask in masks(im).items():
        rr = {c: runs_of(mask, c, top, bottom) for c in cols}
        three = np.array([c for c in cols if len(rr[c]) == 3])
        Y = np.array([sorted(rr[c], reverse=True) for c in three])          # pixel y down: lowest T first
        pts = []
        for j in range(3):                                                 # pass 1, with outlier rejection
            y = Y[:, j]
            med = np.array([np.median(y[max(0, i - 10):i + 11]) for i in range(len(y))])
            ok = np.abs(y - med) < 12
            pts.append(dict(zip(three[ok], y[ok])))
        for j in range(3):                                                 # pass 2: nearest run to each curve
            cs = np.array(sorted(pts[j])); ys = np.array([pts[j][c] for c in cs])
            others = [(np.array(sorted(pts[k])), np.array([pts[k][c] for c in sorted(pts[k])])) for k in range(3) if k != j]
            for c in cols:
                if c in pts[j] or not rr[c]:
                    continue
                pred = np.interp(c, cs, ys)
                cand = [v for v in rr[c] if abs(v - pred) < 10
                        and all(abs(v - np.interp(c, oc, oy)) > abs(v - pred) for oc, oy in others)]
                if cand:
                    pts[j][c] = min(cand, key=lambda v: abs(v - pred))
        for j, G in enumerate(G_ORDER):
            cs = np.array(sorted(pts[j])); yv = np.array([pts[j][c] for c in cs])
            ym = np.array([np.median(yv[max(0, i - 4):i + 5]) for i in range(len(yv))])
            out[(G, Re)] = (xof(cs), Tof(ym))
    return out


def main():
    img, csv = sys.argv[1], sys.argv[2]
    im = np.asarray(Image.open(img).convert('RGB')).astype(int)
    frame, xof, Tof = calibrate(im)
    data = curves(im, frame, xof, Tof)
    xg = np.arange(0.0, 275.01, 0.5)                   # mm from the inlet
    with open(csv, 'w') as f:
        f.write('G,Re,x_mm,T_K\n')
        for (G, Re), (xc, T) in sorted(data.items(), key=lambda kv: (float(kv[0][0]), kv[0][1])):
            xm = 100 + 10 * xc
            xr = xg[xg <= xm.max()]
            Tr = np.interp(xr, xm, T)
            # first segment (inlet to the station at 25 mm): the nine lines start together at the inlet
            # edge, where they overlap; the figure draws it straight, so it is the line through the
            # points read at 12-24 mm, where the lines are apart, carried to the inlet
            seg = (xm >= 12) & (xm <= 24)
            a, b = np.polyfit(xm[seg], T[seg], 1)
            Tr[xr < 25] = a * xr[xr < 25] + b
            for x, v in zip(xr, Tr):
                f.write(f'{G},{Re},{x:.1f},{v:.1f}\n')
            print(f'G {G} Re {Re}: {len(xc)} columns read, x {100 + 10 * xc.min():.0f}-{100 + 10 * xc.max():.0f} mm')
    if len(sys.argv) > 3:                              # overlay for checking
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        left, right, top, bottom = frame
        kx = (xof(right) - xof(left)) / (right - left); ky = (Tof(bottom) - Tof(top)) / (bottom - top)
        fig, ax = plt.subplots(figsize=(im.shape[1] / 100, im.shape[0] / 100), dpi=100)
        ax.imshow(im.astype(np.uint8)); ax.set_axis_off()
        for (G, Re), (xc, T) in data.items():
            Tg = np.interp(xg, 100 + 10 * xc, T)
            px = left + ((xg - 100) / 10 - xof(left)) / kx
            py = top + (Tg - Tof(top)) / ky
            ax.plot(px, py, color='k', lw=0.8, ls=(0, (2, 2)))
        fig.savefig(sys.argv[3], dpi=100, bbox_inches='tight')


if __name__ == '__main__':
    main()
