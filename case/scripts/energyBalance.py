#!/usr/bin/env python3
"""Global energy and mass books of a CHT run, from the function objects of
system/balanceFunctions (latest values). Usage: energyBalance.py <case> [G] (default G 5e6)

  beam in        Q_beam = G * A_window (the beam ray's dAve is the channel axis, so its flux
                 through the window is exactly qRadExt)
  air gain       H_out + H_in            (weightedSum(h) by phi: outflow +, inflow -)
  wall loss      -sum wallHeatFlux       (alpha*snGrad(h) - qr: heat into the fluid minus the
                                          radiation absorbed; the wall has no storage)
  inlet plane    Qr_inlet + Qr_window + Q_beam   (radiation absorbed by the 300 K inlet plane)
  outlet         Qr_outlet               (radiation leaving to the 300 K surroundings)
  closure        Q_beam - (air + walls + inlet plane + outlet)  -> target |closure| < 0.5 % Q_beam
Per join: H_up + H_down = 0 (convective energy continuity); wallHeatFlux(up) + wallHeatFlux(down)
= 0 (face balance: conduction into both fluids = net radiation absorbed by the face)."""
import os, sys, glob
case = sys.argv[1]
G = float(sys.argv[2]) if len(sys.argv) > 2 else 5e6
pp = os.path.join(case, 'postProcessing')
def latest(region, fo, fname):
    dirs = sorted(glob.glob(os.path.join(pp, region, fo, '*')), key=lambda d: float(os.path.basename(d)))
    if not dirs:
        raise FileNotFoundError(f'{region}/{fo}')
    rows = [l.split() for l in open(os.path.join(dirs[-1], fname)) if l.strip() and not l.startswith('#')]
    return rows
def sfv(region, fo):
    r = latest(region, fo, 'surfaceFieldValue.dat')[-1]
    return float(r[0]), float(r[1])
def whf(region):
    rows = latest(region, f'wallHeatFlux_{region}', 'wallHeatFlux.dat')
    t = rows[-1][0]
    return {r[1]: float(r[4]) for r in rows if r[0] == t}
A_window = 0.025*0.025
Qbeam = G*A_window
t, Hin = sfv('inletRegion', 'H_inlet'); _, Hwin = sfv('inletRegion', 'H_opticalWindow'); _, Hout = sfv('exitRegion', 'H_outlet')
air = Hout + Hin + Hwin
_, Qi = sfv('inletRegion', 'Qr_inlet'); _, Qw = sfv('inletRegion', 'Qr_opticalWindow'); _, Qo = sfv('exitRegion', 'Qr_outlet')
inletPlane = Qi + Qw + Qbeam
W = {R: whf(R) for R in ('inletRegion', 'porousRegion', 'exitRegion')}
walls = -(W['inletRegion']['inletWall'] + W['porousRegion']['porousWall'] + W['exitRegion']['exitWall'])
closure = Qbeam - (air + walls + inletPlane + Qo)
print(f'time {t:g}')
print(f'  beam in                   {Qbeam:10.2f} W')
print(f'  air enthalpy gain         {air:10.2f} W   (H_out {Hout:.2f}, H_in {Hin + Hwin:.2f})')
print(f'  wall losses (external)    {walls:10.2f} W   inlet {-W["inletRegion"]["inletWall"]:.2f} / porous {-W["porousRegion"]["porousWall"]:.2f} / exit {-W["exitRegion"]["exitWall"]:.2f}')
print(f'  absorbed by inlet plane   {inletPlane:10.2f} W')
print(f'  lost through outlet       {Qo:10.2f} W')
print(f'  CLOSURE                   {closure:10.2f} W  = {100*closure/Qbeam:+.3f} % of the beam')
_, mi = sfv('inletRegion', 'mdot_inlet'); _, mw = sfv('inletRegion', 'mdot_opticalWindow')
min_ = -(mi + mw)
print(f'  mass in {min_:.8e} kg/s; bulk T out {sfv("exitRegion", "outletBulkT")[1]:.2f} K'
      f' (from H: {300 + air/(min_*1006):.2f} K)')
for up, R1, dn, R2 in [('inletRegion_to_porousRegion', 'inletRegion', 'porousRegion_to_inletRegion', 'porousRegion'),
                       ('porousRegion_to_exitRegion', 'porousRegion', 'exitRegion_to_porousRegion', 'exitRegion')]:
    _, m1 = sfv(R1, f'mdot_{up}'); _, m2 = sfv(R2, f'mdot_{dn}')
    _, h1 = sfv(R1, f'H_{up}'); _, h2 = sfv(R2, f'H_{dn}')
    _, q1 = sfv(R1, f'Qr_{up}'); _, q2 = sfv(R2, f'Qr_{dn}')
    f1, f2 = W[R1][up], W[R2][dn]
    print(f'  join {up}: mdot {m1:+.6e} / {m2:+.6e} (sum {m1+m2:+.1e}); '
          f'H {h1:+.3f} / {h2:+.3f} W (sum {h1+h2:+.3f}); qr absorbed {q1:+.2f} + {q2:+.2f} W; '
          f'face balance {f1 + f2:+.3f} W')
