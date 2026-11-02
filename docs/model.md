# The model in detail

## Modelling choices and the problem each one fixes

The case grew out of a four-region OpenFOAM set-up (upstream air, a 0.5 mm interface slab, porous
zone, downstream air, coupled through mapped patches). Each row is a defect of that set-up and the
fix adopted here.

| Problem | Cause | Fix |
|---|---|---|
| Mass not conserved between regions (0.323 / 0.245 / 0.328 g/s at inlet / porous zone / outlet) | `mappedMixed` U/p between regions does not exchange mass; a mass-conserving partitioned exchange diverged | The flow does not depend on T (constant properties, g = 0), so it is solved **once, in one region** over the whole channel, then split into the three regions and frozen (`frozenFlow true`). `restorePhi.py` rebuilds the join fluxes that `splitMeshRegions` zeroes: every section carries exactly rho U_in A, face by face |
| About 90 % of the beam lost | The beam left the slab through a mapped ray patch, which books no radiative flux | No slab. Both sides of the front (x = 75 mm) and back (x = 125 mm) faces are `greyDiffusiveRadiation`, eps = 0.92: the front face absorbs 92 % of the beam, reflects 8 % and emits into both zones |
| No radiation in the face temperature | `turbulentTemperatureCoupledBaffleMixed` ignores `qr` | `system/radiativeInterfaceT` (codedMixed) on the four join patches solves the face balance k1(Tc1 - Tf) + k2(Tc2 - Tf) + Q_abs - E Tf^4 = 0. Only the incident flux is explicit; emission is implicit (Newton, then the tangent as the mixed condition). Q_abs and E use fvDOM's own ordinate weights, so the face books exactly what fvDOM books in `qr` |
| Air cells next to the front face clipped by the temperature limiter | Air leaves the upstream region through the hot face; with 1.4 mm cells the cell Peclet number was 7-28 and the last cell undershot | Upstream region graded to 0.081 mm at the face (Peclet 0.42 at Re 300, 1.3 at Re 900); this also resolves the real preheating layer (alpha/U about 0.2 mm) |
| Wrong porous-zone settings | Air conductivity in the zone, resistance box 80-120 mm | kappa_con = 24.5 W/mK; d = 3.2468e7 1/m2, f = 1.2536e4 1/m (from K and c_F) on the whole 75-125 mm zone |
| Pressure drop 7.26 Pa instead of 7.5 Pa | The cell-centre velocity wiggles in the 2-3 cells next to a resistance jump (the face fluxes are exact) and the last porous cell was 4.9 mm | Porous zone graded towards both faces (0.61 mm): Delta p = 7.50 / 23.7 / 48.6 Pa at Re 300 / 600 / 900 |
| Negative ray intensities next to the porous faces, a cold strip on one inlet wall at the front face | With an even nTheta the beam ordinate is at theta = 90 +- dTheta/2, so the aligned channel is tilted from the fvDOM polar axis. `greyDiffusiveRadiation` decides from the ordinate direction d whether it leaves a face (fixed emitted value), but fvDOM transports it with dAve. Where the two disagree the emitted value flows out through the face, a sink in every cell along it | **Odd nTheta**: the beam ordinate lies on the equator (d parallel to dAve), the rotation is about the vertical only, and every boundary normal is horizontal or vertical, for which d and dAve always agree in sign. `tools/beamAlign.py` reports such sink rays and `Allrun.split` warns |

## Boundary conditions

| Boundary | Temperature | Radiation |
|---|---|---|
| Inlet plane (inlet and window) | 300 K | grey, eps = 0.92 at 300 K; the window adds the beam G to the beam ordinate (`qRadExt`, `qRadExtDir`) |
| Channel walls | `radiativeWallT`: k (Tc - Tw)/delta + q_abs - c eps sigma Tw^4 = h (Tw - T_inf) + eps sigma (Tw^4 - T_inf^4), h = 5 W/m2K, eps = 0.92, T_inf = 300 K (absorbed radiation, emission into the channel with fvDOM's weights c, loss to the surroundings) | grey, eps = 0.92 |
| Front and back faces of the porous zone (both sides) | `radiativeInterfaceT` (face balance above) | grey, eps = 0.92, opaque |
| Outlet | zero gradient | black at 300 K (open to the surroundings) |

`limitTemperature` (200-5000 K) is only a crash guard: a run in which it clips any cell is invalid.

## Mesh

One block mesh, 209,088 cells, 44 x 44 across (12 cells in each 12.5 mm band next to a wall, 0.63
mm at the wall; 20 uniform cells over the 25 mm window band).

| Region | x (mm) | Axial cells | Spacing |
|---|---|---|---|
| inletRegion | 0-75 | 48 | 3.5 mm at the inlet to 1.0 mm at 70 mm, then to 0.081 mm at the face |
| porousRegion | 75-125 | 30 | 0.61 mm at both faces, 3.45 mm in the middle |
| exitRegion | 125-275 | 30 | 2.3 mm to 9.2 mm |

`python3 case/scripts/genBlockMeshDict.py case/flow/system/blockMeshDict` writes the mesh dictionary;
`case/scripts/design/meshcalc.py` and `outflow_face_1d.py` document the grading and the Peclet
criterion.

## Solution procedure

1. `Allrun.flow <Uin>`: blockMesh, single-region flow with Darcy-Forchheimer on the porous cellZone,
   then a polish with tight pressure tolerances so that the frozen fluxes are divergence free.
2. `Allrun.split <flowCase> <nPhi> <nTheta>`: `beamAlign.py` gives the beam ordinate and the rotation
   that maps +x onto its dAve; `transformPoints -rotateFields` rotates mesh and flow;
   `splitMeshRegions -cellZonesOnly` makes the three regions; `restorePhi.py` restores the join
   fluxes and checks mass conservation face by face.
3. `Allrun.cht <splitCase> <G> <alphaEff>`: energy and fvDOM in the three regions with the flow
   frozen; radiation solved every 20 iterations; under-relaxation of h 0.7 (air), 0.9 (porous zone).
   `WARM=<case>/<time>` starts from the temperature of another run on the same mesh.
4. `Allpost`: an un-rotated copy (`post/`) for post-processing.

With 0.7 in the porous zone, its temperature level (a nearly Neumann problem with high
conductivity) is the slowest mode; at G = 1e6 / Re = 300, where the radiative coupling (~T^3) is
weak, the global residual fell about 25 % per 500 iterations with 0.7 and about 50 % with 0.9. The
converged solution does not depend on it.

## Checks

See the README. The energy books (`system/balanceFunctions`, `system/auditFunctions`,
`energyBalance.py`) close to better than 0.1 % of the beam in every case. The in-run audit (the
`energyAudit global` line, every 100 iterations) books every term and is what `stopWhenConverged.sh`
reads; `energyBalance.py` books the air by its convective enthalpy only, so its closure also holds the
0.3-1.3 W conducted from the air into the 300 K inlet plane.

## Bulk temperature along the channel

`tools/bulkFaces.py` gives the bulk temperature on every face plane from the solver's face mass
fluxes and the face values of its convection scheme (linearUpwind with the cell-limited gradient;
the boundary values on the inlet plane, the porous faces and the outlet), so that every plane is the
convective enthalpy flux of the solver and the outlet value is that of the energy books; Figs 12-13
use it. The cell-layer values of `tools/profiles.py` weight the cells with the
cell-centre velocity and put the zero-gradient outlet value at the centre of the 9.2 mm last cell,
4.6 mm upstream of the outlet, which shows as a kink in the last 10 mm.

Just downstream of the porous zone the bulk temperature drops by 3-15 K within the first exit
cell (2.3 mm), the most at G 1e7 / Re 300. This is in the solution, not in the post-processing: the
walls of the porous zone are held near the temperature of the medium by kappa_con, while the exit
walls exchange heat only with the air and by radiation, and their temperature drops by about 300 K
at x = 125 mm (G 5e6 / Re 300, top wall: 1173 K -> 876 K). The air then meets a much colder wall,
and the new thermal boundary layer takes about 4 W in the first cell (the air-side conduction into
the first exit wall cell is about 10 kW/m2, against about 3 kW/m2 15 mm further on). The walls have
no thickness and no axial conduction, so the wall temperature, and with it this drop, is abrupt.

## Angular discretisation (G = 5e6, Re = 300, same mesh and flow)

| | 6x10 | 12x10 | **6x11 (adopted)** | 8x15 |
|---|---|---|---|---|
| Ordinates | 240 | 480 | 264 | 480 |
| Cost per iteration, 14 cores | 0.31 s | 1.08 s | 0.38 s | 1.12 s |
| Beam ordinate theta / axis elevation out of the horizontal | 81 / 9 deg | 81 / 9 deg | 90 / 0 deg | 90 / 0 deg |
| Ordinates with negative intensity (minimum) | 8 (-7,982 kW/m2/sr) | 18 (-9,036 kW/m2/sr) | 0 | 0 |
| Outlet bulk T | 993.1 K | 993.3 K | 993.6 K | 993.2 K |
| Front / back face maximum T | 2291 / 1217 K | 2290 / 1218 K | 2292 / 1218 K | 2291 / 1218 K |
| Inlet-wall T, departure from the square symmetry (rms) | 11.4 K | 8.4 K | 8.9 K | 4.9 K |

- Bulk, face and porous-zone temperatures do not depend on the ray set (within 2 K), nor do the heat
  flows (within 2.5 W, 0.08 % of the beam). 6x11 and 8x15 are converged to the 0.02 % criterion; the
  even sets end at an in-run audit imbalance of 0.06 % (6x10, 8000 iterations) and 0.03 % (12x10, 4000),
  which leaves up to about 0.35 K on their outlet temperatures.
- Wall temperatures carry a ray effect, a patchy redistribution of the same absorbed heat: about
  8-17 K rms on the inlet walls, larger within 5 mm of the front face, 3-5 K on the exit walls,
  < 0.5 K on the porous walls. Use 8x15 when local wall temperatures next to the front face matter.
- The case is symmetric under the 8 reflections and quarter turns of the square cross-section; the
  ordinate set is not, so the departure from that symmetry measures the angular error
  (`tools/symmetryCheck.py`). With an odd nTheta both mirrors are also symmetries of the ordinate
  set, and only the y <-> z part remains.

## Mesh sensitivity (G = 5e6, Re = 300, ray set 6x11)

`python3 case/scripts/genBlockMeshDict.py <blockMeshDict> <scale>` with scale 0.75 and 1.25 writes the
two other meshes (every cell count scaled, same gradings). Each was run from 300 K to the 0.02 %
criterion. The costs of the coarse and fine meshes were timed back to back with the adopted one on
14 ranks and scaled to its cost on an otherwise idle machine.

| | 33 x 33 x 80 | **44 x 44 x 108 (adopted)** | 55 x 55 x 136 |
|---|---|---|---|
| Cells | 87,120 | 209,088 | 411,400 |
| Upstream cell next to the front face | 0.104 mm | 0.081 mm | 0.066 mm |
| Outlet bulk T | 993.22 K | 993.63 K | 993.93 K |
| Front-face maximum T | 2293.1 K | 2292.0 K | 2292.4 K |
| Peak bulk T (at the front face) | 1631.4 K | 1631.9 K | 1632.1 K |
| Top-wall maximum T | 1428.7 K | 1428.7 K | 1428.9 K |
| Energy closure (energyBalance.py) | +0.022 % | +0.027 % | +0.035 % |
| Cost per iteration, 14 cores | 0.22 s | 0.38 s | 0.86 s |

The adopted mesh is within 0.3 K of the fine one at the outlet and within 0.5 K elsewhere. The
outlet temperature rises by 0.41 and 0.30 K over the two refinements: an apparent order of 0.27,
too low for a meaningful Richardson extrapolation; the differences are close to the level of the
iteration and ray-effect errors.

## Sensitivity to the porous-zone conductivity and to the inlet plane (G = 5e6, Re = 300)

Both runs start from the adopted solution and run to the 0.02 % criterion.

| | adopted | air's conductivity in the porous zone | inlet plane black (eps 1) |
|---|---|---|---|
| Outlet bulk T | 993.6 K | 915.7 K | 993.4 K |
| Air enthalpy gain | 228.2 W | 202.6 W | 228.1 W |
| Wall losses: inlet / porous / exit walls | 913 / 1459 / 263 W | 1373 / 977 / 191 W | 898 / 1457 / 262 W |
| Radiation absorbed by the inlet plane | 243.9 W | 369.2 W | 261.4 W |
| Front-face temperature | 1281-2292 K | 1025-2664 K | 1281-2291 K |
| Back-face temperature | 1133-1218 K | 877-1146 K | 1132-1218 K |

- kappa_con = 24.5 W/mK carries the heat absorbed at the front face into the porous zone and to its
  walls. With the conductivity of air (0.026 W/mK) the face runs about 370 K hotter at the beam
  centre and radiates more back to the inlet plane and the upstream walls: the air gains 11 % less,
  the outlet is 78 K colder, and the back face is three times less uniform.
- The inlet plane, held at the inlet temperature, absorbs 92 % of the radiation that reaches it
  (about 265 W, 8.5 % of the beam), which leaves the channel as through the transparent window of
  the pore-scale model, and reflects 8 %. A black inlet plane, which reflects nothing, changes the
  outlet temperature by 0.3 K.

## OpenFOAM-v2206 notes

- `foamDictionary -set` on `controlDict` writes it back with every `#include` expanded; later edits
  to the included files are then ignored. Edit `controlDict` with `sed`.
- If an edit to a coded function object has no effect, check that `controlDict` still `#include`s
  the file (see the note above); deleting `dynamicCode/<name>` forces a rebuild.
- A top-level `constant/g` is required even though every region has its own.
- fvDOM: `greyDiffusiveRadiationMixed` classifies an ordinate at a boundary by its direction d
  (`(-n & d) > 0`: emitted, fixed value) while transporting it with `dAve & Sf`. Where the two
  disagree in sign the boundary is a sink and intensities go negative (fvDOM does not clip them).
- `reconstructPar` writes wrong `valueFraction`/`refValue` entries for the ray fields `ILambda_*` on
  one wall patch; internal fields and `value` entries are correct. Read boundary-condition
  coefficients from `processor*/<time>`.
- In double precision cos(90 deg) = 6e-17: ordinates exactly tangent to a wall take different
  boundary branches on opposite walls. They carry no flux, but the boundary value enters the limited
  gradient, which leaves mirror differences of 1-3 K in corner cells.
