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
`energyBalance.py`) close to better than 0.1 % of the beam in every case.

## Angular discretisation (G = 5e6, Re = 300, same mesh and flow)

| | 6x10 | 12x10 | **6x11 (adopted)** | 8x15 |
|---|---|---|---|---|
| Ordinates | 240 | 480 | 264 | 480 |
| Cost per iteration, 14 cores | 0.31 s | 1.08 s | 0.38 s | 1.12 s |
| Beam ordinate theta / axis elevation out of the horizontal | 81 / 9 deg | 81 / 9 deg | 90 / 0 deg | 90 / 0 deg |
| Ordinates with negative intensity (minimum) | 8 (-7,982 kW/m2/sr) | 18 (-9,036 kW/m2/sr) | 0 | 0 |
| Outlet bulk T | 993.1 K | 993.3 K | 993.4 K | 993.1 K |
| Front / back face maximum T | 2291 / 1217 K | 2290 / 1218 K | 2292 / 1218 K | 2291 / 1218 K |
| Inlet-wall T, departure from the square symmetry (rms) | 11.4 K | 8.4 K | 8.9 K | 4.9 K |

- Bulk, face and porous-zone temperatures do not depend on the ray set (within 2 K), nor do the heat
  flows (within 2.5 W, 0.08 % of the beam).
- Wall temperatures carry a ray effect, a patchy redistribution of the same absorbed heat: about
  8-17 K rms on the inlet walls, larger within 5 mm of the front face, 3-5 K on the exit walls,
  < 0.5 K on the porous walls. Use 8x15 when local wall temperatures next to the front face matter.
- The case is symmetric under the 8 reflections and quarter turns of the square cross-section; the
  ordinate set is not, so the departure from that symmetry measures the angular error
  (`tools/symmetryCheck.py`). With an odd nTheta both mirrors are also symmetries of the ordinate
  set, and only the y <-> z part remains.

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
