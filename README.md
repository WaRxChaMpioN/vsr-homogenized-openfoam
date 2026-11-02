# vsr-homogenized-openfoam

OpenFOAM case of the **homogenized (macroscopic) model of a volumetric solar receiver** with a TPMS
porous absorber under concentrated, collimated irradiation, used in

> K. Jha, S. S. Rathore, B. Mehta, P. Kumar, *Effective property modelling of radiative heat transfer
> in porous media for concentrated solar energy applications* (submitted, 2026).

The porous absorber is a pseudo-porous zone with Darcy–Forchheimer resistance, an effective
conductivity and an effective absorption coefficient extracted from pore-scale simulations
(K. Jha et al., *Role of radiation in a TPMS porous absorber under concentrated irradiation: a
pore-scale study and effective properties for homogenized models*, under review). Radiation is
solved with the finite-volume discrete ordinates method (fvDOM) of OpenFOAM.

The repository contains the case set-up, the run scripts, the analysis tools and the scripts that
make the paper's figures and tables. It does not contain simulation results (they are regenerated
by the scripts, see *Running*).

## Contents

| Folder | |
|---|---|
| `case/flow` | stage 1: single-region flow case (template) |
| `case/cht` | stage 3: three-region heat-transfer case (template): coded temperature conditions `system/radiativeInterfaceT` (opaque porous faces) and `system/radiativeWallT` (channel walls), energy books `system/balanceFunctions` and `system/auditFunctions`, parameters `system/caseParameters` |
| `case/scripts` | `Allrun.flow`, `Allrun.split`, `Allrun.cht`, `Allpost`, `Allrun.matrix` (all nine G, Re cases), `stopWhenConverged.sh`, `restorePhi.py`, `energyBalance.py`, `genBlockMeshDict.py`, `genFunctions.py`, `rotPoints.py`; `design/`: derivation checks of the boundary-condition maths and the mesh grading |
| `tools` | post-processing and checks (Python) |
| `figures` | the paper's figures and Tables 3-4 from the runs; `data/` holds the pore-scale bulk temperatures digitised from the pore-scale paper (Fig 13) |
| `docs/model.md` | the model in detail: what each modelling choice fixes, the checks, the angular-discretisation study, OpenFOAM-v2206 notes |

## The receiver

| | |
|---|---|
| Channel | 275 mm long, 50 mm x 50 mm cross-section, air in at 300 K |
| Window | 25 mm x 25 mm, centred in the inlet plane; the collimated beam G enters through it along the axis |
| Pseudo-porous zone | 50 mm cube from x = 75 to 125 mm (porosity 0.3, TPMS diamond lattice of alumina) |
| Irradiation G | 1e6, 5e6, 1e7 W/m2 (beam power 625, 3125, 6250 W) |
| Reynolds number | 300, 600, 900 (U_in = 0.1068, 0.2136, 0.3204 m/s) |
| Walls | zero thickness; lose heat to 300 K surroundings by convection (h = 5 W/m2K) and radiation (eps = 0.92) |
| Air | rho 1.225 kg/m3, c_p 1006 J/kgK, k 0.026 W/mK, mu 1.8e-5 kg/ms (constant) |
| Porous zone | kappa_con 24.5 W/mK; K 3.08e-8 m2 and c_F 1.10 (Darcy-Forchheimer d = 3.2468e7 1/m2, f = 1.2536e4 1/m); alpha_eff = 45, 60, 66 1/m for G = 1e6, 5e6, 1e7; no scattering |
| Surfaces | grey, eps = 0.92: walls, both faces of the porous zone, inlet plane (window included); the outlet is open to 300 K surroundings (black) |

## The model in short

- **Frozen flow.** With constant properties and g = 0 the flow does not depend on temperature. It is
  solved once per Re in one region over the whole channel (`chtMultiRegionSimpleFoam`, Darcy-Forchheimer
  on the porous cellZone), then split into the three regions and frozen (`frozenFlow true`). Every
  section carries exactly the inlet mass flow, face by face (`restorePhi.py` rebuilds the fluxes on
  the joins after `splitMeshRegions`).
- **Beam alignment.** The beam is one discrete ordinate. The mesh and the flow are rotated
  (`transformPoints -rotateFields`) so that the channel axis is parallel to that ordinate's transport
  direction dAve (`tools/beamAlign.py`). Use an **odd nTheta** (6x11 here): the beam ordinate then
  lies on the equator and the channel is rotated about the vertical only, so every boundary normal is
  horizontal or vertical and no ordinate is emitted by a boundary while being transported out through
  it. With an even nTheta that inconsistency of fvDOM's boundary condition makes negative intensities
  next to the porous faces (`docs/model.md`; `beamAlign.py` reports these "sink rays" and
  `Allrun.split` warns).
- **Opaque porous faces.** Both faces of the porous zone are grey, opaque surfaces seen from both
  sides: the front face absorbs 92 % of the beam and reflects 8 %. Their temperature
  (`radiativeInterfaceT`) satisfies the face balance k1(Tc1 - Tf) + k2(Tc2 - Tf) + Q_abs - E Tf^4 = 0
  with the absorbed radiation and emission weights taken from fvDOM's own ordinates, emission
  implicit (Newton).
- **Walls** (`radiativeWallT`): k (Tc - Tw)/delta + q_abs = h (Tw - T_inf) + eps sigma (Tw^4 - T_inf^4),
  emission implicit.
- **Mesh.** One block mesh, 209,088 cells: 44 x 44 across (0.63 mm at the walls), axially 48 cells
  in the inlet section (graded to 0.081 mm at the front face, so that the air leaving through the hot
  face keeps a cell Peclet number below 2), 30 in the porous zone (0.61 mm at both faces) and 30 in
  the exit section.
- **Stage 3** (`Allrun.cht`) solves only the energy equation and fvDOM (264 ordinates), with
  under-relaxation 0.7 (air) and 0.9 (porous zone) of h.

## Requirements

- **OpenFOAM-v2206** (openfoam.com) with Scotch and MPI. The coded boundary conditions are compiled at
  run time (`wmake` must work). Developed with the FreeMHD build of v2206; other versions may need
  changes to the coded boundary conditions and function objects.
- **Python 3.8+** with numpy and scipy (tools); matplotlib for the figures.
- For the figures: a LaTeX installation (text is typeset with `text.usetex`) and `pdftoppm` (poppler)
  for PNG previews.

## Running

Source OpenFOAM and set `VSR` to this repository, then work in a run directory outside it (or in
`runs/`, which git ignores):

```sh
source <OpenFOAM-v2206>/etc/bashrc
export VSR=/path/to/vsr-homogenized-openfoam
mkdir -p $VSR/runs && cd $VSR/runs
```

**All nine cases** (flows, splits, G 5e6 / Re 300 from 300 K, then the other eight in two lanes of
10 ranks, each started from its nearest finished neighbour):

```sh
$VSR/case/scripts/Allrun.matrix > log.matrix 2>&1 &
$VSR/case/scripts/stopWhenConverged.sh > log.stopper 2>&1 &   # end each run once converged
```

**One case by hand:**

```sh
# 1. flow, once per Re (Uin 0.1068 | 0.2136 | 0.3204 for Re 300 | 600 | 900)
cp -r $VSR/case/flow flow_Re300 && cp $VSR/case/scripts/Allrun.flow flow_Re300/
(cd flow_Re300 && ./Allrun.flow 0.1068)

# 2. beam alignment and split, once per Re and ray set; must print ALL CHECKS PASSED
mkdir split_Re300_6x11 && cp $VSR/case/scripts/Allrun.split split_Re300_6x11/
(cd split_Re300_6x11 && ./Allrun.split ../flow_Re300 6 11)

# 3. heat transfer, once per (G, Re): G and alpha_eff = 1e6 45 | 5e6 60 | 1e7 66
cp -r $VSR/case/cht G5e6_Re300_6x11 && cp $VSR/case/scripts/Allrun.cht $VSR/case/scripts/Allpost G5e6_Re300_6x11/
(cd G5e6_Re300_6x11 && ./Allrun.cht ../split_Re300_6x11 5e6 60)
#    WARM=<case>/<time> ./Allrun.cht ...   start from the T field of another run on the same mesh
#    SMOKE=1 ./Allrun.cht ...              throw-away test: 5 iterations

# 4. un-rotated copy for post-processing, axial profiles, checks
(cd G5e6_Re300_6x11 && ./Allpost)
python3 $VSR/tools/profiles.py G5e6_Re300_6x11/post <time> profiles.csv
python3 $VSR/tools/negativeRays.py G5e6_Re300_6x11 <time>
```

Set the number of ranks in `system/decomposeParDict` (and `system/<region>/decomposeParDict` in the
heat-transfer case). Wall-clock times on 14 cores (Intel i7-14700): flow 4 min per Re; split under
a minute; heat transfer 50 min from 300 K, 25-45 min from a neighbouring solution; the whole matrix
about 4 h.

## Checks

Every 100 iterations the `energyAudit` function object prints the energy and mass books
(`grep '^energyAudit' log.chtMultiRegionSimpleFoam`); `energyBalance.py` prints them after the run.

| Check | Pass |
|---|---|
| `restorePhi.py` (in `Allrun.split`) | `ALL CHECKS PASSED`: join fluxes equal and opposite face by face, every section = inlet mass flow |
| `energyAudit join ... dMass` | 0 |
| `energyAudit join ... faceBalance` (radiation absorbed by a face = conduction into both sides) | < 1 W |
| `energyAudit join ... dTf` (face temperature seen from both sides) | < 0.1 K |
| `energyAudit wall ... residual` | < 0.1 W |
| `energyAudit global ... residual` (beam = escape + air enthalpy rise + wall loss) | < 0.5 % of the beam (`stopWhenConverged.sh` stops at 0.02 %) |
| `limitTemperature` clipping (`grep 'limited [1-9]'`) | never: any clipping invalidates the energy balance |
| `Solving for Ux` / `p_rgh` in the heat-transfer log | never (flow frozen) |
| `tools/negativeRays.py` | no negative ray intensity |
| `tools/symmetryCheck.py` | y- and z-mirror parts < 0.2 K rms (odd nTheta) |

## Figures and tables of the paper

```sh
export VSR_RUNS=$VSR/runs          # the run directory of Allrun.matrix
cd $VSR/figures
python3 make_figures.py all        # figures/out/figNN_*.pdf/.png, the paper's file names
python3 tables.py                  # Table 4 (and Table 3 if the angular-study runs exist)
python3 make_proof.py              # figures/out/figures_proof.pdf: all figures and tables at the paper's sizes
```

Fig 13 compares with the pore-scale bulk temperatures, digitised from the pore-scale paper's Fig 13
(`figures/data/pore_scale_bulkT_digitised.csv`: outlet values within 0.8 K of that paper's table,
about +-1 K on the gentle parts and +-5-15 K on the steep rise at the front face).

## Tools

| | |
|---|---|
| `beamAlign.py` | fvDOM ordinate set, beam ordinate, rotation for `transformPoints`, boundary sink rays |
| `angularCheck.py` | hemispherical weight error of fvDOM's diffuse emission on each boundary orientation of the rotated channel, for any (nPhi, nTheta) |
| `profiles.py` | axial profiles per cell layer: bulk, mean, min, max and centreline temperature, wall temperature and absorbed radiation |
| `negativeRays.py` | minimum of every ray intensity field |
| `symmetryCheck.py` | angular error from the departure from the square symmetry of the case |
| `fieldDiff.py` | cell-by-cell and face-by-face difference of two runs on the same mesh |
| `paperTargets.py` | G 5e6 / Re 300 run against the values reported in the papers |
| `contours.py`, `foamio.py` | plotting and a minimal ASCII OpenFOAM reader used by the others |

## Citation

If you use this case or its tools, please cite the paper and this repository (`CITATION.cff`; GitHub
shows a *Cite this repository* button):

```bibtex
@article{Jha2026homogenized,
  author  = {Jha, Kaushal and Rathore, Surendra Singh and Mehta, Balkrishna and Kumar, Pradeep},
  title   = {Effective property modelling of radiative heat transfer in porous media for concentrated
             solar energy applications},
  year    = {2026},
  note    = {Submitted}
}

@article{Jha2026porescale,
  author  = {Jha, Kaushal and Rathore, Surendra Singh and Mehta, Balkrishna and Kumar, Pradeep},
  title   = {Role of radiation in a {TPMS} porous absorber under concentrated irradiation: a pore-scale
             study and effective properties for homogenized models},
  year    = {2026},
  note    = {Under review}
}

@software{Jha2026vsrcase,
  author  = {Jha, Kaushal},
  title   = {vsr-homogenized-openfoam: {OpenFOAM} case of a homogenized volumetric solar receiver},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/WaRxChaMpioN/vsr-homogenized-openfoam}
}
```

## License

GNU General Public License v3.0 or later (`LICENSE`), as OpenFOAM.

## Contact

Kaushal Jha, kaushal_jha@mines.edu
