# vsr-homogenized-openfoam

OpenFOAM case of the **homogenized (macroscopic) model of a volumetric solar receiver** with a TPMS
porous absorber under concentrated, collimated irradiation, used in

> K. Jha, S. S. Rathore, B. Mehta, P. Kumar, *Effective property modelling of radiative heat transfer
> in porous media for concentrated solar energy applications* (submitted, 2026).

The porous absorber is a pseudo-porous zone with Darcy-Forchheimer resistance (K and c_F from an
Ergun-type relation), the conductive effective conductivity of the alumina-air matrix (parallel
model) and an effective absorption coefficient extracted from pore-scale simulations (K. Jha et
al., *Role of radiation in a TPMS porous absorber under concentrated irradiation: a pore-scale
study and effective properties for homogenized models*, under review). Radiation is solved with
the finite-volume discrete ordinates method (fvDOM) of OpenFOAM. In comments and outputs,
*paper 1* is that pore-scale study and *paper 2* the paper above.

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
- **Walls** (`radiativeWallT`): k (Tc - Tw)/delta + q_abs - c eps sigma Tw^4 = h (Tw - T_inf) +
  eps sigma (Tw^4 - T_inf^4), with q_abs the absorbed radiation, c eps sigma Tw^4 the emission into
  the channel (fvDOM's own weights) and the right-hand side the loss to the surroundings; emission
  implicit.
- **Mesh.** One block mesh, 209,088 cells: 44 x 44 across (0.63 mm at the walls), axially 48 cells
  in the inlet section (graded to 0.081 mm at the front face, so that the air leaving through the hot
  face keeps a cell Peclet number below 2), 30 in the porous zone (0.61 mm at both faces) and 30 in
  the exit section.
- **Stage 3** (`Allrun.cht`) solves only the energy equation and fvDOM (264 ordinates), with
  under-relaxation 0.7 (air) and 0.9 (porous zone) of h.

## Requirements

- **OpenFOAM-v2206** (openfoam.com, stock release) with Scotch and MPI. The coded boundary conditions
  are compiled at run time (`wmake` must work). Other versions may need changes to the coded
  boundary conditions and function objects.
- **Python 3.8+** with numpy; matplotlib for the figures and for `fieldDiff.py`, `paperTargets.py`
  and `symmetryCheck.py`; scipy for `symmetryCheck.py`.
- **For the figures:** `latex`/`pdflatex` with type1cm, type1ec, cm-super, geometry, underscore,
  amsmath, booktabs, float and graphicx (e.g. texlive-latex-extra, texlive-fonts-recommended,
  cm-super; dvipng is not needed), and `pdftoppm` (poppler), which makes the PNGs. If `latex` is not
  on PATH, set `VSR_TEX_BIN` to its directory.

## Running

Source OpenFOAM and set `VSR` to the **absolute** path of this repository (the scripts change into
each case directory), then work in a run directory outside it or in `runs/`, which git ignores. A
heat-transfer case takes about 3.9 GB on disk, all nine cases about 36 GB.

```sh
source <OpenFOAM-v2206>/etc/bashrc
export VSR=/absolute/path/to/vsr-homogenized-openfoam      # or: cd <clone> && export VSR=$(pwd)
mkdir -p $VSR/runs && cd $VSR/runs
```

**All nine cases** (flows, splits, G 5e6 / Re 300 from 300 K, then the other eight in two lanes,
each started from its nearest finished neighbour). `Allrun.matrix` uses `NPBASE` ranks (default 14)
for the flows and the first case, then two jobs of `NP` ranks at once (default 10); it overwrites
the decomposition of every case it makes. The jobs take hours, so detach them from the terminal:

```sh
NPBASE=14 NP=10 nohup $VSR/case/scripts/Allrun.matrix > log.matrix 2>&1 &
nohup $VSR/case/scripts/stopWhenConverged.sh > log.stopper 2>&1 &   # end each run once converged
```

Other settings of `Allrun.matrix`: `END` (6000, iterations for a new Re), `ENDG` (8000, for a new G),
`ENDCOLD` (8000, the first case) and `NPHI`/`NTHETA` (6/11). `stopWhenConverged.sh` ends a run when
the global residual stays below `TOL` % of the beam (default 0.02; every case below
ended under it); it stops by itself when the matrix ends.

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
#    from a neighbouring solution (WARM is relative to the new case directory), e.g.
#      sed -i 's/^endTime .*/endTime 4000;/' G1e6_Re300_6x11/system/controlDict
#      (cd G1e6_Re300_6x11 && WARM=../G5e6_Re300_6x11/8000 ./Allrun.cht ../split_Re300_6x11 1e6 45)
#    SMOKE=1 ./Allrun.cht ...              throw-away test: 5 iterations

# 4. un-rotated copy for post-processing, axial profiles (where the figure scripts read them), checks
(cd G5e6_Re300_6x11 && ./Allpost)
mkdir -p G5e6_Re300_6x11/figures
python3 $VSR/tools/profiles.py G5e6_Re300_6x11/post <time> G5e6_Re300_6x11/figures/profiles_<time>.csv
python3 $VSR/tools/bulkFaces.py G5e6_Re300_6x11/post <time> G5e6_Re300_6x11/figures/bulk_faces_<time>.csv
python3 $VSR/tools/negativeRays.py G5e6_Re300_6x11 <time>
```

The templates decompose into 14 ranks; for another number set it in every decomposeParDict:
`sed -i "s/^numberOfSubdomains .*/numberOfSubdomains 8;/" <case>/system/decomposeParDict <case>/system/*/decomposeParDict`.
Wall-clock times on 14 ranks (Intel i7-14700): flow 4 min per Re; split under a minute; heat
transfer 0.38 s per iteration, i.e. 50 min for 8000 iterations from 300 K and 25 min for 4000 from a
neighbouring solution. In `Allrun.matrix` (two jobs of 10 ranks) 35-50 min per case and about 90 min
for G 1e6 / Re 300; the whole matrix about 4.5 h.

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
| `tools/symmetryCheck.py 6x11=<case>/post:<time>` | y- and z-mirror parts < 0.2 K rms (odd nTheta) |
| `log.drift` (`Allrun.matrix`; or `tools/drift.py <case> <t0> <t1>`) | well under 1 K over the last write interval |

## Results

Outlet bulk temperature of the nine cases (ray set 6x11) against the pore-scale model; deviation
= (T_pore - T)/T_pore; efficiency = m_dot c_p (T_out - 300 K)/(G A_window), the share of the beam
carried away by the air. Every case passes the checks above (no clipping, no negative intensity,
energy closure better than 0.1 % of the beam, temperature drift < 0.6 K over the last 1000
iterations).

| G (W/m2) | Re | Outlet bulk T, this case | Pore scale | Deviation | Efficiency (this case / pore scale) | Energy closure |
|---|---|---|---|---|---|---|
| 1e6 | 300 | 699 K | 712 K | +1.8 % | 21.0 % / 21.7 % | +0.007 % |
| 1e6 | 600 | 688 K | 679 K | -1.2 % | 40.8 % / 39.9 % | -0.004 % |
| 1e6 | 900 | 651 K | 639 K | -1.8 % | 55.4 % / 53.6 % | -0.005 % |
| 5e6 | 300 | 993 K | 1019 K | +2.5 % | 7.3 % / 7.6 % | +0.067 % |
| 5e6 | 600 | 1042 K | 1038 K | -0.4 % | 15.6 % / 15.5 % | +0.015 % |
| 5e6 | 900 | 1053 K | 1041 K | -1.1 % | 23.8 % / 23.4 % | +0.010 % |
| 1e7 | 300 | 1135 K | 1171 K | +3.1 % | 4.4 % / 4.6 % | +0.020 % |
| 1e7 | 600 | 1201 K | 1196 K | -0.4 % | 9.5 % / 9.4 % | +0.012 % |
| 1e7 | 900 | 1227 K | 1205 K | -1.9 % | 14.6 % / 14.3 % | +0.009 % |

Pressure drop across the porous zone: 7.5, 23.7 and 48.6 Pa at Re = 300, 600 and 900.

## Figures and tables of the paper

```sh
export VSR_RUNS=$VSR/runs          # the run directory of Allrun.matrix
cd $VSR/figures
python3 make_figures.py all        # figures/out/figNN_*.pdf/.png, the paper's file names
python3 tables.py                  # Table 4, and Table 3 if the angular-study runs exist (below)
python3 make_proof.py              # figures/out/figures_proof.pdf: all figures and tables at the paper's sizes
```

Table 3 (angular sensitivity) needs the runs `G5e6_Re300_6x10`, `G5e6_Re300_12x10` and
`G5e6_Re300_8x15` next to `G5e6_Re300_6x11`: repeat steps 2-4 of *One case by hand* at G 5e6 /
Re 300 with `nPhi nTheta` = `6 10`, `12 10` and `8 15` (the two even sets print the sink-ray warning,
which is what the table shows). `tables.py` reads the outlet and face temperatures from their
function objects and the wall maximum from `post/`, and says which runs are missing.

Fig 13 compares with the pore-scale bulk temperatures, digitised from the pore-scale paper's Fig 13
(`figures/data/pore_scale_bulkT_digitised.csv`: outlet values within 0.8 K of that paper's table,
about +-1 K on the gentle parts and +-5-15 K on the steep rise at the front face).

## Tools

| | |
|---|---|
| `beamAlign.py` | fvDOM ordinate set, beam ordinate, rotation for `transformPoints`, boundary sink rays |
| `angularCheck.py` | hemispherical weight error of fvDOM's diffuse emission on each boundary orientation of the rotated channel, for any (nPhi, nTheta) |
| `profiles.py` | axial profiles per cell layer: bulk, mean, min, max and centreline temperature, wall temperature and absorbed radiation |
| `bulkFaces.py` | bulk temperature on every face plane from the face mass fluxes, inlet to outlet (Figs 12-13); its outlet value is that of the energy balance |
| `negativeRays.py <case> <time>` | minimum of every ray intensity field (the case itself, not `post/`) |
| `symmetryCheck.py [--figure f.png] label=<case>/post:<time> ...` | angular error from the departure from the square symmetry of the case |
| `fieldDiff.py <caseA>/post:<t> <caseB>/post:<t>` | cell-by-cell and face-by-face difference of two runs on the same mesh |
| `drift.py <case> <t0> <t1>` | change of T between two write times (convergence) |
| `paperTargets.py` | G 5e6 / Re 300 run against the values reported in the papers |
| `contours.py`, `foamio.py` | plotting and a minimal ASCII OpenFOAM reader used by the others |

## Citation

If you use this case or its tools, please cite the paper and this repository (`CITATION.cff`; GitHub
shows a *Cite this repository* button):

```bibtex
@unpublished{Jha2026homogenized,
  author  = {Jha, Kaushal and Rathore, Surendra Singh and Mehta, Balkrishna and Kumar, Pradeep},
  title   = {Effective property modelling of radiative heat transfer in porous media for concentrated
             solar energy applications},
  year    = {2026},
  note    = {Submitted}
}

@unpublished{Jha2026porescale,
  author  = {Jha, Kaushal and Rathore, Surendra Singh and Mehta, Balkrishna and Kumar, Pradeep},
  title   = {Role of radiation in a {TPMS} porous absorber under concentrated irradiation: a pore-scale
             study and effective properties for homogenized models},
  year    = {2026},
  note    = {Under review}
}

@misc{Jha2026vsrcase,
  author       = {Jha, Kaushal},
  title        = {vsr-homogenized-openfoam: {OpenFOAM} case of a homogenized volumetric solar receiver
                  with a {TPMS} porous absorber},
  year         = {2026},
  note         = {Version 1.0.0},
  howpublished = {\url{https://github.com/WaRxChaMpioN/vsr-homogenized-openfoam}}
}
```

## License

GNU General Public License v3.0 or later (`LICENSE`), as OpenFOAM.

## Contact

Kaushal Jha, kaushal_jha@mines.edu
