# forceBalanceCAFoam

**A force-balance-based dynamic contact angle model for multiple droplets in OpenFOAM**

forceBalanceCAFoam extends the OpenFOAM-7 `interFoam` solver and its `dynamicAlphaContactAngle` boundary condition. The dynamic contact angle is driven by a **droplet-level force balance** between aerodynamic drag and adhesion, instead of the empirical velocity scale `u_theta` of the standard boundary condition. Each droplet in the domain is identified, its forces are integrated over its own interface, and the contact angle is set either within the contact angle hysteresis (CAH) window (pinned droplet) or beyond it (sliding droplet).

Main features:

- **Parallel droplet identification**: a flood-fill algorithm gives every droplet a unique, globally consistent label, also when a droplet spans several processors.
- **Per-droplet force integration**: a PLIC-like reconstruction of the interface gives the interfacial area in each cell, over which pressure and viscous contributions are integrated. The adhesion force follows ElSherbini & Jacobi (2006).
- **Explicit CAH treatment**: the contact angle is interpolated inside the hysteresis window while adhesion prevails, and follows the Kistler (advancing) and Nichita et al. (receding) correlations once the droplet slides.
- **Solver coupling**: droplets in adhesion are pinned (if they were at rest) or decelerated (if they were moving); coalescence is handled naturally by the labelling.

The code was developed for water removal in the gas channels of PEM fuel cells, but applies to any wall-bounded two-phase flow with droplets.

---

## Repository structure

```
forceBalanceCAFoam/
├── src/forceBalanceCA/                     Boundary condition library (libforceBalanceCA)
│   ├── forceBalanceContactAngleFvPatchScalarField.{C,H}
│   ├── readProperties.H                    Reads transportProperties and g
│   ├── findDrop.H                          Parallel flood-fill droplet identification (vofID field)
│   ├── dropG.H, dropVel.H                  Centre of mass, mass and mean velocity of each droplet
│   ├── forcesOnDroplet.H                   Interface reconstruction, drag and adhesion forces
│   ├── lateralScaleFactor.H                Blending between advancing and receding fronts
│   ├── contactAngleHysteresis.H            Contact angle inside the CAH window (pinned droplets)
│   ├── contactAngleCorrelations.H          Dynamic contact angle correlations (sliding droplets)
│   ├── printGpositionVelocity.H, printForces.H   Per-droplet output files
│   └── utilities/                          Plane search, plane–cell intersection, polygon area
├── applications/solvers/forceBalanceCAFoam/   interFoam (OpenFOAM-7) + adh_check.H
├── Allwmake, Allwclean
├── LICENSE / Licence.txt                   GNU GPL v3
└── README.md
```

## Requirements

| Item | Version used for development and testing |
|---|---|
| OpenFOAM | 7 (OpenFOAM Foundation, openfoam.org) |
| Operating system | Linux, Ubuntu 20.04 |
| Compiler / build | GCC as configured by OpenFOAM, `wmake` |
| Parallel runs | MPI (as configured by OpenFOAM) |

Other OpenFOAM versions or the ESI distribution (openfoam.com) have not been tested.

## Installation

1. Source the OpenFOAM-7 environment, e.g.

   ```bash
   source /opt/openfoam7/etc/bashrc
   ```

2. Clone the repository and compile:

   ```bash
   git clone https://github.com/Christian-ante/forceBalanceCAFoam.git
   cd forceBalanceCAFoam
   ./Allwmake
   ```

   This builds:
   - the library `libforceBalanceCA.so` in `$FOAM_USER_LIBBIN`;
   - the solver `forceBalanceCAFoam` in `$FOAM_USER_APPBIN`.

3. To remove the build products: `./Allwclean`.

## Usage

### 1. Load the library — `system/controlDict`

```
libs ("libforceBalanceCA.so");
```

### 2. Boundary condition — `0/alpha.water`

Apply the boundary condition to the wall(s) on which droplets are deposited:

```
wall
{
    type            forceBalanceContactAngle;
    theta0          50;         // static (equilibrium) contact angle [deg]
    thetaA          55;         // limiting advancing contact angle [deg]
    thetaR          45;         // limiting receding contact angle [deg]
    kk              1;          // geometric shape factor of the adhesion force (2D cases)
    limit           gradient;   // as in the standard alphaContactAngle conditions
    value           uniform 0;
}
```

The values above are those of the PEMFC substrate used in the paper; replace them with the angles measured for your surface. No velocity scale (`uTheta`) is needed.

### 3. Material properties — `constant/transportProperties`

The boundary condition reads the properties from `transportProperties` and **requires the two phases to be named `water` and `air`** (volume fraction field `alpha.water`):

```
phases (water air);

water
{
    transportModel  Newtonian;
    nu              1e-06;
    rho             1000;
}

air
{
    transportModel  Newtonian;
    nu              1.48e-05;
    rho             1;
}

sigma           0.07;
```

It reads `rho` and `nu` of both phases and `sigma`, written as plain scalars as in the standard OpenFOAM-7 interFoam tutorials (no extra entries are needed). The dynamic viscosities are computed as `mu = rho*nu`: the air viscosity is used for the viscous drag, the water viscosity for the capillary number. The gravity vector is read from `constant/g`.

### 4. Run

```bash
forceBalanceCAFoam                              # serial
decomposePar
mpirun -np <N> forceBalanceCAFoam -parallel     # parallel
```

Run the solver **from the case directory**: the boundary condition and the solver exchange information through text files written in the working directory (see below).

## Assumptions and limitations

- **Mesh**: uniform Cartesian mesh with **cubic cells** (square cells in 2D). The interface reconstruction is not valid on non-orthogonal or graded meshes.
- **Flow direction**: drag, gravity component and sliding are evaluated along the **+x direction**. Orient the case so that the airflow (or the downhill direction for inclined walls) is +x.
- **Rigid sliding**: the contact line velocity is approximated by the droplet centre-of-mass velocity; a droplet either slides as a whole or stays pinned.
- **2D cases**: the adhesion correlation is three-dimensional; in 2D the shape factor `kk` accounts for the unresolved transverse direction and must be chosen by the user.
- **Droplet IDs**: when droplets merge, the number of droplets decreases and the indices in the output files shift, so an index does not always identify the same droplet over time.
- **Single wall patch**: apply the boundary condition to one patch only. Each patch keeps its own state, but the exchange files read by the solver are shared; a warning is printed if more than one patch uses `forceBalanceContactAngle` (merge the walls, e.g. with `createPatch`).
- Hard-coded numerical parameters (in the source): liquid threshold `alpha > 0.02` (boundary condition) and `alpha > 0.01` (solver); reconstruction tolerance `0.01` on the cell volume fraction; flow quantities sampled in the air, 3 cells away from the interface along its normal.
- **Drag evaluation**: the pressure contribution uses the dynamic pressure `0.5*rho_air*|U|^2` of the sampled air cell (not the static pressure field); the viscous contribution uses `mu_air*(grad(U) + grad(U)^T)` of the same cell.

## Output

In addition to the standard fields, the following are written:

| File / field | Content |
|---|---|
| `vofID` (field, time directories) | Droplet label of each cell (0 = no liquid) |
| `CenterOfMass_<i>.txt` | Time, centre-of-mass position and mean velocity of droplet *i* |
| `ForcesOnDroplet_<i>.txt` | Time, pressure and viscous force, gravity (x), total drag (x), adhesion force, drag-to-adhesion ratio |
| `dropNumber.txt`, `droplet_<i>_label.txt` | Number of droplets and their labels (read by the solver) |
| `Adhesion_Check_<i>.txt`, `dropDeceleration_<i>.txt` | Present only while droplet *i* is in adhesion; deceleration applied by the solver |
| `WARNING_2D/3D_wrongIntersectionNumber_*.txt` | Diagnostic files written if the plane–cell intersection is inconsistent |

The `.txt` files are written in the directory from which the solver is launched.

## Tutorial

*A ready-to-run example case will be added in `tutorials/`.*

## How to cite

If you use forceBalanceCAFoam, please cite the software article (currently under submission to *SoftwareX*):

> C. Antetomaso, S. S. Merola, M. Andersson, H. Garg, *forceBalanceCAFoam: a force-balance-based dynamic contact angle model for multiple droplets in OpenFOAM*, SoftwareX (submitted).

A DOI for the archived release will be added here.

## License

forceBalanceCAFoam is derived from OpenFOAM and is released under the **GNU General Public License v3.0** (see `LICENSE`). OpenFOAM® is a registered trademark of OpenCFD Ltd; this software is not approved or endorsed by OpenCFD Ltd or the OpenFOAM Foundation.

## Acknowledgements

Partial financial support from the European Union, NextGenerationEU, National Sustainable Mobility Center – MOST, CN00000023, Italian Ministry of University and Research Decree n. 1033 – 17/06/2022, Spoke 12, CUP B43C22000440001.

## Contact

Christian Antetomaso — CNR STEMS, Naples, Italy — christian.antetomaso@stems.cnr.it
