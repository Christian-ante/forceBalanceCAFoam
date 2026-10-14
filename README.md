# forceBalanceCAFoam

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23367987.svg)](https://doi.org/10.5281/zenodo.23367987)

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
├── src/forceBalanceCA/                       Library libforceBalanceCA
│   ├── dropletSet/dropletSet.{C,H}           Droplet data shared by boundary condition and solver:
│   │                                         parallel flood-fill (vofID), mass, centre of mass,
│   │                                         velocity, interface reconstruction and drag
│   ├── forceBalanceContactAngleFvPatchScalarField.{C,H}
│   │                                         Boundary condition: adhesion on the wall, drag-to-adhesion
│   │                                         ratio, hysteresis window and dynamic contact angle correlations
│   └── utilities/                            Plane search, plane–cell intersection, polygon area
├── applications/solvers/forceBalanceCAFoam/  interFoam (OpenFOAM-7) + droplet pinning/deceleration
├── tutorials/twoDroplets/                    Example case: two droplets in a channel
├── Allwmake, Allwclean
├── LICENSE / Licence.txt                     GNU GPL v3
└── README.md
```

### How the parts interact

At each evaluation of the contact angle, the boundary condition updates the `dropletSet`, an object stored in the mesh database (`objectRegistry`) and therefore shared with the solver. The `dropletSet` identifies the droplets and computes their kinematics and drag; the boundary condition adds the adhesion force on its wall, computes the drag-to-adhesion ratio, sets the contact angle and stores the adhesion state of each droplet. At the end of the time step the solver reads this state from the `dropletSet` and pins or decelerates the droplets in adhesion. No files are used to exchange data.

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

The solver is already linked to the library; the `libs` entry is needed so that the other OpenFOAM utilities (e.g. `decomposePar`, `foamToVTK`) recognise the boundary condition.

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

### 4. Numerical parameters (optional) — `constant/dropletProperties`

If the file or an entry is missing, the default value below is used, so existing cases run unchanged:

```
alphaMin            0.01;   // volume fraction above which a cell belongs to a droplet
alphaMaxInterface   0.99;   // volume fraction below which a liquid cell is an interface cell
sampleDistance      3;      // distance (cells) of the drag sampling point from the interface
planeTolerance      0.01;   // tolerance on the cell volume fraction in the plane search
planeMaxIter        100;    // maximum bisection iterations of the plane search
planePoints         10;     // sampling points per cell edge in the plane search
planeDRange         10;     // half-width of the initial interval of the plane intercept
```

The values used are printed in the log at start-up.

### 5. Run

```bash
forceBalanceCAFoam                              # serial
decomposePar
mpirun -np <N> forceBalanceCAFoam -parallel     # parallel
```

Run the solver from the case directory: the per-droplet output files are written in the working directory (see below).

## Assumptions and limitations

- **Mesh**: uniform Cartesian mesh with **cubic cells** (square cells in 2D). The interface reconstruction is not valid on non-orthogonal or graded meshes.
- **Flow direction**: drag, gravity component and sliding are evaluated along the **+x direction**. Orient the case so that the airflow (or the downhill direction for inclined walls) is +x.
- **Rigid sliding**: the contact line velocity is approximated by the droplet centre-of-mass velocity; a droplet either slides as a whole or stays pinned.
- **2D cases**: the adhesion correlation is three-dimensional; in 2D the shape factor `kk` accounts for the unresolved transverse direction and must be chosen by the user.
- **Droplet IDs**: when droplets merge, the number of droplets decreases and the indices in the output files shift, so an index does not always identify the same droplet over time.
- **Single wall patch**: apply the boundary condition to one patch only. Adhesion force and state are stored once per droplet, so with several patches the last one evaluated prevails; a warning is printed if more than one patch uses `forceBalanceContactAngle` (merge the walls, e.g. with `createPatch`).
- **Static mesh**: mesh refinement or topology changes are not supported (the run stops with an error).
- **Numerical parameters**: see `constant/dropletProperties` above; drag and pressure are sampled in the air, `sampleDistance` cells away from the interface along its normal.
- **Drag evaluation**: the pressure contribution uses the dynamic pressure `0.5*rho_air*|U|^2` of the sampled air cell (not the static pressure field); the viscous contribution uses `mu_air*(grad(U) + grad(U)^T)` of the same cell.

## Output

In addition to the standard fields, the following are written:

| File / field | Content |
|---|---|
| `vofID` (field, time directories) | Droplet label of each cell (0 = no liquid) |
| `CenterOfMass_<i>.txt` | Time, centre-of-mass position and mean velocity of droplet *i* |
| `ForcesOnDroplet_<i>.txt` | Time, pressure and viscous force, gravity (x), total drag (x), adhesion force, drag-to-adhesion ratio |
| `WARNING_2D/3D_wrongIntersectionNumber_*.txt` | Diagnostic files written if the plane–cell intersection is inconsistent |

The `.txt` files are written in the directory from which the solver is launched, at time 0 and at the output times.

## Tutorial

### `tutorials/twoDroplets`

Two sessile water droplets on the lower wall of a straight channel exposed to an air flow (B600 upstream, B400 downstream, 1 mm apart), as in the PEMFC configuration of the paper but in a shorter channel:

| | |
|---|---|
| Domain | 2 × 1 × 1 mm (x × y × z), uniform cubic cells of 20 μm (100 × 50 × 50 = 250 000 cells) |
| Droplets | Spherical caps with footprint diameter 0.6 mm (B600, x = 0.5 mm) and 0.4 mm (B400, x = 1.5 mm) |
| Wall | `forceBalanceContactAngle` on `lowerWall`: θ0 = 50°, θA = 55°, θR = 45° |
| Air flow | 15 m/s at the inlet, starting after 0.5 ms of relaxation |
| Run | 2 ms on 6 processors |

All the case parameters (geometry, droplets, angles, inlet velocity, times) are collected in `system/input_par`.

```bash
cd tutorials/twoDroplets
./Allrun        # blockMesh, setFields, decomposePar, forceBalanceCAFoam in parallel, reconstructPar
./Allclean      # back to the initial state
```

To change the number of processors, edit `numberOfSubdomains` and `n` in `system/decomposeParDict`. The droplet data are written in `CenterOfMass_<i>.txt` and `ForcesOnDroplet_<i>.txt`; the droplet labels are in the field `vofID`.

## How to cite

If you use forceBalanceCAFoam, please cite the software article (currently under submission to *SoftwareX*):

> C. Antetomaso, S. S. Merola, M. Andersson, H. Garg, *forceBalanceCAFoam: a force-balance-based dynamic contact angle model for multiple droplets in OpenFOAM*, SoftwareX (submitted).

and the archived software:

> C. Antetomaso, S. S. Merola, M. Andersson, H. Garg, *forceBalanceCAFoam* (v1.0.0), Zenodo, 2026. https://doi.org/10.5281/zenodo.23367988

The DOI https://doi.org/10.5281/zenodo.23367987 always resolves to the latest version. The citation metadata are also in `CITATION.cff` (GitHub: *Cite this repository*).

## License

forceBalanceCAFoam is derived from OpenFOAM and is released under the **GNU General Public License v3.0** (see `LICENSE`). OpenFOAM® is a registered trademark of OpenCFD Ltd; this software is not approved or endorsed by OpenCFD Ltd or the OpenFOAM Foundation.

## Acknowledgements

Partial financial support from the European Union, NextGenerationEU, National Sustainable Mobility Center – MOST, CN00000023, Italian Ministry of University and Research Decree n. 1033 – 17/06/2022, Spoke 12, CUP B43C22000440001.

## Contact

Christian Antetomaso — CNR STEMS, Naples, Italy — christian.antetomaso@stems.cnr.it
