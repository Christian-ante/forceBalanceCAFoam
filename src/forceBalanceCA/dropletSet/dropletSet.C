/*---------------------------------------------------------------------------*\
  forceBalanceCAFoam: force-balance dynamic contact angle for OpenFOAM-7
-------------------------------------------------------------------------------
License
    This file is part of forceBalanceCAFoam, derived from OpenFOAM.

    forceBalanceCAFoam is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dropletSet.H"
#include "fvcGrad.H"
#include "globalIndex.H"
#include "IOdictionary.H"
#include "uniformDimensionedFields.H"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>

// Geometric utilities of the interface reconstruction (global namespace)
#include "../utilities/vector2stringHelper.H"
#include "../utilities/interfaceArea/findInterfacePlane.H"
#include "../utilities/interfaceArea/findIntersection.H"
#include "../utilities/interfaceArea/intersectionArea.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dropletSet, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dropletSet::dropletSet(const fvMesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    sigma_(0),
    rhoWater_(0),
    rhoAir_(0),
    muWater_(0),
    muAir_(0),
    g_(Zero),
    alphaMin_(0.01),
    alphaMaxInterface_(0.99),
    sampleDistance_(3),
    planeTolerance_(0.01),
    planeMaxIter_(100),
    planePoints_(10),
    planeDRange_(10),
    vofID_
    (
        IOobject
        (
            "vofID",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE        // written at output times
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    nInterface_
    (
        IOobject
        (
            "nInterface",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false                       // not registered
        ),
        mesh,
        dimensionedVector("zero", dimless, Zero)
    ),
    droplets_(),
    labelToIndex_(),
    lastCentreWriteTime_(-1.0),
    lastForceWriteTime_(-1.0)
{
    readProperties();
}


Foam::dropletSet& Foam::dropletSet::New(const fvMesh& mesh)
{
    if (mesh.foundObject<dropletSet>(typeName))
    {
        return mesh.lookupObjectRef<dropletSet>(typeName);
    }

    dropletSet* ptr = new dropletSet(mesh);
    ptr->store();       // ownership passed to the mesh registry

    return *ptr;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::dropletSet::readProperties()
{
    // Fluid properties: constant/transportProperties
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const scalar def = -999;    // Value returned if an entry cannot be read

    sigma_ = transportProperties.lookupOrDefault<scalar>("sigma", def);
    rhoAir_ =
        transportProperties.subDict("air").lookupOrDefault<scalar>("rho", def);
    rhoWater_ =
        transportProperties.subDict("water").lookupOrDefault<scalar>("rho", def);
    const scalar nuAir =
        transportProperties.subDict("air").lookupOrDefault<scalar>("nu", def);
    const scalar nuWater =
        transportProperties.subDict("water").lookupOrDefault<scalar>("nu", def);

    if
    (
        sigma_ == def || rhoAir_ == def || rhoWater_ == def
     || nuAir == def || nuWater == def
    )
    {
        FatalErrorInFunction
            << "Cannot read material properties (water/air rho and nu, sigma)!"
            << nl << "Please check your transportProperties file."
            << exit(FatalError);
    }

    muWater_ = rhoWater_*nuWater;   // Used for the capillary number
    muAir_ = rhoAir_*nuAir;         // Used for the viscous drag (air side)

    // Gravity: registered by the solver, otherwise read from constant/g
    if (mesh_.foundObject<uniformDimensionedVectorField>("g"))
    {
        g_ = mesh_.lookupObject<uniformDimensionedVectorField>("g").value();
    }
    else
    {
        uniformDimensionedVectorField g
        (
            IOobject
            (
                "g",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        g_ = g.value();
    }

    // Numerical parameters: constant/dropletProperties (optional)
    IOdictionary dropletProperties
    (
        IOobject
        (
            "dropletProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ),
        dictionary()
    );

    alphaMin_ = dropletProperties.lookupOrDefault<scalar>("alphaMin", 0.01);
    alphaMaxInterface_ =
        dropletProperties.lookupOrDefault<scalar>("alphaMaxInterface", 0.99);
    sampleDistance_ =
        dropletProperties.lookupOrDefault<label>("sampleDistance", 3);
    planeTolerance_ =
        dropletProperties.lookupOrDefault<scalar>("planeTolerance", 0.01);
    planeMaxIter_ = dropletProperties.lookupOrDefault<label>("planeMaxIter", 100);
    planePoints_ = dropletProperties.lookupOrDefault<label>("planePoints", 10);
    planeDRange_ = dropletProperties.lookupOrDefault<scalar>("planeDRange", 10);

    Info<< "dropletSet: alphaMin " << alphaMin_
        << ", alphaMaxInterface " << alphaMaxInterface_
        << ", sampleDistance " << sampleDistance_
        << ", planeTolerance " << planeTolerance_
        << ", planeMaxIter " << planeMaxIter_
        << ", planePoints " << planePoints_
        << ", planeDRange " << planeDRange_ << endl;
}


void Foam::dropletSet::identifyDroplets(const volScalarField& alpha)
{
    volScalarField& vofID = vofID_;
    vofID = dimensionedScalar("zero", dimless, 0);

    // Global cell index, to avoid conflicts between processors
    globalIndex globalNumbering(mesh_.nCells());

    // Flood-fill on the local sub-domain: each connected liquid region takes
    // the global index of its first (lowest) liquid cell
    forAll(mesh_.cells(), cellI)
    {
        if (alpha[cellI] > alphaMin_ && vofID[cellI] == 0.0)
        {
            const label globalCellID = globalNumbering.toGlobal(cellI);

            vofID[cellI] = globalCellID;

            DynamicList<label> cells(1);
            cells.append(mesh_.cellCells()[cellI]);

            while (!cells.empty())
            {
                DynamicList<label> neighbours(0);

                forAll(cells, i)
                {
                    const label cellI_ = cells[i];

                    if (alpha[cellI_] > alphaMin_ && vofID[cellI_] == 0.0)
                    {
                        vofID[cellI_] = globalCellID;
                        neighbours.append(mesh_.cellCells()[cellI_]);
                    }
                }
                cells.clear();
                cells.transfer(neighbours);
            }
        }
    }
    vofID.correctBoundaryConditions();

    // Parallel: propagate the lowest label across processor boundaries
    if (Pstream::parRun())
    {
        bool corrCheck = false;

        for (int iter = 0; iter < Pstream::nProcs() - 1; iter++)
        {
            forAll(mesh_.boundary(), patchI)
            {
                const fvPatch& pPatch = mesh_.boundary()[patchI];

                if (pPatch.coupled())
                {
                    const scalarField nbrID
                    (
                        vofID.boundaryField()[patchI].patchNeighbourField()
                    );

                    forAll(pPatch, faceI)
                    {
                        const label cellI = pPatch.faceCells()[faceI];
                        const scalar ID_own = vofID[cellI];
                        const scalar ID_nei = nbrID[faceI];

                        if (ID_own > ID_nei && ID_nei > 0)
                        {
                            vofID[cellI] = ID_nei;
                            corrCheck = true;

                            DynamicList<label> cellsCorr(1);
                            cellsCorr.append(mesh_.cellCells()[cellI]);

                            while (!cellsCorr.empty())
                            {
                                DynamicList<label> neighboursCorr(1);

                                forAll(cellsCorr, i)
                                {
                                    if (vofID[cellsCorr[i]] > vofID[cellI])
                                    {
                                        vofID[cellsCorr[i]] = vofID[cellI];
                                        neighboursCorr.append
                                        (
                                            mesh_.cellCells()[cellsCorr[i]]
                                        );
                                    }
                                }
                                cellsCorr.clear();
                                cellsCorr.transfer(neighboursCorr);
                            }
                        }
                    }
                }
            }
            vofID.correctBoundaryConditions();

            reduce(corrCheck, orOp<bool>());
            if (!corrCheck)
            {
                break;
            }
            corrCheck = false;
        }
    }

    // List of the droplet labels, identical on all processors
    List<label> labels;

    if (Pstream::parRun())
    {
        DynamicList<label> localLabels;
        forAll(vofID, cellI)
        {
            const label id = vofID[cellI];
            if (id > 0 && findIndex(localLabels, id) == -1)
            {
                localLabels.append(id);
            }
        }

        List<List<label>> gatheredLabels(Pstream::nProcs());
        gatheredLabels[Pstream::myProcNo()] = localLabels;
        Pstream::gatherList(gatheredLabels);

        DynamicList<label> globalLabels;
        if (Pstream::master())
        {
            forAll(gatheredLabels, procI)
            {
                forAll(gatheredLabels[procI], i)
                {
                    const label id = gatheredLabels[procI][i];
                    if (findIndex(globalLabels, id) == -1)
                    {
                        globalLabels.append(id);
                    }
                }
            }
        }
        labels = globalLabels;
        Pstream::scatter(labels);
    }
    else
    {
        DynamicList<label> localLabels;
        forAll(vofID, cellI)
        {
            const label id = vofID[cellI];
            if (id > 0 && findIndex(localLabels, id) == -1)
            {
                localLabels.append(id);
            }
        }
        labels = localLabels;
    }

    // Droplet list and label -> index map
    droplets_.assign(labels.size(), droplet());
    labelToIndex_.clear();
    forAll(labels, i)
    {
        droplets_[i].id = labels[i];
        labelToIndex_.insert(labels[i], i);
    }
}


void Foam::dropletSet::computeKinematics
(
    const volScalarField& alpha,
    const volVectorField& U
)
{
    const label nDrops = size();

    List<label> nCells(nDrops, Zero);
    List<vector> sumCentre(nDrops, Zero);
    List<vector> sumVelocity(nDrops, Zero);
    List<scalar> mass(nDrops, Zero);

    // Single pass over the cells
    forAll(mesh_.cells(), celli)
    {
        if (alpha[celli] > alphaMin_)
        {
            const label i = index(vofID_[celli]);
            if (i >= 0)
            {
                nCells[i]++;
                sumCentre[i] += mesh_.C()[celli];
                mass[i] += rhoWater_*mesh_.V()[celli]*alpha[celli];
                sumVelocity[i] += U[celli];
            }
        }
    }

    for (label i = 0; i < nDrops; i++)
    {
        if (Pstream::parRun())
        {
            reduce(nCells[i], sumOp<label>());
            reduce(sumCentre[i], sumOp<vector>());
            reduce(mass[i], sumOp<scalar>());
            reduce(sumVelocity[i], sumOp<vector>());
        }

        droplet& d = droplets_[i];
        d.mass = mass[i];
        d.centre = nCells[i] > 0 ? sumCentre[i]/nCells[i] : vector(Zero);
        d.velocity = nCells[i] > 0 ? sumVelocity[i]/nCells[i] : vector(Zero);
    }
}


void Foam::dropletSet::computeDrag
(
    const volScalarField& alpha,
    const volVectorField& U
)
{
    using Point = std::vector<double>;
    using Edge = std::pair<Point, Point>;
    using EdgeList = std::vector<Edge>;
    using PointCloud = std::vector<Point>;
    using VectorList = std::vector<std::vector<double>>;

    // Interface normal (points into the liquid) and velocity gradient
    const volVectorField gradAlpha(fvc::grad(alpha));
    const dimensionedScalar smallDim("smallDim", gradAlpha.dimensions(), SMALL);
    nInterface_ = gradAlpha/(mag(gradAlpha) + smallDim);
    const volVectorField& nInter = nInterface_;

    const volTensorField gradU(fvc::grad(U));

    // Build the cell decomposition used by findCell on all processors
    // (collective in parallel: must not depend on the local droplets)
    (void)mesh_.tetBasePtIs();

    const label nDrops = size();
    List<vector> pressureForce(nDrops, Zero);
    List<vector> viscousForce(nDrops, Zero);
    List<scalar> interfaceArea(nDrops, Zero);

    // Point clouds used by the plane search, rebuilt only if the cell size
    // changes (uniform mesh: built once)
    scalar lCloud2D = -1, lCloud3D = -1;
    PointCloud cloud2D, cloud3D;

    const int numPoints = planePoints_;

    // Single pass over the interface cells
    forAll(mesh_.cells(), celli)
    {
        if (!(alpha[celli] > alphaMin_ && alpha[celli] < alphaMaxInterface_))
        {
            continue;
        }

        const label i = index(vofID_[celli]);
        if (i < 0)
        {
            continue;
        }

        double dA = 0;  // Interface area in the cell

        // Plane a*x + b*y + c*z = d approximating the interface
        const double a = nInter[celli][0];
        const double b = nInter[celli][1];
        const double c = nInter[celli][2];

        // NB: the cell must be a cube (also in 2D)
        const double l = cbrt(mesh_.V()[celli]);

        double D_min = -planeDRange_;
        double D_max = planeDRange_;

        if (std::abs(c) < 0.05)         // 2D: y = A*x + D
        {
            if (std::abs(a) > 0.95 || std::abs(b) > 0.95)
            {
                dA = l*l;               // Plane parallel to a cell face
            }
            else
            {
                const double A = -a/b;

                EdgeList edges =
                {
                    {{0, 0}, {l, 0}},
                    {{0, 0}, {0, l}},
                    {{l, 0}, {l, l}},
                    {{0, l}, {l, l}}
                };

                if (l != lCloud2D)
                {
                    cloud2D.clear();
                    cloud2D.reserve(numPoints*numPoints);
                    const double step = l/(numPoints - 1);
                    for (int ii = 0; ii < numPoints; ++ii)
                    {
                        for (int jj = 0; jj < numPoints; ++jj)
                        {
                            cloud2D.push_back({ii*step, jj*step});
                        }
                    }
                    lCloud2D = l;
                }

                const double D = findInterfacePlane
                (
                    l, A, D_min, D_max, planeMaxIter_, planeTolerance_,
                    alpha[celli], cloud2D
                );
                VectorList intersections = findIntersection(A, D, edges);

                if (intersections.size() == 2)
                {
                    dA = intersectionArea(intersections, l);
                }
                else
                {
                    std::stringstream fileName;
                    fileName << "WARNING_2D_wrongIntersectionNumber="
                        << intersections.size() << "_time="
                        << mesh_.time().value() << "[s].txt";
                    std::ofstream outFileW(fileName.str(), std::ios::out);
                    outFileW.close();

                    if (intersections.size() > 2)
                    {
                        dA = intersectionArea(intersections, l);
                    }
                }
            }
        }
        else                            // 3D: z = A*x + B*y + D
        {
            const double A = -a/c;
            const double B = -b/c;

            EdgeList edges =
            {
                {{0, 0, 0}, {l, 0, 0}},
                {{0, 0, 0}, {0, l, 0}},
                {{0, 0, 0}, {0, 0, l}},
                {{l, 0, 0}, {l, l, 0}},
                {{l, 0, 0}, {l, 0, l}},
                {{0, l, 0}, {l, l, 0}},
                {{0, l, 0}, {0, l, l}},
                {{0, 0, l}, {l, 0, l}},
                {{0, 0, l}, {0, l, l}},
                {{l, l, 0}, {l, l, l}},
                {{l, 0, l}, {l, l, l}},
                {{0, l, l}, {l, l, l}}
            };

            if (l != lCloud3D)
            {
                cloud3D.clear();
                cloud3D.reserve(numPoints*numPoints*numPoints);
                const double step = l/(numPoints - 1);
                for (int ii = 0; ii < numPoints; ++ii)
                {
                    for (int jj = 0; jj < numPoints; ++jj)
                    {
                        for (int kk = 0; kk < numPoints; ++kk)
                        {
                            cloud3D.push_back({ii*step, jj*step, kk*step});
                        }
                    }
                }
                lCloud3D = l;
            }

            const double D = findInterfacePlane
            (
                l, A, B, D_min, D_max, planeMaxIter_, planeTolerance_,
                alpha[celli], cloud3D
            );
            VectorList intersections = findIntersection(A, B, D, edges);

            if (intersections.size() < 3 || intersections.size() > 6)
            {
                std::stringstream fileName;
                fileName << "WARNING_3D_wrongIntersectionNumber="
                    << intersections.size() << "_time="
                    << mesh_.time().value() << "[s].txt";
                std::ofstream outFileW(fileName.str(), std::ios::out);
                outFileW.close();
            }

            if (intersections.size() >= 3)
            {
                dA = intersectionArea(intersections);
            }
        }

        // Pressure and viscous contributions, sampled in the air at
        // sampleDistance cells from the interface along its normal
        const vector samplePoint =
            mesh_.C()[celli] - scalar(sampleDistance_)*nInter[celli]*l;

        label sampleCell = mesh_.findCell(samplePoint);    // Octree search
        if (sampleCell < 0)     // Outside the local sub-domain or the mesh
        {
            sampleCell = mesh_.findNearestCell(samplePoint);
        }

        // Dynamic pressure of the approaching air, 1/2*rho_air*|U|^2
        const scalar pDyn = 0.5*rhoAir_*magSqr(U[sampleCell]);
        pressureForce[i] += nInter[celli]*pDyn*dA;

        // Viscous stress mu_air*(grad(U) + grad(U)^T) on the outward normal
        // of the droplet (nInter points into the droplet)
        const tensor viscousStress(muAir_*twoSymm(gradU[sampleCell]));
        viscousForce[i] += viscousStress & (-nInter[celli]*dA);

        interfaceArea[i] += dA;
    }

    for (label i = 0; i < nDrops; i++)
    {
        vector dragForce = pressureForce[i] + viscousForce[i];

        if (Pstream::parRun())
        {
            reduce(pressureForce[i], sumOp<vector>());
            reduce(viscousForce[i], sumOp<vector>());
            reduce(dragForce, sumOp<vector>());
            reduce(interfaceArea[i], sumOp<scalar>());
        }

        droplet& d = droplets_[i];
        d.pressureForce = pressureForce[i];
        d.viscousForce = viscousForce[i];
        d.dragForce = dragForce;
        d.interfaceArea = interfaceArea[i];
        d.gravityForceX = d.mass*g_.x();
    }
}


void Foam::dropletSet::writeCentreOfMass()
{
    if (!Pstream::master())
    {
        return;
    }

    const scalar currentTime = mesh_.time().value();

    if
    (
        (mesh_.time().outputTime() || currentTime == 0)
     && currentTime != lastCentreWriteTime_
    )
    {
        lastCentreWriteTime_ = currentTime;

        for (label i = 0; i < size(); i++)
        {
            const std::string fileName =
                "CenterOfMass_" + std::to_string(i + 1) + ".txt";
            std::ofstream outFile(fileName, std::ios::out | std::ios::app);
            if (outFile.is_open())
            {
                if (outFile.tellp() == 0)
                {
                    outFile << "Time [s]\tCenter Of Mass Position [m]"
                        "\tCenter Of Mass Velocity [m/s]\n";
                }
                outFile << std::fixed << std::setprecision(6)
                    << currentTime << "\t"
                    << vectorToString(droplets_[i].centre) << "\t"
                    << vectorToString(droplets_[i].velocity) << "\n";
                outFile.close();
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dropletSet::update()
{
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>("alpha.water");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    // The droplet fields are sized on the initial mesh (the interface
    // reconstruction requires a fixed, uniform Cartesian mesh anyway)
    if (vofID_.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "The number of cells has changed (mesh refinement or topology"
            << " change):" << nl
            << "    this is not supported by forceBalanceCAFoam."
            << exit(FatalError);
    }

    identifyDroplets(alpha);
    computeKinematics(alpha, U);
    writeCentreOfMass();
    computeDrag(alpha, U);
}


void Foam::dropletSet::setAdhesion
(
    const label i,
    const scalar Fadh,
    const scalar S
)
{
    droplets_[i].adhesionForce = Fadh;
    droplets_[i].ratio = S;
}


void Foam::dropletSet::setAdhesionState
(
    const label i,
    const bool inAdhesion,
    const scalar deceleration
)
{
    droplets_[i].inAdhesion = inAdhesion;
    droplets_[i].deceleration = deceleration;
}


void Foam::dropletSet::writeForces()
{
    if (!Pstream::master())
    {
        return;
    }

    const scalar currentTime = mesh_.time().value();

    if
    (
        (mesh_.time().outputTime() || currentTime == 0)
     && currentTime != lastForceWriteTime_
    )
    {
        lastForceWriteTime_ = currentTime;

        for (label i = 0; i < size(); i++)
        {
            const droplet& d = droplets_[i];

            const std::string fileName =
                "ForcesOnDroplet_" + std::to_string(i + 1) + ".txt";
            std::ofstream outFile(fileName, std::ios::out | std::ios::app);
            if (outFile.is_open())
            {
                if (outFile.tellp() == 0)
                {
                    outFile << "Time [s]\tPressure Force [N]\tViscous Force [N]"
                        "\tG Force X [N]\tTotal Drag X-dir [N]"
                        "\tAdhesion X-dir [N]\tRatio[-]\n";
                }
                outFile << std::fixed << std::setprecision(8)
                    << currentTime << "\t"
                    << vectorToString(d.pressureForce) << "\t"
                    << vectorToString(d.viscousForce) << "\t"
                    << d.gravityForceX << "\t"
                    << d.dragForce[0] << "\t"
                    << d.adhesionForce << "\t"
                    << d.ratio << "\n";
                outFile.close();
            }
        }
    }
}


void Foam::dropletSet::correctVelocity
(
    volVectorField& U,
    const volScalarField& alpha
) const
{
    const scalar deltaT = mesh_.time().deltaTValue();

    for (label i = 0; i < size(); i++)
    {
        const droplet& d = droplets_[i];

        if (!d.inAdhesion)
        {
            continue;
        }

        if (d.deceleration == 0)
        {
            // The droplet was at rest: pin it (cells next to the walls)
            forAll(mesh_.boundaryMesh(), patchi)
            {
                const polyPatch& pp = mesh_.boundaryMesh()[patchi];

                if (pp.name() == "defaultFaces")    // 2D front/back
                {
                    continue;
                }

                const labelUList& faceCells = pp.faceCells();

                forAll(faceCells, facei)
                {
                    const label celli = faceCells[facei];

                    if (alpha[celli] > alphaMin_ && vofID_[celli] == d.id)
                    {
                        U[celli] = Zero;
                    }
                }
            }
        }
        else
        {
            // The droplet was moving: decelerate it, without reversing it
            forAll(U, celli)
            {
                if (alpha[celli] > alphaMin_ && vofID_[celli] == d.id)
                {
                    U[celli][0] =
                        max(0.0, U[celli][0] - d.deceleration*deltaT);
                }
            }
        }
    }
}


// ************************************************************************* //
