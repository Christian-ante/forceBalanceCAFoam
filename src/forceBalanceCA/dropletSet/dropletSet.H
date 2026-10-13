/*---------------------------------------------------------------------------*\
  forceBalanceCAFoam: force-balance dynamic contact angle for OpenFOAM-7
-------------------------------------------------------------------------------
License
    This file is part of forceBalanceCAFoam, derived from OpenFOAM.

    forceBalanceCAFoam is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

Class
    Foam::dropletSet

Description
    Set of the liquid droplets present in the domain, shared between the
    forceBalanceContactAngle boundary condition and the forceBalanceCAFoam
    solver through the mesh object registry (no exchange files).

    At each update the set:
      - identifies the droplets with a parallel flood-fill algorithm
        (field vofID: one global label per droplet, 0 = no liquid);
      - computes mass, centre of mass and mean velocity of each droplet;
      - reconstructs the interface (PLIC-like, cubic cells) and integrates
        the pressure and viscous drag acting on each droplet.

    The adhesion force, the drag-to-adhesion ratio and the adhesion state are
    set by the boundary condition (they depend on the wall patch); the solver
    then calls correctVelocity() to pin or decelerate droplets in adhesion.

    Optional numerical parameters are read from constant/dropletProperties;
    if the file (or an entry) is missing, the default value is used.

SourceFiles
    dropletSet.C

\*---------------------------------------------------------------------------*/

#ifndef dropletSet_H
#define dropletSet_H

#include "fvMesh.H"
#include "volFields.H"
#include "regIOobject.H"
#include "Map.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Struct droplet Declaration
\*---------------------------------------------------------------------------*/

//- Data of a single droplet
struct droplet
{
    //- Droplet label (value of vofID in the droplet cells)
    label id = 0;

    //- Droplet mass [kg]
    scalar mass = 0;

    //- Centre of mass (mean of the liquid cell centres) [m]
    vector centre = Zero;

    //- Mean velocity of the liquid cells [m/s]
    vector velocity = Zero;

    //- Pressure (dynamic pressure) contribution to the drag [N]
    vector pressureForce = Zero;

    //- Viscous contribution to the drag [N]
    vector viscousForce = Zero;

    //- Total drag (pressure + viscous) [N]
    vector dragForce = Zero;

    //- Reconstructed interface area [m2]
    scalar interfaceArea = 0;

    //- Weight component along x [N]
    scalar gravityForceX = 0;

    // Set by the contact angle boundary condition

        //- Maximum adhesion force [N]
        scalar adhesionForce = 0;

        //- Drag-to-adhesion ratio S [-]
        scalar ratio = 0;

        //- True if the droplet is in adhesion (S <= 1)
        bool inAdhesion = false;

        //- Deceleration applied by the solver while in adhesion [m/s2]
        //  (0 = the droplet was at rest: pin it)
        scalar deceleration = 0;
};


/*---------------------------------------------------------------------------*\
                         Class dropletSet Declaration
\*---------------------------------------------------------------------------*/

class dropletSet
:
    public regIOobject
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        // Fluid properties (constant/transportProperties, constant/g)

            scalar sigma_;
            scalar rhoWater_;
            scalar rhoAir_;
            scalar muWater_;
            scalar muAir_;
            vector g_;

        // Numerical parameters (constant/dropletProperties, optional)

            //- Volume fraction above which a cell belongs to a droplet
            scalar alphaMin_;

            //- Volume fraction below which a liquid cell is an interface cell
            scalar alphaMaxInterface_;

            //- Distance (in cells) of the drag sampling point from the
            //  interface, along the interface normal
            label sampleDistance_;

            //- Tolerance on the cell volume fraction for the plane search
            scalar planeTolerance_;

            //- Maximum number of bisection iterations for the plane search
            label planeMaxIter_;

            //- Number of sampling points per cell edge for the plane search
            label planePoints_;

            //- Half-width of the initial interval of the plane intercept D
            scalar planeDRange_;

        //- Droplet label of each cell (0 = no liquid)
        volScalarField vofID_;

        //- Unit normal to the interface, grad(alpha)/|grad(alpha)|
        //  (points into the liquid)
        volVectorField nInterface_;

        //- Droplets, in the order in which they are identified
        std::vector<droplet> droplets_;

        //- Droplet label -> index in droplets_
        Map<label> labelToIndex_;

        //- Last time written to the CenterOfMass_<i>.txt files
        scalar lastCentreWriteTime_;

        //- Last time written to the ForcesOnDroplet_<i>.txt files
        scalar lastForceWriteTime_;


    // Private Member Functions

        //- Read the fluid properties and the numerical parameters
        void readProperties();

        //- Label the droplets (parallel flood-fill) and fill droplets_
        void identifyDroplets(const volScalarField& alpha);

        //- Mass, centre of mass and mean velocity of each droplet
        void computeKinematics
        (
            const volScalarField& alpha,
            const volVectorField& U
        );

        //- Interface reconstruction and drag of each droplet
        void computeDrag
        (
            const volScalarField& alpha,
            const volVectorField& U
        );

        //- Append centre of mass and velocity to CenterOfMass_<i>.txt
        void writeCentreOfMass();

        //- Disallow default bitwise copy construction and assignment
        dropletSet(const dropletSet&) = delete;
        void operator=(const dropletSet&) = delete;


public:

    //- Runtime type information
    TypeName("dropletSet");


    // Constructors

        //- Construct for the given mesh
        explicit dropletSet(const fvMesh& mesh);

        //- Return the dropletSet stored in the mesh registry, creating it
        //  on first use
        static dropletSet& New(const fvMesh& mesh);


    //- Destructor
    virtual ~dropletSet() = default;


    // Member Functions

        // Update

            //- Identify the droplets and compute kinematics and drag
            //  with the current alpha.water and U fields
            void update();

            //- Store adhesion force and drag-to-adhesion ratio of droplet i
            void setAdhesion(const label i, const scalar Fadh, const scalar S);

            //- Store the adhesion state of droplet i
            void setAdhesionState
            (
                const label i,
                const bool inAdhesion,
                const scalar deceleration
            );

            //- Append the forces to ForcesOnDroplet_<i>.txt (output times)
            void writeForces();

            //- Pin (if at rest) or decelerate (if moving) the droplets in
            //  adhesion; called by the solver at the end of the time step
            void correctVelocity
            (
                volVectorField& U,
                const volScalarField& alpha
            ) const;


        // Access

            //- Number of droplets
            label size() const
            {
                return label(droplets_.size());
            }

            //- Droplet i
            const droplet& operator[](const label i) const
            {
                return droplets_[i];
            }

            //- Index of the droplet with the given label (-1 if none)
            label index(const scalar id) const
            {
                return labelToIndex_.found(label(id))
                  ? labelToIndex_[label(id)]
                  : -1;
            }

            const volScalarField& vofID() const
            {
                return vofID_;
            }

            const volVectorField& interfaceNormal() const
            {
                return nInterface_;
            }

            scalar alphaMin() const
            {
                return alphaMin_;
            }

            scalar sigma() const
            {
                return sigma_;
            }

            scalar muWater() const
            {
                return muWater_;
            }


        // I/O

            //- Nothing is written through the registry (vofID is written
            //  as a field at output times)
            virtual bool writeData(Ostream&) const
            {
                return true;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
