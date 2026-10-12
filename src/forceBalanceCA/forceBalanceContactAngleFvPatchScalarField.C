/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "forceBalanceContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include <vector>
#include <iomanip>
#include <fstream>
#include "globalIndex.H"

#include "./utilities/vector2stringHelper.H"
#include "./utilities/interfaceArea/findInterfacePlane.H"
#include "./utilities/interfaceArea/findIntersection.H"
#include "./utilities/interfaceArea/intersectionArea.H"

// * * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::forceBalanceContactAngleFvPatchScalarField::
forceBalanceContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    theta0_(0.0),
    thetaA_(0.0),
    thetaR_(0.0),
    kk(0.0)
{}


Foam::forceBalanceContactAngleFvPatchScalarField::
forceBalanceContactAngleFvPatchScalarField
(
    const forceBalanceContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
    kk(gcpsf.kk)
{}


Foam::forceBalanceContactAngleFvPatchScalarField::
forceBalanceContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleFvPatchScalarField(p, iF, dict),
    theta0_(readScalar(dict.lookup("theta0"))),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    thetaR_(readScalar(dict.lookup("thetaR"))),
    kk(readScalar(dict.lookup("kk")))
{
    evaluate();
}


Foam::forceBalanceContactAngleFvPatchScalarField::
forceBalanceContactAngleFvPatchScalarField
(
    const forceBalanceContactAngleFvPatchScalarField& gcpsf
)
:
    alphaContactAngleFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
    kk(gcpsf.kk)
{}


Foam::forceBalanceContactAngleFvPatchScalarField::
forceBalanceContactAngleFvPatchScalarField
(
    const forceBalanceContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_),
    kk(gcpsf.kk)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::forceBalanceContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
	#include "readProperties.H"						// Read transportProperties file and check for reading error
	#include "findDrop.H"							// Find and label different droplets in the domain
    #include "dropG.H"								// Evaluate center of mass of water droplets
	#include "dropVel.H"							// Evaluate mean velocity of water droplets
	#include "printGpositionVelocity.H"			// Print center of mass position and velocity over time in a txt file (for each droplet)
	#include "forcesOnDroplet.H"					// Evaluate sum of forces acting on each droplet
	#include "printForces.H"						// Print forces acting on each droplet over time in a txt file (for each droplet)
	#include "lateralScaleFactor.H"					// Scale the value of contact angle for lateral wall of each droplet
	#include "contactAngleHysteresis.H"				// Evaluate Contact Angle Hysteresis (linear with drag/adhesion ratio)
	#include "contactAngleCorrelations.H"			// Evaluate dynamic contact angle (based on experimental correlations)
}


void Foam::forceBalanceContactAngleFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleFvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;
    os.writeKeyword("kk") << kk << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        forceBalanceContactAngleFvPatchScalarField
    );
}
