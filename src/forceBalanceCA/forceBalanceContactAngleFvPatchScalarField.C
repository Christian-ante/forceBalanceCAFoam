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
#include "dropletSet.H"
#include <cmath>
#include <vector>

// * * * * * * * * * * * * * * Local constants * * * * * * * * * * * * * * //

namespace
{
    constexpr double deg2rad = M_PI/180.0;     // Degrees -> radians
    constexpr double rad2deg = 180.0/M_PI;     // Radians -> degrees
}

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
    kk(0.0),
    thetaDynOld_(p.size(), 0.0),
    dropAlreadySliding_()
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
    kk(gcpsf.kk),
    thetaDynOld_(p.size(), gcpsf.theta0_),
    dropAlreadySliding_(gcpsf.dropAlreadySliding_)
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
    kk(readScalar(dict.lookup("kk"))),
    thetaDynOld_(p.size(), theta0_),
    dropAlreadySliding_()
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
    kk(gcpsf.kk),
    thetaDynOld_(gcpsf.thetaDynOld_),
    dropAlreadySliding_(gcpsf.dropAlreadySliding_)
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
    kk(gcpsf.kk),
    thetaDynOld_(gcpsf.thetaDynOld_),
    dropAlreadySliding_(gcpsf.dropAlreadySliding_)
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::forceBalanceContactAngleFvPatchScalarField::checkSinglePatch() const
{
    static bool checked = false;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Not during construction: the boundary field may be incomplete
    if (checked || mesh.time().timeIndex() <= mesh.time().startTimeIndex())
    {
        return;
    }

    const volScalarField& alpha =
        mesh.lookupObject<volScalarField>("alpha.water");

    label nPatches = 0;
    forAll(alpha.boundaryField(), patchi)
    {
        if (alpha.boundaryField()[patchi].type() == typeName)
        {
            nPatches++;
        }
    }

    if (nPatches > 1)
    {
        WarningInFunction
            << typeName << " is applied to " << nPatches << " patches." << nl
            << "    Adhesion forces and droplet states are stored once per"
            << " droplet, so the last" << nl
            << "    patch evaluated overwrites the others: merge the walls into"
            << " a single patch" << nl
            << "    (e.g. with createPatch)." << endl;
    }

    checked = true;
}


void Foam::forceBalanceContactAngleFvPatchScalarField::computeAdhesion
(
    dropletSet& droplets
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const volScalarField& alpha =
        mesh.lookupObject<volScalarField>("alpha.water");
    const volScalarField& vofID = droplets.vofID();
    const labelUList& faceCells = patch().faceCells();

    // Wetted area of each droplet on this patch (cubic cells)
    List<scalar> wettedArea(droplets.size(), Zero);
    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];
        if (alpha[celli] > droplets.alphaMin())
        {
            const label i = droplets.index(vofID[celli]);
            if (i >= 0)
            {
                const scalar l = cbrt(mesh.V()[celli]);
                wettedArea[i] += l*l;
            }
        }
    }

    for (label i = 0; i < droplets.size(); i++)
    {
        reduce(wettedArea[i], sumOp<scalar>());

        // ElSherbini & Jacobi (2006): maximum adhesion force, with the
        // diameter of the circle having the same area as the wetted area
        const scalar Dadh = std::sqrt(4*wettedArea[i]/M_PI);
        const scalar Fadh =
            kk*24.0*droplets.sigma()*Dadh
           *(cos(thetaR_*deg2rad) - cos(thetaA_*deg2rad))/std::pow(M_PI, 3);

        // Drag-to-adhesion ratio, including the weight along x
        const droplet& d = droplets[i];
        const scalar S =
            Fadh != 0
          ? std::abs((d.dragForce[0] + d.gravityForceX)/Fadh)
          : 999;        // Droplet not in contact with the wall

        droplets.setAdhesion(i, Fadh, S);
    }
}


void Foam::forceBalanceContactAngleFvPatchScalarField::lateralScaleFactor
(
    const dropletSet& droplets,
    const labelList& faceDroplet,
    scalarField& factor,
    scalarField& uwall
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const volVectorField& nInter = droplets.interfaceNormal();
    const labelUList& faceCells = patch().faceCells();

    if (mesh.time().value() != 0)       // Skip the first time step
    {
        for (label i = 0; i < droplets.size(); i++)
        {
            // Max/min contact angle applied to the droplet at the previous
            // evaluation
            double maxTheta = 0, minTheta = 180;

            forAll(faceCells, facei)
            {
                if (faceDroplet[facei] == i)
                {
                    if (thetaDynOld_[facei] > maxTheta)
                    {
                        maxTheta = thetaDynOld_[facei];
                    }
                    else if (thetaDynOld_[facei] < minTheta)
                    {
                        minTheta = thetaDynOld_[facei];
                    }
                }
            }
            reduce(maxTheta, maxOp<double>());
            reduce(minTheta, minOp<double>());

            // Max (front) and min (rear) x-projection of the interface normal
            const double nxMax = sin(maxTheta*deg2rad);
            const double nxMin = sin(minTheta*deg2rad);

            forAll(faceCells, facei)
            {
                if (faceDroplet[facei] == i)
                {
                    const scalar nx = nInter[faceCells[facei]][0];

                    if (nx <= 0.0)      // Advancing front
                    {
                        factor[facei] = std::abs(nx/nxMax);
                        uwall[facei] = -droplets[i].velocity[0];
                    }
                    else                // Receding front
                    {
                        factor[facei] = std::abs(nx/nxMin);
                        uwall[facei] = droplets[i].velocity[0];
                    }

                    if (factor[facei] > 1)      // Limit the factor to 1
                    {
                        factor[facei] = 1;
                    }
                }
            }
        }
    }

    // Reset for the next evaluation: the droplet faces are overwritten when
    // their contact angle is evaluated
    thetaDynOld_ = theta0_;
}


void Foam::forceBalanceContactAngleFvPatchScalarField::contactAngleHysteresis
(
    dropletSet& droplets,
    const labelList& faceDroplet,
    const scalarField& factor,
    const scalarField& uwall,
    scalarField& thetaDyn
) const
{
    for (label i = 0; i < droplets.size(); i++)
    {
        const droplet& d = droplets[i];

        if (d.ratio > 1)
        {
            droplets.setAdhesionState(i, false, 0);
            continue;
        }

        // Adhesion: pin the droplet if it was at rest, otherwise decelerate it
        const scalar deceleration =
            dropAlreadySliding_[i]
          ? std::abs(d.dragForce[0] - d.adhesionForce)/d.mass
          : 0;
        droplets.setAdhesionState(i, true, deceleration);

        forAll(thetaDyn, facei)
        {
            if (faceDroplet[facei] == i)
            {
                if (uwall[facei] >= 0)  // Receding front
                {
                    thetaDyn[facei] =
                        theta0_ - d.ratio*(theta0_ - thetaR_)*factor[facei];
                }
                else                    // Advancing front
                {
                    thetaDyn[facei] =
                        theta0_ + d.ratio*(thetaA_ - theta0_)*factor[facei];
                }

                thetaDynOld_[facei] = thetaDyn[facei];
            }
        }
    }
}


void Foam::forceBalanceContactAngleFvPatchScalarField::contactAngleCorrelations
(
    const dropletSet& droplets,
    const labelList& faceDroplet,
    const scalarField& factor,
    const scalarField& uwall,
    scalarField& thetaDyn
) const
{
    const scalarField Ca(droplets.muWater()*uwall/droplets.sigma());

    for (label i = 0; i < droplets.size(); i++)
    {
        if (droplets[i].ratio <= 1)
        {
            continue;
        }

        forAll(thetaDyn, facei)
        {
            if (faceDroplet[facei] == i)
            {
                if (uwall[facei] >= 0)
                {
                    // Receding front: Nichita et al. (2010)
                    const double thetaRad = thetaR_*deg2rad;
                    const double innerValue =
                        std::cbrt(std::pow(thetaRad, 3) - 72.0*Ca[facei]);
                    const double delta = theta0_ - innerValue*rad2deg;
                    thetaDyn[facei] = theta0_ - delta*factor[facei];
                }
                else
                {
                    // Advancing front: Kistler (1993), Hoffman function
                    const double thetaRad = thetaA_*deg2rad;
                    const double fHoffInv = std::pow(thetaRad, 3)/72;
                    const double xHoff = Ca[facei] + fHoffInv;
                    const double temp = acos
                    (
                        1.0 - 2.0*tanh
                        (
                            5.16*pow
                            (
                                xHoff/(1.0 + 1.31*pow(xHoff, 0.99)),
                                0.706
                            )
                        )
                    );
                    const double delta = temp*rad2deg - theta0_;
                    thetaDyn[facei] = theta0_ + delta*factor[facei];
                }

                thetaDynOld_[facei] = thetaDyn[facei];
            }
        }

        dropAlreadySliding_[i] = true;      // The droplet is sliding
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::forceBalanceContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    checkSinglePatch();

    // Shared droplet data: identification, kinematics and drag
    dropletSet& droplets = dropletSet::New(mesh);
    droplets.update();

    // Adhesion on this patch and drag-to-adhesion ratio
    computeAdhesion(droplets);
    droplets.writeForces();

    // Per-patch state sized on the patch and on the number of droplets
    if (thetaDynOld_.size() != patch().size())
    {
        thetaDynOld_.setSize(patch().size(), theta0_);
    }
    if (label(dropAlreadySliding_.size()) < droplets.size())
    {
        dropAlreadySliding_.resize(droplets.size(), false);
    }

    // Droplet wetting each face of the patch (-1 = none)
    const labelUList& faceCells = patch().faceCells();
    labelList faceDroplet(patch().size(), -1);
    forAll(faceCells, facei)
    {
        faceDroplet[facei] = droplets.index(droplets.vofID()[faceCells[facei]]);
    }

    scalarField factor(patch().size(), 0.0);
    scalarField uwall(patch().size(), 0.0);
    lateralScaleFactor(droplets, faceDroplet, factor, uwall);

    tmp<scalarField> tthetaDyn(new scalarField(patch().size(), theta0_));
    scalarField& thetaDyn = tthetaDyn.ref();

    contactAngleHysteresis(droplets, faceDroplet, factor, uwall, thetaDyn);
    contactAngleCorrelations(droplets, faceDroplet, factor, uwall, thetaDyn);

    return tthetaDyn;
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
