/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "ShihQuadraticKE.H"
#include "bound.H"
#include "wallFvPatch.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "makeMomentumTransportModel.H"

makeMomentumTransportModelTypes
(
    geometricOneField,
    geometricOneField,
    incompressibleMomentumTransportModel
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ShihQuadraticKE, 0);
addToRunTimeSelectionTable
(
    RASincompressibleMomentumTransportModel,
    ShihQuadraticKE,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void ShihQuadraticKE::boundEpsilon()
{
    epsilon_ = max(epsilon_, Cmu_*sqr(k_)/(nutMaxCoeff_*nu()));
}


void ShihQuadraticKE::correctNut()
{
    correctNonlinearStress(fvc::grad(U_));
}


void ShihQuadraticKE::correctNonlinearStress(const volTensorField& gradU)
{
    volSymmTensorField S(symm(gradU));
    volTensorField W(skew(gradU));

    volScalarField sBar((k_/epsilon_)*sqrt(2.0)*mag(S));
    volScalarField wBar((k_/epsilon_)*sqrt(2.0)*mag(W));

    volScalarField Cmu((2.0/3.0)/(Cmu1_ + sBar + Cmu2_*wBar));

    boundEpsilon();
    nut_ = Cmu*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    nonlinearStress_ =
        k_*sqr(k_/epsilon_)/(Cbeta_ + pow3(sBar))
       *(
           Cbeta1_*dev(innerSqr(S))
         + Cbeta2_*twoSymm(S&W)
         + Cbeta3_*dev(symm(W&W))
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ShihQuadraticKE::ShihQuadraticKE
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    nonlinearEddyViscosity<incompressible::RASModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    Ceps1_("Ceps1", coeffDict(), 1.44),
    Ceps2_("Ceps2", coeffDict(), 1.92),
    sigmak_("sigmak", coeffDict(), 1.0),
    sigmaEps_("sigmaEps", coeffDict(), 1.3),
    Cmu_("Cmu", coeffDict(), 0.09),
    Cmu1_("Cmu1", coeffDict(), 1.25),
    Cmu2_("Cmu2", coeffDict(), 0.9),
    Cbeta_("Cbeta", coeffDict(), 1000.0),
    Cbeta1_("Cbeta1", coeffDict(), 3.0),
    Cbeta2_("Cbeta2", coeffDict(), 15.0),
    Cbeta3_("Cbeta3", coeffDict(), -19.0),

    k_
    (
        IOobject
        (
            this->groupName("k"),
            runTime_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            this->groupName("epsilon"),
            runTime_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    bound(k_, kMin_);
    boundEpsilon();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool ShihQuadraticKE::read()
{
    if (nonlinearEddyViscosity<incompressible::RASModel>::read())
    {
        Ceps1_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        Cmu_.readIfPresent(coeffDict());
        Cmu1_.readIfPresent(coeffDict());
        Cmu2_.readIfPresent(coeffDict());
        Cbeta_.readIfPresent(coeffDict());
        Cbeta1_.readIfPresent(coeffDict());
        Cbeta2_.readIfPresent(coeffDict());
        Cbeta3_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void ShihQuadraticKE::correct()
{
    if (!turbulence_)
    {
        return;
    }

    nonlinearEddyViscosity<incompressible::RASModel>::correct();

    tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    volScalarField G
    (
        GName(),
        (nut_*twoSymm(gradU) - nonlinearStress_) && gradU
    );


    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        Ceps1_*G*epsilon_/k_
      - fvm::Sp(Ceps2_*epsilon_/k_, epsilon_)
    );

    epsEqn.ref().relax();
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    boundEpsilon();


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
      ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity and non-linear stress
    correctNonlinearStress(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
