/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "shockFluid.H"
#include "fvMeshStitcher.H"
#include "localEulerDdtScheme.H"
#include "hydrostaticInitialisation.H"
#include "fvcMeshPhi.H"
#include "fvcVolumeIntegrate.H"
#include "fvcReconstruct.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(shockFluid, 0);
    addToRunTimeSelectionTable(solver, shockFluid, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::shockFluid::correctCoNum(const surfaceScalarField& amaxSf)
{
    const scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

    CoNum_ =
        0.5*gMax(sumAmaxSf/mesh.V().primitiveField())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5
       *(gSum(sumAmaxSf)/gSum(mesh.V().primitiveField()))
       *runTime.deltaTValue();

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::solvers::shockFluid::clearTemporaryFields()
{
    rho_pos.clear();
    rho_neg.clear();

    rhoU_pos.clear();
    rhoU_neg.clear();

    U_pos.clear();
    U_neg.clear();

    p_pos.clear();
    p_neg.clear();

    a_pos.clear();
    a_neg.clear();

    aSf.clear();

    aphiv_pos.clear();
    aphiv_neg.clear();

    devTau.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::shockFluid::shockFluid(fvMesh& mesh)
:
    fluidSolver(mesh),

    thermoPtr_(psiThermo::New(mesh)),

    thermo_(thermoPtr_()),

    p_(thermo_.p()),

    rho_
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo_.renameRho()
    ),

    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    phi_
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho_*U_) & mesh.Sf()
    ),

    K("K", 0.5*magSqr(U_)),

    inviscid
    (
        max(thermo_.mu().primitiveField()) > 0
      ? false
      : true
    ),

    momentumTransport
    (
        inviscid
      ? autoPtr<compressibleMomentumTransportModel>(nullptr)
      : compressible::momentumTransportModel::New
        (
            rho_,
            U_,
            phi_,
            thermo_
        )
    ),

    thermophysicalTransport
    (
        inviscid
      ? autoPtr<fluidThermoThermophysicalTransportModel>(nullptr)
      : fluidThermoThermophysicalTransportModel::New
        (
            momentumTransport(),
            thermo_
        )
    ),

    fluxScheme
    (
        mesh.schemes().dict().lookupOrDefault<word>("fluxScheme", "Kurganov")
    ),

    thermo(thermo_),
    p(p_),
    rho(rho_),
    U(U_),
    phi(phi_)
{
    thermo.validate(type(), "e");

    if (momentumTransport.valid())
    {
        momentumTransport->validate();
        mesh.schemes().setFluxRequired(U.name());
    }

    fluxPredictor();

    if (transient())
    {
        const surfaceScalarField amaxSf
        (
            max(mag(aphiv_pos()), mag(aphiv_neg()))
        );

        correctCoNum(amaxSf);
    }
    else if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::shockFluid::~shockFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockFluid::preSolve()
{
    {
        const surfaceScalarField amaxSf
        (
            max(mag(aphiv_pos()), mag(aphiv_neg()))
        );

        if (transient())
        {
            correctCoNum(amaxSf);
        }
        else if (LTS)
        {
            setRDeltaT(amaxSf);
        }
    }

    fvModels().preUpdateMesh();

    if (mesh.topoChanging() || mesh.stitcher().stitches())
    {
        pos.clear();
        neg.clear();

        clearTemporaryFields();
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::shockFluid::prePredictor()
{
    fluxPredictor();
    correctDensity();
}


void Foam::solvers::shockFluid::momentumTransportPredictor()
{
    if (!inviscid)
    {
        momentumTransport->predict();
    }
}


void Foam::solvers::shockFluid::thermophysicalTransportPredictor()
{
    if (!inviscid)
    {
        thermophysicalTransport->predict();
    }
}


void Foam::solvers::shockFluid::momentumTransportCorrector()
{
    if (!inviscid)
    {
        momentumTransport->correct();
    }
}


void Foam::solvers::shockFluid::thermophysicalTransportCorrector()
{
    if (!inviscid)
    {
        thermophysicalTransport->correct();
    }
}


void Foam::solvers::shockFluid::postSolve()
{}


// ************************************************************************* //
