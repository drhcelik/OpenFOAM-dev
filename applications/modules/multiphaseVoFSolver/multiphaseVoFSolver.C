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

#include "multiphaseVoFSolver.H"
#include "localEulerDdtScheme.H"
#include "fvcAverage.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(multiphaseVoFSolver, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::multiphaseVoFSolver::correctCoNum()
{
    VoFSolver::correctCoNum();

    const scalarField sumPhi
    (
        mixture.nearInterface()().primitiveField()
       *fvc::surfaceSum(mag(phi))().primitiveField()
    );

    alphaCoNum =
        0.5*gMax(sumPhi/mesh.V().primitiveField())*runTime.deltaTValue();

    const scalar meanAlphaCoNum =
        0.5
       *(gSum(sumPhi)/gSum(mesh.V().primitiveField()))
       *runTime.deltaTValue();

    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
        << " max: " << alphaCoNum << endl;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::multiphaseVoFSolver::read()
{
    VoFSolver::read();

    const dictionary& alphaControls = mesh.solution().solverDict("alpha");

    cAlpha = alphaControls.lookup<scalar>("cAlpha");

    nAlphaSubCyclesPtr =
        Function1<scalar>::New
        (
            alphaControls.found("nAlphaSubCycles")
          ? "nAlphaSubCycles"
          : "nSubCycles",
            dimless,
            dimless,
            alphaControls
        );

    MULEScontrols.read(alphaControls);

    return true;
}


void Foam::solvers::multiphaseVoFSolver::correctInterface()
{}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::multiphaseVoFSolver::surfaceTensionForce() const
{
    return mixture.surfaceTensionForce(U);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::multiphaseVoFSolver::multiphaseVoFSolver
(
    fvMesh& mesh,
    autoPtr<multiphaseVoFMixture> mixturePtr
)
:
    VoFSolver(mesh, autoPtr<VoFMixture>(mixturePtr.ptr())),

    mixture(refCast<multiphaseVoFMixture>(VoFSolver::mixture_)),

    phases(mixture.phases())
{
    read();

    if (transient())
    {
        correctCoNum();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::multiphaseVoFSolver::~multiphaseVoFSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multiphaseVoFSolver::preSolve()
{
    VoFSolver::preSolve();
}


void Foam::solvers::multiphaseVoFSolver::prePredictor()
{
    alphaPredictor();
}


// ************************************************************************* //
