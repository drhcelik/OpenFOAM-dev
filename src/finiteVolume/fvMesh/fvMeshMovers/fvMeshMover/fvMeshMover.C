/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "fvMeshMover.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshMover, 0);
    defineRunTimeSelectionTable(fvMeshMover, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMover::fvMeshMover(fvMesh& mesh)
:
    mesh_(mesh)
{}


Foam::fvMeshMover::velocityMotionCorrection::velocityMotionCorrection
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    velocityFields_(dict.lookupOrDefault("velocityFields", wordList()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMover::~fvMeshMover()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshMover::velocityMotionCorrection::update() const
{
    forAll(velocityFields_, i)
    {
        if (mesh_.foundObject<volVectorField>(velocityFields_[i]))
        {
            mesh_.lookupObjectRef<volVectorField>
            (
                velocityFields_[i]
            ).correctBoundaryConditions();
        }
    }
}

// ************************************************************************* //
