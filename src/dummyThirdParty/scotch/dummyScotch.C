/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "scotch.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

static const char* notImplementedMessage =
"You are trying to use scotch but do not have the scotch library loaded."
"\nThis message is from the dummy scotch stub library instead.\n"
"\n"
"Please install scotch and make sure that libscotch.so is in your "
"LD_LIBRARY_PATH.\n"
"The scotch library can then be built in "
"$FOAM_SRC/parallel/decompose/decompositionMethods/scotch\n";


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(scotch, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotch,
        decomposer
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotch,
        distributor
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decompositionMethods::scotch::check
(
    const int retVal,
    const char* str
)
{}


Foam::label Foam::decompositionMethods::scotch::decompose
(
    const fileName& meshPath,
    const List<label>& adjncy,
    const List<label>& xadj,
    const scalarField& cWeights,

    List<label>& finalDecomp
)
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::scotch::scotch
(
    const dictionary& decompositionDict,
    const dictionary& methodDict
)
:
    decompositionMethod(decompositionDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::scotch::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::decompositionMethods::scotch::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
)
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::decompositionMethods::scotch::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


// ************************************************************************* //
