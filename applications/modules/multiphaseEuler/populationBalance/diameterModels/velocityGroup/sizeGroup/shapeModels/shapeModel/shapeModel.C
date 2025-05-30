/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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

#include "shapeModel.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(shapeModel, 0);
    defineRunTimeSelectionTable(shapeModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModel::shapeModel(const sizeGroup& group)
:
    group_(group)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::shapeModel>
Foam::diameterModels::shapeModel::New
(
    const dictionary& dict,
    const sizeGroup& group,
    const dictionary& groupDict
)
{
    word shapeModelType(dict.lookup("shapeModel"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(shapeModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown shapeModel type "
            << shapeModelType << endl << endl
            << "Valid shapeModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, group, groupDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::shapeModel::~shapeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::diameterModels::shapeModel&
Foam::diameterModels::shapeModel::model(const sizeGroup& group)
{
    return group.group().popBal().sizeGroups()[group.i()].shape();
}


const Foam::diameterModels::sizeGroup&
Foam::diameterModels::shapeModel::group() const
{
    return group_;
}


void Foam::diameterModels::shapeModel::correct()
{}


void Foam::diameterModels::shapeModel::addCoalescence
(
    const volScalarField::Internal& Su,
    const sizeGroup& fj,
    const sizeGroup& fk
)
{}


void Foam::diameterModels::shapeModel::addBreakup
(
    const volScalarField::Internal& Su,
    const sizeGroup& fj
)
{}


void Foam::diameterModels::shapeModel::reset()
{}


// ************************************************************************* //
