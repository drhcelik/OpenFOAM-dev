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

#include "extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::extrudeModel> Foam::extrudeModel::New
(
    const dictionary& dict
)
{
    const dictionary& modelDict =
        dict.isDict("extrudeModel")
      ? dict.subDict("extrudeModel")
      : dict;

    const word modelType
    (
        modelDict.lookupBackwardsCompatible({"type", "extrudeModel"})
    );

    Info<< "Selecting extrudeModel " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown extrudeModel type "
            << modelType << nl << nl
            << "Valid extrudeModel types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return autoPtr<extrudeModel>(cstrIter()(modelDict));
}


// ************************************************************************* //
