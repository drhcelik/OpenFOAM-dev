/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

Class
    Foam::chemistryReductionMethods::DRG

Description
    Implementation of the Directed Relation Graph (DRG) method

SourceFiles
    DRG.C

\*---------------------------------------------------------------------------*/

#ifndef DRG_H
#define DRG_H

#include "chemistryReductionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryReductionMethods
{

/*---------------------------------------------------------------------------*\
                             Class DRG Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class DRG
:
    public chemistryReductionMethod<ThermoType>
{
    // Private Data

        //- List of label for the search initiating set
        labelList searchInitSet_;

public:

    //- Runtime type information
    TypeName("DRG");


    // Constructors

        //- Construct from components
        DRG
        (
            const dictionary& dict,
            chemistryModel<ThermoType>& chemistry
        );


    // Destructor
    virtual ~DRG();


    // Member Functions

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            List<label>& ctos,
            DynamicList<label>& stoc,
            const label li
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryReductionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DRG.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
