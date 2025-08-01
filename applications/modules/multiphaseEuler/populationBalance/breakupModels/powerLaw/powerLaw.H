/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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
    Foam::populationBalance::breakupModels::powerLaw

Description
    Power-law breakup model

SourceFiles
    powerLaw.C

\*---------------------------------------------------------------------------*/

#ifndef powerLaw_H
#define powerLaw_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                          Class powerLaw Declaration
\*---------------------------------------------------------------------------*/

class powerLaw
:
    public breakupModel
{
    // Private Data

        //- Power
        scalar power_;


public:

    //- Runtime type information
    TypeName("powerLaw");


    // Constructor

        powerLaw
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~powerLaw()
    {}


    // Member Functions

        //- Set total breakupRate
        virtual void setBreakupRate
        (
            volScalarField::Internal& breakupRate,
            const label i
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupModels
} // End namespace populationBalance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
