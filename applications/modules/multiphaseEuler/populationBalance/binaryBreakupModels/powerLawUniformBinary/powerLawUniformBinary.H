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
    Foam::populationBalance::binaryBreakupModels::powerLawUniformBinary

Description
    Powerlaw kernel with a uniform daughter size distribution.

SourceFiles
    powerLawUniformBinary.C

\*---------------------------------------------------------------------------*/

#ifndef powerLawUniformBinary_H
#define powerLawUniformBinary_H

#include "binaryBreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace binaryBreakupModels
{

/*---------------------------------------------------------------------------*\
                    Class powerLawUniformBinary Declaration
\*---------------------------------------------------------------------------*/

class powerLawUniformBinary
:
    public binaryBreakupModel
{
    // Private Data

        //- Power
        scalar power_;


public:

    //- Runtime type information
    TypeName("powerLawUniformBinary");

    // Constructor

        powerLawUniformBinary
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~powerLawUniformBinary()
    {}


    // Member Functions

        //- Add to binary breakupRate
        virtual void addToBinaryBreakupRate
        (
            volScalarField::Internal& binaryBreakupRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryBreakupModels
} // End namespace populationBalance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
