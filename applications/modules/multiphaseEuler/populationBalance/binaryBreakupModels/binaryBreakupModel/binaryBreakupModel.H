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
    Foam::populationBalance::binaryBreakupModel

Description
    Base class for binary breakup models that provide a breakup rate between a
    size class pair directly, i.e. without explicitly stating the daughter size
    distribution function.

See also
    Foam::populationBalanceModel

SourceFiles
    binaryBreakupModel.C

\*---------------------------------------------------------------------------*/

#ifndef binaryBreakupModel_H
#define binaryBreakupModel_H

#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{

/*---------------------------------------------------------------------------*\
                     Class binaryBreakupModel Declaration
\*---------------------------------------------------------------------------*/

class binaryBreakupModel
{
protected:

    // Protected Data

        //- Reference to the populationBalanceModel
        const populationBalanceModel& popBal_;


public:

    //- Runtime type information
    TypeName("binaryBreakupModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            binaryBreakupModel,
            dictionary,
            (
                const populationBalanceModel& popBal,
                const dictionary& dict
            ),
            (popBal, dict)
        );


    //- Class used for the read-construction of
    //  PtrLists of binary breakup models
    class iNew
    {
        const populationBalanceModel& popBal_;

    public:

        iNew(const populationBalanceModel& popBal)
        :
            popBal_(popBal)
        {}

        autoPtr<binaryBreakupModel> operator()(Istream& is) const
        {
            word type(is);
            dictionary dict(is);
            return binaryBreakupModel::New(type, popBal_, dict);
        }
    };


    // Constructor

        binaryBreakupModel
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );

        autoPtr<binaryBreakupModel> clone() const
        {
            NotImplemented;
            return autoPtr<binaryBreakupModel>(nullptr);
        }


    // Selector

        static autoPtr<binaryBreakupModel> New
        (
            const word& type,
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~binaryBreakupModel()
    {}


    // Member Functions

        //- Precompute diameter independent expressions
        virtual void precompute();

        //- Add to binary breakupRate
        virtual void addToBinaryBreakupRate
        (
            volScalarField::Internal& binaryBreakupRate,
            const label i,
            const label j
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
