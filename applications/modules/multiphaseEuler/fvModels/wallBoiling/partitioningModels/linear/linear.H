/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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
    Foam::wallBoilingModels::partitioningModels::linear

Description
    Linear wall heat flux partitioning model.

    Proposed threshold liquid fractions:
      - alphaLiquid1 0.1
      - alphaLiquid0 0.05

    \verbatim
        Ioilev, A., Samigulin, M., Ustinenko (2007).
        Advances in the modeling of cladding heat transfer
        and critical heat flux in boiling water reactor fuel assemblies.
        In Proc. 12th International Topical Meeting on
        Nuclear Reactor Thermal Hydraulics (NURETH-12),
        Pittsburgh, Pennsylvania, USA.
    \endverbatim

SourceFiles
    linear.C

\*---------------------------------------------------------------------------*/

#ifndef linear_H
#define linear_H

#include "partitioningModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{

/*---------------------------------------------------------------------------*\
                           Class linear Declaration
\*---------------------------------------------------------------------------*/

class linear
:
    public partitioningModel
{
    // Private Data

        //- Lower threshold liquid phase fraction
        const scalar alphaLiquid0_;

        //- Upper threshold liquid phase fraction
        const scalar alphaLiquid1_;


    // Private Member Functions

        //- Calculate and return the wet fraction
        template<class ScalarFieldType>
        tmp<ScalarFieldType> calculate
        (
            const ScalarFieldType& alphaLiquid
        ) const;


public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from a dictionary
        linear(const dictionary& dict);

        //- Copy construct
        linear(const linear& model);

        //- Construct and return a clone
        virtual autoPtr<partitioningModel> clone() const
        {
            return autoPtr<partitioningModel>(new linear(*this));
        }


    //- Destructor
    virtual ~linear();


    // Member Functions

        //- Return the wet fraction
        virtual tmp<scalarField> wetFraction
        (
            const scalarField& alphaLiquid
        ) const;

        //- Return the wet fraction
        virtual tmp<volScalarField::Internal> wetFraction
        (
            const volScalarField::Internal& alphaLiquid
        ) const;

        //- Return the wet fraction
        virtual tmp<volScalarField> wetFraction
        (
            const volScalarField& alphaLiquid
        ) const;

        //- Write to stream
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace partitioningModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
