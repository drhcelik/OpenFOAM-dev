/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
    Foam::InertPhaseModel

Description
    Class which represents an inert phase, with no reactions. Returns zero
    reaction rate and heat.

SourceFiles
    InertPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef InertPhaseModel_H
#define InertPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class InertPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class InertPhaseModel
:
    public BasePhaseModel
{
public:

    // Constructors

        InertPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~InertPhaseModel();


    // Member Functions

        //- Return the consumption rate for a specie
        virtual tmp<volScalarField::Internal> R(const label speciei) const;

        //- Return consumption rate matrix for a specie
        virtual tmp<fvScalarMatrix> R(volScalarField& Yi) const;

        //- Return the heat release rate
        virtual tmp<volScalarField> Qdot() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "InertPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
