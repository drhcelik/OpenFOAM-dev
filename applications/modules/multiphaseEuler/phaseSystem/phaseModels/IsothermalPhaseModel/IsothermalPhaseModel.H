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
    Foam::IsothermalPhaseModel

Description
    Class which represents a phase for which the temperature remains
    constant. Returns an empty energy equation and updates the energy
    corresponding to pressure changes only when correctThermo is called.

SourceFiles
    IsothermalPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef IsothermalPhaseModel_H
#define IsothermalPhaseModel_H

#include "ThermophysicalTransportPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class IsothermalPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class IsothermalPhaseModel
:
    public BasePhaseModel
{
public:

    // Constructors

        IsothermalPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~IsothermalPhaseModel();


    // Member Functions

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Return whether the phase is isothermal
        virtual bool isothermal() const;

        //- Return the effective thermal conductivity on a patch
        virtual tmp<scalarField> kappaEff(const label patchi) const;

        //- Return the enthalpy equation
        virtual tmp<fvScalarMatrix> heEqn();
};


/*---------------------------------------------------------------------------*\
        Typedef IsothermalThermophysicalTransportPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
using IsothermalThermophysicalTransportPhaseModel =
    IsothermalPhaseModel<ThermophysicalTransportPhaseModel<BasePhaseModel>>;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IsothermalPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
