/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "Kusters.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace breakupModels
{
    defineTypeNameAndDebug(Kusters, 0);
    addToRunTimeSelectionTable
    (
        breakupModel,
        Kusters,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::breakupModels::Kusters::Kusters
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    breakupModel(popBal, dict),
    B_("B", dimensionSet(0, 3, -3, 0, 0), dict),
    dP_("dP", dimLength, dict),
    kc_("kc", dimless, dict, 1),
    Df_("Df", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalance::breakupModels::Kusters::setBreakupRate
(
    volScalarField::Internal& breakupRate,
    const label i
)
{
    using Foam::constant::mathematical::pi;

    tmp<volScalarField> tdi = popBal_.d(i);
    const volScalarField::Internal& di = tdi();

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volScalarField::Internal& epsilonc = tepsilonc();
    tmp<volScalarField> tnu(popBal_.continuousPhase().fluidThermo().nu());
    const volScalarField::Internal nuc = tnu();

    breakupRate =
        sqrt(4*epsilonc/(15*pi*nuc))
       *exp(- B_/(dP_*0.5*pow(pow(di/dP_, Df_)/kc_, 1/Df_))/epsilonc);
}


// ************************************************************************* //
