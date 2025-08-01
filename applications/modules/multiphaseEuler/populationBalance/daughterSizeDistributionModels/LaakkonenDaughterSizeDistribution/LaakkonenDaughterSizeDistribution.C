/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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

#include "LaakkonenDaughterSizeDistribution.H"
#include "addToRunTimeSelectionTable.H"
#include "breakupModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace daughterSizeDistributionModels
{
    defineTypeNameAndDebug(LaakkonenDaughterSizeDistribution, 0);
    addToRunTimeSelectionTable
    (
        daughterSizeDistributionModel,
        LaakkonenDaughterSizeDistribution,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::daughterSizeDistributionModels::
LaakkonenDaughterSizeDistribution::LaakkonenDaughterSizeDistribution
(
    const breakupModel& breakup,
    const dictionary& dict
)
:
    daughterSizeDistributionModel(breakup, dict),
    C4_("C4", dimless, dict, 18.25)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalance::daughterSizeDistributionModels::
LaakkonenDaughterSizeDistribution::~LaakkonenDaughterSizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::populationBalance::daughterSizeDistributionModels::
LaakkonenDaughterSizeDistribution::antiderivative
(
    const dimensionedScalar& xk,
    const dimensionedScalar& v,
    const dimensionedScalar& bndr,
    const dimensionedScalar range
) const
{
    return
        (4.0/3.0 + C4_/3)
       *(
            pow(xk, -C4_ - 3)*pow(xk - v, C4_)*(v - xk)
           *(
                (C4_ + 1)*(C4_ + 2)*(C4_ + 3)*pow3(v)
              - (C4_ + 1)*(C4_ + 2)*(bndr*(C4_ + 4) - 3*xk)*sqr(v)
              - 2*v*xk*(C4_ + 1)*(bndr*(C4_ + 4) - 3*xk)
              - 2*bndr*C4_*sqr(xk)
              + 6*pow3(xk)
              - 8*bndr*sqr(xk)
            )
        )/(2*range*(C4_ + 4));
}


Foam::dimensionedScalar
Foam::populationBalance::daughterSizeDistributionModels::
LaakkonenDaughterSizeDistribution::calcNik
(
    const label i,
    const label k
) const
{
    const populationBalanceModel& popBal = breakup_.popBal();

    const dimensionedScalar& x0 = popBal.v(0);
    const dimensionedScalar xi = popBal.v(i) - x0;
    const dimensionedScalar xk = popBal.v(k) - x0;

    if (i == 0)
    {
        const dimensionedScalar xii = popBal.v(i+1) - x0;

        if (k == 0)
        {
            return 1;
        }

        return
            antiderivative(xk, xi, xii, xii - xi)
          - antiderivative(xk, xii, xii, xii - xi);
    }
    else if (i == k)
    {
        const dimensionedScalar x = popBal.v(i-1) - x0;

        return
            antiderivative(xk, xi, x, xi - x)
          - antiderivative(xk, x, x, xi - x);
    }
    else
    {
        const dimensionedScalar x = popBal.v(i-1) - x0;
        const dimensionedScalar xii = popBal.v(i+1) - x0;

        return
            antiderivative(xk, xi, xii, xii - xi)
          - antiderivative(xk, xii, xii, xii - xi)
          + antiderivative(xk, xi, x, xi - x)
          - antiderivative(xk, x, x, xi - x);
    }
}


// ************************************************************************* //
