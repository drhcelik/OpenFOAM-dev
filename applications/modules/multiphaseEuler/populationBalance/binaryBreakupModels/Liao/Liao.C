/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "Liao.H"
#include "fvcGrad.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{
    defineTypeNameAndDebug(Liao, 0);
    addToRunTimeSelectionTable
    (
        binaryBreakupModel,
        Liao,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::binaryBreakupModels::Liao::Liao
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binaryBreakupModel(popBal, dict),
    LiaoBase(popBal, dict),
    BTurb_("BTurb", dimless, dict, 1),
    BShear_("BShear", dimless, dict, 1),
    BEddy_("BEddy", dimless, dict, 1),
    BFric_("BFric", dimless, dict, 0.25),
    turbulence_(dict.lookup("turbulence")),
    laminarShear_(dict.lookup("laminarShear")),
    turbulentShear_(dict.lookup("turbulentShear")),
    interfacialFriction_(dict.lookup("interfacialFriction"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::binaryBreakupModels::Liao::precompute()
{
    LiaoBase::precompute();
}


void Foam::diameterModels::binaryBreakupModels::Liao::addToBinaryBreakupRate
(
    volScalarField::Internal& binaryBreakupRate,
    const label i,
    const label j
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(fi.phase()));
    const volScalarField::Internal& sigma = tsigma();

    tmp<volScalarField> tmuc(popBal_.continuousPhase().fluidThermo().mu());
    const volScalarField::Internal& muc = tmuc();

    const dimensionedScalar dk(cbrt(pow3(fj.dSph()) - pow3(fi.dSph())));

    const volScalarField::Internal tauCrit1
    (
        6*sigma/fj.dSph()*(sqr(fi.dSph()/fj.dSph()) + sqr(dk/fj.dSph()) - 1)
    );

    const volScalarField::Internal tauCrit2
    (
        sigma/min(dk, fi.dSph())
    );

    const volScalarField::Internal tauCrit(max(tauCrit1, tauCrit2));

    if (turbulence_)
    {
        tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
        const volScalarField::Internal& epsilonc = tepsilonc();

        const volScalarField::Internal tauTurb
        (
            pos(fj.dSph() - kolmogorovLengthScale_)*BTurb_*rhoc
           *sqr(cbrt(epsilonc*fj.dSph()))
        );

        binaryBreakupRate +=
            pos(tauTurb - tauCrit)
           /fj.dSph()
           *sqrt(mag(tauTurb - tauCrit)/rhoc)
           /fj.x();
    }

    if (laminarShear_)
    {
        const volScalarField::Internal tauShear
        (
            BShear_*muc*shearStrainRate_
        );

        binaryBreakupRate +=
            pos(tauShear - tauCrit)
           /fj.dSph()
           *sqrt(mag(tauShear - tauCrit)/rhoc)
           /fj.x();
    }

    if (turbulentShear_)
    {
        const volScalarField::Internal tauEddy
        (
            pos0(kolmogorovLengthScale_ - fj.dSph())
           *BEddy_
           *muc
           *eddyStrainRate_
        );

        binaryBreakupRate +=
            pos(tauEddy - tauCrit)
           /fj.dSph()
           *sqrt(mag(tauEddy - tauCrit)/rhoc)/fj.x();
    }

    if (interfacialFriction_)
    {
        const volScalarField::Internal tauFric
        (
            BFric_*0.5*rhoc*sqr(uTerminal_[j])*Cd_[j]
        );

        binaryBreakupRate +=
            pos(tauFric - tauCrit)
           /fj.dSph()
           *sqrt(mag(tauFric - tauCrit)/rhoc)/fj.x();
    }
}


// ************************************************************************* //
