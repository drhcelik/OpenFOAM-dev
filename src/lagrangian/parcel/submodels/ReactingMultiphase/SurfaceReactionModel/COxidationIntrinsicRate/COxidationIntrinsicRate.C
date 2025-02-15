/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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

#include "COxidationIntrinsicRate.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationIntrinsicRate<CloudType>::COxidationIntrinsicRate
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceReactionModel<CloudType>(dict, owner, typeName),
    Sb_(this->coeffDict().template lookup<scalar>("Sb")),
    C1_(this->coeffDict().template lookup<scalar>("C1")),
    rMean_(this->coeffDict().template lookup<scalar>("rMean")),
    theta_(this->coeffDict().template lookup<scalar>("theta")),
    Ai_(this->coeffDict().template lookup<scalar>("Ai")),
    Ei_(this->coeffDict().template lookup<scalar>("Ei")),
    Ag_(this->coeffDict().template lookup<scalar>("Ag")),
    tau_(this->coeffDict().lookupOrDefault("tau", sqrt(2.0))),
    CsLocalId_(-1),
    O2GlobalId_(owner.composition().carrierId("O2")),
    CO2GlobalId_(owner.composition().carrierId("CO2")),
    WC_(0.0),
    WO2_(0.0),
    HcCO2_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Set local copies of thermo properties
    WO2_ = owner.composition().carrier().WiValue(O2GlobalId_);
    const scalar WCO2 = owner.composition().carrier().WiValue(CO2GlobalId_);
    WC_ = WCO2 - WO2_;

    HcCO2_ = owner.composition().carrier().hfiValue(CO2GlobalId_);

    if (Sb_ < 0)
    {
        FatalErrorInFunction
            << "Stoichiometry of reaction, Sb, must be greater than zero" << nl
            << exit(FatalError);
    }

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;
}


template<class CloudType>
Foam::COxidationIntrinsicRate<CloudType>::COxidationIntrinsicRate
(
    const COxidationIntrinsicRate<CloudType>& srm
)
:
    SurfaceReactionModel<CloudType>(srm),
    Sb_(srm.Sb_),
    C1_(srm.C1_),
    rMean_(srm.rMean_),
    theta_(srm.theta_),
    Ai_(srm.Ai_),
    Ei_(srm.Ei_),
    Ag_(srm.Ag_),
    tau_(srm.tau_),
    CsLocalId_(srm.CsLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    WC_(srm.WC_),
    WO2_(srm.WO2_),
    HcCO2_(srm.HcCO2_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationIntrinsicRate<CloudType>::
~COxidationIntrinsicRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::COxidationIntrinsicRate<CloudType>::calculate
(
    const scalar dt,
    const label celli,
    const scalar d,
    const scalar T,
    const scalar Tc,
    const scalar pc,
    const scalar rhoc,
    const scalar mass,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    const scalarField& YMixture,
    const scalar N,
    scalarField& dMassGas,
    scalarField& dMassLiquid,
    scalarField& dMassSolid,
    scalarField& dMassSRCarrier
) const
{
    // Fraction of remaining combustible material
    const label idSolid = this->owner().composition().idSolid();
    const scalar Ychar = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion until combustible fraction is consumed
    if (Ychar < small)
    {
        return 0.0;
    }

    const parcelThermo& thermo = this->owner().thermo();
    const fluidMulticomponentThermo& carrierThermo =
        this->owner().composition().carrier();

    // Local mass fraction of O2 in the carrier phase []
    const scalar YO2 = carrierThermo.Y(O2GlobalId_)[celli];

    // Quick exit if oxidant not present
    if (YO2 < rootVSmall)
    {
        return 0.0;
    }

    // Diffusion rate coefficient [m^2/s]
    const scalar D0 = C1_/d*pow(0.5*(T + Tc), 0.75);

    // Apparent density of pyrolysis char [kg/m^3]
    const scalar rhop = 6.0*mass/(constant::mathematical::pi*pow3(d));

    // Knusden diffusion coefficient [m^2/s]
    const scalar Dkn = 97.0*rMean_*sqrt(T/WO2_);

    // Effective diffusion [m^2/s]
    const scalar De = theta_/sqr(tau_)/(1.0/Dkn + 1/D0);

    // Cell carrier phase O2 species density [kg/m^3]
    const scalar rhoO2 = rhoc*YO2;

    // Partial pressure O2 [Pa]
    const scalar ppO2 = rhoO2/WO2_*RR*Tc;

    // Intrinsic reactivity [1/s]
    const scalar ki = Ai_*exp(-Ei_/RR/T);

    // Thiele modulus []
    const scalar phi =
        max(0.5*d*sqrt(Sb_*rhop*Ag_*ki*ppO2/(De*rhoO2)), rootVSmall);

    // Effectiveness factor []
    const scalar eta = max(3.0/sqr(phi)*(phi/tanh(phi) - 1.0), 0.0);

    // Chemical rate [kmol/m2/s]
    const scalar R = eta*d/6.0*rhop*Ag_*ki;

    // Particle surface area [m^2]
    const scalar Ap = constant::mathematical::pi*sqr(d);

    // Change in C mass [kg]
    scalar dmC = Ap*rhoc*RR*Tc*YO2/WO2_*D0*R/(D0 + R)*dt;

    // Limit mass transfer by availability of C
    dmC = min(mass*Ychar, dmC);

    // Molar consumption [kmol]
    const scalar dOmega = dmC/WC_;

    // Change in O2 mass [kg]
    const scalar dmO2 = dOmega*Sb_*WO2_;

    // Mass of newly created CO2 [kg]
    const scalar dmCO2 = dOmega*(WC_ + Sb_*WO2_);

    // Update local particle C mass
    dMassSolid[CsLocalId_] += dOmega*WC_;

    // Update carrier O2 and CO2 mass
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[CO2GlobalId_] += dmCO2;

    const scalar hsC = thermo.solids().properties()[CsLocalId_].hs(T);

    // carrier sensible enthalpy exchange handled via change in mass

    // Heat of reaction [J]
    return dmC*hsC - dmCO2*HcCO2_;
}


// ************************************************************************* //
