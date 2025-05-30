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

\*---------------------------------------------------------------------------*/

#include "InertPhaseModel.H"
#include "phaseSystem.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::InertPhaseModel<BasePhaseModel>::InertPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::InertPhaseModel<BasePhaseModel>::~InertPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField::Internal>
Foam::InertPhaseModel<BasePhaseModel>::R(const label speciei) const
{
    return
        volScalarField::Internal::New
        (
            IOobject::groupName("R_" + this->Y()[speciei].name(), this->name()),
            this->mesh(),
            dimensionedScalar(dimDensity/dimTime, 0)
        );
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::InertPhaseModel<BasePhaseModel>::R(volScalarField& Yi) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix(Yi, dimMass/dimTime)
    );
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::InertPhaseModel<BasePhaseModel>::Qdot() const
{
    return volScalarField::New
    (
        IOobject::groupName("Qdot", this->name()),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimTime/dimVolume, 0)
    );
}


// ************************************************************************* //
