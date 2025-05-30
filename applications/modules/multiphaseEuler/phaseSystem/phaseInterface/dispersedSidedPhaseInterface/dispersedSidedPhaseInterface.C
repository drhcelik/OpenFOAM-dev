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

#include "dispersedSidedPhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
        dispersedSidedPhaseInterface,
        separatorsToTypeName
        ({
            dispersedPhaseInterface::separator(),
            sidedPhaseInterface::separator()
        }).c_str(),
        0
    );
    addToRunTimeSelectionTable
    (
        phaseInterface,
        dispersedSidedPhaseInterface,
        word
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::dispersedSidedPhaseInterface::same
(
    const phaseInterface& interface,
    bool strict
) const
{
    return
        (!strict || isType<dispersedSidedPhaseInterface>(interface))
     && dispersedPhaseInterface::same(interface, false)
     && sidedPhaseInterface::same(interface, false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersedSidedPhaseInterface::dispersedSidedPhaseInterface
(
    const phaseModel& dispersed,
    const phaseModel& continuous,
    const phaseModel& phase
)
:
    phaseInterface(dispersed, continuous),
    dispersedPhaseInterface(dispersed, continuous),
    sidedPhaseInterface(phase, phaseInterface::otherPhase(phase))
{}


Foam::dispersedSidedPhaseInterface::dispersedSidedPhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    dispersedPhaseInterface(fluid, name),
    sidedPhaseInterface(fluid, name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dispersedSidedPhaseInterface::~dispersedSidedPhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::dispersedSidedPhaseInterface::name() const
{
    return
        dispersedPhaseInterface::name()
      + '_'
      + sidedPhaseInterface::separator()
      + '_'
      + phase().name();
}


// ************************************************************************* //
