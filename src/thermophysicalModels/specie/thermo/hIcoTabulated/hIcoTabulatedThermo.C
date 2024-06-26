/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2024 OpenFOAM Foundation
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

#include "hIcoTabulatedThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hIcoTabulatedThermo<EquationOfState>::hIcoTabulatedThermo
(
    const word& name,
    const dictionary& dict
)
:
    EquationOfState(name, dict),
    hf_
    (
        dict
       .subDict("thermodynamics")
       .lookupBackwardsCompatible<scalar>({"hf", "Hf"})
    ),
    sf_
    (
        dict
       .subDict("thermodynamics")
       .lookupBackwardsCompatible<scalar>({"sf", "Sf"})
    ),
    Cp_
    (
        "Cp",
        {dimTemperature, dimSpecificHeatCapacity},
        dict.subDict("thermodynamics").subDict("Cp")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hIcoTabulatedThermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("hf", hf_);
    dict.add("sf", sf_);
    dict.add("Cp", Cp_.values());
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hIcoTabulatedThermo<EquationOfState>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
