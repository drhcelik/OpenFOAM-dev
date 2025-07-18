/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2025 OpenFOAM Foundation
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

#include "rotorDisk.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalar Foam::fv::rotorDisk::rhoRef() const
{
    return rhoRef_;
}


inline Foam::scalar Foam::fv::rotorDisk::omega() const
{
    return omega_;
}


inline const Foam::fvCellZone& Foam::fv::rotorDisk::zone() const
{
    return zone_;
}


inline const Foam::List<Foam::point>& Foam::fv::rotorDisk::x() const
{
    return x_;
}


inline const Foam::coordinateSystems::cylindrical&
Foam::fv::rotorDisk::coordSys() const
{
    return coordSys_;
}


// ************************************************************************* //
