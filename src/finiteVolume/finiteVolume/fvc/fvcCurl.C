/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "fvcCurl.H"
#include "fvcGrad.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>>
curl
(
    const VolField<Type>& vf
)
{
    word nameCurlVf = "curl(" + vf.name() + ')';

    // Gausses theorem curl
    // tmp<VolField<Type>> tcurlVf =
    // fvc::surfaceIntegrateExtrapolate(vf.mesh().Sf() ^ fvc::interpolate(vf));

    // Calculate curl as the Hodge dual of the skew-symmetric part of grad
    tmp<VolField<Type>> tcurlVf =
        2.0*(*skew(fvc::grad(vf, nameCurlVf)));

    tcurlVf.ref().rename(nameCurlVf);

    return tcurlVf;
}


template<class Type>
tmp<VolField<Type>>
curl
(
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolField<Type>> Curl(fvc::curl(tvf()));
    tvf.clear();
    return Curl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
