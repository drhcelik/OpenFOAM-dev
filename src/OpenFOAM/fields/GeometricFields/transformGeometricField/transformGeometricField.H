/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

InClass
    Foam::transformGeometricField

Description
    Spatial transformation functions for FieldFields.

SourceFiles
    transformGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef transformGeometricField_H
#define transformGeometricField_H

#include "transform.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
void transform
(
    GeometricField<Type, GeoMesh>&,
    const GeometricField<tensor, GeoMesh>&,
    const GeometricField<Type, GeoMesh>&
);

template<class Type, class GeoMesh>
tmp<GeometricField<Type, GeoMesh>> transform
(
    const GeometricField<tensor, GeoMesh>&,
    const GeometricField<Type, GeoMesh>&
);

template<class Type, class GeoMesh>
tmp<GeometricField<Type, GeoMesh>> transform
(
    const GeometricField<tensor, GeoMesh>&,
    const tmp<GeometricField<Type, GeoMesh>>&
);

template<class Type, class GeoMesh>
tmp<GeometricField<Type, GeoMesh>> transform
(
    const tmp<GeometricField<tensor, GeoMesh>>&,
    const GeometricField<Type, GeoMesh>&
);

template<class Type, class GeoMesh>
tmp<GeometricField<Type, GeoMesh>> transform
(
    const tmp<GeometricField<tensor, GeoMesh>>&,
    const tmp<GeometricField<Type, GeoMesh>>&
);


template<class Type, class GeoMesh>
void transform
(
    GeometricField<Type, GeoMesh>&,
    const dimensionedTensor&,
    const GeometricField<Type, GeoMesh>&
);

template<class Type, class GeoMesh>
tmp<GeometricField<Type, GeoMesh>> transform
(
    const dimensionedTensor&,
    const GeometricField<Type, GeoMesh>&
);

template<class Type, class GeoMesh>
tmp<GeometricField<Type, GeoMesh>> transform
(
    const dimensionedTensor&,
    const tmp<GeometricField<Type, GeoMesh>>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "transformGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
