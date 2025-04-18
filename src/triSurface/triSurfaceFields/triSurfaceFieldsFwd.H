/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#ifndef triSurfaceFieldsFwd_H
#define triSurfaceFieldsFwd_H

#include "DimensionedFieldFwd.H"
#include "fieldTypes.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfaceGeoMesh;

typedef DimensionedField<label, triSurfaceGeoMesh, Field> triSurfaceLabelField;
typedef DimensionedField<scalar, triSurfaceGeoMesh, Field>
    triSurfaceScalarField;
typedef DimensionedField<vector, triSurfaceGeoMesh, Field>
    triSurfaceVectorField;
typedef DimensionedField<sphericalTensor, triSurfaceGeoMesh, Field>
    triSurfaceSphericalTensorField;
typedef DimensionedField<symmTensor, triSurfaceGeoMesh, Field>
    triSurfaceSymmTensorField;
typedef DimensionedField<tensor, triSurfaceGeoMesh, Field>
    triSurfaceTensorField;

class triSurfacePointGeoMesh;

typedef DimensionedField<label, triSurfacePointGeoMesh, Field>
    triSurfacePointLabelField;
typedef DimensionedField<scalar, triSurfacePointGeoMesh, Field>
    triSurfacePointScalarField;
typedef DimensionedField<vector, triSurfacePointGeoMesh, Field>
    triSurfacePointVectorField;
typedef DimensionedField<sphericalTensor, triSurfacePointGeoMesh, Field>
    triSurfacePointSphericalTensorField;
typedef DimensionedField<symmTensor, triSurfacePointGeoMesh, Field>
    triSurfacePointSymmTensorField;
typedef DimensionedField<tensor, triSurfacePointGeoMesh, Field>
    triSurfacePointTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
