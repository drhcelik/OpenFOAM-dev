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

Class
    Foam::nonConformalErrorPointPatchField

Description
    Constraint pointPatchField for nonConformalError patches.

SourceFiles
    nonConformalErrorPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalErrorPointPatchField_H
#define nonConformalErrorPointPatchField_H

#include "pointPatchField.H"
#include "nonConformalErrorPointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonConformalErrorPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonConformalErrorPointPatchField
:
    public pointPatchField<Type>
{

public:

    //- Runtime type information
    TypeName(nonConformalErrorPointPatch::typeName_());


    // Constructors

        //- Construct from patch and nonConformalError field
        nonConformalErrorPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, nonConformalError field and dictionary
        nonConformalErrorPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        nonConformalErrorPointPatchField
        (
            const nonConformalErrorPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting nonConformalError field reference
        nonConformalErrorPointPatchField
        (
            const nonConformalErrorPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting nonConformalError field
        //  reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new nonConformalErrorPointPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member Operators

        //- Inherit assignment
        using pointPatchField<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonConformalErrorPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
