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

Class
    Foam::dynamicPressureFvPatchScalarField

Description
    This boundary condition provides a dynamic pressure condition. It subtracts
    a kinetic energy term from a reference pressure to get a value which is
    fixed on the boundary. It forms the base class for the total and
    entrainment pressure conditions.

See also
    Foam::totalPressureFvPatchField

SourceFiles
    dynamicPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicPressureFvPatchScalarField_H
#define dynamicPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dynamicPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class dynamicPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected Data

        //- Name of the density field used to normalise the mass flux
        const word rhoName_;

        //- Name of the compressibility field used to calculate the wave speed
        const word psiName_;

        //- Heat capacity ratio
        const scalar gamma_;

        //- Reference pressure
        scalarField p0_;


    // Protected Member Functions

        using fixedValueFvPatchScalarField::updateCoeffs;

        //- Update the coefficients associated with the patch field
        //  using the given external static pressure, and difference between
        //  the external kinetic energy and patch kinetic energy fields
        void updateCoeffs
        (
            const scalarField& p0p,
            const scalarField& K0mKp
        );


public:

    //- Runtime type information
    TypeName("dynamicPressure");


    // Constructors

        //- Construct from patch, internal field and dictionary
        dynamicPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given dynamicPressureFvPatchScalarField
        //  onto a new patch
        dynamicPressureFvPatchScalarField
        (
            const dynamicPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        dynamicPressureFvPatchScalarField
        (
            const dynamicPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        dynamicPressureFvPatchScalarField
        (
            const dynamicPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
