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
    Foam::atmBoundaryLayerInletEpsilonFvPatchScalarField

Description
    This boundary condition specifies an inlet value for the turbulence
    dissipation, \f$\epsilon\f$, appropriate for atmospheric boundary layers.

    See Foam::atmBoundaryLayer for details.

    Example of the boundary condition specification:
    \verbatim
    ground
    {
        type            atmBoundaryLayerInletEpsilon;
        flowDir         (1 0 0);
        zDir            (0 0 1);
        Uref            10.0;
        Zref            20.0;
        z0              uniform 0.1;
        zGround         uniform 0.0;
    }
    \endverbatim

See also
    Foam::atmBoundaryLayer,
    Foam::atmBoundaryLayerInletVelocityFvPatchVectorField,
    Foam::atmBoundaryLayerInletKFvPatchScalarField

SourceFiles
    atmBoundaryLayerInletEpsilonFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmBoundaryLayerInletEpsilonFvPatchScalarField_H
#define atmBoundaryLayerInletEpsilonFvPatchScalarField_H

#include "fvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "atmBoundaryLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class atmBoundaryLayerInletEpsilonFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmBoundaryLayerInletEpsilonFvPatchScalarField
:
    public inletOutletFvPatchScalarField,
    public atmBoundaryLayer
{

public:

    //- Runtime type information
    TypeName("atmBoundaryLayerInletEpsilon");


    // Constructors

        //- Construct from patch, internal field and dictionary
        atmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  atmBoundaryLayerInletEpsilonFvPatchScalarField onto a new patch
        atmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const atmBoundaryLayerInletEpsilonFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        atmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const atmBoundaryLayerInletEpsilonFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new atmBoundaryLayerInletEpsilonFvPatchScalarField(*this, iF)
            );
        }


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
        using inletOutletFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
