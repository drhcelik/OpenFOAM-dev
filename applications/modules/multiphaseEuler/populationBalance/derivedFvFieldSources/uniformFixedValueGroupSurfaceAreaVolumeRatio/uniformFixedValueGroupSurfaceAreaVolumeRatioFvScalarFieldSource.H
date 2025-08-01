/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource

Description
    This source condition creates a uniform fixed value source term in the
    group surface-area-volume-ratio (kappa) equations of the fractal shape
    model of a population balance model.

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type            uniformFixedValueGroupSurfaceAreaVolumeRatio;
        libs            ("libmultiphaseEulerFvModels.so");
        uniformValue    1.5e10;
    }
    \endverbatim

SourceFiles
    uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource_H
#define uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource_H

#include "fvFieldSources.H"
#include "groupPropertyFvScalarField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
Class uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldS... Declaration
\*---------------------------------------------------------------------------*/

class uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource
:
    public fvScalarFieldSource,
    public groupPropertyFvScalarField
{
private:

    // Private Data

        //- Uniform value
        autoPtr<Function1<scalar>> uniformValue_;


public:

    //- Runtime type information
    TypeName("uniformFixedValueGroupSurfaceAreaVolumeRatio");


    // Constructors

        //- Construct from internal field and dictionary
        uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource
        (
            const
            uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvScalarFieldSource>
            (
                new
                uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~uniformFixedValueGroupSurfaceAreaVolumeRatioFvScalarFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<DimensionedField<scalar, volMesh>> sourceValue
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the source value
        virtual tmp<scalarField> sourceValue
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;

        //- Return the internal coefficient
        virtual tmp<DimensionedField<scalar, volMesh>> internalCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
