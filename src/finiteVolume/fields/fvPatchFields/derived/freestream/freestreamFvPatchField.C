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

\*---------------------------------------------------------------------------*/

#include "freestreamFvPatchField.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::freestreamFvPatchField<Type>::freestreamFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchField<Type>(p, iF)
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    freestreamValue() =
        Field<Type>("freestreamValue", iF.dimensions(), dict, p.size());

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", iF.dimensions(), dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(freestreamValue());
    }
}


template<class Type>
Foam::freestreamFvPatchField<Type>::freestreamFvPatchField
(
    const freestreamFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    inletOutletFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::freestreamFvPatchField<Type>::freestreamFvPatchField
(
    const freestreamFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    inletOutletFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::freestreamFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (this->phiName_ != "phi")
    {
        writeEntry(os, "phi", this->phiName_);
    }
    writeEntry(os, "freestreamValue", freestreamValue());
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
