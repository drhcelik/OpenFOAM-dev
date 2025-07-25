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

#include "totalFlowRateAdvectiveDiffusiveFvPatchScalarField.H"
#include "surfaceFields.H"
#include "fluidThermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    massFluxFraction_
    (
        dict.lookupOrDefault<scalar>("massFluxFraction", dimless, 1.0)
    )
{
    refValue() = 1.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            Field<scalar>("value", iF.dimensions(), dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(refValue());
    }
}

Foam::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const totalFlowRateAdvectiveDiffusiveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    massFluxFraction_(ptf.massFluxFraction_)
{}


Foam::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const totalFlowRateAdvectiveDiffusiveFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    massFluxFraction_(tppsf.massFluxFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);
}


void Foam::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);
}


void Foam::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const label patchi = patch().index();

    const fluidThermophysicalTransportModel& ttm =
        db().lookupType<fluidThermophysicalTransportModel>
        (
            internalField().group()
        );

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const scalarField alphap
    (
        ttm.kappaEff(patchi)/ttm.thermo().Cpv().boundaryField()[patchi]
    );

    refValue() = massFluxFraction_;
    refGrad() = 0.0;

    valueFraction() =
        1.0
        /
        (
            1.0 +
            alphap*patch().deltaCoeffs()*patch().magSf()/max(mag(phip), small)
        );

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar phi = gSum(-phip*(*this));

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " mass flux[Kg/s]:" << phi
            << endl;
    }
}


void Foam::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "phi", phiName_);
    writeEntry(os, "rho", rhoName_);
    writeEntry(os, "massFluxFraction", massFluxFraction_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
    );

}

// ************************************************************************* //
