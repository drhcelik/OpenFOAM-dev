/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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

#include "externalCoupledTemperatureMixedFvPatchScalarField.H"
#include "fluidThermophysicalTransportModel.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::writeHeader
(
    OFstream& os
) const
{
    os  << "# Values: magSf value qDot htc" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    externalCoupledMixedFvPatchScalarField(p, iF, dict)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    externalCoupledMixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ecmpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalCoupledMixedFvPatchScalarField(ecmpf, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::externalCoupledTemperatureMixedFvPatchScalarField::
~externalCoupledTemperatureMixedFvPatchScalarField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::transferData
(
    OFstream& os
) const
{
    if (log())
    {
        Info<< type() << ": " << this->patch().name()
            << ": writing data to " << os.name()
            << endl;
    }

    const label patchi = patch().index();

    // heat flux [W/m^2]
    scalarField qDot(this->patch().size(), 0.0);

    const fluidThermophysicalTransportModel& ttm =
        db().lookupType<fluidThermophysicalTransportModel>
        (
            internalField().group()
        );

    // patch temperature [K]
    const scalarField Tp(*this);

    qDot = ttm.kappaEff(patchi)*snGrad();

    // near wall cell temperature [K]
    const scalarField Tc(patchInternalField());

    // heat transfer coefficient [W/m^2/K]
    const scalarField htc(qDot/(Tp - Tc + rootVSmall));

    if (Pstream::parRun())
    {
        int tag = Pstream::msgType() + 1;

        List<Field<scalar>> magSfs(Pstream::nProcs());
        magSfs[Pstream::myProcNo()].setSize(this->patch().size());
        magSfs[Pstream::myProcNo()] = this->patch().magSf();
        Pstream::gatherList(magSfs, tag);

        List<Field<scalar>> values(Pstream::nProcs());
        values[Pstream::myProcNo()].setSize(this->patch().size());
        values[Pstream::myProcNo()] = Tp;
        Pstream::gatherList(values, tag);

        List<Field<scalar>> qDots(Pstream::nProcs());
        qDots[Pstream::myProcNo()].setSize(this->patch().size());
        qDots[Pstream::myProcNo()] = qDot;
        Pstream::gatherList(qDots, tag);

        List<Field<scalar>> htcs(Pstream::nProcs());
        htcs[Pstream::myProcNo()].setSize(this->patch().size());
        htcs[Pstream::myProcNo()] = htc;
        Pstream::gatherList(htcs, tag);

        if (Pstream::master())
        {
            forAll(values, proci)
            {
                const Field<scalar>& magSf = magSfs[proci];
                const Field<scalar>& value = values[proci];
                const Field<scalar>& qDot = qDots[proci];
                const Field<scalar>& htc = htcs[proci];

                forAll(magSf, facei)
                {
                    os  << magSf[facei] << token::SPACE
                        << value[facei] << token::SPACE
                        << qDot[facei] << token::SPACE
                        << htc[facei] << token::SPACE
                        << nl;
                }
            }

            os.flush();
        }
    }
    else
    {
        const Field<scalar>& magSf(this->patch().magSf());

        forAll(patch(), facei)
        {
            os  << magSf[facei] << token::SPACE
                << Tp[facei] << token::SPACE
                << qDot[facei] << token::SPACE
                << htc[facei] << token::SPACE
                << nl;
        }

        os.flush();
    }
}


void Foam::externalCoupledTemperatureMixedFvPatchScalarField::evaluate
(
    const Pstream::commsTypes comms
)
{
    externalCoupledMixedFvPatchScalarField::evaluate(comms);
}


void Foam::externalCoupledTemperatureMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    externalCoupledMixedFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalCoupledTemperatureMixedFvPatchScalarField
    );
}


// ************************************************************************* //
