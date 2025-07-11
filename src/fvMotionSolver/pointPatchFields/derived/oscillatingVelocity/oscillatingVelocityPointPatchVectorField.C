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

#include "oscillatingVelocityPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oscillatingVelocityPointPatchVectorField::
oscillatingVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict),
    amplitude_(dict.lookup<vector>("amplitude", dimLength)),
    omega_(dict.lookup<scalar>("omega", unitRadians/dimTime))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dimLength, dict, p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


oscillatingVelocityPointPatchVectorField::
oscillatingVelocityPointPatchVectorField
(
    const oscillatingVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    p0_(mapper(ptf.p0_))
{}


oscillatingVelocityPointPatchVectorField::
oscillatingVelocityPointPatchVectorField
(
    const oscillatingVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oscillatingVelocityPointPatchVectorField::map
(
    const pointPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    const oscillatingVelocityPointPatchVectorField& oVptf =
        refCast<const oscillatingVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchVectorField::map(oVptf, mapper);

    mapper(p0_, oVptf.p0_);
}


void oscillatingVelocityPointPatchVectorField::reset
(
    const pointPatchVectorField& ptf
)
{
    const oscillatingVelocityPointPatchVectorField& oVptf =
        refCast<const oscillatingVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchVectorField::reset(oVptf);

    p0_.reset(oVptf.p0_);
}


void oscillatingVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();
    const pointPatch& p = this->patch();

    Field<vector>::operator=
    (
        (p0_ + amplitude_*sin(omega_*t.value()) - p.localPoints())
       /t.deltaTValue()
    );

    fixedValuePointPatchVectorField::updateCoeffs();
}


void oscillatingVelocityPointPatchVectorField::write(Ostream& os) const
{
    pointPatchVectorField::write(os);
    writeEntry(os, "amplitude", amplitude_);
    writeEntry(os, "omega", omega_);
    writeEntry(os, "p0", p0_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    oscillatingVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
