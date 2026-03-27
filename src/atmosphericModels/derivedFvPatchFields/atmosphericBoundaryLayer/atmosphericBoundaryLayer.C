/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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

#include "atmosphericBoundaryLayer.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(atmosphericBoundaryLayer, 0);
}

const Foam::word Foam::atmosphericBoundaryLayer::dictName
(
    atmosphericBoundaryLayer::typeName + "Properties"
);

const Foam::scalar Foam::atmosphericBoundaryLayer::kappaDefault_ = 0.41;
const Foam::scalar Foam::atmosphericBoundaryLayer::CmuDefault_ = 0.09;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::atmosphericBoundaryLayer::z0
(
    const vectorField& C
) const
{
    return tmp<scalarField>(new scalarField(C.size(), z0_->value(0, 0)));
}


Foam::tmp<Foam::scalarField> Foam::atmosphericBoundaryLayer::Ustar
(
    const scalarField& z0
) const
{
    return kappa_*Uref_/(log((Zref_ + z0)/z0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmosphericBoundaryLayer::atmosphericBoundaryLayer
(
    const objectRegistry& db
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            db.time().constant(),
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    flowDir_(normalised(lookup<vector>("flowDir", dimless))),
    zDir_(normalised(lookup<vector>("zDir", dimless))),
    kappa_(lookupOrDefault<scalar>("kappa", dimless, kappaDefault_)),
    Cmu_(lookupOrDefault<scalar>("Cmu", dimless, CmuDefault_)),
    Uref_(lookup<scalar>("Uref", dimVelocity)),
    Zref_(lookup<scalar>("Zref", dimLength)),
    z0_
    (
        Function2<scalar>::New
        (
            "z0",
            dimLength,
            dimLength,
            dimLength,
            *this
        )
    ),
    zGround_
    (
        Function2<scalar>::New
        (
            "zGround",
            dimLength,
            dimLength,
            dimLength,
            *this
        )
    ),
    offset_(found("Ulower")),
    Ulower_(lookupOrDefault<scalar>("Ulower", dimVelocity, 0)),
    kLower_(lookupOrDefault<scalar>("kLower", dimEnergy/dimMass, 0)),
    epsilonLower_
    (
        lookupOrDefault<scalar>
        (
            "epsilonLower",
            dimEnergy/dimMass/dimTime,
            0
        )
    )
{
    if (mag(flowDir_) < small || mag(zDir_) < small)
    {
        FatalErrorInFunction
            << "magnitude of n or z must be greater than zero"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

const Foam::atmosphericBoundaryLayer& Foam::atmosphericBoundaryLayer::New
(
    const objectRegistry& db
)
{
    if (db.foundObject<atmosphericBoundaryLayer>(dictName))
    {
        return db.lookupObject<atmosphericBoundaryLayer>(dictName);
    }
    else
    {
        atmosphericBoundaryLayer* ptr = new atmosphericBoundaryLayer(db);

        ptr->store();

        return *ptr;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::atmosphericBoundaryLayer::U
(
    const vectorField& p
) const
{
    const scalarField z0(p.size(), z0_->value(0, 0));
    const scalarField zGround(p.size(), zGround_->value(0, 0));

    const scalarField Un
    (
        (Ustar(z0)/kappa_)
       *log(max((zDir_ & p) - zGround + z0, z0)/z0)
    );

    if (offset_)
    {
        return flowDir_*Un + flowDir_*Ulower_;
    }
    else
    {
        return flowDir_*Un;
    }
}


Foam::tmp<Foam::scalarField> Foam::atmosphericBoundaryLayer::k
(
    const vectorField& p
) const
{
    const scalarField z0(p.size(), z0_->value(0, 0));
    const scalarField zGround(p.size(), zGround_->value(0, 0));

    tmp<scalarField> tk
    (
        sqr(Ustar(z0))/sqrt(Cmu_)
    );

    if (offset_)
    {
        const scalarField z((zDir_ & p) - zGround);
        tk.ref() = pos0(z)*tk() + neg(z)*kLower_;
    }

    return tk;
}


Foam::tmp<Foam::scalarField> Foam::atmosphericBoundaryLayer::epsilon
(
    const vectorField& p
) const
{
    const scalarField z0(p.size(), z0_->value(0, 0));
    const scalarField zGround(p.size(), zGround_->value(0, 0));

    tmp<scalarField> tepsilon
    (
        pow3(Ustar(z0))/(kappa_*((zDir_ & p) - zGround + z0))
    );

    if (offset_)
    {
        const scalarField z((zDir_ & p) - zGround);
        tepsilon.ref() = pos0(z)*tepsilon() + neg(z)*epsilonLower_;
    }

    return tepsilon;
}


void Foam::atmosphericBoundaryLayer::write(Ostream& os) const
{
    writeEntry(os, "flowDir", flowDir_);
    writeEntry(os, "zDir", zDir_);
    writeEntry(os, "kappa", kappa_);
    writeEntry(os, "Cmu", Cmu_);
    writeEntry(os, "Uref", Uref_);
    writeEntry(os, "Zref", Zref_);

    if (offset_)
    {
        writeEntry(os, "Ulower", Ulower_);
        writeEntry(os, "kLower", kLower_);
        writeEntry(os, "epsilonLower", epsilonLower_);
    }

    writeEntry
    (
        os,
        dimLength,
        dimLength,
        dimLength,
        z0_()
    );

    writeEntry
    (
        os,
        dimLength,
        dimLength,
        dimLength,
        zGround_()
    );
}


// ************************************************************************* //
