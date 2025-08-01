/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2025 OpenFOAM Foundation
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

#include "limitTemperature.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitTemperature, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        limitTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::limitTemperature::readCoeffs(const dictionary& dict)
{
    Tmin_ = dict.lookup<scalar>("min");
    Tmax_ = dict.lookup<scalar>("max");
    fieldName_ = dict.lookupOrDefault<word>("field", word::null);
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitTemperature::limitTemperature
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict),
    zone_(mesh, coeffs(dict)),
    Tmin_(-vGreat),
    Tmax_(vGreat),
    fieldName_(word::null),
    phaseName_(word::null)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::limitTemperature::constrainedFields() const
{
    if (fieldName_ != word::null)
    {
        return wordList(1, IOobject::groupName(fieldName_, phaseName_));
    }
    else
    {
        const basicThermo& thermo =
            mesh().lookupObject<basicThermo>
            (
                IOobject::groupName(physicalProperties::typeName, phaseName_)
            );

        return wordList(1, thermo.he().name());
    }
}


bool Foam::fv::limitTemperature::constrain(volScalarField& he) const
{
    const labelList& cells = zone_.zone();

    if (he.dimensions() == dimTemperature)
    {
        scalarField& Tc = he.primitiveFieldRef();

        forAll(cells, i)
        {
            const label celli = cells[i];
            Tc[celli] = max(min(Tc[celli], Tmax_), Tmin_);
        }

        // Handle boundaries in the case of 'all'
        if (zone_.all())
        {
            volScalarField::Boundary& Tbf = he.boundaryFieldRef();

            forAll(Tbf, patchi)
            {
                fvPatchScalarField& Tp = Tbf[patchi];

                if (!Tp.fixesValue())
                {
                    forAll(Tp, facei)
                    {
                        Tp[facei] = max(min(Tp[facei], Tmax_), Tmin_);
                    }
                }
            }
        }
    }
    else
    {
        const basicThermo& thermo =
            mesh().lookupObject<basicThermo>
            (
                IOobject::groupName(physicalProperties::typeName, phaseName_)
            );

        const scalarField Tmin(cells.size(), Tmin_);
        const scalarField Tmax(cells.size(), Tmax_);

        const scalarField heMin(thermo.he(Tmin, cells));
        const scalarField heMax(thermo.he(Tmax, cells));

        scalarField& hec = he.primitiveFieldRef();

        forAll(cells, i)
        {
            const label celli = cells[i];
            hec[celli] = max(min(hec[celli], heMax[i]), heMin[i]);
        }

        // Handle boundaries in the case of 'all'
        if (zone_.all())
        {
            volScalarField::Boundary& bf = he.boundaryFieldRef();

            forAll(bf, patchi)
            {
                fvPatchScalarField& hep = bf[patchi];

                if (!hep.fixesValue())
                {
                    const scalarField Tminp(hep.size(), Tmin_);
                    const scalarField Tmaxp(hep.size(), Tmax_);

                    const scalarField heMinp(thermo.he(Tminp, patchi));
                    const scalarField heMaxp(thermo.he(Tmaxp, patchi));

                    forAll(hep, facei)
                    {
                        hep[facei] =
                            max(min(hep[facei], heMaxp[facei]), heMinp[facei]);
                    }
                }
            }
        }
    }

    return cells.size();
}


bool Foam::fv::limitTemperature::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::limitTemperature::topoChange(const polyTopoChangeMap& map)
{
    zone_.topoChange(map);
}


void Foam::fv::limitTemperature::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
}


void Foam::fv::limitTemperature::distribute(const polyDistributionMap& map)
{
    zone_.distribute(map);
}


bool Foam::fv::limitTemperature::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        zone_.read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
