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

#include "volFieldValue.H"
#include "fvMesh.H"
#include "volFields.H"
#include "writeFile.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(volFieldValue, 0);
    addToRunTimeSelectionTable(fieldValue, volFieldValue, dictionary);
    addToRunTimeSelectionTable(functionObject, volFieldValue, dictionary);
}
}
}

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volFieldValue::operationType,
    12
> Foam::functionObjects::fieldValues::volFieldValue::operationTypeNames_
{
    "none",
    "sum",
    "sumMag",
    "average",
    "volAverage",
    "volIntegrate",
    "min",
    "max",
    "minMag",
    "maxMag",
    "CoV",
    "UI"
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldValues::volFieldValue::
writeFileHeaderLocation()
{
    switch (operation_)
    {
        case operationType::minMag:
        case operationType::maxMag:
            file() << tab << "location" << tab << "cell";
            if (Pstream::parRun())
            {
                file() << tab << "processor";
            }
            break;
        default:
            break;
    }
}


template<>
void Foam::functionObjects::fieldValues::volFieldValue::
writeFileHeaderLocation<Foam::scalar>()
{
    switch (operation_)
    {
        case operationType::min:
        case operationType::max:
        case operationType::minMag:
        case operationType::maxMag:
            file() << tab << "location" << tab << "cell";
            if (Pstream::parRun())
            {
                file() << tab << "processor";
            }
            break;
        default:
            break;
    }
}


void Foam::functionObjects::fieldValues::volFieldValue::writeFileHeader
(
    const label i
)
{
    if (operation_ != operationType::none)
    {
        zone_.writeFileHeader(*this, file());

        writeCommented(file(), "Time");

        if (writeNCells_) file() << tab << "Cells";
        if (writeVolume_) file() << tab << "Volume";

        forAll(fields_, fieldi)
        {
            file() << tab << operationTypeNames_[operation_] << "(";

            forAll(weightFieldNames_, i)
            {
                file() << weightFieldNames_[i] << ',';
            }

            file() << fields_[fieldi] << ")";

            if (writeLocation_)
            {
                #define writeFileHeaderLocationFieldType(fieldType, none)      \
                    if (validField<fieldType>(fields_[fieldi]))                \
                    {                                                          \
                        writeFileHeaderLocation<fieldType>();                  \
                    }
                FOR_ALL_FIELD_TYPES(writeFileHeaderLocationFieldType)
                #undef writeHeaderLocationFieldType
            }
        }

        file() << endl;
    }
}


bool Foam::functionObjects::fieldValues::volFieldValue::processValues
(
    const Field<scalar>& values,
    const scalarField& weights,
    const scalarField& V,
    Result<scalar>& result
) const
{
    switch (operation_)
    {
        case operationType::min:
        {
            compareScalars(values, vGreat, result, lessOp<scalar>());
            return true;
        }
        case operationType::minMag:
        {
            compareScalars(mag(values), vGreat, result, lessOp<scalar>());
            return true;
        }
        case operationType::max:
        {
            compareScalars(values, -vGreat, result, greaterOp<scalar>());
            return true;
        }
        case operationType::maxMag:
        {
            compareScalars(mag(values), -vGreat, result, greaterOp<scalar>());
            return true;
        }
        default:
        {
            // Fall through to same-type operations
            return processValuesTypeType(values, weights, V, result);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    zone_(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    scaleFactor_(1),
    writeNCells_(dict.lookupOrDefault("writeNumberOfCells", false)),
    writeVolume_(dict.lookupOrDefault("writeVolume", false)),
    writeLocation_(dict.lookupOrDefault("writeLocation", false))
{
    read(dict);
}


Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    zone_(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    scaleFactor_(1),
    writeNCells_(dict.lookupOrDefault("writeNumberOfCells", false)),
    writeVolume_(dict.lookupOrDefault("writeVolume", false)),
    writeLocation_(dict.lookupOrDefault("writeLocation", false))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volFieldValue::~volFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::volFieldValue::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);

    if (dict.found("weightFields"))
    {
        dict.lookup("weightFields") >> weightFieldNames_;
    }
    else if (dict.found("weightField"))
    {
        weightFieldNames_.setSize(1);

        dict.lookup("weightField") >> weightFieldNames_[0];
    }

    dict.readIfPresent("scaleFactor", scaleFactor_);

    // Report configuration
    Info<< type() << ' ' << name() << " read:" << nl;
    Info<< "    number of cells = " << zone_.nGlobalCells() << nl;
    if (zone_.nGlobalCells())
    {
        Info<< "    volume = " << zone_.V() << nl;
    }
    Info<< "    operation = " << operationTypeNames_[operation_] << nl;
    if (weightFieldNames_.size() == 1)
    {
        Info<< "    weight field = " << weightFieldNames_[0] << nl;
    }
    if (weightFieldNames_.size() > 1)
    {
        Info<< "    weight fields =";
        forAll(weightFieldNames_, i) Info<< ' ' << weightFieldNames_[i];
        Info<< nl;
    }
    if (scaleFactor_ != scalar(1))
    {
        Info<< "    scale factor = " << scaleFactor_;
    }
    Info<< endl;

    return true;
}


bool Foam::functionObjects::fieldValues::volFieldValue::write()
{
    // Look to see if any fields exist. Use the flag to suppress output later.
    bool anyFields = false;
    forAll(fields_, i)
    {
        #define validFieldType(fieldType, none)                          \
            anyFields = anyFields || validField<fieldType>(fields_[i]);
        FOR_ALL_FIELD_TYPES(validFieldType);
        #undef validFieldType
    }
    if (!anyFields && fields_.size() > 1) // (error for 1 will happen below)
    {
        cannotFindObjects(fields_);
    }

    // Initialise the file, write the header, etc...
    if (anyFields && operation_ != operationType::none)
    {
        fieldValue::write();
    }

    // Write the time
    if (anyFields && operation_ != operationType::none && Pstream::master())
    {
        writeTime(file());
    }

    // Write the number of faces and/or the area if necessary
    if (anyFields && operation_ != operationType::none && Pstream::master())
    {
        if (writeNCells_)
        {
            file() << tab << zone_.nGlobalCells();
        }
        if (writeVolume_)
        {
            file() << tab << zone_.V();
        }
    }
    if (writeNCells_)
    {
        Log << "    number of cells = " << zone_.nGlobalCells() << endl;
    }
    if (writeVolume_)
    {
        Log << "    volume = " << zone_.V() << endl;
    }

    // Construct the weight field and the volumes
    scalarField weights(zone_.nCells(), 1);
    forAll(weightFieldNames_, i)
    {
        weights *= getFieldValues<scalar>(weightFieldNames_[i]);
    }
    const scalarField V(filterField(fieldValue::mesh_.V()));

    // Process the fields
    forAll(fields_, i)
    {
        const word& fieldName = fields_[i];
        bool ok = false;

        #define writeValuesFieldType(fieldType, none)                          \
            ok = ok || writeValues<fieldType>(fieldName, weights, V);
        FOR_ALL_FIELD_TYPES(writeValuesFieldType)
        #undef writeValuesFieldType

        if (!ok)
        {
            cannotFindObject(fieldName);
        }
    }

    // Finalise
    if (anyFields && operation_ != operationType::none && Pstream::master())
    {
        file() << endl;
    }
    Log << endl;

    return true;
}


void Foam::functionObjects::fieldValues::volFieldValue::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &this->mesh())
    {
        fieldValue::movePoints(mesh);
        zone_.movePoints();
    }
}


void Foam::functionObjects::fieldValues::volFieldValue::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        fieldValue::topoChange(map);
        zone_.topoChange(map);
    }
}


void Foam::functionObjects::fieldValues::volFieldValue::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        fieldValue::mapMesh(map);
        zone_.mapMesh(map);
    }
}


void Foam::functionObjects::fieldValues::volFieldValue::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        fieldValue::distribute(map);
        zone_.distribute(map);
    }
}


// ************************************************************************* //
