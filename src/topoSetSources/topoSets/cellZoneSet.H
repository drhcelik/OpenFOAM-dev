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
    Foam::cellZoneSet

Description
    Like cellSet but -reads data from cellZone -updates cellZone when writing.

SourceFiles
    cellZone.C

\*---------------------------------------------------------------------------*/

#ifndef cellZoneSet_H
#define cellZoneSet_H

#include "cellSet.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cellZoneSet Declaration
\*---------------------------------------------------------------------------*/

class cellZoneSet
:
    public cellSet
{
    // Private Data

        const polyMesh& mesh_;

        labelList addressing_;


public:

    //- Runtime type information
    TypeName("cellZoneSet");


    // Constructors

        //- Construct from objectRegistry and name
        cellZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            readOption r=MUST_READ,
            writeOption w=NO_WRITE
        );

        //- Construct from additional size of labelHashSet
        cellZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            const label,
            writeOption w=NO_WRITE
        );

        //- Construct from existing set
        cellZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet&,
            writeOption w=NO_WRITE
        );



    //- Destructor
    virtual ~cellZoneSet();


    // Member Functions

        const labelList& addressing() const
        {
            return addressing_;
        }

        labelList& addressing()
        {
            return addressing_;
        }

        //- Sort addressing and make cellSet part consistent with addressing
        void updateSet();

        //- Invert contents. (insert all members 0..maxLen-1 which were not in
        //  set)
        virtual void invert(const label maxLen);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const topoSet& set);

        //- Add elements present in set.
        virtual void addSet(const topoSet& set);

        //- Delete elements present in set.
        virtual void deleteSet(const topoSet& set);

        //- Sync cellSet across coupled patches; update cellZone from cellSet
        virtual void sync(const polyMesh& mesh);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const;

        //- Write cellZone
        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool write = true
        ) const;

        //- Update any stored data for new labels
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
