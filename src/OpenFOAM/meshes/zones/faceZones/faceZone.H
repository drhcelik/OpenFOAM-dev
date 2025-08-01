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
    Foam::faceZone

Description
    Named list of face indices representing a sub-set of the mesh faces

    with optional flipMap to provide orientation relative to the orientation
    of the corresponding mesh faces

    Used by mesh-manipulation tools, field initialisation and for sources and
    functionObjects that apply to sub-sets of the faces.

See also
    faceZoneList
    Zone
    zoneGenerator

SourceFiles
    faceZone.C

\*---------------------------------------------------------------------------*/

#ifndef faceZone_H
#define faceZone_H

#include "Zone.H"
#include "boolList.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faceZoneList;

/*---------------------------------------------------------------------------*\
                          Class faceZone Declaration
\*---------------------------------------------------------------------------*/

class faceZone
:
    public Zone<faceZone, faceZoneList>
{
    // Private Data

        //- The name associated with the zone-labels dictionary entry
        static const word labelsName_;

        //- faceZone is oriented if the flipMap has been set
        bool oriented_;

        //- Flip map for all faces in the zone.  Set to true if the
        //  face needs to be flipped to achieve the correct orientation.
        boolList flipMap_;


        // Demand-driven private data

            //- Primitive patch made out of correctly flipped faces
            mutable primitiveFacePatch* patchPtr_;

            //- Master cell layer
            mutable labelList* masterCellsPtr_;

            //- Slave cell layer
            mutable labelList* slaveCellsPtr_;

            //- Global edge addressing
            mutable labelList* mePtr_;


    // Private Member Functions

        //- Build primitive patch
        void calcFaceZonePatch() const;

        //- Return map of local face indices
        const Map<label>& faceLookupMap() const;

        //- Calculate master and slave face layer
        void calcCellLayers() const;

        //- Check addressing
        void checkAddressing() const;

        //- Reset the indices and flipMap from the given map
        void reset(const Map<bool>& indices);

        //- Reset the indices from the given set
        void reset(const labelHashSet& indices);


public:

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        static const char* const labelsName;

        //- Runtime type information
        ClassName("faceZone");


    // Constructors

        //- Construct from components
        faceZone
        (
            const word& name,
            const labelUList& addr,
            const boolList& fm,
            const faceZoneList& mz,
            const bool moveUpdate = false,
            const bool topoUpdate = false
        );

        //- Construct from components without flipMap
        faceZone
        (
            const word& name,
            const labelUList& addr,
            const faceZoneList& mz,
            const bool moveUpdate = false,
            const bool topoUpdate = false
        );

        //- Construct from components, moving contents
        faceZone
        (
            const word& name,
            labelList&& addr,
            boolList&& fm,
            const faceZoneList&,
            const bool moveUpdate = false,
            const bool topoUpdate = false
        );

        //- Construct from components without flipMap, moving contents
        faceZone
        (
            const word& name,
            labelList&& addr,
            const faceZoneList&,
            const bool moveUpdate = false,
            const bool topoUpdate = false
        );

        //- Construct from dictionary
        faceZone
        (
            const word& name,
            const dictionary&,
            const faceZoneList&
        );

        //- Construct given the original zone and resetting the
        //  face list and mesh zones information
        faceZone
        (
            const faceZone&,
            const word& name,
            const labelUList& addr,
            const boolList& fm,
            const faceZoneList&
        );

        //- Construct given the original zone and resetting the
        //  face list and mesh zones information without flipMap
        faceZone
        (
            const faceZone&,
            const word& name,
            const labelUList& addr,
            const faceZoneList&
        );

        //- Construct given the original zone, resetting the
        //  face list and mesh zones information
        faceZone
        (
            const faceZone&,
            labelList&& addr,
            boolList&& fm,
            const faceZoneList&
        );

        //- Construct given the original zone, resetting the
        //  face list and mesh zones information without flipMap
        faceZone
        (
            const faceZone&,
            labelList&& addr,
            const faceZoneList&
        );

        //- Disallow default bitwise copy construction
        faceZone(const faceZone&) = delete;

        //- Construct and return a clone
        autoPtr<faceZone> clone() const;

        //- Construct and return a clone with a new name
        autoPtr<faceZone> clone(const word& name) const;

        //- Construct and return a clone, resetting the mesh zones
        autoPtr<faceZone> clone(const faceZoneList& mz) const;

        //- Construct and return a clone, resetting the face list
        //  and mesh zones
        autoPtr<faceZone> clone
        (
            const labelUList& addr,
            const boolList& fm,
            const faceZoneList& mz
        ) const;

        //- Construct and return a clone, resetting the face list
        //  and mesh zones
        autoPtr<faceZone> clone
        (
            const labelUList& addr,
            const faceZoneList& mz
        ) const;


    //- Destructor
    ~faceZone();


    // Member Functions

        //- Return true if the faceZone is oriented, i.e. the flipMap is set
        bool oriented() const
        {
            return oriented_;
        }

        //- Return face flip map
        const boolList& flipMap() const
        {
            if (!oriented_)
            {
                FatalErrorInFunction
                    << "flipMap requested for the unoriented faceZone "
                    << name()
                    << exit(FatalError);
            }

            return flipMap_;
        }

        //- Return reference to primitive patch
        const primitiveFacePatch& patch() const;


        // Addressing into mesh

            //- Return labels of master cells (cells next to the master face
            //  zone in the prescribed direction)
            const labelList& masterCells() const;

            //- Return labels of slave cells
            const labelList& slaveCells() const;

            //- Return global edge index for local edges
            const labelList& meshEdges() const;


        //- Clear addressing
        void clearAddressing();

        //- Reset addressing and flip map (clearing demand-driven data)
        void resetAddressing(const labelUList&, const boolList&);

        //- Reset addressing (clearing demand-driven data)
        void resetAddressing(const labelUList&);

        //- Check zone definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Check whether all procs have faces synchronised. Return
        //  true if in error.
        bool checkParallelSync(const bool report = false) const;

        //- Insert given indices and corresponding face flips into zone
        void insert(const Map<bool>& newIndices);

        //- Insert given indices into zone
        void insert(const labelHashSet& newIndices);

        //- Swap two faceZones
        void swap(faceZone&);

        //- Correct patch after moving points
        void movePoints(const pointField&);

        //- Update zone using the given map
        void topoChange(const polyTopoChangeMap&);

        //- Write dictionary
        void writeDict(Ostream&) const;


    // Member Operators

        //- Assignment to zone, clearing demand-driven data
        void operator=(const faceZone&);

        //- Move assignment to zone, clearing demand-driven data
        void operator=(faceZone&&);

        //- Assignment operator to indices
        void operator=(const labelUList&);

        //- Move assignment of indices
        void operator=(labelList&&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
