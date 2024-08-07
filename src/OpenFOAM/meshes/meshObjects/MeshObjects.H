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

Class
    Foam::DeletableMeshObject

Description
    MeshObject types:

    - DeletableMeshObject:
        mesh object to be deleted after any mesh change
    - MoveableMeshObject:
        mesh object to be updated after mesh motion otherwise deleted
    - DistributeableMeshObject
        mesh object to be updated after mesh redistribution or motion
        otherwise deleted
    - TopoChangeableMeshObject:
        mesh object to be updated after mesh topology change,
        mesh-to-mesh mapping, redistribution or motion otherwise deleted
    - RepatchableMeshObject:
        mesh object to be updated on patch or topology change,
        mesh-to-mesh mapping, redistribution or motion otherwise deleted

\*---------------------------------------------------------------------------*/

#ifndef MeshObjects_H
#define MeshObjects_H

#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;
class meshObjects;

/*---------------------------------------------------------------------------*\
                    Class DeletableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class DeletableMeshObject
{
    //- Reference to the regIOobject of the base-class
    regIOobject& io_;

    //- The regIOobject reference is used by the meshObjects functions
    friend class meshObjects;

public:

    template<class Type>
    DeletableMeshObject(Type& mo)
    :
        io_(mo)
    {}

    //- Virtual destructor to make class polymorphic
    virtual ~DeletableMeshObject() = default;
};


/*---------------------------------------------------------------------------*\
                     Class MoveableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class MoveableMeshObject
:
    public DeletableMeshObject<Mesh>
{
public:

    template<class Type>
    MoveableMeshObject(Type& mo)
    :
        DeletableMeshObject<Mesh>(mo)
    {}

    //- Update for mesh motion
    virtual bool movePoints() = 0;
};


/*---------------------------------------------------------------------------*\
                   Class DistributeableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class DistributeableMeshObject
:
    public MoveableMeshObject<Mesh>
{
public:

    template<class Type>
    DistributeableMeshObject(Type& mo)
    :
        MoveableMeshObject<Mesh>(mo)
    {}

    //- Redistribute or update using the given distribution map
    virtual void distribute(const polyDistributionMap& map) = 0;
};


/*---------------------------------------------------------------------------*\
                    Class TopoChangeableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class TopoChangeableMeshObject
:
    public DistributeableMeshObject<Mesh>
{
public:

    template<class Type>
    TopoChangeableMeshObject(Type& mo)
    :
        DistributeableMeshObject<Mesh>(mo)
    {}

    //- Update topology using the given map
    virtual void topoChange(const polyTopoChangeMap& map) = 0;

    //- Update from another mesh using the given map
    virtual void mapMesh(const polyMeshMap& map) = 0;
};


/*---------------------------------------------------------------------------*\
                       Class RepatchableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class RepatchableMeshObject
:
    public TopoChangeableMeshObject<Mesh>
{
public:

    template<class Type>
    RepatchableMeshObject(Type& mo)
    :
        TopoChangeableMeshObject<Mesh>(mo)
    {}

    //- Reordered/removed trailing patches. If validBoundary call is parallel
    //  synced and all add the same patch with same settings
    virtual void reorderPatches
    (
        const labelUList& newToOld,
        const bool validBoundary
    ) = 0;

    //- Inserted patch at patchi
    virtual void addPatch(const label patchi) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
