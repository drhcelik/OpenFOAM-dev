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

#include "polyMesh.H"
#include "primitiveMesh.H"
#include "globalMeshData.H"
#include "meshObjects.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMesh::removeBoundary()
{
    if (debug)
    {
        InfoInFunction << "Removing boundary patches." << endl;
    }

    // Remove the point zones
    boundary_.clear();
    boundary_.setSize(0);

    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::printAllocated() const
{
    primitiveMesh::printAllocated();

    Pout<< "polyMesh allocated :" << endl;

    if (tetBasePtIsPtr_.valid())
    {
        Pout<< "    Tet base points" << endl;
    }
}


void Foam::polyMesh::clearGeom()
{
    if (debug)
    {
        InfoInFunction << "Clearing geometric data" << endl;
    }

    // Clear all geometric mesh objects
    meshObjects::clear<pointMesh, DeletableMeshObject>(*this);
    meshObjects::clear<polyMesh, DeletableMeshObject>(*this);

    primitiveMesh::clearGeom();

    boundary_.clearGeom();

    // Reset valid directions (could change with rotation)
    geometricD_ = Zero;
    solutionD_ = Zero;
}


void Foam::polyMesh::clearAddressing(const bool isMeshUpdate)
{
    if (debug)
    {
        InfoInFunction
            << "Clearing topology  isMeshUpdate:" << isMeshUpdate << endl;
    }

    if (isMeshUpdate)
    {
        // Part of a mesh update. Keep meshObjects that have an topoChange
        // callback
        meshObjects::clearUpto
        <
            pointMesh,
            DeletableMeshObject,
            TopoChangeableMeshObject
        >
        (
            *this
        );
        meshObjects::clearUpto
        <
            polyMesh,
            DeletableMeshObject,
            TopoChangeableMeshObject
        >
        (
            *this
        );
    }
    else
    {
        meshObjects::clear<pointMesh, DeletableMeshObject>(*this);
        meshObjects::clear<polyMesh, DeletableMeshObject>(*this);
    }

    primitiveMesh::clearAddressing();

    // parallelData depends on the processorPatch ordering so force
    // recalculation
    globalMeshDataPtr_.clear();

    // Reset valid directions
    geometricD_ = Zero;
    solutionD_ = Zero;

    // Update zones
    pointZones_.clearAddressing();
    faceZones_.clearAddressing();
    cellZones_.clearAddressing();

    // Remove the stored tet base points
    tetBasePtIsPtr_.clear();
}


void Foam::polyMesh::clearPrimitives()
{
    resetMotion();

    points_.setSize(0);
    faces_.setSize(0);
    owner_.setSize(0);
    neighbour_.setSize(0);

    clearedPrimitives_ = true;
}


void Foam::polyMesh::clearOut()
{
    clearGeom();
    clearAddressing();
}


void Foam::polyMesh::clearTetBasePtIs()
{
    if (debug)
    {
        InfoInFunction << "Clearing tet base points" << endl;
    }

    tetBasePtIsPtr_.clear();
}


// ************************************************************************* //
