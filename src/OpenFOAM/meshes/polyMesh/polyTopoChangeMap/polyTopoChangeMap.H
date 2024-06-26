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
    Foam::polyTopoChangeMap

Description
    Class containing mesh-to-mesh mapping information after a change
    in polyMesh topology.

    General:
        - pointMap/faceMap/cellMap: \n
          from current mesh back to previous mesh.
          (so to 'pull' the information onto the current mesh)
        - reversePointMap/faceMap/cellMap: \n
          from previous mesh to current. (so to 'push' information)

    In the topology change points/faces/cells
    - can be unchanged. (faces might be renumbered though)
    - can be removed (into nothing)
    - can be removed into/merged with existing same entity
      (so point merged with other point, face with other face, cell with
       other cell. Note that probably only cell with cell is relevant)
    - can be added from existing same 'master' entity
      (so point from point, face from face and cell from cell)
    - can be appended: added 'out of nothing'.

    All this information is necessary to correctly map fields.

    \par points

    - unchanged:
        - pointMap[pointi] contains old point label
        - reversePointMap[oldPointi] contains new point label
    - removed:
        - reversePointMap[oldPointi] contains -1
    - merged into point:
        - reversePointMap[oldPointi] contains <-1 : -newPointi-2
        - pointMap[pointi] contains the old master point label
        - pointsFromPoints gives for pointi all the old point labels
          (including the old master point!)
    - added-from-same:
        - pointMap[pointi] contains the old master point label
    - appended:
        - pointMap[pointi] contains -1

    \par faces

    - unchanged:
        - faceMap[facei] contains old face label
        - reverseFaceMap[oldFacei] contains new face label
    - removed:
        - reverseFaceMap[oldFacei] contains -1
    - merged into face:
        - reverseFaceMap[oldFacei] contains <-1 : -newFacei-2
        - faceMap[facei] contains the old master face label
        - facesFromFaces gives for facei all the old face labels
          (including the old master face!)
    - added-from-same:
        - faceMap[facei] contains the old master face label
    - appended:
        - faceMap[facei] contains -1

    \par cells

    - unchanged:
        - cellMap[celli] contains old cell label
        - reverseCellMap[oldCelli] contains new cell label
    - removed:
        - reverseCellMap[oldCelli] contains -1
    - merged into cell:
        - reverseCellMap[oldCelli] contains <-1 : -newCelli-2
        - cellMap[celli] contains the old master cell label
        - cellsFromCells gives for celli all the old cell labels
          (including the old master cell!)
    - added-from-same:
        - cellMap[celli] contains the old master cell label
    - appended:
        - cellMap[celli] contains -1

SourceFiles
    polyTopoChangeMap.C

\*---------------------------------------------------------------------------*/

#ifndef polyTopoChangeMap_H
#define polyTopoChangeMap_H

#include "labelList.H"
#include "objectMap.H"
#include "pointField.H"
#include "HashSet.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                         Class polyTopoChangeMap Declaration
\*---------------------------------------------------------------------------*/

class polyTopoChangeMap
{
    // Private Data

        //- Reference to polyMesh
        const polyMesh& mesh_;

        //- Number of old live points
        const label nOldPoints_;

        //- Number of old live faces
        const label nOldFaces_;

        //- Number of old live cells
        const label nOldCells_;

        //- Old point map.
        //  Contains the old point label for all new points.
        //  - for preserved points this is the old point label.
        //  - for added points this is the master point ID
        //  - for points added with no master, this is -1
        //  Size of the list equals the size of new points
        const labelList pointMap_;

        //- Points resulting from merging points
        const List<objectMap> pointsFromPointsMap_;

        //- Old face map.
        //  Contains a list of old face labels for every new face.
        //  Size of the list equals the number of new faces
        //  - for preserved faces this is the old face label.
        //  - for faces added from faces this is the master face ID
        //  - for faces added with no master, this is -1
        //  - for faces added from points or edges, this is -1
        const labelList faceMap_;

        //- Faces resulting from merging faces
        const List<objectMap> facesFromFacesMap_;

        //- Old cell map.
        //  Contains old cell label for all preserved cells.
        //  Size of the list equals the number or preserved cells
        const labelList cellMap_;

        //- Cells resulting from merging cells
        const List<objectMap> cellsFromCellsMap_;

        //- Reverse point map
        const labelList reversePointMap_;

        //- Reverse face map
        const labelList reverseFaceMap_;

        //- Reverse cell map
        const labelList reverseCellMap_;

        //- Map of flipped face flux faces
        const labelHashSet flipFaceFlux_;

        //- Patch mesh point renumbering
        const labelListList patchPointMap_;

        //- List of the old patch sizes
        const labelList oldPatchSizes_;

        //- List of the old patch start labels
        const labelList oldPatchStarts_;

        //- List of numbers of mesh points per old patch
        const labelList oldPatchNMeshPoints_;

        //- Optional old cell volumes (for mapping)
        autoPtr<scalarField> oldCellVolumesPtr_;


public:

    // Constructors

        //- Construct from mesh only, no topology change
        polyTopoChangeMap(const polyMesh& mesh);

        //- Construct from components reusing storage
        polyTopoChangeMap
        (
            const polyMesh& mesh,
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            labelList&& pointMap,
            List<objectMap>&& pointsFromPoints,
            labelList&& faceMap,
            List<objectMap>&& facesFromFaces,
            labelList&& cellMap,
            List<objectMap>&& cellsFromCells,
            labelList&& reversePointMap,
            labelList&& reverseFaceMap,
            labelList&& reverseCellMap,
            labelHashSet&& flipFaceFlux,
            labelListList&& patchPointMap,
            labelList&& oldPatchSizes,
            labelList&& oldPatchStarts,
            labelList&& oldPatchNMeshPoints,
            autoPtr<scalarField>&& oldCellVolumesPtr
        );

        //- Disallow default bitwise copy construction
        polyTopoChangeMap(const polyTopoChangeMap&) = delete;


    // Member Functions

        // Access

            //- Return polyMesh
            const polyMesh& mesh() const
            {
                return mesh_;
            }

            //- Number of old points
            label nOldPoints() const
            {
                return nOldPoints_;
            }

            //- Number of old internal faces
            label nOldInternalFaces() const
            {
                return oldPatchStarts_[0];
            }

            //- Number of old faces
            label nOldFaces() const
            {
                return nOldFaces_;
            }

            //- Number of old cells
            label nOldCells() const
            {
                return nOldCells_;
            }

            //- Old point map.
            //  Contains the old point label for all new points.
            //  For preserved points this is the old point label.
            //  For added points this is the master point ID
            const labelList& pointMap() const
            {
                return pointMap_;
            }

            //- Points originating from points
            const List<objectMap>& pointsFromPointsMap() const
            {
                return pointsFromPointsMap_;
            }

            //- Old face map.
            //  Contains a list of old face labels for every new face.
            //  Warning: this map contains invalid entries for new faces
            const labelList& faceMap() const
            {
                return faceMap_;
            }

            //- Faces originating from faces
            const List<objectMap>& facesFromFacesMap() const
            {
                return facesFromFacesMap_;
            }

            //- Old cell map.
            //  Contains old cell label for all preserved cells.
            const labelList& cellMap() const
            {
                return cellMap_;
            }

            //- Cells originating from cells
            const List<objectMap>& cellsFromCellsMap() const
            {
                return cellsFromCellsMap_;
            }


            // Reverse maps

                //- Reverse point map
                //  Contains new point label for all old and added points
                const labelList& reversePointMap() const
                {
                    return reversePointMap_;
                }

                //- If point is removed return point (on new mesh) it merged
                //  into
                label mergedPoint(const label oldPointi) const
                {
                    label i = reversePointMap_[oldPointi];

                    if (i == -1)
                    {
                        return i;
                    }
                    else if (i < -1)
                    {
                        return -i-2;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "old point label " << oldPointi
                            << " has reverseMap " << i << endl
                            << "Only call mergedPoint for removed points."
                            << abort(FatalError);
                        return -1;
                    }
                }

                //- Reverse face map
                //  Contains new face label for all old and added faces
                const labelList& reverseFaceMap() const
                {
                    return reverseFaceMap_;
                }

                //- If face is removed return face (on new mesh) it merged into
                label mergedFace(const label oldFacei) const
                {
                    label i = reverseFaceMap_[oldFacei];

                    if (i == -1)
                    {
                        return i;
                    }
                    else if (i < -1)
                    {
                        return -i-2;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "old face label " << oldFacei
                            << " has reverseMap " << i << endl
                            << "Only call mergedFace for removed faces."
                            << abort(FatalError);
                        return -1;
                    }
                }

                //- Reverse cell map
                //  Contains new cell label for all old and added cells
                const labelList& reverseCellMap() const
                {
                    return reverseCellMap_;
                }

                //- If cell is removed return cell (on new mesh) it merged into
                label mergedCell(const label oldCelli) const
                {
                    label i = reverseCellMap_[oldCelli];

                    if (i == -1)
                    {
                        return i;
                    }
                    else if (i < -1)
                    {
                        return -i-2;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "old cell label " << oldCelli
                            << " has reverseMap " << i << endl
                            << "Only call mergedCell for removed cells."
                            << abort(FatalError);
                        return -1;
                    }
                }

                //- Map of flipped face flux faces
                const labelHashSet& flipFaceFlux() const
                {
                    return flipFaceFlux_;
                }

                //- Patch point renumbering
                //  For every preserved point on a patch give the old position.
                //  For added points, the index is set to -1
                const labelListList& patchPointMap() const
                {
                    return patchPointMap_;
                }


            //- Return list of the old patch sizes
            const labelList& oldPatchSizes() const
            {
                return oldPatchSizes_;
            }

            //- Return list of the old patch start labels
            const labelList& oldPatchStarts() const
            {
                return oldPatchStarts_;
            }

            //- Return numbers of mesh points per old patch
            const labelList& oldPatchNMeshPoints() const
            {
                return oldPatchNMeshPoints_;
            }


            // Geometric mapping data

                bool hasOldCellVolumes() const
                {
                    return oldCellVolumesPtr_.valid();
                }

                const scalarField& oldCellVolumes() const
                {
                    return oldCellVolumesPtr_();
                }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polyTopoChangeMap&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
