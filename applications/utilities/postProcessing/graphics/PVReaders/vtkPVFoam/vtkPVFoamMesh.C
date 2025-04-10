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

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "fvMeshSubset.H"
#include "LagrangianMesh.H"
#include "vtkPVFoamReader.h"
#include "uindirectPrimitivePatch.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::convertMeshVolume
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeVolume_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    // Resize for decomposed polyhedra
    regionPolyDecomp_.setSize(range.size());

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    // Convert the internalMesh
    // This looks like more than one part, but it isn't
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word partName = "internalMesh";

        if (!partStatus_[partId])
        {
            continue;
        }

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            mesh,
            regionPolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshlagrangian
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangelagrangian_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word cloudName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        vtkPolyData* vtkmesh = lagrangianVTKMesh(mesh, cloudName);

        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, cloudName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshLagrangian
(
    vtkMultiBlockDataSet* output,
    int& blockNo,
    PtrList<LagrangianMesh>& LmeshPtrs
)
{
    arrayRange& range = arrayRangeLagrangian_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    LmeshPtrs.clear();

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word LagrangianName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        autoPtr<LagrangianMesh> LmeshPtr;

        vtkPolyData* vtkmesh =
            LagrangianVTKMesh(mesh, LagrangianName, LmeshPtr);

        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, LagrangianName);

            vtkmesh->Delete();

            LmeshPtrs.append(LmeshPtr);

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPatches
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangePatches_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        if (!partStatus_[partId])
        {
            continue;
        }

        const word patchName = getPartName(partId);

        labelHashSet
            patchIds(patches.patchSet(List<wordRe>(1, wordRe(patchName))));

        if (debug)
        {
            InfoInFunction
                << "Creating VTK mesh for patches [" << patchIds <<"] "
                << patchName << endl;
        }

        vtkPolyData* vtkmesh = nullptr;
        if (patchIds.size() == 1)
        {
            vtkmesh = patchVTKMesh(patchName, patches[patchIds.begin().key()]);
        }
        else
        {
            // Patch group. Collect patch faces.
            label sz = 0;
            forAllConstIter(labelHashSet, patchIds, iter)
            {
                sz += patches[iter.key()].size();
            }
            labelList meshFaceLabels(sz);
            sz = 0;
            forAllConstIter(labelHashSet, patchIds, iter)
            {
                const polyPatch& pp = patches[iter.key()];
                forAll(pp, i)
                {
                    meshFaceLabels[sz++] = pp.start()+i;
                }
            }
            UIndirectList<face> fcs(mesh.faces(), meshFaceLabels);
            uindirectPrimitivePatch pp(fcs, mesh.points());

            vtkmesh = patchVTKMesh(patchName, pp);
        }


        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, patchName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshCellZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeCellZones_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    // Resize for decomposed polyhedra
    zonePolyDecomp_.setSize(range.size());

    if (range.empty())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    const cellZoneList& zMesh = mesh.cellZones();
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word zoneName = getPartName(partId);
        const label  zoneId = zMesh.findIndex(zoneName);

        if (!partStatus_[partId] || zoneId < 0)
        {
            continue;
        }

        if (debug)
        {
            InfoInFunction
                << "Creating VTK mesh for cellZone[" << zoneId << "] "
                << zoneName << endl;
        }

        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(zMesh[zoneId]);

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            subsetter.subMesh(),
            zonePolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            // superCells + addPointCellLabels must contain global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                zonePolyDecomp_[datasetNo].superCells()
            );
            inplaceRenumber
            (
                subsetter.cellMap(),
                zonePolyDecomp_[datasetNo].addPointCellLabels()
            );

            // Copy pointMap as well, otherwise pointFields fail
            zonePolyDecomp_[datasetNo].pointMap() = subsetter.pointMap();

            AddToBlock(output, vtkmesh, range, datasetNo, zoneName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshCellSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeCellSets_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    // Resize for decomposed polyhedra
    csetPolyDecomp_.setSize(range.size());

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word partName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        if (debug)
        {
            InfoInFunction
                << "Creating VTK mesh for cellSet=" << partName << endl;
        }

        const cellSet cSet(mesh, partName);
        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(cSet);

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            subsetter.subMesh(),
            csetPolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            // superCells + addPointCellLabels must contain global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                csetPolyDecomp_[datasetNo].superCells()
            );
            inplaceRenumber
            (
                subsetter.cellMap(),
                csetPolyDecomp_[datasetNo].addPointCellLabels()
            );

            // Copy pointMap as well, otherwise pointFields fail
            csetPolyDecomp_[datasetNo].pointMap() = subsetter.pointMap();

            AddToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeFaceZones_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (range.empty())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    const faceZoneList& zMesh = mesh.faceZones();
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word zoneName = getPartName(partId);
        const label  zoneId = zMesh.findIndex(zoneName);

        if (!partStatus_[partId] || zoneId < 0)
        {
            continue;
        }

        if (debug)
        {
            InfoInFunction
                << "Creating VTKmesh for faceZone[" << zoneId << "] "
                << zoneName << endl;
        }

        vtkPolyData* vtkmesh = patchVTKMesh(zoneName, zMesh[zoneId].patch());

        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, zoneName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeFaceSets_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word partName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        if (debug)
        {
            InfoInFunction
                << "Creating VTK mesh for faceSet=" << partName << endl;
        }

        const faceSet fSet(mesh, partName);

        vtkPolyData* vtkmesh = faceSetVTKMesh(mesh, fSet);
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPointZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangePointZones_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    if (range.size())
    {
        const pointZoneList& zMesh = mesh.pointZones();
        for (int partId = range.start(); partId < range.end(); ++partId)
        {
            word zoneName = getPartName(partId);
            label zoneId = zMesh.findIndex(zoneName);

            if (!partStatus_[partId] || zoneId < 0)
            {
                continue;
            }

            vtkPolyData* vtkmesh = pointZoneVTKMesh(mesh, zMesh[zoneId]);
            if (vtkmesh)
            {
                AddToBlock(output, vtkmesh, range, datasetNo, zoneName);
                vtkmesh->Delete();

                partDataset_[partId] = datasetNo++;
            }
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}



void Foam::vtkPVFoam::convertMeshPointSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangePointSets_;
    range.block(blockNo);      // Set output block
    label datasetNo = 0;       // Restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        InfoInFunction << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        word partName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        if (debug)
        {
            InfoInFunction
                << "Creating VTK mesh for pointSet=" << partName << endl;
        }

        const pointSet pSet(mesh, partName);

        vtkPolyData* vtkmesh = pointSetVTKMesh(mesh, pSet);
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // Anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        printMemory();
    }
}


// ************************************************************************* //
