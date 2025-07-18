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

InClass
    vtkPVFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoamSurfaceField_H
#define vtkPVFoamSurfaceField_H

#include "vtkPVFoam.H"
#include "vtkPVFoamReader.h"
#include "vtkOpenFOAMTupleRemap.H"

// OpenFOAM includes
#include "domainDecomposition.H"
#include "emptyFvsPatchField.H"
#include "faceSet.H"
#include "surfaceFields.H"
#include "fvFieldReconstructor.H"

// VTK includes
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPVFoam::convertSurfaceField
(
    const VolField<Type>& tf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const label datasetNo,
    const fvMesh& mesh,
    const labelList& faceLabels
)
{
    const label nComp = pTraits<Type>::nComponents;
    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray* cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(faceLabels.size());
    cellData->SetNumberOfComponents(nComp);
    cellData->Allocate(nComp*faceLabels.size());
    cellData->SetName(tf.name().c_str());

    DebugInFunction
        << "Converting Surface field: " << tf.name()
        << " size=" << tf.size() << " (" << faceLabels.size()
        << "), nComp=" << nComp << endl;

    float vec[nComp];

    // For interior faces: average owner/neighbour
    // For boundary faces: owner
    forAll(faceLabels, i)
    {
        const label facei = faceLabels[i];
        if (facei < nInternalFaces)
        {
            const Type t = 0.5*(tf[faceOwner[facei]] + tf[faceNeigh[facei]]);

            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
        }
        else
        {
            const label patchi = mesh.boundaryMesh().whichPatch(facei);
            const label pFacei = mesh.boundaryMesh()[patchi].whichFace(facei);
            const Type& t = tf.boundaryField()[patchi][pFacei];

            for (direction d=0; d<nComp; ++d)
            {
                vec[d] = component(t, d);
            }
        }
        vtkOpenFOAMTupleRemap<Type>(vec);

        cellData->InsertTuple(i, vec);
    }


    vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, range, datasetNo)
    )   ->GetCellData()
        ->AddArray(cellData);

    cellData->Delete();
}


template<class Type>
void Foam::vtkPVFoam::convertSurfaceField
(
    const SurfaceField<Type>& tf,
    vtkMultiBlockDataSet* output,
    const arrayRange& range,
    const label datasetNo,
    const fvMesh& mesh,
    const labelList& faceLabels
)
{
    const label nComp = pTraits<Type>::nComponents;
    const label nInternalFaces = mesh.nInternalFaces();

    vtkFloatArray* cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(faceLabels.size());
    cellData->SetNumberOfComponents(nComp);
    cellData->Allocate(nComp*faceLabels.size());
    cellData->SetName(tf.name().c_str());

    DebugInFunction
        << "Converting Surface field: " << tf.name()
        << " size=" << tf.size() << " (" << faceLabels.size()
        << "), nComp=" << nComp << endl;

    // To avoid whichPatch first flatten the field
    Field<Type> flatFld(mesh.nFaces(), Zero);
    SubField<Type>(flatFld, nInternalFaces) = tf.internalField();
    forAll(tf.boundaryField(), patchi)
    {
        const fvsPatchField<Type>& fvs = tf.boundaryField()[patchi];

        SubField<Type>
        (
            flatFld,
            fvs.size(),
            fvs.patch().start()
        ) = fvs;
    }

    forAll(faceLabels, i)
    {
        const label facei = faceLabels[i];

        float vec[nComp];
        for (direction d=0; d<nComp; ++d)
        {
            vec[d] = component(flatFld[facei], d);
        }
        vtkOpenFOAMTupleRemap<Type>(vec);

        cellData->InsertTuple(i, vec);
    }

    vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, range, datasetNo)
    )   ->GetCellData()
        ->AddArray(cellData);

    cellData->Delete();
}


template<class Type>
void Foam::vtkPVFoam::convertSurfaceFields
(
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output
)
{
    forAllConstIter(IOobjectList, objects, iter)
    {
        // Restrict to GeometricField<Type, ...>
        if (iter()->headerClassName() != SurfaceField<Type>::typeName)
        {
            continue;
        }

        // Load the field
        tmp<SurfaceField<Type>> ttf;

        FatalIOError.throwExceptions();

        try
        {
            if (reader_->GetDecomposedCase())
            {
                if (!fvReconstructorPtr_.valid())
                {
                    fvReconstructorPtr_.set
                    (
                        new fvFieldReconstructor
                        (
                            procMeshesPtr_->completeMesh(),
                            procMeshesPtr_->procMeshes(),
                            procMeshesPtr_->procFaceAddressing(),
                            procMeshesPtr_->procCellAddressing(),
                            procMeshesPtr_->procFaceAddressingBf()
                        )
                    );
                }

                ttf =
                    fvReconstructorPtr_
                  ->reconstructFvSurfaceField<Type>(*iter());
            }
            else
            {
                ttf =
                    new SurfaceField<Type>
                    (
                        *iter(),
                        procMeshesPtr_->completeMesh()
                    );
            }
        }
        catch (IOerror& err)
        {
            Warning<< err << endl;
            continue;
        }

        FatalIOError.dontThrowExceptions();

        const SurfaceField<Type>& tf = ttf();

        // Convert patches - if activated
        for
        (
            int partId = arrayRangePatches_.start();
            partId < arrayRangePatches_.end();
            ++partId
        )
        {
            const word patchName = getPartName(partId);
            const label datasetNo = partDataset_[partId];
            const label patchId = tf.mesh().boundaryMesh().findIndex(patchName);

            if (!partStatus_[partId] || datasetNo < 0 || patchId < 0)
            {
                continue;
            }

            const fvsPatchField<Type>& ptf = tf.boundaryField()[patchId];

            if (!isType<emptyFvsPatchField<Type>>(ptf))
            {
                convertPatchField
                (
                    tf.name(),
                    ptf,
                    output,
                    arrayRangePatches_,
                    datasetNo
                );
            }
        }

        // Convert face zones - if activated
        for
        (
            int partId = arrayRangeFaceZones_.start();
            partId < arrayRangeFaceZones_.end();
            ++partId
        )
        {
            const word zoneName = getPartName(partId);
            const label datasetNo = partDataset_[partId];

            if (!partStatus_[partId] || datasetNo < 0)
            {
                continue;
            }

            const faceZoneList& zMesh = tf.mesh().faceZones();
            const label zoneId = zMesh.findIndex(zoneName);

            if (zoneId < 0)
            {
                continue;
            }

            convertSurfaceField
            (
                tf,
                output,
                arrayRangeFaceZones_,
                datasetNo,
                tf.mesh(),
                zMesh[zoneId]
            );
        }

        // Convert face sets - if activated
        for
        (
            int partId = arrayRangeFaceSets_.start();
            partId < arrayRangeFaceSets_.end();
            ++partId
        )
        {
            const word selectName = getPartName(partId);
            const label datasetNo = partDataset_[partId];

            if (!partStatus_[partId] || datasetNo < 0)
            {
                continue;
            }

            const autoPtr<faceSet> fSetPtr =
                reader_->GetDecomposedCase()
              ? procMeshesPtr_->reconstructSet<faceSet>(selectName)
              : autoPtr<faceSet>(new faceSet(tf.mesh(), selectName));

            convertSurfaceField
            (
                tf,
                output,
                arrayRangeFaceSets_,
                datasetNo,
                tf.mesh(),
                fSetPtr().toc()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
