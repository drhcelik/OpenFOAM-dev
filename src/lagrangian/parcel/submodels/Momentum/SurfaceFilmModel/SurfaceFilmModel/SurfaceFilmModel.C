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

#include "SurfaceFilmModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    g_(owner.g()),
    ejectedParcelType_(0),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    deltaFilmPatch_(0),
    nParcelsTransferred_(0),
    nParcelsInjected_(0)
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    g_(owner.g()),
    ejectedParcelType_
    (
        this->coeffDict().lookupOrDefault("ejectedParcelType", -1)
    ),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    deltaFilmPatch_(),
    nParcelsTransferred_(0),
    nParcelsInjected_(0)
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const SurfaceFilmModel<CloudType>& sfm
)
:
    CloudSubModelBase<CloudType>(sfm),
    g_(sfm.g_),
    ejectedParcelType_(sfm.ejectedParcelType_),
    massParcelPatch_(sfm.massParcelPatch_),
    diameterParcelPatch_(sfm.diameterParcelPatch_),
    deltaFilmPatch_(sfm.deltaFilmPatch_),
    nParcelsTransferred_(sfm.nParcelsTransferred_),
    nParcelsInjected_(sfm.nParcelsInjected_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::~SurfaceFilmModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
template<class TrackCloudType>
void Foam::SurfaceFilmModel<CloudType>::inject(TrackCloudType& cloud)
{
    const meshSearch& searchEngine = meshSearch::New(this->owner().mesh());

    const labelList& filmPatches = this->filmPatches();

    forAll(filmPatches, filmi)
    {
        const label filmPatchi = filmPatches[filmi];

        const fvMesh& mesh = this->owner().mesh();
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        const labelList& injectorCellsPatch = pbm[filmPatchi].faceCells();

        // Get and cache the properties of the droplets ejected from the film
        cacheFilmFields(filmi);

        const vectorField& Cf = mesh.C().boundaryField()[filmPatchi];
        const vectorField& Sf = mesh.Sf().boundaryField()[filmPatchi];
        const scalarField& magSf = mesh.magSf().boundaryField()[filmPatchi];

        label nLocateBoundaryHits = 0;

        if (massParcelPatch_.size())
        {
            forAll(injectorCellsPatch, j)
            {
                if (massParcelPatch_[j] > 0)
                {
                    const label celli = injectorCellsPatch[j];

                    const scalar offset = max
                    (
                        diameterParcelPatch_[j],
                        deltaFilmPatch_[j]
                    );

                    const point pos = Cf[j] - 1.1*offset*Sf[j]/magSf[j];

                    // Create a new parcel
                    parcelType* pPtr =
                        new parcelType
                        (
                            searchEngine,
                            pos,
                            celli,
                            nLocateBoundaryHits
                        );

                    // Check/set new parcel thermo properties
                    cloud.setParcelThermoProperties(*pPtr);

                    setParcelProperties(*pPtr, j);

                    if (pPtr->nParticle() > 0.001)
                    {
                        // Check new parcel properties
                        cloud.checkParcelProperties(*pPtr, -1);

                        // Add the new parcel to the cloud
                        cloud.addParticle(pPtr);

                        nParcelsInjected_++;
                    }
                    else
                    {
                        // TODO: cache mass and re-distribute?
                        delete pPtr;
                    }
                }
            }
        }

        reduce(nLocateBoundaryHits, sumOp<label>());
        if (nLocateBoundaryHits != 0)
        {
            WarningInFunction
                << "Injection by surface film model for cloud "
                << this->owner().name()
                << " on patch " << pbm[filmPatchi].name()
                << " did not accurately locate " << nLocateBoundaryHits
                << " particles" << endl;
        }
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::info(Ostream& os)
{
    label nTrans0 =
        this->template getModelProperty<label>("nParcelsTransferred");

    label nInject0 =
        this->template getModelProperty<label>("nParcelsInjected");

    label nTransTotal =
        nTrans0 + returnReduce(nParcelsTransferred_, sumOp<label>());

    label nInjectTotal =
        nInject0 + returnReduce(nParcelsInjected_, sumOp<label>());

    os  << "    Parcels absorbed into film      = " << nTransTotal << nl
        << "    New film detached parcels       = " << nInjectTotal << endl;

    if (this->writeTime())
    {
        this->setModelProperty("nParcelsTransferred", nTransTotal);
        this->setModelProperty("nParcelsInjected", nInjectTotal);
        nParcelsTransferred_ = 0;
        nParcelsInjected_ = 0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SurfaceFilmModelNew.C"

// ************************************************************************* //
