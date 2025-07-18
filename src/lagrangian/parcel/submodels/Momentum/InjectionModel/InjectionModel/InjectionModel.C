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

#include "InjectionModel.H"
#include "meshTools.H"
#include "meshSearch.H"
#include "volFields.H"
#include "Scale.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::readMassTotal
(
    const dictionary& dict,
    CloudType& owner
)
{
    if (dict.found("nParticle"))
    {
        if (dict.found("massTotal"))
        {
            IOWarningInFunction(dict)
                << "If nParticle is specified then the massTotal "
                << "setting has no effect " << endl;
        }

        return NaN;
    }

    if (owner.solution().steadyState())
    {
        FatalErrorInFunction
            << "The " << type() << " injection model is not compatible with "
            << "steady state solution"
            << exit(FatalError);

        return NaN;
    }

    return dict.lookup<scalar>("massTotal", dimMass);
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::readDuration
(
    const dictionary& dict,
    CloudType& owner
)
{
    if (owner.solution().steadyState())
    {
        return vGreat;
    }

    return dict.lookup<scalar>("duration", owner.db().time().userUnits());
}


template<class CloudType>
Foam::autoPtr<Foam::Function1<Foam::scalar>>
Foam::InjectionModel<CloudType>::readMassFlowRate
(
    const dictionary& dict,
    CloudType& owner,
    const scalar duration
)
{
    const bool haveMassFlowRate = dict.found("massFlowRate");
    const bool haveMassTotal = dict.found("massTotal");

    if (dict.found("nParticle"))
    {
        if (haveMassFlowRate || haveMassTotal)
        {
            IOWarningInFunction(dict)
                << "If nParticle is specified then massFlowRate and massTotal "
                << "settings have no effect " << endl;
        }

        return
            autoPtr<Function1<scalar>>
            (
                new Function1s::Constant<scalar>("NaN", NaN)
            );
    }

    if (owner.solution().steadyState() && haveMassTotal)
    {
        FatalIOErrorInFunction(dict)
            << "Cannot specify the massTotal of a steady injection. Use "
            << "massFlowRate instead." << exit(FatalIOError);
    }

    if (haveMassFlowRate && haveMassTotal)
    {
        FatalIOErrorInFunction(dict)
            << "Cannot specify both massFlowRate and massTotal. Use one or "
            << "the other." << exit(FatalIOError);
    }

    if (owner.solution().steadyState() || haveMassFlowRate)
    {
        return
            Function1<scalar>::New
            (
                "massFlowRate",
                this->owner().db().time().userUnits(),
                dimMass/dimTime,
                dict
            );
    }

    const scalar massTotal = dict.lookup<scalar>("massTotal", dimMass);

    if (!dict.found("flowRateProfile"))
    {
        return
            autoPtr<Function1<scalar>>
            (
                new Function1s::Constant<scalar>
                (
                    "massFlowRate",
                    massTotal/duration
                )
            );
    }

    autoPtr<Function1<scalar>> flowRateProfile
    (
        Function1<scalar>::New
        (
            "flowRateProfile",
            this->owner().db().time().userUnits(),
            dimless,
            dict
        )
    );

    const scalar sumFlowRateProfile = flowRateProfile->integral(0, duration);

    return
        autoPtr<Function1<scalar>>
        (
            new Function1s::Scale<scalar>
            (
                "massFlowRate",
                Function1s::Constant<scalar>("m", massTotal/sumFlowRateProfile),
                Function1s::Constant<scalar>("one", scalar(1)),
                flowRateProfile()
            )
        );
}


template<class CloudType>
Foam::autoPtr<Foam::Function1<Foam::scalar>>
Foam::InjectionModel<CloudType>::readParcelsPerSecond
(
    const dictionary& dict,
    CloudType& owner
)
{
    return
        Function1<scalar>::New
        (
            "parcelsPerSecond",
            this->owner().db().time().userUnits(),
            dimless/dimTime,
            dict
        );
}


template<class CloudType>
Foam::label Foam::InjectionModel<CloudType>::index() const
{
    forAll(this->owner().injectors(), i)
    {
        if (this->owner().injectors()(i) == this)
        {
            return i;
        }
    }

    return -1;
}


template<class CloudType>
bool Foam::InjectionModel<CloudType>::findCellAtPosition
(
    const meshSearch& searchEngine,
    const point& position,
    barycentric& coordinates,
    label& celli,
    label& tetFacei,
    label& tetPti,
    bool errorOnNotFound
)
{
    // Subroutine for finding the cell
    auto findProcAndCell = [this,&searchEngine](const point& pos)
    {
        // Find the containing cell
        label celli = searchEngine.findCell(pos);

        // Synchronise so only a single processor finds this position
        label proci = celli >= 0 ? Pstream::myProcNo() : -1;
        reduce(proci, maxOp<label>());
        if (proci != Pstream::myProcNo())
        {
            celli = -1;
        }

        return labelPair(proci, celli);
    };

    point pos = position;

    // Try and find the cell at the given position
    const labelPair procAndCelli = findProcAndCell(pos);
    label proci = procAndCelli.first();
    celli = procAndCelli.second();

    // Didn't find it. The point may be awkwardly on an edge or face. Try
    // again, but move the point into its nearest cell a little bit.
    if (proci == -1)
    {
        pos += small*(this->owner().mesh().C()[celli] - pos);
        const labelPair procAndCelli = findProcAndCell(pos);
        proci = procAndCelli.first();
        celli = procAndCelli.second();
    }

    // Didn't find it. Error or return false.
    if (proci == -1)
    {
        if (errorOnNotFound)
        {
            FatalErrorInFunction
                << "Cannot find parcel injection cell. "
                << "Parcel position = " << position << nl
                << exit(FatalError);
        }

        return false;
    }

    // Found it. Construct the barycentric coordinates.
    if (proci == Pstream::myProcNo())
    {
        label nLocateBoundaryHits = 0;
        particle p(searchEngine, pos, celli, nLocateBoundaryHits);

        if (nLocateBoundaryHits != 0)
        {
            WarningInFunction
                << "Injection model " << this->modelName()
                << " for cloud " << this->owner().name()
                << " did not accurately locate the position "
                << pos << " within the mesh" << endl;
        }

        coordinates = p.coordinates();
        celli = p.cell();
        tetFacei = p.tetFace();
        tetPti = p.tetPt();
    }

    return true;
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::constrainPosition
(
    typename CloudType::parcelType::trackingData& td,
    typename CloudType::parcelType& parcel
)
{
    const vector d = parcel.deviationFromMeshCentre(td.mesh);

    if (d == vector::zero)
    {
        return;
    }

    const label facei = parcel.face();

    // If the parcel is not on a face, then just track it to the mesh centre
    if (facei == -1)
    {
        parcel.track(td.mesh, - d, 0);
    }

    // If the parcel is on a face, then track in two steps, going slightly into
    // the current cell. This prevents a boundary hit from ending the track
    // prematurely.
    if (facei != -1)
    {
        const vector pc =
            td.mesh.cellCentres()[parcel.cell()] - parcel.position(td.mesh);

        parcel.track(td.mesh, - d/2 + rootSmall*pc, 0);
        parcel.track(td.mesh, - d/2 - rootSmall*pc, 0);
    }

    // Restore any face-association that got changed during tracking
    parcel.face() = facei;
}


template<class CloudType>
Foam::label Foam::InjectionModel<CloudType>::sizeSampleQ() const
{
    switch (uniformParcelSize_)
    {
        case uniformParcelSize::nParticle:
            return 0;
        case uniformParcelSize::surfaceArea:
            return 2;
        case uniformParcelSize::volume:
            return 3;
    }

    return -labelMax;
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::setNumberOfParticles
(
    PtrList<parcelType>& parcelPtrs,
    const scalar mass
) const
{
    auto size = [&](const parcelType& p)
    {
        switch (uniformParcelSize_)
        {
            case uniformParcelSize::nParticle:
                return scalar(1);
            case uniformParcelSize::surfaceArea:
                return p.areaS();
            case uniformParcelSize::volume:
                return p.volume();
        }
        return NaN;
    };

    // Determine the total mass and size of all created particles
    scalar sumMassBySize = 0;
    forAll(parcelPtrs, parceli)
    {
        if (parcelPtrs.set(parceli))
        {
            const parcelType& p = parcelPtrs[parceli];
            sumMassBySize += p.mass()/size(p);
        }
    }

    reduce(sumMassBySize, sumOp<scalar>());

    // Set the numbers of particles on each parcel
    forAll(parcelPtrs, parceli)
    {
        if (parcelPtrs.set(parceli))
        {
            parcelType& p = parcelPtrs[parceli];
            p.nParticle() = mass/size(p)/sumMassBySize;
        }
    }

    // Check that the constraints are correct
    if (debug)
    {
        scalar massN = 0, minSizeN = vGreat, maxSizeN = -vGreat;
        forAll(parcelPtrs, parceli)
        {
            if (parcelPtrs.set(parceli))
            {
                const parcelType& p = parcelPtrs[parceli];
                massN += p.nParticle()*p.mass();
                minSizeN = min(minSizeN, p.nParticle()*size(p));
                maxSizeN = max(minSizeN, p.nParticle()*size(p));
            }
        }

        reduce(massN, sumOp<scalar>());

        if (mag(massN - mass) > rootSmall*(massN + mass)/2)
        {
            FatalErrorInFunction
                << "Parcels do not have the required mass"
                << exit(FatalError);
        }

        reduce(minSizeN, minOp<scalar>());
        reduce(maxSizeN, maxOp<scalar>());

        if (maxSizeN - minSizeN > rootSmall*(maxSizeN + minSizeN)/2)
        {
            FatalErrorInFunction
                << "Parcel sizes are not uniform"
                << exit(FatalError);
        }
    }
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::preInject
(
    typename parcelType::trackingData& td
)
{}


template<class CloudType>
void Foam::InjectionModel<CloudType>::postInject
(
    const label nParcelsAdded,
    const scalar massAdded,
    typename parcelType::trackingData& td
)
{
    const label allNParcelsAdded = returnReduce(nParcelsAdded, sumOp<label>());

    if (allNParcelsAdded > 0)
    {
        Info<< nl
            << "Cloud: " << this->owner().name()
            << " injector: " << this->modelName() << nl
            << "    Added " << allNParcelsAdded << " new parcels" << nl << endl;
    }

    // Increment total number of parcels added
    nParcelsInjected_ += allNParcelsAdded;

    // Increment total mass injected
    massInjected_ += returnReduce(massAdded, sumOp<scalar>());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    SOI_(0),
    massInjected_(this->template getModelProperty<scalar>("massInjected")),
    nParcelsInjected_
    (
        this->template getModelProperty<scalar>("nParcelsInjected")
    ),
    massDeferred_
    (
        this->template getModelProperty<scalar>("massDeferred")
    ),
    nParcelsDeferred_
    (
        this->template getModelProperty<scalar>("nParcelsDeferred")
    ),
    nParticleFixed_(-vGreat),
    uniformParcelSize_(uniformParcelSize::nParticle)
{}


template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName,
    const word& modelType
)
:
    CloudSubModelBase<CloudType>(modelName, owner, dict, typeName, modelType),
    SOI_(0),
    massInjected_(this->template getModelProperty<scalar>("massInjected")),
    nParcelsInjected_
    (
        this->template getModelProperty<scalar>("nParcelsInjected")
    ),
    massDeferred_
    (
        this->template getModelProperty<scalar>("massDeferred")
    ),
    nParcelsDeferred_
    (
        this->template getModelProperty<scalar>("nParcelsDeferred")
    ),
    nParticleFixed_(dict.lookupOrDefault<scalar>("nParticle", -vGreat)),
    uniformParcelSize_
    (
        uniformParcelSizeNames_
        [
            !dict.found("parcelBasisType") && nParticleFixed_ > 0
          ? dict.lookupOrDefault<word>
            (
                "uniformParcelSize",
                uniformParcelSizeNames_[uniformParcelSize::nParticle]
            )
          : dict.lookup<word>("uniformParcelSize")
        ]
    )
{
    // Provide some info. Also serves to initialise mesh dimensions. This may
    // be needed for parallel runs due to lazy evaluation of valid mesh
    // dimensions.
    Info<< "    Constructing " << owner.mesh().nGeometricD() << "-D injection"
        << endl;

    if
    (
        nParticleFixed_ > 0
     && uniformParcelSize_ != uniformParcelSize::nParticle
    )
    {
        FatalIOErrorInFunction(dict)
            << "If nParticle is specified then the uniformParcelSize must be "
            << uniformParcelSizeNames_[uniformParcelSize::nParticle]
            << exit(FatalIOError);
    }

    if (owner.solution().transient())
    {
        SOI_ = dict.lookup<scalar>("SOI", owner.db().time().userUnits());
    }
}


template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel
(
    const InjectionModel<CloudType>& im
)
:
    CloudSubModelBase<CloudType>(im),
    SOI_(im.SOI_),
    massInjected_(im.massInjected_),
    nParcelsInjected_(im.nParcelsInjected_),
    massDeferred_(im.massDeferred_),
    nParcelsDeferred_(im.nParcelsDeferred_),
    nParticleFixed_(im.nParticleFixed_),
    uniformParcelSize_(im.uniformParcelSize_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::~InjectionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::InjectionModel<CloudType>::topoChange()
{}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::averageParcelMass()
{
    const scalar deltaT =
        this->owner().solution().transient() ? timeEnd() - timeStart() : 1;

    return massToInject(0, deltaT)/nParcelsToInject(0, deltaT);
}


template<class CloudType>
template<class TrackCloudType>
void Foam::InjectionModel<CloudType>::inject
(
    TrackCloudType& cloud,
    typename CloudType::parcelType::trackingData& td
)
{
    const polyMesh& mesh = this->owner().mesh();

    const meshSearch& searchEngine = meshSearch::New(mesh);

    const scalar time1 = this->owner().db().time().value();
    const scalar time0 =
        this->owner().db().time().value()
      - this->owner().db().time().deltaTValue();

    preInject(td);

    // Reset counters
    label nParcelsAdded = 0;
    scalar massAdded = 0;

    // Get amounts to inject
    label nParcels;
    scalar mass = NaN;
    {
        // Injection has started. Get amounts between times.
        const scalar t0 = time0 - SOI_, t1 = time1 - SOI_;

        // Get the number of parcels to inject, round it down to the nearest
        // integer, and then store the excess to apply at a later time
        const scalar nParcelsNoRound =
            nParcelsToInject(t0, t1) + nParcelsDeferred_;
        nParcels = floor(nParcelsNoRound);
        nParcelsDeferred_ = nParcelsNoRound - nParcels;

        // If we need mass...
        if (nParticleFixed_ < 0)
        {
            // Get the mass to inject, and if there are no parcels then store
            // this to apply at a later time
            mass = massToInject(t0, t1) + massDeferred_;
            massDeferred_ = nParcels > 0 ? 0 : mass;

            // Special case. If there is no mass, then don't create anything.
            // This is a hack so that mass flow rates that switch on and off
            // still work tolerably with a constant number rate. Really,
            // though, cases like these should be specified with a
            // correspondingly varying number rate.
            if (mass == 0)
            {
                nParcels = 0;
                nParcelsDeferred_ = 0;
            }
        }
    }

    // Do injection
    if (nParcels > 0)
    {
        // Duration of injection period during this timestep
        const scalar deltaT =
            max
            (
                scalar(0),
                min
                (
                    td.trackTime(),
                    min
                    (
                        time1 - SOI_,
                        timeEnd() - time0
                    )
                )
            );

        // Pad injection time if injection starts during this timestep
        const scalar padTime = max(scalar(0), SOI_ - time0);

        // Create new parcels linearly across carrier phase timestep
        PtrList<parcelType> parcelPtrs(nParcels);
        forAll(parcelPtrs, parceli)
        {
            // Calculate the pseudo time of injection for parcel 'parceli'
            scalar timeInj = time0 + padTime + deltaT*parceli/nParcels;

            // Determine the injection coordinates and owner cell,
            // tetFace and tetPt
            barycentric coordinates = barycentric::uniform(NaN);
            label celli = -1, tetFacei = -1, tetPti = -1, facei = -1;
            setPositionAndCell
            (
                searchEngine,
                parceli,
                nParcels,
                timeInj,
                coordinates,
                celli,
                tetFacei,
                tetPti,
                facei
            );

            if (celli > -1)
            {
                // Lagrangian timestep
                const scalar dt = timeInj - time0;

                // Create a new parcel
                parcelPtrs.set
                (
                    parceli,
                    new parcelType
                    (
                        mesh,
                        coordinates,
                        celli,
                        tetFacei,
                        tetPti,
                        facei
                    )
                );
                parcelType& p = parcelPtrs[parceli];

                // Correct the position for reduced-dimension cases
                constrainPosition(td, p);

                // Check/set new parcel thermo properties
                cloud.setParcelThermoProperties(p);

                // Assign new parcel properties in injection model
                setProperties(parceli, nParcels, timeInj, td, p);

                // Check/set new parcel injection properties
                cloud.checkParcelProperties(p, index());

                // Apply correction to velocity for 2-D cases
                meshTools::constrainDirection
                (
                    mesh,
                    mesh.solutionD(),
                    p.U()
                );

                // Modify the step fraction so that the particles are
                // injected continually through the time-step
                p.stepFraction() = dt/td.trackTime();

                // Set the number of particles. If not fixed, this will set
                // a junk value, which will get corrected below.
                p.nParticle() = nParticleFixed_;
            }
        }

        // Set the number of particles so that the introduced mass is correct
        // and the uniform size is as specified
        if (nParticleFixed_ < 0)
        {
            setNumberOfParticles(parcelPtrs, mass);
        }

        // Add the new parcels
        forAll(parcelPtrs, parceli)
        {
            if (parcelPtrs.set(parceli))
            {
                parcelType& p = parcelPtrs[parceli];
                nParcelsAdded ++;
                massAdded += p.nParticle()*p.mass();
                cloud.addParticle(parcelPtrs.set(parceli, nullptr).ptr());
            }
        }
    }

    postInject(nParcelsAdded, massAdded, td);
}


template<class CloudType>
template<class TrackCloudType>
void Foam::InjectionModel<CloudType>::injectSteadyState
(
    TrackCloudType& cloud,
    typename CloudType::parcelType::trackingData& td
)
{
    const polyMesh& mesh = this->owner().mesh();

    const meshSearch& searchEngine = meshSearch::New(mesh);

    preInject(td);

    // Reset counters
    label nParcelsAdded = 0;
    scalar massAdded = 0;

    // Get amounts to inject based on first second of injection
    const label nParcels = floor(nParcelsToInject(0, 1));
    const scalar mass = nParticleFixed_ < 0 ? massToInject(0, 1) : NaN;

    // Do injection
    if (nParcels > 0)
    {
        PtrList<parcelType> parcelPtrs(nParcels);
        forAll(parcelPtrs, parceli)
        {
            // Determine the injection coordinates and owner cell,
            // tetFace and tetPt
            barycentric coordinates = barycentric::uniform(NaN);
            label celli = -1, tetFacei = -1, tetPti = -1, facei = -1;
            setPositionAndCell
            (
                searchEngine,
                parceli,
                nParcels,
                0,
                coordinates,
                celli,
                tetFacei,
                tetPti,
                facei
            );

            if (celli > -1)
            {
                // Create a new parcel
                parcelPtrs.set
                (
                    parceli,
                    new parcelType
                    (
                        mesh,
                        coordinates,
                        celli,
                        tetFacei,
                        tetPti,
                        facei
                    )
                );
                parcelType& p = parcelPtrs[parceli];

                // Correct the position for reduced-dimension cases
                constrainPosition(td, p);

                // Check/set new parcel thermo properties
                cloud.setParcelThermoProperties(p);

                // Assign new parcel properties in injection model
                setProperties(parceli, nParcels, 0, td, p);

                // Check/set new parcel injection properties
                cloud.checkParcelProperties(p, index());

                // Apply correction to velocity for 2-D cases
                meshTools::constrainDirection(mesh, mesh.solutionD(), p.U());

                // Initial step fraction
                p.stepFraction() = 0;

                // Set the number of particles. If not fixed, this will set
                // a junk value, which will get corrected below.
                p.nParticle() = nParticleFixed_;
            }
        }

        // Set the number of particles so that the introduced mass is correct
        // and the uniform size is as specified
        if (nParticleFixed_ < 0)
        {
            setNumberOfParticles(parcelPtrs, mass);
        }

        // Add the new parcels
        forAll(parcelPtrs, parceli)
        {
            if (parcelPtrs.set(parceli))
            {
                parcelType& p = parcelPtrs[parceli];
                nParcelsAdded ++;
                massAdded += p.nParticle()*p.mass();
                cloud.addParticle(parcelPtrs.set(parceli, nullptr).ptr());
            }
        }
    }

    postInject(nParcelsAdded, massAdded, td);
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::info(Ostream& os)
{
    os  << "    " << this->modelName() << ":" << nl
        << "        number of parcels added     = " << nParcelsInjected_ << nl
        << "        mass introduced             = " << massInjected_ << nl;

    if (this->writeTime())
    {
        this->setModelProperty("massInjected", massInjected_);
        this->setModelProperty("nParcelsInjected", nParcelsInjected_);
        this->setModelProperty("massDeferred", massDeferred_);
        this->setModelProperty("nParcelsDeferred", nParcelsDeferred_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InjectionModelNew.C"

// ************************************************************************* //
