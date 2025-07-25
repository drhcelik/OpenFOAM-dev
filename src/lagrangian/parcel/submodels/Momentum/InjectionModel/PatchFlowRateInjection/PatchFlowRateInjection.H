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
    Foam::PatchFlowRateInjection

Description
    Patch injection in which the quantity of particulate material introduced is
    specified as a fraction of the carrier flow rate.

    User specifies:
      - Name of patch
      - Injection duration
      - The particulate/carrier volume or mass ratio
      - The number concentration of parcels

    Properties:
      - Initial parcel velocity given by local flow velocity
      - Parcel diameters obtained by distribution model
      - Parcels injected randomly across the patch

SourceFiles
    PatchFlowRateInjection.C

\*---------------------------------------------------------------------------*/

#ifndef PatchFlowRateInjection_H
#define PatchFlowRateInjection_H

#include "InjectionModel.H"
#include "patchInjectionBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distribution;

/*---------------------------------------------------------------------------*\
                   Class PatchFlowRateInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchFlowRateInjection
:
    public InjectionModel<CloudType>,
    public patchInjectionBase
{
    // Private Data

        //- Name of carrier (mass or volume) flux field
        const word phiName_;

        //- Name of carrier density field
        const word rhoName_;

        //- Injection duration [s]
        const scalar duration_;

        //- The ratio of particulate volume to carrier volume [-]
        autoPtr<Function1<scalar>> volumeRatio_;

        //- The ratio of particulate mass to carrier mass [-]
        autoPtr<Function1<scalar>> massRatio_;

        //- Parcels to introduce per unit volume of particulates [m^-3]
        const scalar parcelConcentration_;

        //- Parcel size distribution model
        const autoPtr<distribution> sizeDistribution_;


public:

    //- Runtime type information
    TypeName("patchFlowRateInjection");


    // Constructors

        //- Construct from dictionary
        PatchFlowRateInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        PatchFlowRateInjection(const PatchFlowRateInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new PatchFlowRateInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PatchFlowRateInjection();


    // Member Functions

        //- Inherit topoChange from patchInjectionBase
        using patchInjectionBase::topoChange;

        //- Set injector locations when mesh is updated
        virtual void topoChange();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

        //- Return the total volumetric flow rate across the patch [m^3/s]
        scalar volumetricFlowRate() const;

        //- Return the total mass flow rate across the patch [kg/s]
        scalar massFlowRate() const;

        //- Number of parcels to introduce relative to SOI
        virtual scalar nParcelsToInject(const scalar time0, const scalar time1);

        //- Parcel mass to introduce relative to SOI
        virtual scalar massToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Inherit setPositionAndCell from patchInjectionBase
            using patchInjectionBase::setPositionAndCell;

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const meshSearch& searchEngine,
                const label parcelI,
                const label nParcels,
                const scalar time,
                barycentric& coordinates,
                label& celli,
                label& tetFacei,
                label& tetPti,
                label& facei
            );

            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType::trackingData& td,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchFlowRateInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
