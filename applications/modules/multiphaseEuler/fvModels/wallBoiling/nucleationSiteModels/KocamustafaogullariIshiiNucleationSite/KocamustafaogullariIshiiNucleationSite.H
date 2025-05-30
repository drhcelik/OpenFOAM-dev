/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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
    Foam::wallBoilingModels::nucleationSiteModels::
    KocamustafaogullariIshiiNucleationSite

Description
    A correlation for nucleation site density.

    Reference:
    \verbatim
        Kocamustafaogullari, G., & Ishii, M. (1983).
        Interfacial area and nucleation site density in boiling systems.
        International Journal of Heat and Mass Transfer, 26(9), 1377-1387.
    \endverbatim

SourceFiles
    KocamustafaogullariIshiiNucleationSite.C

\*---------------------------------------------------------------------------*/

#ifndef KocamustafaogullariIshiiNucleationSite_H
#define KocamustafaogullariIshiiNucleationSite_H

#include "nucleationSiteModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace nucleationSiteModels
{

/*---------------------------------------------------------------------------*\
           Class KocamustafaogullariIshiiNucleationSite Declaration
\*---------------------------------------------------------------------------*/

class KocamustafaogullariIshiiNucleationSite
:
    public nucleationSiteModel
{
    // Private data:

        //- Nucleation density coefficient
        const dimensionedScalar Cn_;


    // Private Member Functions

        //- Calculate and return the nucleation site density
        template<class ScalarFieldType>
        tmp<ScalarFieldType> calculate
        (
            const ScalarFieldType& Tsatw,
            const ScalarFieldType& L,
            const ScalarFieldType& dDep,
            const ScalarFieldType& Tw,
            const ScalarFieldType& rhoLiquid,
            const ScalarFieldType& rhoVapour,
            const ScalarFieldType& sigma
        ) const;


public:

    //- Runtime type information
    TypeName("KocamustafaogullariIshii");

    // Constructors

        //- Construct from a dictionary
        KocamustafaogullariIshiiNucleationSite(const dictionary& dict);

        //- Copy construct
        KocamustafaogullariIshiiNucleationSite
        (
            const KocamustafaogullariIshiiNucleationSite& model
        );

        //- Construct and return a clone
        virtual autoPtr<nucleationSiteModel> clone() const
        {
            return autoPtr<nucleationSiteModel>
            (
                new KocamustafaogullariIshiiNucleationSite(*this)
            );
        }


    //- Destructor
    virtual ~KocamustafaogullariIshiiNucleationSite();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> nucleationSiteDensity
        (
            const phaseModel& liquid,
            const phaseModel& vapour,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L,
            const scalarField& dDep,
            const scalarField& fDep
        ) const;

        //- Calculate and return the nucleation-site density
        virtual tmp<volScalarField::Internal> nucleationSiteDensity
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const phaseModel& solid,
            const volScalarField::Internal& Tf,
            const volScalarField::Internal& Tsatw,
            const volScalarField::Internal& L,
            const volScalarField::Internal& dDep,
            const volScalarField::Internal& fDep
        ) const;

        //- Calculate and return the nucleation-site density
        virtual tmp<volScalarField> nucleationSiteDensity
        (
            const phaseModel& liquid,
            const phaseModel& vapour,
            const phaseModel& solid,
            const volScalarField& Tf,
            const volScalarField& Tsatw,
            const volScalarField& L,
            const volScalarField& dDep,
            const volScalarField& fDep
        ) const;

        //- Write to stream
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationSiteModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
