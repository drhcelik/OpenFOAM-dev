/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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
    Foam::wallBoilingModels::nucleationSiteModel

Description
    Base class for nucleation site density models

SourceFiles
    nucleationSiteModel.C
    nucleationSiteModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef nucleationSiteModel_H
#define nucleationSiteModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{

/*---------------------------------------------------------------------------*\
                     Class nucleationSiteModel Declaration
\*---------------------------------------------------------------------------*/

class nucleationSiteModel
{
public:

    //- Runtime type information
    TypeName("nucleationSiteModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        nucleationSiteModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    // Constructors

        //- Construct null
        nucleationSiteModel();

        //- Copy construct
        nucleationSiteModel(const nucleationSiteModel& model);

        //- Construct and return a clone
        virtual autoPtr<nucleationSiteModel> clone() const = 0;


    // Selectors

        //- Select null constructed
        static autoPtr<nucleationSiteModel> New(const dictionary& dict);


    //- Destructor
    virtual ~nucleationSiteModel();


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
        ) const = 0;

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
        ) const = 0;

        //- Calculate and return the nucleation-site density
        virtual tmp<volScalarField> nucleationSiteDensity
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const phaseModel& solid,
            const volScalarField& Tf,
            const volScalarField& Tsatw,
            const volScalarField& L,
            const volScalarField& dDep,
            const volScalarField& fDep
        ) const = 0;

        //- Write to stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nucleationSiteModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
