/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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
    Foam::RASModels::phasePressureModel

Description
    Particle-particle phase-pressure RAS model

    The derivative of the phase-pressure with respect to the phase-fraction
    is evaluated as

        g0*min(exp(preAlphaExp*(alpha - alphaMax)), expMax)

    The default model coefficients correspond to the following:
    \verbatim
        phasePressureCoeffs
        {
            preAlphaExp     500;
            expMax          1000;
            g0              1000;
        }
    \endverbatim

SourceFiles
    phasePressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef phasePressureModel_H
#define phasePressureModel_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                     Class phasePressureModel Declaration
\*---------------------------------------------------------------------------*/

class phasePressureModel
:
    public eddyViscosity<RASModel<phaseCompressible::momentumTransportModel>>
{
    // Private Data

        const phaseModel& phase_;

        // Phase pressure coefficients

            //- Pre-exponential factor
            scalar preAlphaExp_;

            //- Maximum limit of the exponential
            scalar expMax_;

            //- g0
            dimensionedScalar g0_;


    // Private Member Functions

        void correctNut()
        {}

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        tmp<volScalarField> pPrime() const;


public:

    //- Runtime type information
    TypeName("phasePressure");


    // Constructors

        //- Construct from components
        phasePressureModel
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        phasePressureModel(const phasePressureModel&) = delete;


    //- Destructor
    virtual ~phasePressureModel();


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            return this->nut();
        }

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const
        {
            return this->nut(patchi);
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const;

        //- Return the stress tensor (0) [m^2/s^2]
        virtual tmp<volSymmTensorField> R() const;

        //- Return the face-phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<surfaceScalarField> pPrimef() const;

        //- Return the effective stress
        virtual tmp<surfaceVectorField> devTau() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const;

        //- Solve the kinetic theory equations and correct the viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phasePressureModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
