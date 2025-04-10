/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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
    Foam::fvMeshStitchers::stationary

Description
    Mesh stitcher for stationary meshes

SourceFiles
    stationary_fvMeshStitcher.C

\*---------------------------------------------------------------------------*/

#ifndef stationary_fvMeshStitcher_H
#define stationary_fvMeshStitcher_H

#include "fvMeshStitcher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshStitchers
{

/*---------------------------------------------------------------------------*\
                Class fvMeshStitchers::stationary Declaration
\*---------------------------------------------------------------------------*/

class stationary
:
    public fvMeshStitcher
{
    // Private Member Functions

        //- Return whether or not this stitcher supports a changing mesh
        virtual bool changing() const
        {
            return false;
        }

        //- Correct the mesh fluxes following a conform operation
        virtual void conformCorrectMeshPhi(surfaceScalarField& phi);

        //- Initialise correction of the mesh fluxes
        virtual void createNonConformalCorrectMeshPhiGeometry
        (
            SurfaceFieldBoundary<label>& polyFacesBf,
            surfaceVectorField& SfSf,
            surfaceVectorField& CfSf
        );

        //- Correct the mesh fluxes following an unconform operation
        virtual void unconformCorrectMeshPhi
        (
            const SurfaceFieldBoundary<label>& polyFacesBf,
            surfaceVectorField& SfSf,
            surfaceVectorField& CfSf,
            surfaceScalarField& phi
        );


public:

    //- Runtime type information
    TypeName("stationary");


    // Constructors

        //- Construct from fvMesh
        explicit stationary(fvMesh&);


    //- Destructor
    virtual ~stationary();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshStitchers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
