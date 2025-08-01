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
    Foam::blockFaces::projectFace

Description
    Projects the given set of face points onto the selected surface of the
    geometry provided as a searchableSurfaces object.

SourceFiles
    projectFace.C

\*---------------------------------------------------------------------------*/

#ifndef projectFace_blockFace_H
#define projectFace_blockFace_H

#include "blockFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockFaces
{

/*---------------------------------------------------------------------------*\
                         Class projectFace Declaration
\*---------------------------------------------------------------------------*/

class projectFace
:
    public blockFace
{
    // Private Data

        //- The surface onto which the points are projected
        const searchableSurface& surface_;


    // Private Member Functions

        const searchableSurface& lookupSurface
        (
            const searchableSurfaceList& geometry,
            Istream& is
        ) const;

        //- Convert i,j to single index
        label index
        (
            const labelPair& n,
            const labelPair& coord
        ) const;

        //- Calculate lambdas (but unnormalised)
        void calcLambdas
        (
            const labelPair& n,
            const pointField& points,
            scalarField& lambdaI,
            scalarField& lambdaJ
        ) const;


public:

    //- Runtime type information
    TypeName("project");


    // Constructors

        //- Construct from Istream setting pointsList
        projectFace
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaceList& geometry,
            Istream&
        );

        //- Disallow default bitwise copy construction
        projectFace(const projectFace&) = delete;


    //- Destructor
    virtual ~projectFace()
    {}


    // Member Functions

        //- Project the given points onto the surface
        virtual void project
        (
            const blockDescriptor&,
            const label blockFacei,
            pointField& points
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const projectFace&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockFaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
