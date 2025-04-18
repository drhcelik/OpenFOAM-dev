/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::coupledFvPatch

Description
    An abstract base class for patches that couple regions of the
    computational domain e.g. cyclic and processor-processor links.

SourceFiles
    coupledFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFvPatch_H
#define coupledFvPatch_H

#include "fvPatch.H"
#include "lduInterface.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coupledFvPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledFvPatch
:
    public lduInterface,
    public fvPatch
{
    // Private Data

        //- Reference to the coupled polyPatch
        const coupledPolyPatch& coupledPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights
        (
            scalarField& w,
            const vectorField& nbrSf,
            const vectorField& nbrDelta
        ) const;


public:

    //- Runtime type information
    TypeName(coupledPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        coupledFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            coupledPolyPatch_(refCast<const coupledPolyPatch>(patch))
        {}


    //- Destructor
    virtual ~coupledFvPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return coupledPolyPatch_.coupled();
            }

            //- Does this side own the patch ?
            virtual bool owner() const
            {
                return coupledPolyPatch_.owner();
            }

            //- Does the coupled side own the patch ?
            virtual bool neighbour() const
            {
                return !owner();
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const = 0;

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return fvPatch::faceCells();
            }

            //- Make patch weighting factors
            virtual void makeWeights(scalarField& w) const = 0;

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const = 0;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const = 0;

            //- Inherit initInternalFieldTransfer from lduInterface
            using lduInterface::initInternalFieldTransfer;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                labelUList& iF
            ) const
            {}

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
