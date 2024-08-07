/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::nonConformalCyclicPolyPatch

Description
    Non-conformal cyclic poly patch. As nonConformalCoupledPolyPatch, but the
    neighbouring patch is local and known and is made available by this class.

See also
    Foam::nonConformalCoupledPolyPatch

SourceFiles
    nonConformalCyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCyclicPolyPatch_H
#define nonConformalCyclicPolyPatch_H

#include "cyclicPolyPatch.H"
#include "nonConformalCoupledPolyPatch.H"
#include "intersectionPatchToPatch.H"
#include "raysPatchToPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class nonConformalCyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCyclicPolyPatch
:
    public cyclicPolyPatch,
    public nonConformalCoupledPolyPatch
{
public:

    // Public types

        //- Enumeration for the condition that triggers re-calculation of the
        //  intersection following mesh motion.
        //
        //  'always' will re-calculate the intersection, 'never' will not, and
        //  'detect' will compare the new and old points and re-calculate if
        //  there is any difference.
        //
        //  The comparison that 'detect' does has an expense associated with it,
        //  so it should only really be used if the patch is both moving and
        //  static at different times. If the patch is always in motion or is
        //  fully stationary, then it is more efficient to specify 'always' or
        //  'never'.
        //
        //  'always' is the default, as it is safest to assume that motion
        //  always affects the intersection. Most non-conformal interfaces are
        //  used for sliding interfaces that are in motion for the entirety of
        //  the simulation, so 'always' is also typically more efficient than
        //  'detect'.
        enum class moveUpdate
        {
            always,
            detect,
            never
        };

        //- Names of the move-update conditions
        static const NamedEnum<moveUpdate, 3> moveUpdateNames_;


protected:

    // Protected data

        //- Is the intersection engine up to date? Zero if out of date. One if
        //  out of date (or potentially not out of date) because of motion. Two
        //  if up to date.
        mutable label intersectionIsValid_;

        //- Patch-to-patch intersection engine
        mutable patchToPatches::intersection intersection_;

        //- Is the rays engine up to date? Zero if out of date. One if out of
        //  date (or potentially not out of date) because of motion. Two if up
        //  to date.
        mutable label raysIsValid_;

        //- Patch-to-patch rays engine
        mutable patchToPatches::rays rays_;

        //- The condition that triggers re-calculation following motion
        const moveUpdate moveUpdate_;


    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers& pBufs, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initTopoChange(PstreamBuffers&);

        //- Clear geometry
        virtual void clearGeom();

        //- Reset the patch name
        virtual void rename(const wordList& newNames);

        //- Reset the patch index
        virtual void reorder(const labelUList& newToOldIndex);


public:

    //- Runtime type information
    TypeName("nonConformalCyclic");


    // Constructors

        //- Construct from components
        nonConformalCyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from components
        nonConformalCyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& nbrPatchName,
            const word& origPatchName,
            const cyclicTransform& transform=cyclicTransform(true)
        );

        //- Construct from dictionary
        nonConformalCyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        nonConformalCyclicPolyPatch
        (
            const nonConformalCyclicPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        nonConformalCyclicPolyPatch
        (
            const nonConformalCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName,
            const word& origPatchName
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new nonConformalCyclicPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new nonConformalCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName(),
                    origPatchName()
                )
            );
        }


    //- Destructor
    virtual ~nonConformalCyclicPolyPatch();


    // Member Functions

        //- Neighbour patch
        const nonConformalCyclicPolyPatch& nbrPatch() const;

        //- Inherit the cyclic owner method
        using cyclicPolyPatch::owner;

        //- Inherit the cyclic neighbour method
        using cyclicPolyPatch::neighbour;

        //- Inherit the cyclic transform method
        using cyclicPolyPatch::transform;

        //- Is this patch coupled? Returns false. For NCC patches the poly
        //  mesh is considered non-coupled whilst the finite volume mesh is
        //  considered coupled.
        virtual bool coupled() const;

        //- Access the intersection engine
        const patchToPatches::intersection& intersection() const;

        //- Access the rays engine
        const patchToPatches::rays& rays() const;

        //- Compute a ray intersection across the coupling
        remote ray
        (
            const scalar fraction,
            const label origFacei,
            const vector& p,
            const vector& n,
            point& nbrP
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
