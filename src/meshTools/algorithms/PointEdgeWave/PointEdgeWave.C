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

#include "PointEdgeWave.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamCombineReduceOps.H"
#include "debug.H"
#include "typeInfo.H"
#include "globalMeshData.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type, class TrackingData>
Foam::scalar Foam::PointEdgeWave<Type, TrackingData>::propagationTol_ = 0.01;

template<class Type, class TrackingData>
int Foam::PointEdgeWave<Type, TrackingData>::defaultTrackingData_ = -1;


namespace Foam
{
    //- Reduction class. If x and y are not equal assign value.
    template<class Type, class TrackingData>
    class combineEqOp
    {
        TrackingData& td_;

        public:
            combineEqOp(TrackingData& td)
            :
                td_(td)
            {}

        void operator()(Type& x, const Type& y) const
        {
            if (!x.valid(td_) && y.valid(td_))
            {
                x = y;
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class TrackingData>
void Foam::PointEdgeWave<Type, TrackingData>::transform
(
    const polyPatch& patch,
    const labelList& patchPointLabels,
    const transformer& transform,
    List<Type>& pointInfo
) const
{
    forAll(pointInfo, i)
    {
        pointInfo[i].transform(patch, patchPointLabels[i], transform, td_);
    }
}


template<class Type, class TrackingData>
bool Foam::PointEdgeWave<Type, TrackingData>::updatePoint
(
    const label pointi,
    const label neighbourEdgeI,
    const Type& neighbourInfo,
    Type& pointInfo
)
{
    // Update info for pointi, at position pt, with information from
    // neighbouring edge.
    // Updates:
    //      - changedPoint_, changedPoints_, nChangedPoints_,
    //      - statistics: nEvals_, nUnvisitedPoints_

    nEvals_++;

    bool wasValid = pointInfo.valid(td_);

    bool propagate =
        pointInfo.updatePoint
        (
            mesh_,
            pointi,
            neighbourEdgeI,
            neighbourInfo,
            propagationTol_,
            td_
        );

    if (propagate)
    {
        if (!changedPoint_[pointi])
        {
            changedPoint_[pointi] = true;
            changedPoints_[nChangedPoints_++] = pointi;
        }
    }

    if (!wasValid && pointInfo.valid(td_))
    {
        --nUnvisitedPoints_;
    }

    return propagate;
}


template<class Type, class TrackingData>
bool Foam::PointEdgeWave<Type, TrackingData>::updatePoint
(
    const label pointi,
    const Type& neighbourInfo,
    Type& pointInfo
)
{
    // Update info for pointi, at position pt, with information from
    // same point.
    // Updates:
    //      - changedPoint_, changedPoints_, nChangedPoints_,
    //      - statistics: nEvals_, nUnvisitedPoints_

    nEvals_++;

    bool wasValid = pointInfo.valid(td_);

    bool propagate =
        pointInfo.updatePoint
        (
            mesh_,
            pointi,
            neighbourInfo,
            propagationTol_,
            td_
        );

    if (propagate)
    {
        if (!changedPoint_[pointi])
        {
            changedPoint_[pointi] = true;
            changedPoints_[nChangedPoints_++] = pointi;
        }
    }

    if (!wasValid && pointInfo.valid(td_))
    {
        --nUnvisitedPoints_;
    }

    return propagate;
}


template<class Type, class TrackingData>
bool Foam::PointEdgeWave<Type, TrackingData>::updateEdge
(
    const label edgeI,
    const label neighbourPointi,
    const Type& neighbourInfo,
    Type& edgeInfo
)
{
    // Update info for edgeI, at position pt, with information from
    // neighbouring point.
    // Updates:
    //      - changedEdge_, changedEdges_, nChangedEdges_,
    //      - statistics: nEvals_, nUnvisitedEdge_

    nEvals_++;

    bool wasValid = edgeInfo.valid(td_);

    bool propagate =
        edgeInfo.updateEdge
        (
            mesh_,
            edgeI,
            neighbourPointi,
            neighbourInfo,
            propagationTol_,
            td_
        );

    if (propagate)
    {
        if (!changedEdge_[edgeI])
        {
            changedEdge_[edgeI] = true;
            changedEdges_[nChangedEdges_++] = edgeI;
        }
    }

    if (!wasValid && edgeInfo.valid(td_))
    {
        --nUnvisitedEdges_;
    }

    return propagate;
}


template<class Type, class TrackingData>
template<class PatchType>
Foam::label Foam::PointEdgeWave<Type, TrackingData>::countPatchType() const
{
    label nPatches = 0;

    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<PatchType>(mesh_.boundaryMesh()[patchi]))
        {
            nPatches++;
        }
    }
    return nPatches;
}


template<class Type, class TrackingData>
void Foam::PointEdgeWave<Type, TrackingData>::handleProcPatches()
{
    // 1. Send all point info on processor patches.

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    DynamicList<Type> patchInfo;
    DynamicList<label> thisPoints;
    DynamicList<label> nbrPoints;

    forAll(mesh_.globalData().processorPatches(), i)
    {
        label patchi = mesh_.globalData().processorPatches()[i];
        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchi]);

        patchInfo.clear();
        patchInfo.reserve(procPatch.nPoints());
        thisPoints.clear();
        thisPoints.reserve(procPatch.nPoints());
        nbrPoints.clear();
        nbrPoints.reserve(procPatch.nPoints());

        // Get all changed points in reverse order
        const labelList& neighbPoints = procPatch.nbrPoints();
        forAll(neighbPoints, thisPointi)
        {
            label meshPointi = procPatch.meshPoints()[thisPointi];
            if (changedPoint_[meshPointi])
            {
                patchInfo.append(allPointInfo_[meshPointi]);
                thisPoints.append(thisPointi);
                nbrPoints.append(neighbPoints[thisPointi]);
            }
        }

        // if (debug)
        //{
        //    Pout<< "Processor patch " << patchi << ' ' << procPatch.name()
        //        << " communicating with " << procPatch.neighbProcNo()
        //        << "  Sending:" << patchInfo.size() << endl;
        //}

        UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);
        toNeighbour << nbrPoints << patchInfo;
    }


    pBufs.finishedSends();

    //
    // 2. Receive all point info on processor patches.
    //

    forAll(mesh_.globalData().processorPatches(), i)
    {
        label patchi = mesh_.globalData().processorPatches()[i];
        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchi]);

        List<Type> patchInfo;
        labelList patchPoints;

        {
            UIPstream fromNeighbour(procPatch.neighbProcNo(), pBufs);
            fromNeighbour >> patchPoints >> patchInfo;
        }

        // if (debug)
        //{
        //    Pout<< "Processor patch " << patchi << ' ' << procPatch.name()
        //        << " communicating with " << procPatch.neighbProcNo()
        //        << "  Received:" << patchInfo.size() << endl;
        //}

        // Transform info across the interface
        transform(procPatch, patchPoints, procPatch.transform(), patchInfo);

        // Merge received info
        const labelList& meshPoints = procPatch.meshPoints();
        forAll(patchInfo, i)
        {
            label meshPointi = meshPoints[patchPoints[i]];

            if (!allPointInfo_[meshPointi].equal(patchInfo[i], td_))
            {
                updatePoint
                (
                    meshPointi,
                    patchInfo[i],
                    allPointInfo_[meshPointi]
                );
            }
        }
    }

    // Collocated points should be handled by face based transfer
    // (since that is how connectivity is worked out)
    // They are also explicitly equalised in handleCollocatedPoints to
    // guarantee identical values.
}


template<class Type, class TrackingData>
void Foam::PointEdgeWave<Type, TrackingData>::handleCyclicPatches()
{
    // 1. Send all point info on cyclic patches.

    DynamicList<Type> nbrInfo;
    DynamicList<label> nbrPoints;
    DynamicList<label> thisPoints;

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];

        if (isA<cyclicPolyPatch>(patch))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patch);

            nbrInfo.clear();
            nbrInfo.reserve(cycPatch.nPoints());
            nbrPoints.clear();
            nbrPoints.reserve(cycPatch.nPoints());
            thisPoints.clear();
            thisPoints.reserve(cycPatch.nPoints());

            // Collect nbrPatch points that have changed
            {
                const cyclicPolyPatch& nbrPatch = cycPatch.nbrPatch();
                const edgeList& pairs = cycPatch.coupledPoints();
                const labelList& meshPoints = nbrPatch.meshPoints();

                forAll(pairs, pairI)
                {
                    label thisPointi = pairs[pairI][0];
                    label nbrPointi = pairs[pairI][1];
                    label meshPointi = meshPoints[nbrPointi];

                    if (changedPoint_[meshPointi])
                    {
                        nbrInfo.append(allPointInfo_[meshPointi]);
                        nbrPoints.append(nbrPointi);
                        thisPoints.append(thisPointi);
                    }
                }
            }

            // if (debug)
            //{
            //    Pout<< "Cyclic patch " << patchi << ' ' << patch.name()
            //        << "  Changed : " << nbrInfo.size()
            //        << endl;
            //}

            // Transform info across the interface
            transform(cycPatch, thisPoints, cycPatch.transform(), nbrInfo);

            // Merge received info
            const labelList& meshPoints = cycPatch.meshPoints();
            forAll(nbrInfo, i)
            {
                label meshPointi = meshPoints[thisPoints[i]];

                if (!allPointInfo_[meshPointi].equal(nbrInfo[i], td_))
                {
                    updatePoint
                    (
                        meshPointi,
                        nbrInfo[i],
                        allPointInfo_[meshPointi]
                    );
                }
            }
        }
    }
}


template<class Type, class TrackingData>
Foam::label Foam::PointEdgeWave<Type, TrackingData>::handleCollocatedPoints()
{
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh_.globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const distributionMap& slavesMap = gmd.globalPointSlavesMap();
    const labelListList& slaves = gmd.globalPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, pointi)
    {
        elems[pointi] = allPointInfo_[meshPoints[pointi]];
    }

    // Pull slave data onto master (which might or might not have any
    // initialised points). No need to update transformed slots.
    slavesMap.distribute(elems, false);

    // Combine master data with slave data
    combineEqOp<Type, TrackingData> cop(td_);

    forAll(slaves, pointi)
    {
        Type& elem = elems[pointi];

        const labelList& slavePoints = slaves[pointi];

        // Combine master with untransformed slave data
        forAll(slavePoints, j)
        {
            cop(elem, elems[slavePoints[j]]);
        }

        // Copy result back to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elem;
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, pointi)
    {
        if (elems[pointi].valid(td_))
        {
            label meshPointi = meshPoints[pointi];

            Type& elem = allPointInfo_[meshPointi];

            bool wasValid = elem.valid(td_);

            // Like updatePoint but bypass Type::updatePoint with its tolerance
            // checking
            // if (!elem.valid(td_) || !elem.equal(elems[pointi], td_))
            if (!elem.equal(elems[pointi], td_))
            {
                nEvals_++;
                elem = elems[pointi];

                // See if element now valid
                if (!wasValid && elem.valid(td_))
                {
                    --nUnvisitedPoints_;
                }

                // Update database of changed points
                if (!changedPoint_[meshPointi])
                {
                    changedPoint_[meshPointi] = true;
                    changedPoints_[nChangedPoints_++] = meshPointi;
                }
            }
        }
    }

    // Sum nChangedPoints over all procs
    label totNChanged = nChangedPoints_;

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class TrackingData>
Foam::PointEdgeWave<Type, TrackingData>::PointEdgeWave
(
    const polyMesh& mesh,
    const labelList& startPoints,
    const List<Type>& startPointsInfo,

    UList<Type>& allPointInfo,
    UList<Type>& allEdgeInfo,

    const label maxIter,
    TrackingData& td
)
:
    mesh_(mesh),
    allPointInfo_(allPointInfo),
    allEdgeInfo_(allEdgeInfo),
    td_(td),
    changedPoint_(mesh_.nPoints(), false),
    changedPoints_(mesh_.nPoints()),
    nChangedPoints_(0),
    changedEdge_(mesh_.nEdges(), false),
    changedEdges_(mesh_.nEdges()),
    nChangedEdges_(0),
    nCyclicPatches_(countPatchType<cyclicPolyPatch>()),
    nEvals_(0),
    nUnvisitedPoints_(mesh_.nPoints()),
    nUnvisitedEdges_(mesh_.nEdges())
{
    if (allPointInfo_.size() != mesh_.nPoints())
    {
        FatalErrorInFunction
            << "size of pointInfo work array is not equal to the number"
            << " of points in the mesh" << endl
            << "    pointInfo   :" << allPointInfo_.size() << endl
            << "    mesh.nPoints:" << mesh_.nPoints()
            << exit(FatalError);
    }
    if (allEdgeInfo_.size() != mesh_.nEdges())
    {
        FatalErrorInFunction
            << "size of edgeInfo work array is not equal to the number"
            << " of edges in the mesh" << endl
            << "    edgeInfo   :" << allEdgeInfo_.size() << endl
            << "    mesh.nEdges:" << mesh_.nEdges()
            << exit(FatalError);
    }


    // Set from initial changed points data
    setPointInfo(startPoints, startPointsInfo);

    if (debug)
    {
        Info<< typeName << ": Seed points               : "
            << returnReduce(nChangedPoints_, sumOp<label>()) << endl;
    }

    // Iterate until nothing changes
    label iter = iterate(maxIter);

    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter." << endl
            << "    maxIter:" << maxIter << endl
            << "    nChangedPoints:" << nChangedPoints_ << endl
            << "    nChangedEdges:" << nChangedEdges_ << endl
            << exit(FatalError);
    }
}


template<class Type, class TrackingData>
Foam::PointEdgeWave<Type, TrackingData>::PointEdgeWave
(
    const polyMesh& mesh,
    UList<Type>& allPointInfo,
    UList<Type>& allEdgeInfo,
    TrackingData& td
)
:
    mesh_(mesh),
    allPointInfo_(allPointInfo),
    allEdgeInfo_(allEdgeInfo),
    td_(td),
    changedPoint_(mesh_.nPoints(), false),
    changedPoints_(mesh_.nPoints()),
    nChangedPoints_(0),
    changedEdge_(mesh_.nEdges(), false),
    changedEdges_(mesh_.nEdges()),
    nChangedEdges_(0),
    nCyclicPatches_(countPatchType<cyclicPolyPatch>()),
    nEvals_(0),
    nUnvisitedPoints_(mesh_.nPoints()),
    nUnvisitedEdges_(mesh_.nEdges())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class TrackingData>
Foam::PointEdgeWave<Type, TrackingData>::~PointEdgeWave()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type, class TrackingData>
Foam::label Foam::PointEdgeWave<Type, TrackingData>::getUnsetPoints() const
{
    return nUnvisitedPoints_;
}


template<class Type, class TrackingData>
Foam::label Foam::PointEdgeWave<Type, TrackingData>::getUnsetEdges() const
{
    return nUnvisitedEdges_;
}


// Copy point information into member data
template<class Type, class TrackingData>
void Foam::PointEdgeWave<Type, TrackingData>::setPointInfo
(
    const labelList& startPoints,
    const List<Type>& startPointsInfo
)
{
    forAll(startPoints, startPointi)
    {
        label pointi = startPoints[startPointi];

        bool wasValid = allPointInfo_[pointi].valid(td_);

        // Copy info for pointi
        allPointInfo_[pointi] = startPointsInfo[startPointi];

        // Maintain count of unset points
        if (!wasValid && allPointInfo_[pointi].valid(td_))
        {
            --nUnvisitedPoints_;
        }

        // Mark pointi as changed, both on list and on point itself.

        if (!changedPoint_[pointi])
        {
            changedPoint_[pointi] = true;
            changedPoints_[nChangedPoints_++] = pointi;
        }
    }

    // Sync
    handleCollocatedPoints();
}


template<class Type, class TrackingData>
Foam::label Foam::PointEdgeWave<Type, TrackingData>::edgeToPoint()
{
    for
    (
        label changedEdgeI = 0;
        changedEdgeI < nChangedEdges_;
        changedEdgeI++
    )
    {
        label edgeI = changedEdges_[changedEdgeI];

        if (!changedEdge_[edgeI])
        {
            FatalErrorInFunction
                << "edge " << edgeI
                << " not marked as having been changed" << nl
                << "This might be caused by multiple occurrences of the same"
                << " seed point." << abort(FatalError);
        }


        const Type& neighbourWallInfo = allEdgeInfo_[edgeI];

        // Evaluate all connected points (= edge endpoints)
        const edge& e = mesh_.edges()[edgeI];

        forAll(e, eI)
        {
            Type& currentWallInfo = allPointInfo_[e[eI]];

            if (!currentWallInfo.equal(neighbourWallInfo, td_))
            {
                updatePoint
                (
                    e[eI],
                    edgeI,
                    neighbourWallInfo,
                    currentWallInfo
                );
            }
        }

        // Reset status of edge
        changedEdge_[edgeI] = false;
    }

    // Handled all changed edges by now
    nChangedEdges_ = 0;

    if (nCyclicPatches_ > 0)
    {
        // Transfer changed points across cyclic halves
        handleCyclicPatches();
    }
    if (Pstream::parRun())
    {
        // Transfer changed points from neighbouring processors.
        handleProcPatches();
    }

    // Sum nChangedPoints over all procs
    label totNChanged = nChangedPoints_;

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


template<class Type, class TrackingData>
Foam::label Foam::PointEdgeWave<Type, TrackingData>::pointToEdge()
{
    const labelListList& pointEdges = mesh_.pointEdges();

    for
    (
        label changedPointi = 0;
        changedPointi < nChangedPoints_;
        changedPointi++
    )
    {
        label pointi = changedPoints_[changedPointi];

        if (!changedPoint_[pointi])
        {
            FatalErrorInFunction
                << "Point " << pointi
                << " not marked as having been changed" << nl
                << "This might be caused by multiple occurrences of the same"
                << " seed point." << abort(FatalError);
        }

        const Type& neighbourWallInfo = allPointInfo_[pointi];

        // Evaluate all connected edges

        const labelList& edgeLabels = pointEdges[pointi];
        forAll(edgeLabels, edgeLabelI)
        {
            label edgeI = edgeLabels[edgeLabelI];

            Type& currentWallInfo = allEdgeInfo_[edgeI];

            if (!currentWallInfo.equal(neighbourWallInfo, td_))
            {
                updateEdge
                (
                    edgeI,
                    pointi,
                    neighbourWallInfo,
                    currentWallInfo
                );
            }
        }

        // Reset status of point
        changedPoint_[pointi] = false;
    }

    // Handled all changed points by now
    nChangedPoints_ = 0;

    // if (debug)
    //{
    //    Pout<< "Changed edges             : " << nChangedEdges_ << endl;
    //}

    // Sum nChangedPoints over all procs
    label totNChanged = nChangedEdges_;

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


template<class Type, class TrackingData>
Foam::label Foam::PointEdgeWave<Type, TrackingData>::iterate
(
    const label maxIter
)
{
    if (nCyclicPatches_ > 0)
    {
        // Transfer changed points across cyclic halves
        handleCyclicPatches();
    }
    if (Pstream::parRun())
    {
        // Transfer changed points from neighbouring processors.
        handleProcPatches();
    }

    nEvals_ = 0;

    label iter = 0;

    while (iter < maxIter)
    {
        while (iter < maxIter)
        {
            if (debug)
            {
                Info<< typeName << ": Iteration " << iter << endl;
            }

            label nEdges = pointToEdge();

            if (debug)
            {
                Info<< typeName << ": Total changed edges       : "
                    << nEdges << endl;
            }

            if (nEdges == 0)
            {
                break;
            }

            label nPoints = edgeToPoint();

            if (debug)
            {
                Info<< typeName << ": Total changed points      : "
                    << nPoints << nl
                    << typeName << ": Total evaluations         : "
                    << returnReduce(nEvals_, sumOp<label>()) << nl
                    << typeName << ": Remaining unvisited points: "
                    << returnReduce(nUnvisitedPoints_, sumOp<label>()) << nl
                    << typeName << ": Remaining unvisited edges : "
                    << returnReduce(nUnvisitedEdges_, sumOp<label>()) << nl
                    << endl;
            }

            if (nPoints == 0)
            {
                break;
            }

            iter++;
        }


        // Enforce collocated points are exactly equal. This might still mean
        // non-collocated points are not equal though. WIP.
        label nPoints = handleCollocatedPoints();
        if (debug)
        {
            Info<< typeName << ": Collocated point sync     : "
                << nPoints << nl << endl;
        }

        if (nPoints == 0)
        {
            break;
        }
    }

    return iter;
}


// ************************************************************************* //
