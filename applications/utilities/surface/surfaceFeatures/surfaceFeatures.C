/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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

Description
    Identifies features in a surface geometry and writes them to file,
    based on control parameters specified by the user.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurface_searchableSurface.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "surfaceFeatures.H"
#include "triSurfaceFields.H"
#include "vtkWritePolyData.H"
#include "systemDict.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    autoPtr<surfaceFeatures> extractFromFile
    (
        const fileName& featureEdgeFile,
        const triSurface& surf,
        const Switch& geometricTestOnly
    )
    {
        edgeMesh eMesh(featureEdgeFile);

        // Sometimes duplicate edges are present. Remove them.
        eMesh.mergeEdges();

        Info<< nl << "Reading existing feature edges from file "
            << featureEdgeFile << nl
            << "Selecting edges purely based on geometric tests: "
            << geometricTestOnly.asText() << endl;

        return  autoPtr<surfaceFeatures>
        (
            new surfaceFeatures
            (
                surf,
                eMesh.points(),
                eMesh.edges(),
                1e-6,
                geometricTestOnly
            )
        );
    }


    autoPtr<surfaceFeatures> extractFromSurface
    (
        const triSurface& surf,
        const Switch& geometricTestOnly,
        const scalar includedAngle
    )
    {
        Info<< nl
            << "Constructing feature set from included angle "
            << includedAngle << nl
            << "Selecting edges purely based on geometric tests: "
            << geometricTestOnly.asText() << endl;

        return  autoPtr<surfaceFeatures>
        (
            new surfaceFeatures
            (
                surf,
                includedAngle,
                0,
                0,
                geometricTestOnly
            )
        );
    }


    autoPtr<surfaceFeatures> surfaceFeatureSet
    (
        const fileName& surfaceFileName,
        const triSurface& surf,
        const dictionary& dict,
        const scalar includedAngle
    )
    {
        const Switch geometricTestOnly = dict.lookupOrDefault<Switch>
        (
            "geometricTestOnly",
            "no"
        );

        if (dict.found("files"))
        {
            HashTable<fileName, fileName> fileNames(dict.lookup("files"));

            if (fileNames.found(surfaceFileName))
            {
                return extractFromFile
                (
                    fileNames[surfaceFileName],
                    surf,
                    geometricTestOnly
                );
            }
            else
            {
                return extractFromSurface
                (
                    surf,
                    geometricTestOnly,
                    includedAngle
                );
            }
        }
        else
        {
            return extractFromSurface
            (
                surf,
                geometricTestOnly,
                includedAngle
            );
        }
    }


    void extractFeatures
    (
        const fileName& surfaceFileName,
        const Time& runTime,
        const dictionary& dict
    )
    {
        const fileName sFeatFileName = surfaceFileName.lessExt().name();

        Info<< "Surface            : " << surfaceFileName << nl << endl;

        const Switch writeVTK =
            dict.lookupOrDefault<Switch>("writeVTK", "off");
        const Switch writeObj =
            dict.lookupOrDefault<Switch>("writeObj", "off");
        const Switch verboseObj =
            dict.lookupOrDefault<Switch>("verboseObj", "off");

        const Switch curvature =
            dict.lookupOrDefault<Switch>("curvature", "off");
        const Switch featureProximity =
            dict.lookupOrDefault<Switch>("featureProximity", "off");


        Info<< nl << "Feature line extraction is only valid on closed manifold "
            << "surfaces." << endl;

        // Read
        // ~~~~

        triSurface surf
        (
            runTime.path()
           /runTime.constant()
           /searchableSurface::geometryDir(runTime)
           /surfaceFileName
        );

        Info<< "Statistics:" << endl;
        surf.writeStats(Info);
        Info<< endl;

        const faceList faces(surf.faces());

        // Either construct features from surface & featureAngle or read set.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const scalar includedAngle =
            dict.lookup<scalar>("includedAngle", unitDegrees);

        autoPtr<surfaceFeatures> set
        (
            surfaceFeatureSet
            (
                surfaceFileName,
                surf,
                dict,
                includedAngle
            )
        );

        // Trim set
        // ~~~~~~~~

        if (dict.isDict("trimFeatures"))
        {
            dictionary trimDict = dict.subDict("trimFeatures");

            scalar minLen =
                trimDict.lookupOrAddDefault<scalar>("minLen", -great);

            label minElem = trimDict.lookupOrAddDefault<label>("minElem", 0);

            // Trim away small groups of features
            if (minElem > 0 || minLen > 0)
            {
                Info<< "Removing features of length < "
                    << minLen << endl;
                Info<< "Removing features with number of edges < "
                    << minElem << endl;

                set().trimFeatures(minLen, minElem, includedAngle);
            }
        }


        // Subset
        // ~~~~~~

        // Convert to marked edges, points
        List<surfaceFeatures::edgeStatus> edgeStat(set().toStatus());

        if (dict.isDict("subsetFeatures"))
        {
            const dictionary& subsetDict = dict.subDict
            (
                "subsetFeatures"
            );

            if (subsetDict.found("insideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("insideBox")());

                Info<< "Selecting edges inside bb " << bb;
                if (writeObj)
                {
                    Info << " see insideBox.obj";
                    bb.writeOBJ("insideBox.obj");
                }
                Info<< endl;

                selectBox(surf, bb, true, edgeStat);
            }
            else if (subsetDict.found("outsideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("outsideBox")());

                Info<< "Removing all edges inside bb " << bb;
                if (writeObj)
                {
                    Info<< " see outsideBox.obj" << endl;
                    bb.writeOBJ("outsideBox.obj");
                }
                Info<< endl;

                selectBox(surf, bb, false, edgeStat);
            }

            const Switch nonManifoldEdges =
                subsetDict.lookupOrDefault<Switch>("nonManifoldEdges", "yes");

            if (!nonManifoldEdges)
            {
                Info<< "Removing all non-manifold edges"
                    << " (edges with > 2 connected faces) unless they"
                    << " cross multiple regions" << endl;

                selectManifoldEdges(surf, 1e-5, includedAngle, edgeStat);
            }

            const Switch openEdges =
                subsetDict.lookupOrDefault<Switch>("openEdges", "yes");

            if (!openEdges)
            {
                Info<< "Removing all open edges"
                    << " (edges with 1 connected face)" << endl;

                forAll(edgeStat, edgei)
                {
                    if (surf.edgeFaces()[edgei].size() == 1)
                    {
                        edgeStat[edgei] = surfaceFeatures::NONE;
                    }
                }
            }

            if (subsetDict.found("plane"))
            {
                const plane cutPlane(subsetDict.subDict("plane"));

                selectCutEdges(surf, cutPlane, edgeStat);

                Info<< "Only edges that intersect the plane with normal "
                    << cutPlane.normal()
                    << " and base point " << cutPlane.refPoint()
                    << " will be included as feature edges."<< endl;
            }
        }


        surfaceFeatures newSet(surf);
        newSet.setFromStatus(edgeStat, includedAngle);

        Info<< nl
            << "Initial feature set:" << nl
            << "    feature points : " << newSet.featurePoints().size() << nl
            << "    feature edges  : " << newSet.featureEdges().size() << nl
            << "    of which" << nl
            << "        region edges   : " << newSet.nRegionEdges() << nl
            << "        external edges : " << newSet.nExternalEdges() << nl
            << "        internal edges : " << newSet.nInternalEdges() << nl
            << endl;

        boolList surfBaffleRegions(surf.patches().size(), false);

        wordList surfBaffleNames;
        dict.readIfPresent("baffles", surfBaffleNames);

        forAll(surf.patches(), pI)
        {
            const word& name = surf.patches()[pI].name();

            if (findIndex(surfBaffleNames, name) != -1)
            {
                Info<< "Adding baffle region " << name << endl;
                surfBaffleRegions[pI] = true;
            }
        }

        // Extracting and writing a extendedFeatureEdgeMesh
        extendedFeatureEdgeMesh feMesh
        (
            newSet,
            runTime,
            sFeatFileName + ".extendedFeatureEdgeMesh",
            surfBaffleRegions
        );


        if (dict.isDict("addFeatures"))
        {
            const word addFeName = dict.subDict("addFeatures")["name"];
            Info<< "Adding (without merging) features from " << addFeName
                << nl << endl;

            extendedFeatureEdgeMesh addFeMesh
            (
                IOobject
                (
                    addFeName,
                    runTime.time().constant(),
                    "extendedFeatureEdgeMesh",
                    runTime.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            Info<< "Read " << addFeMesh.name() << nl;
            addFeMesh.writeStats(Info);

            feMesh.add(addFeMesh);
        }


        Info<< nl
            << "Final feature set:" << nl;
        feMesh.writeStats(Info);

        Info<< nl << "Writing extendedFeatureEdgeMesh to "
            << feMesh.relativeObjectPath() << endl;

        mkDir(feMesh.path());

        if (writeObj)
        {
            feMesh.writeObj
            (
                feMesh.path()/surfaceFileName.lessExt().name(),
                verboseObj
            );
        }

        feMesh.write();

        // Write a featureEdgeMesh for backwards compatibility
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                surfaceFileName.lessExt().name() + ".eMesh",
                runTime.constant(),
                searchableSurface::geometryDir(runTime),
                runTime,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            feMesh.points(),
            feMesh.edges()
        );

        Info<< nl << "Writing featureEdgeMesh to "
            << bfeMesh.relativeObjectPath() << endl;

        bfeMesh.regIOobject::write();

        // Data to write out in VTK format
        wordList writeVTKFieldNames;
        boolList writeVTKFieldIsPointValues;
        #define DeclareWriteVTKFieldTypeValues(Type, nullArg) \
            PtrList<const Field<Type>> writeVTKField##Type##Values;
        FOR_ALL_FIELD_TYPES(DeclareWriteVTKFieldTypeValues);
        #undef DeclareWriteVTKFieldTypeValues

        // Find distance between close features
        if (dict.isDict("closeness"))
        {
            Info<< nl << "Extracting internal and external closeness of "
                << "surface." << endl;

            const dictionary& closenessDict = dict.subDict("closeness");

            const Switch faceCloseness =
                closenessDict.lookupOrDefault<Switch>("faceCloseness", "off");
            const Switch pointCloseness =
                closenessDict.lookupOrDefault<Switch>("pointCloseness", "off");

            const scalar internalAngleTolerance
            (
                closenessDict.lookupOrDefault<scalar>
                (
                    "internalAngleTolerance",
                    unitDegrees,
                    80
                )
            );

            const scalar externalAngleTolerance
            (
                closenessDict.lookupOrDefault<scalar>
                (
                    "externalAngleTolerance",
                    unitDegrees,
                    80
                )
            );

            // Searchable triSurface
            const searchableSurfaces::triSurface searchSurf
            (
                IOobject
                (
                    sFeatFileName + ".closeness",
                    runTime.constant(),
                    searchableSurface::geometryDir(runTime),
                    runTime
                ),
                surf
            );

            if (faceCloseness)
            {
                Pair<tmp<triSurfaceScalarField>> closenessFields
                (
                    searchSurf.extractCloseness
                    (
                        internalAngleTolerance,
                        externalAngleTolerance
                    )
                );

                Info<< "    writing "
                    << closenessFields.first()->name() << endl;
                closenessFields.first()->write();

                Info<< "    writing "
                    << closenessFields.second()->name() << endl;
                closenessFields.second()->write();

                if (writeVTK)
                {
                    writeVTKFieldNames.append("internalCloseness");
                    writeVTKFieldIsPointValues.append(false);
                    writeVTKFieldscalarValues.append
                    (
                        new scalarField(closenessFields.first())
                    );

                    writeVTKFieldNames.append("externalCloseness");
                    writeVTKFieldIsPointValues.append(false);
                    writeVTKFieldscalarValues.append
                    (
                        new scalarField(closenessFields.second())
                    );
                }
            }

            if (pointCloseness)
            {
                Pair<tmp<triSurfacePointScalarField >> closenessFields
                (
                    searchSurf.extractPointCloseness
                    (
                        internalAngleTolerance,
                        externalAngleTolerance
                    )
                );

                Info<< "    writing "
                    << closenessFields.first()->name() << endl;
                closenessFields.first()->write();

                Info<< "    writing "
                    << closenessFields.second()->name() << endl;
                closenessFields.second()->write();

                if (writeVTK)
                {
                    const faceList faces(searchSurf.faces());
                    const Map<label>& meshPointMap = searchSurf.meshPointMap();

                    const triSurfacePointScalarField&
                        internalClosenessPointField = closenessFields.first();

                    const triSurfacePointScalarField&
                        externalClosenessPointField = closenessFields.second();

                    scalarField internalCloseness(searchSurf.nPoints(), great);
                    scalarField externalCloseness(searchSurf.nPoints(), great);

                    forAll(meshPointMap, pi)
                    {
                        internalCloseness[pi] =
                            internalClosenessPointField[meshPointMap[pi]];

                        externalCloseness[pi] =
                            externalClosenessPointField[meshPointMap[pi]];
                    }

                    writeVTKFieldNames.append("internalPointCloseness");
                    writeVTKFieldIsPointValues.append(true);
                    writeVTKFieldscalarValues.append
                    (
                        new scalarField(internalCloseness)
                    );

                    writeVTKFieldNames.append("externalPointCloseness");
                    writeVTKFieldIsPointValues.append(true);
                    writeVTKFieldscalarValues.append
                    (
                        new scalarField(externalCloseness)
                    );
                }
            }
        }


        if (curvature)
        {
            Info<< nl << "Extracting curvature of surface at the points."
                << endl;

            triSurfacePointScalarField k
            (
                IOobject
                (
                    sFeatFileName + ".curvature",
                    runTime.constant(),
                    searchableSurface::geometryDir(runTime),
                    runTime
                ),
                surf,
                dimLength,
                surf.curvature()
            );

            k.write();

            if (writeVTK)
            {
                writeVTKFieldNames.append("curvature");
                writeVTKFieldIsPointValues.append(true);
                writeVTKFieldscalarValues.append(new scalarField(k));
            }
        }


        if (featureProximity)
        {
            Info<< nl << "Extracting proximity of close feature points and "
                << "edges to the surface" << endl;

            const scalar searchDistance =
                dict.lookup<scalar>("maxFeatureProximity");

            scalarField featureProximity(surf.size(), searchDistance);

            forAll(surf, fi)
            {
                const triPointRef& tri = surf[fi].tri(surf.points());

                const Tuple2<point, scalar> circle = tri.circumCircle();
                const point& c = circle.first();
                const scalar rSqr =
                    min(sqr(4*circle.second()), sqr(searchDistance));

                pointIndexHitList hitList;

                feMesh.allNearestFeatureEdges(c, rSqr, hitList);
                featureProximity[fi] = min
                (
                    feMesh.minDisconnectedDist(hitList),
                    featureProximity[fi]
                );

                feMesh.allNearestFeaturePoints(c, rSqr, hitList);
                featureProximity[fi] = min
                (
                    minDist(hitList),
                    featureProximity[fi]
                );
            }

            triSurfaceScalarField featureProximityField
            (
                IOobject
                (
                    sFeatFileName + ".featureProximity",
                    runTime.constant(),
                    searchableSurface::geometryDir(runTime),
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimLength,
                featureProximity
            );

            featureProximityField.write();

            if (writeVTK)
            {
                writeVTKFieldNames.append("featureProximity");
                writeVTKFieldIsPointValues.append(false);
                writeVTKFieldscalarValues.append
                (
                    new scalarField(featureProximity)
                );
            }
        }

        if (writeVTK)
        {
            #define WriteVTKResizeFieldTypeValues(Type, nullArg) \
                writeVTKField##Type##Values.resize(writeVTKFieldNames.size());
            FOR_ALL_FIELD_TYPES(WriteVTKResizeFieldTypeValues)
            #undef WriteVTKResizeFieldTypeValues

            vtkWritePolyData::write
            (
                runTime.path()
               /runTime.constant()
               /searchableSurface::geometryDir(runTime)
               /sFeatFileName + "Features.vtk",
                sFeatFileName,
                runTime.writeFormat() == IOstream::BINARY,
                surf.points(),
                labelList(),
                labelListList(),
                faces,
                writeVTKFieldNames,
                writeVTKFieldIsPointValues,
                UPtrList<const Field<label>>(writeVTKFieldNames.size())
                #define WriteVTKFieldTypeValuesParameter(Type, nullArg) \
                    , UPtrList<const Field<Type>>(writeVTKField##Type##Values)
                FOR_ALL_FIELD_TYPES(WriteVTKFieldTypeValuesParameter)
                #undef WriteVTKFieldTypeValuesParameter
            );
        }

        Info<< endl;
    }


    void extractFeatures
    (
        const fileNameList& surfaceFileNames,
        const Time& runTime,
        const dictionary& dict
    )
    {
        forAll(surfaceFileNames, i)
        {
            extractFeatures(surfaceFileNames[i], runTime, dict);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "extract and write surface features to file"
    );
    argList::noParallel();

    #include "addDictOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    const dictionary dict(systemDict("surfaceFeaturesDict", args, runTime));

    if (dict.found("surfaces"))
    {
        extractFeatures
        (
            fileNameList(dict.lookup("surfaces")),
            runTime,
            dict
        );
    }
    else
    {
        forAllConstIter(dictionary, dict, iter)
        {
            if (!iter().isDict())
            {
                continue;
            }

            extractFeatures
            (
                fileNameList(iter().dict().lookup("surfaces")),
                runTime,
                iter().dict()
            );
        }
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}



// ************************************************************************* //
