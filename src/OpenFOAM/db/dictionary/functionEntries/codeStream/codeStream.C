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

#include "codeStream.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "Time.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(codeStream, 0);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        dictionaryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        codeStream,
        execute,
        primitiveEntryIstream
    );
}
}


const Foam::wordList Foam::functionEntries::codeStream::codeKeys
(
    {"code", "codeInclude", "localCode"}
);

const Foam::wordList Foam::functionEntries::codeStream::codeDictVars
(
    {"dict", word::null, word::null}
);

const Foam::word Foam::functionEntries::codeStream::codeTemplateC =
    "codeStreamTemplate.C";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionEntries::codeStream::masterOnlyRead
(
    const dictionary& dict
)
{
    const dictionary& topDict = dict.topDict();

    if (debug)
    {
        Pout<< "codeStream : dictionary:" << dict.name()
            << " master-only-reading:" << topDict.global()
            << endl;
    }

    return topDict.global();
}


Foam::functionEntries::codeStream::streamingFunctionType
Foam::functionEntries::codeStream::getFunction
(
    const dictionary& contextDict,
    const dictionary& codeDict
)
{
    // Get code, codeInclude, ...
    const dynamicCodeContext context
    (
        contextDict,
        codeDict,
        codeKeys,
        codeDictVars
    );

    // codeName: codeStream + _<sha1>
    // codeDir : _<sha1>
    const std::string sha1Str(context.sha1().str(true));
    dynamicCode dynCode("codeStream" + sha1Str, sha1Str);

    // Load library if not already loaded
    // Version information is encoded in the libPath (encoded with the SHA1)
    const fileName libPath = dynCode.libPath();

    // See if library is loaded
    void* lib = libs.findLibrary(libPath);

    if (debug && !lib)
    {
        Info<< "Using #codeStream with " << libPath << endl;
    }

    // Nothing loaded
    // avoid compilation if possible by loading an existing library
    if (!lib)
    {
        // Cached access to dl libs.
        // Guarantees clean up upon destruction of Time.
        if (libs.open(libPath, false))
        {
            lib = libs.findLibrary(libPath);
        }
        else
        {
            // Uncached opening of libPath. Do not complain if cannot be loaded
            lib = dlOpen(libPath, false);
        }
    }

    // Create library if required
    if (!lib)
    {
        const bool create =
            Pstream::master()
         || (regIOobject::fileModificationSkew <= 0);   // not NFS

        if (create)
        {
            if (!dynCode.upToDate(context))
            {
                // Filter with this context
                dynCode.reset(context);

                // Compile filtered C template
                dynCode.addCompileFile(codeTemplateC);

                // Define Make/options
                dynCode.setMakeOptions
                (
                    "EXE_INC = -g \\\n"
                  + context.options()
                  + "\n\nLIB_LIBS = \\\n"
                  + "    -lOpenFOAM \\\n"
                  + context.libs()
                );

                if (!dynCode.copyOrCreateFiles(true))
                {
                    FatalIOErrorInFunction
                    (
                        contextDict
                    )   << "Failed writing files for" << nl
                        << dynCode.libRelPath() << nl
                        << exit(FatalIOError);
                }
            }

            if (!dynCode.wmakeLibso())
            {
                FatalIOErrorInFunction
                (
                    contextDict
                )   << "Failed wmake " << dynCode.libRelPath() << nl
                    << exit(FatalIOError);
            }
        }

        // Only block if not master only reading of a global dictionary
        if
        (
           !masterOnlyRead(contextDict)
         && regIOobject::fileModificationSkew > 0
        )
        {
            // Determine and communicate the master file size. Scattering
            // blocks the other processes until the master has finished
            // compiling.
            off_t masterSize = Pstream::master() ? fileSize(libPath) : -1;
            Pstream::scatter(masterSize);

            // Determine the local file size. This may be incorrect if NFS is
            // taking its time, in which case we wait and try again.
            off_t mySize = Pstream::master() ? masterSize : fileSize(libPath);

            if (debug)
            {
                Pout<< endl<< "on processor " << Pstream::myProcNo()
                    << " have masterSize:" << masterSize
                    << " and localSize:" << mySize
                    << endl;
            }

            if (mySize < masterSize)
            {
                if (debug)
                {
                    Pout<< "Local file " << libPath
                        << " not of same size (" << mySize
                        << ") as master ("
                        << masterSize << "). Waiting for "
                        << regIOobject::fileModificationSkew
                        << " seconds." << endl;
                }
                Foam::sleep(regIOobject::fileModificationSkew);

                // Recheck local size
                mySize = Foam::fileSize(libPath);

                if (mySize < masterSize)
                {
                    FatalIOErrorInFunction
                    (
                        contextDict
                    )   << "Cannot read (NFS mounted) library " << nl
                        << libPath << nl
                        << "on processor " << Pstream::myProcNo()
                        << " detected size " << mySize
                        << " whereas master size is " << masterSize
                        << " bytes." << nl
                        << "If your case is not NFS mounted"
                        << " (so distributed) set fileModificationSkew"
                        << " to 0"
                        << exit(FatalIOError);
                }
            }

            if (debug)
            {
                Pout<< endl<< "on processor " << Pstream::myProcNo()
                    << " after waiting: have masterSize:" << masterSize
                    << " and localSize:" << mySize
                    << endl;
            }
        }

        if (libs.open(libPath, false))
        {
            if (debug)
            {
                Pout<< "Opening cached dictionary:" << libPath << endl;
            }

            lib = libs.findLibrary(libPath);
        }
        else
        {
            // Uncached opening of libPath
            if (debug)
            {
                Pout<< "Opening uncached dictionary:" << libPath << endl;
            }

            lib = dlOpen(libPath, true);
        }
    }

    if (!lib)
    {
        FatalIOErrorInFunction
        (
            contextDict
        )   << "Failed loading library " << libPath << nl
            << "Did you add all libraries to the 'libs' entry"
            << " in system/controlDict?"
            << exit(FatalIOError);
    }

    bool haveLib = lib;
    if (!masterOnlyRead(contextDict))
    {
        reduce(haveLib, andOp<bool>());
    }

    if (!haveLib)
    {
        FatalIOErrorInFunction
        (
            contextDict
        )   << "Failed loading library " << libPath
            << " on some processors."
            << exit(FatalIOError);
    }


    // Find the function handle in the library
    const streamingFunctionType function =
        reinterpret_cast<streamingFunctionType>
        (
            dlSym(lib, dynCode.codeName())
        );


    if (!function)
    {
        FatalIOErrorInFunction
        (
            contextDict
        )   << "Failed looking up symbol " << dynCode.codeName()
            << " in library " << lib << exit(FatalIOError);
    }

    return function;
}


Foam::string Foam::functionEntries::codeStream::run
(
    const dictionary& contextDict,
    Istream& is
)
{
    if (debug)
    {
        Info<< "Using #codeStream at line " << is.lineNumber()
            << " in file " <<  contextDict.name() << endl;
    }

    dynamicCode::checkSecurity
    (
        "functionEntries::codeStream::execute(..)",
        contextDict
    );

    // Construct codeDict for codeStream
    // Parent dictionary provided for string expansion and variable substitution
    const dictionary codeDict("#codeStream", contextDict, is);

    const streamingFunctionType function = getFunction(contextDict, codeDict);

    // Use function to write stream
    OStringStream os(is.format());
    (*function)(os, contextDict);

    // Return the string containing the results of the code execution
    return os.str();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeStream::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    return insert(contextDict, run(contextDict, is));
}


bool Foam::functionEntries::codeStream::execute
(
    const dictionary& contextDict,
    primitiveEntry& thisEntry,
    Istream& is
)
{
    return insert(contextDict, thisEntry, run(contextDict, is));
}


// ************************************************************************* //
