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

Class
    Foam::dictionary

Description
    A list of keywords followed by any number of values (e.g. words and
    numbers) or sub-dictionaries

    The keywords can represent patterns which are matched using Posix regular
    expressions. The general order for searching is as follows:
    - exact match
    - pattern match (in reverse order)
    - optional recursion into the enclosing (parent) dictionaries

    The dictionary class is the base class for IOdictionary.
    It also serves as a bootstrap dictionary for the objectRegistry data
    dictionaries since, unlike the IOdictionary class, it does not use an
    objectRegistry itself to work.

SourceFiles
    dictionaryStatics.C
    dictionary.C
    dictionaryIO.C
    dictionaryTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef dictionary_H
#define dictionary_H

#include "className.H"
#include "entry.H"
#include "IDLList.H"
#include "DLList.H"
#include "ITstream.H"
#include "Pair.H"
#include "fieldTypes.H"
#include "Tuple3.H"
#include <tuple>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class dictionary;
class primitiveEntry;
class regExp;
class SHA1Digest;
class unitConversion;

Istream& operator>>(Istream&, dictionary&);
Ostream& operator<<(Ostream&, const dictionary&);

/*---------------------------------------------------------------------------*\
                       Class dictionaryName Declaration
\*---------------------------------------------------------------------------*/

class dictionaryName
{
    // Private Data

        fileName name_;


public:

    // Constructors

        //- Construct dictionaryName null
        dictionaryName()
        {}

        //- Construct dictionaryName as copy of the given fileName
        dictionaryName(const fileName& name)
        :
            name_(name)
        {}

        //- Move constructor
        dictionaryName(dictionaryName&& name)
        :
            name_(move(name.name_))
        {}


    // Member Functions

        //- Return the dictionary name
        const fileName& name() const
        {
            return name_;
        }

        //- Return the dictionary name
        fileName& name()
        {
            return name_;
        }

        //- Return the local dictionary name (final part of scoped name)
        const word dictName() const
        {
            const word scopedName = name_.name();

            string::size_type i = scopedName.rfind('/');

            if (i == scopedName.npos)
            {
                return scopedName;
            }
            else
            {
                return scopedName.substr(i + 1, scopedName.npos);
            }
        }


    // Member Operators

        void operator=(const dictionaryName& name)
        {
            name_ = name.name_;
        }

        void operator=(dictionaryName&& name)
        {
            name_ = move(name.name_);
        }
};


/*---------------------------------------------------------------------------*\
                         Class dictionary Declaration
\*---------------------------------------------------------------------------*/

class dictionary
:
    public dictionaryName,
    public IDLList<entry>
{
    // Private Data

        //- HashTable of the entries held on the DL-list for quick lookup
        HashTable<entry*> hashedEntries_;

        //- Parent dictionary
        const dictionary& parent_;

        //- Current stream/file pointer
        mutable const Istream* filePtr_;

        //- Entries of matching patterns
        DLList<entry*> patternEntries_;

        //- Patterns as precompiled regular expressions
        DLList<autoPtr<regExp>> patternRegexps_;


    // Private Member Functions

        //- Find and return an entry data stream pointer if present
        //  otherwise return nullptr.
        //  Allows scoping using '/' with special handling for '!' and '..'.
        const entry* lookupScopedSubEntryPtr
        (
            const word&,
            bool recursive,
            bool patternMatch
        ) const;

        //- Search patterns table for exact match or regular expression match
        bool findInPatterns
        (
            const bool patternMatch,
            const word& Keyword,
            DLList<entry*>::const_iterator& wcLink,
            DLList<autoPtr<regExp>>::const_iterator& reLink
        ) const;

        //- Search patterns table for exact match or regular expression match
        bool findInPatterns
        (
            const bool patternMatch,
            const word& Keyword,
            DLList<entry*>::iterator& wcLink,
            DLList<autoPtr<regExp>>::iterator& reLink
        );

        //- Check that no unit conversions are being performed
        void assertNoConvertUnits
        (
            const char* typeName,
            const word& keyword,
            const unitConversion& defaultUnits,
            ITstream& is
        ) const;

        //- Read a value, check its dimensions and convert its units
        template<class T>
        T readTypeAndConvertUnits
        (
            const word& keyword,
            const unitConversion& defaultUnits,
            ITstream& is
        ) const;

        //- Read a value from the token stream
        template<class T>
        T readType
        (
            const word& keyword,
            const unitConversion& defaultUnits,
            ITstream& is
        ) const;

        //- Read a value from the token stream
        template<class T>
        T readType(const word& keyword, ITstream& is) const;

        //- Assign multiple entries, overwriting any existing entries
        template<class ... Entries, size_t ... Indices>
        void set
        (
            const std::tuple<const Entries& ...>&,
            const std::integer_sequence<size_t, Indices ...>&
        );

        //- Assign multiple entries, overwriting any existing entries
        template<class ... Entries>
        void set(const std::tuple<const Entries& ...>&);


    // Private Classes

        class includedDictionary;


public:

    //- Declare friendship with the entry class for IO
    friend class entry;


    // Declare name of the class and its debug switch
    ClassName("dictionary");


    // Public static data

        //- Null dictionary
        static const dictionary null;

        //- If true write optional keywords and values
        //  if not present in dictionary
        static int writeOptionalEntries;


    // Static Member Functions

        //- Construct an entries tuple from which to make a dictionary
        template<class ... Entries>
        static std::tuple<const Entries& ...> entries(const Entries& ...);


    // Constructors

        //- Construct top-level dictionary null
        dictionary();

        //- Construct top-level empty dictionary with given name
        dictionary(const fileName& name);

        //- Construct an empty sub-dictionary with given name and parent
        dictionary(const fileName& name, const dictionary& parentDict);

        //- Construct given the name, parent dictionary and Istream,
        //  reading entries until lastEntry or EOF
        dictionary
        (
            const fileName& name,
            const dictionary& parentDict,
            Istream&
        );

        //- Construct top-level dictionary from Istream,
        //  reading entries until EOF, optionally keeping the header
        dictionary(Istream&, const bool keepHeader=false);

        //- Construct as copy given the parent dictionary
        dictionary(const dictionary& parentDict, const dictionary&);

        //- Construct top-level dictionary as copy
        dictionary(const dictionary&);

        //- Construct top-level dictionary as copy from pointer to dictionary.
        //  A null pointer is treated like an empty dictionary.
        dictionary(const dictionary*);

        //- Construct top-level dictionary with given entries
        template<class ... Entries>
        dictionary(const std::tuple<const Entries& ...>&);

        //- Construct top-level dictionary with given name and entries
        template<class ... Entries>
        dictionary
        (
            const fileName& name,
            const std::tuple<const Entries& ...>&
        );

        //- Construct dictionary with given name, parent and entries
        template<class ... Entries>
        dictionary
        (
            const fileName& name,
            const dictionary& parentDict,
            const std::tuple<const Entries& ...>&
        );

        //- Construct dictionary as copy and add a list of entries
        template<class ... Entries>
        dictionary
        (
            const dictionary& dict,
            const std::tuple<const Entries& ...>&
        );

        //- Construct and return clone
        autoPtr<dictionary> clone() const;

        //- Construct top-level dictionary on freestore from Istream
        static autoPtr<dictionary> New(Istream&);


    //- Destructor
    virtual ~dictionary();


    // Member Functions

        //- Return the parent dictionary
        const dictionary& parent() const
        {
            return parent_;
        }

        //- Return whether this dictionary is null
        bool isNull() const
        {
            return this == &null;
        }

        //- Return the top of the tree
        const dictionary& topDict() const;

        //- Return the scoped keyword with which this dictionary can be
        //  accessed from the top dictionary in the tree
        word topDictKeyword() const;

        //- Return line number of first token in dictionary
        virtual label startLineNumber() const;

        //- Return line number of last token in dictionary
        virtual label endLineNumber() const;

        //- Return the SHA1 digest of the dictionary contents
        SHA1Digest digest() const;

        //- Return the dictionary as a list of tokens
        tokenList tokens() const;


        // Search and lookup

            //- Search dictionary for given keyword
            //  If recursive, search parent dictionaries
            //  If patternMatch, use regular expressions
            bool found
            (
                const word&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find and return an entry data stream pointer if present
            //  otherwise return nullptr.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions
            const entry* lookupEntryPtr
            (
                const word&,
                bool recursive,
                bool patternMatch
            ) const;

            //- Find and return an entry data stream pointer for manipulation
            //  if present otherwise return nullptr.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            entry* lookupEntryPtr
            (
                const word&,
                bool recursive,
                bool patternMatch
            );

            //- Find and return an entry data stream if present, trying a list
            //  of keywords in sequence, otherwise return nullptr.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions
            const entry* lookupEntryPtrBackwardsCompatible
            (
                const wordList&,
                bool recursive,
                bool patternMatch
            ) const;

            //- Find and return an entry data stream if present otherwise error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            const entry& lookupEntry
            (
                const word&,
                bool recursive,
                bool patternMatch
            ) const;

            //- Find and return an entry data stream if present, trying a list
            //  of keywords in sequence, otherwise error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions
            const entry& lookupEntryBackwardsCompatible
            (
                const wordList&,
                bool recursive,
                bool patternMatch
            ) const;

            //- Find and return an entry data stream
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            ITstream& lookup
            (
                const word&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find and return an entry data stream, trying a list of keywords
            //  in sequence
            //  if not found throw a fatal error relating to the first keyword
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            ITstream& lookupBackwardsCompatible
            (
                const wordList&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find and return a T, if not found throw a fatal error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookup
            (
                const word&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find and return a T, with dimension checking and unit
            //  conversions, and if not found throw a fatal error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookup
            (
                const word&,
                const unitConversion&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find and return a T, trying a list of keywords in sequence,
            //  and if not found throw a fatal error relating to the first
            //  (preferred) keyword.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupBackwardsCompatible
            (
                const wordList&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find and return a T, with dimension checking and unit
            //  conversions, trying a list of keywords in sequence, and if not
            //  found throw a fatal error relating to the first (preferred)
            //  keyword.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupBackwardsCompatible
            (
                const wordList&,
                const unitConversion&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find and return a T, if not found return the given default
            //  value.
            template<class T>
            T lookupOrDefault
            (
                const word&,
                const T&,
                const bool writeDefault = writeOptionalEntries > 0
            ) const;

            //- Find and return a T with dimension checking and unit
            //  conversions, and if not found return the given default value.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupOrDefault
            (
                const word&,
                const unitConversion&,
                const T&,
                const bool writeDefault = writeOptionalEntries > 0
            ) const;

            //- Find and return a T, trying a list of keywords in sequence,
            //  and if not found throw a fatal error relating to the first
            //  (preferred) keyword
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupOrDefaultBackwardsCompatible
            (
                const wordList&,
                const T&
            ) const;

            //- Find and return a T, with dimension checking and unit
            //  conversions, trying a list of keywords in sequence, and if not
            //  found throw a fatal error relating to the first (preferred)
            //  keyword
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupOrDefaultBackwardsCompatible
            (
                const wordList&,
                const unitConversion&,
                const T&
            ) const;

            //- Find and return a T, if not found return the given
            //  default value, and add to dictionary.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupOrAddDefault
            (
                const word&,
                const T&
            );

            //- Find an entry if present, and assign to T.
            //  Returns true if the entry was found.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            bool readIfPresent
            (
                const word&,
                T&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find an entry if present, and assign to T, with dimension
            //  checking and unit conversions.
            //  Returns true if the entry was found.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            bool readIfPresent
            (
                const word&,
                const unitConversion&,
                T&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find and return an entry data stream pointer if present,
            //  otherwise return nullptr.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            //  Allows scoping using '/' with special handling for '!' and '..'.
            const entry* lookupScopedEntryPtr
            (
                const word&,
                bool recursive,
                bool patternMatch
            ) const;

            //- Find and return a T,
            //  if not found throw a fatal error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            //  Allows scoping using '/' with special handling for '!' and '..'.
            template<class T>
            T lookupScoped
            (
                const word&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find return the reference to the compound T,
            //  if not found or not a compound throw a fatal error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            //  Allows scoping using '/' with special handling for '!' and '..'.
            template<class T>
            const T& lookupCompoundScoped
            (
                const word& keyword,
                bool recursive,
                bool patternMatch
            ) const;

            //- Check if entry is a sub-dictionary
            bool isDict(const word&) const;

            //- Find and return a sub-dictionary pointer if present
            //  otherwise return nullptr.
            const dictionary* subDictPtr(const word&) const;

            //- Find and return a sub-dictionary pointer if present
            //  otherwise return nullptr.
            dictionary* subDictPtr(const word&);

            //- Find and return a sub-dictionary
            const dictionary& subDict(const word&) const;

            //- Find and return a sub-dictionary for manipulation
            dictionary& subDict(const word&);

            //- Find and return a sub-dictionary, trying a list of keywords in
            //  sequence, otherwise error.
            const dictionary& subDictBackwardsCompatible(const wordList&) const;

            //- Find and return a sub-dictionary
            //  or empty dictionary if the sub-dictionary does not exist
            const dictionary& subOrEmptyDict
            (
                const word&,
                const bool mustRead = false
            ) const;

            //- Find and return a sub-dictionary if found
            //  otherwise return this dictionary
            const dictionary& optionalSubDict(const word&) const;

            //- Find and return a sub-dictionary by scoped lookup
            //  i.e. the keyword may contain scope characters.
            //  If the keyword is null this dictionary is returned
            const dictionary& scopedDict(const word&) const;

            //- Find and return a sub-dictionary by scoped lookup
            //  i.e. the keyword may contain scope characters.
            //  If the keyword is null this dictionary is returned
            dictionary& scopedDict(const word&);

            //- Return the table of contents
            wordList toc() const;

            //- Return the sorted table of contents
            wordList sortedToc() const;

            //- Return the list of available keys or patterns
            List<keyType> keys(bool patterns=false) const;


        // Editing

            //- Substitute the given keyword prepended by '$' with the
            //  corresponding sub-dictionary entries
            bool substituteKeyword(const word& keyword);

            //- Add a new entry
            //  With the merge option, dictionaries are interwoven and
            //  primitive entries are overwritten
            bool add(entry*, bool mergeEntry=false);

            //- Add an entry
            //  With the merge option, dictionaries are interwoven and
            //  primitive entries are overwritten
            void add(const entry&, bool mergeEntry=false);

            //- Add a word entry
            //  optionally overwrite an existing entry
            void add(const keyType&, const word&, bool overwrite=false);

            //- Add a string entry
            //  optionally overwrite an existing entry
            void add(const keyType&, const string&, bool overwrite=false);

            //- Add a label entry
            //  optionally overwrite an existing entry
            void add(const keyType&, const label, bool overwrite=false);

            //- Add a scalar entry
            //  optionally overwrite an existing entry
            void add(const keyType&, const scalar, bool overwrite=false);

            //- Add a dictionary entry
            //  optionally merge with an existing sub-dictionary
            void add
            (
                const keyType&,
                const dictionary&,
                bool mergeEntry=false
            );

            //- Add a T entry
            //  optionally overwrite an existing entry
            template<class T>
            void add(const keyType&, const T&, bool overwrite=false);

            //- Assign a new entry, overwrite any existing entry
            void set(entry*);

            //- Assign a new entry, overwrite any existing entry
            void set(const entry&);

            //- Assign a dictionary entry, overwrite any existing entry
            void set(const keyType&, const dictionary&);

            //- Assign a dictionary entry, overwrite any existing entry
            template<class ... Entries>
            void set(const keyType&, const std::tuple<const Entries& ...>&);

            //- Assign a T entry, overwrite any existing entry
            template<class T>
            void set(const keyType&, const T&);

            //- Assign multiple entries, overwriting any existing entries
            template<class ... Entries>
            void set(const entry& e, const Entries& ...);

            //- Assign multiple T entries, overwriting any existing entries
            template<class T, class ... Entries>
            void set(const keyType&, const T&, const Entries& ...);

            //- Remove an entry specified by keyword
            bool remove(const word&);

            //- Remove entries specified by keywords
            void remove(const wordList&);

            //- Change the keyword for an entry,
            //  optionally forcing overwrite of an existing entry
            bool changeKeyword
            (
                const keyType& oldKeyword,
                const keyType& newKeyword,
                bool forceOverwrite=false
            );

            //- Merge entries from the given dictionary.
            //  Also merge sub-dictionaries as required.
            bool merge(const dictionary&);

            //- Clear the dictionary
            void clear();

            //- Transfer the contents of the argument and annul the argument.
            void transfer(dictionary&);


        // Read

            //- Read dictionary from Istream, optionally keeping the header
            bool read(Istream&, const bool keepHeader=false);

            //- Return true if the dictionary global,
            //  i.e. the same on all processors.
            //  Defaults to false, must be overridden by global IO dictionaries
            virtual bool global() const;


        // Write

            //- Write dictionary, normally with sub-dictionary formatting
            void write(Ostream&, const bool subDict=true) const;


    // Member Operators

        //- Find and return entry
        ITstream& operator[](const word&) const;

        void operator=(const dictionary&);

        //- Include entries from the given dictionary.
        //  Warn, but do not overwrite existing entries.
        void operator+=(const dictionary&);

        //- Conditionally include entries from the given dictionary.
        //  Do not overwrite existing entries.
        void operator|=(const dictionary&);

        //- Unconditionally include entries from the given dictionary.
        //  Overwrite existing entries.
        void operator<<=(const dictionary&);


    // IOstream Operators

        //- Read dictionary from Istream
        friend Istream& operator>>(Istream&, dictionary&);

        //- Write dictionary to Ostream
        friend Ostream& operator<<(Ostream&, const dictionary&);
};


// Private Classes

class dictionary::includedDictionary
:
    public dictionary
{
    // Private Data

        //- Global IO status inherited from the parent dictionary
        bool global_;


public:

    // Constructors

        //- Construct an included dictionary for the given parent
        //  setting the "global" status dictionary without setting the parent
        includedDictionary
        (
            const fileName& fName,
            const dictionary& parentDict
        );


    //- Destructor
    virtual ~includedDictionary()
    {}


    // Member Functions

        //- Return true if the dictionary global,
        //  i.e. the same on all processors.
        //  Inherited from the parent dictionary into which this is included
        virtual bool global() const
        {
            return global_;
        }
};


// Template Specialisations

//- Specialise readType for types for which unit conversions can be performed
#define DECLARE_SPECIALISED_READ_TYPE(T, nullArg)                              \
                                                                               \
    template<>                                                                 \
    T Foam::dictionary::readType                                               \
    (                                                                          \
        const word& keyword,                                                   \
        const unitConversion& defaultUnits,                                    \
        ITstream& is                                                           \
    ) const;                                                                   \
                                                                               \
    template<>                                                                 \
    T Foam::dictionary::readType                                               \
    (                                                                          \
        const word& keyword,                                                   \
        ITstream& is                                                           \
    ) const;

#define DECLARE_SPECIALISED_READ_PAIR_TYPE(T, nullArg)                         \
    DECLARE_SPECIALISED_READ_TYPE(Pair<T>, nullArg)

#define DECLARE_SPECIALISED_READ_LIST_TYPE(T, nullArg)                         \
    DECLARE_SPECIALISED_READ_TYPE(List<T>, nullArg)

FOR_ALL_FIELD_TYPES(DECLARE_SPECIALISED_READ_TYPE);
FOR_ALL_FIELD_TYPES(DECLARE_SPECIALISED_READ_PAIR_TYPE);
FOR_ALL_FIELD_TYPES(DECLARE_SPECIALISED_READ_LIST_TYPE);

#undef DECLARE_SPECIALISED_READ_TYPE
#undef DECLARE_SPECIALISED_READ_PAIR_TYPE
#undef DECLARE_SPECIALISED_READ_LIST_TYPE


// Global Operators

//- Combine dictionaries.
//  Starting from the entries in dict1 and then including those from dict2.
//  Warn, but do not overwrite the entries from dict1.
dictionary operator+(const dictionary& dict1, const dictionary& dict2);

//- Combine dictionaries.
//  Starting from the entries in dict1 and then including those from dict2.
//  Do not overwrite the entries from dict1.
dictionary operator|(const dictionary& dict1, const dictionary& dict2);


// Global Functions

//- Parse dictionary substitution argument list
void dictArgList
(
    const Tuple2<string, label>& argString,
    word& configName,
    List<Tuple2<wordRe, label>>& args,
    List<Tuple3<word, string, label>>& namedArgs
);

//- Parse dictionary substitution argument list
void dictArgList
(
    const Tuple2<string, label>& argString,
    List<Tuple2<wordRe, label>>& args,
    List<Tuple3<word, string, label>>& namedArgs
);

//- Extracts dict name and keyword
Pair<word> dictAndKeyword(const word& scopedName);

//- Return the list of configuration files in
//  user/group/shipped directories.
//  The search scheme allows for version-specific and
//  version-independent files using the following hierarchy:
//  - \b user settings:
//    - ~/.OpenFOAM/\<VERSION\>/caseDicts/functions
//    - ~/.OpenFOAM/caseDicts/functions
//  - \b group (site) settings (when $WM_PROJECT_SITE is set):
//    - $WM_PROJECT_SITE/\<VERSION\>/etc/caseDicts/functions
//    - $WM_PROJECT_SITE/etc/caseDicts/functions
//  - \b group (site) settings (when $WM_PROJECT_SITE is not set):
//    - $WM_PROJECT_INST_DIR/site/\<VERSION\>/etc/
//          caseDicts/functions
//    - $WM_PROJECT_INST_DIR/site/etc/caseDicts/functions
//  - \b other (shipped) settings:
//    - $WM_PROJECT_DIR/etc/caseDicts/functions
wordList listAllConfigFiles
(
    const fileName& configFilesPath
);

//- Search for configuration file for given region
//  and if not present also search the case directory as well as the
//  user/group/shipped directories.
//  The search scheme allows for version-specific and
//  version-independent files using the following hierarchy:
//  - \b user settings:
//    - ~/.OpenFOAM/\<VERSION\>/caseDicts/functions
//    - ~/.OpenFOAM/caseDicts/functions
//  - \b group (site) settings (when $WM_PROJECT_SITE is set):
//    - $WM_PROJECT_SITE/\<VERSION\>/etc/caseDicts/functions
//    - $WM_PROJECT_SITE/etc/caseDicts/functions
//  - \b group (site) settings (when $WM_PROJECT_SITE is not set):
//    - $WM_PROJECT_INST_DIR/site/\<VERSION\>/etc/
//          caseDicts/functions
//    - $WM_PROJECT_INST_DIR/site/etc/caseDicts/functions
//  - \b other (shipped) settings:
//    - $WM_PROJECT_DIR/etc/caseDicts/functions
//
//  \return The path of the configuration file if found
//  otherwise null
fileName findConfigFile
(
    const word& configName,
    const fileName& configFilesPath,
    const word& configFilesDir,
    const word& region = word::null
);

//- Expand arg within the dict context and return
string expandArg
(
    const string& arg,
    dictionary& dict,
    const label lineNumber
);

//- Add the keyword value pair to dict
//  setting the given lineNumber for the entry
void addArgEntry
(
    dictionary& dict,
    const word& keyword,
    const string& value,
    const label lineNumber
);

//- Read the specified configuration file
//  parsing the optional arguments included in the string
//  'argString', inserting 'field' or 'fields' entries as required
//  and merging the resulting configuration dictionary into
//  'parentDict'.
//
//  Parses the optional arguments:
//      'Q(U)' -> configFileName = Q; args = (U)
//             -> field U;
//
//  Supports named arguments:
//      'patchAverage(patch=inlet, p,U)'
//  or
//      'patchAverage(patch=inlet, field=(p U))'
//       -> configFileName = patchAverage;
//          args = (patch=inlet, p,U)
//       -> patch inlet;
//          fields (p U);
bool readConfigFile
(
    const word& configType,
    const Tuple2<string, label>& argString,
    dictionary& parentDict,
    const fileName& configFilesPath,
    const word& configFilesDir,
    const word& region = word::null
);

//- Write a dictionary entry
void writeEntry(Ostream& os, const dictionary& dict);

//- Helper function to write the keyword and entry
template<class EntryType>
void writeEntry(Ostream& os, const word& entryName, const EntryType& value);

//- Helper function to write the keyword and entry
template<class EntryType>
void writeEntry
(
    Ostream& os,
    const word& entryName,
    const unitConversion& defaultUnits,
    const EntryType& value
);

//- Helper function to write the keyword and entry only if the
//  values are not equal. The value is then output as value2
template<class EntryType>
void writeEntryIfDifferent
(
    Ostream& os,
    const word& entryName,
    const EntryType& value1,
    const EntryType& value2
);

//- Helper function to write the keyword and entry only if the
//  values are not equal. The value is then output as value2
template<class EntryType>
void writeEntryIfDifferent
(
    Ostream& os,
    const word& entryName,
    const unitConversion& defaultUnits,
    const EntryType& value1,
    const EntryType& value2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dictionaryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
