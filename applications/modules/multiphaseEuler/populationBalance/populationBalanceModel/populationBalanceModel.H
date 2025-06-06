/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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
    Foam::diameterModels::populationBalanceModel

Description
    Model for tracking the evolution of a dispersed phase size distribution due
    to coalescence (synonymous with coagulation, aggregation, agglomeration)
    and breakup events as well as density or phase changes. Provides an
    approximate solution of the population balance equation by means of a class
    method. The underlying theory is described in the article of Lehnigk et al.
    (2021).

    The size distribution, expressed through a volume-based number density
    function, is discretised using the fixot pivot technique of Kumar and
    Ramkrishna (1996). Thereby, the population balance equation is transformed
    into a series of transport equations for the particle (bubble, droplet)
    number concentrations in separate size classes that are coupled through
    their source terms. The discretisation is based on representative particle
    volumes, which are provided by the user through the corresponding sphere
    equivalent diameters.

    Since the representative volumes are fixed a priori and the total dispersed
    phase volume already available from solving the phase continuity equation,
    the model only determines the evolution of the individual size class
    fractions

    \f[
        f_{i,\varphi} = \frac{\alpha_{i,\varphi}}{\alpha_{\varphi}}\,,
    \f]

    where \f$\alpha_{i,\varphi}\f$ is the volume fraction of the size class and
    \f$\alpha_{\varphi}\f$ the total phase fraction of phase \f$\varphi\f$.

    The source terms are formulated such that the first and second moment of
    the distribution, i.e. the total particle number and volume, are conserved
    irrespective of the discretisation of the size domain. The treatment of
    particle breakup depends on the selected breakup submodels. For models
    which provide a total breakup frequency and a separate daughter size
    distribution function, the formulation provided Kumar and Ramkrishna (1996)
    is utilised, which is applicable both for binary and multiple breakup
    events. Currently, only field-independent daughter size distribution models
    are allowed. In case of binary breakup models that provide the breakup
    frequency between a size class pair directly, the formulation of Liao et
    al. (2018) is adopted, which is computationally more efficient compared to
    first extracting the field-dependent daughter size distribution and then
    consuming it in the formulation of Kumar and Ramkrishna. The source terms
    describing a drift of the size distribution through particle growth or
    shrinkage are derived using upwind differencing, thus ensuring conservation
    of the total particle number and volume. Note that due to the volume-based
    implementation, both density as well as phase change lead to a drift of the
    size distribution function. Further, users can specify multiple submodels
    for each mechanism, whose contributions are added up.

    The model also allows to distribute the size classes over multiple
    representative phases with identical physical properties that collectively
    define the dispersed phase. Thereby, size class fields can be transported
    with different velocity fields in order to account for the size dependency
    of the particle motion. A possible mass transfer between representative
    phases by means of coalescence, breakup and drift is taken into account.
    Similarly, the spatial evolution of secondary particle properties such as
    the particle surface area can be tracked.

    The key variable during a simulation is the Sauter diameter, which is
    computed from the size class fractions of the corresponding phase. The
    underlying size distribution can be extracted from the simulation using the
    functionObject 'sizeDistribution'. Integral and mean properties of a size
    distribution can be computed with the functionObject 'moments'.

    Verification cases for the population balance modelling functionality are
    provided in test/multiphaseEuler/populationBalance.

    References:
    \verbatim
        Lehnigk, R., Bainbridge, W., Liao, Y., Lucas, D., Niemi, T.,
        Peltola, J., & Schlegel, F. (2021).
        An open‐source population balance modeling framework for the simulation
        of polydisperse multiphase flows.
        AIChE Journal, 68(3), e17539.
    \endverbatim

    \verbatim
        Coalescence and breakup term formulation:
        Kumar, S., & Ramkrishna, D. (1996).
        On the solution of population balance equations by discretization-I. A
        fixed pivot technique.
        Chemical Engineering Science, 51(8), 1311-1332.
    \endverbatim

    \verbatim
        Binary breakup term formulation:
        Liao, Y., Oertel, R., Kriebitzsch, S., Schlegel, F., & Lucas, D. (2018).
        A discrete population balance equation for binary breakage.
        International Journal for Numerical Methods in Fluids, 87(4), 202-215.
    \endverbatim

Usage
    Excerpt from an exemplary phaseProperties dictionary:
    \verbatim
    populationBalances (bubbles);

    air
    {
        type            pureIsothermalPhaseModel;

        diameterModel   velocityGroup;

        velocityGroupCoeffs
        {
            populationBalance    bubbles;

            shapeModel           spherical;

            sizeGroups
            (
                f1 { dSph 1e-3; }
                f2 { dSph 2e-3; }
                f3 { dSph 3e-3; }
                f4 { dSph 4e-3; }
                f5 { dSph 5e-3; }
                ...
            );
        }

        residualAlpha   1e-6;
    }

    ...

    populationBalanceCoeffs
    {
        bubbles
        {
            continuousPhase water;

            coalescenceModels
            (
                LehrMilliesMewes{}
            );

            binaryBreakupModels
            (
                LehrMilliesMewes{}
            );

            breakupModels
            ();
        }
    }
    \endverbatim

See also
    Foam::PopulationBalancePhaseSystem
    Foam::diameterModels::sizeGroup
    Foam::diameterModels::velocityGroup
    Foam::diameterModels::SecondaryPropertyModel
    Foam::diameterModels::coalescenceModel
    Foam::diameterModels::breakupModel
    Foam::diameterModels::daughterSizeDistributionModel
    Foam::diameterModels::binaryBreakupModel
    Foam::functionObjects::sizeDistribution
    Foam::functionObjects::moments

SourceFiles
    populationBalanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceModel_H
#define populationBalanceModel_H

#include "sizeGroup.H"
#include "phaseSystem.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distribution;
class phaseCompressibleMomentumTransportModel;

namespace diameterModels
{

class coalescenceModel;
class breakupModel;
class binaryBreakupModel;

/*---------------------------------------------------------------------------*\
                   Class populationBalanceModel Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceModel
:
    public regIOobject
{
public:

    // Public Type Definitions

        //- Table of interfacial mass transfer rates
        typedef
            HashPtrTable
            <
                volScalarField::Internal,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            >
            dmdtfTable;


    // Public Classes

        //- Class to accumulate population balance sub-class pointers
        class groups
        :
            public regIOobject
        {
            // Private Data

                //- Velocity group pointers
                HashTable<const velocityGroup*> velocityGroupPtrs_;

                //- Size group pointers
                UPtrList<sizeGroup> sizeGroups_;


            // Private Constructor

                //- Construct with a name on a database
                groups(const word& popBalName, const objectRegistry& db)
                :
                    regIOobject
                    (
                        IOobject
                        (
                            name(popBalName),
                            db.time().constant(),
                            db
                        )
                    )
                {}


            // Private Member Functions

                //- Return the name of the pointer cache
                static word name(const word& popBalName)
                {
                    return popBalName + ":groups";
                }


        public:

            // Constructors

                //- Lookup in the registry or construct new
                static groups& New
                (
                    const word& popBalName,
                    const objectRegistry& db
                )
                {
                    if (db.foundObject<groups>(name(popBalName)))
                    {
                        return db.lookupObjectRef<groups>(name(popBalName));
                    }

                    groups* ps = new groups(popBalName, db);

                    ps->store();

                    return *ps;
                }


            // Member Functions

                //- Return the current number of size groups
                label nSizeGroups()
                {
                    return sizeGroups_.size();
                }

                //- Insert a velocity group into the table
                void insert(velocityGroup& group)
                {
                    velocityGroupPtrs_.insert(group.phase().name(), &group);

                    const label i0 = nSizeGroups();

                    sizeGroups_.resize(i0 + group.sizeGroups().size());

                    forAll(group.sizeGroups(), i)
                    {
                        sizeGroups_.set(i0 + i, &group.sizeGroups()[i]);

                        if
                        (
                            i0 + i != 0
                         && sizeGroups_[i0 + i - 1].x().value()
                         >= sizeGroups_[i0 + i].x().value()
                        )
                        {
                            FatalErrorInFunction
                                << "Size groups must be specified in order of "
                                << "their representative size"
                                << exit(FatalError);
                        }
                    }
                }

                //- Retrieve the pointers
                static void retrieve
                (
                    const populationBalanceModel& popBal,
                    HashTable<const velocityGroup*>& velGroupPtrs,
                    UPtrList<sizeGroup>& szGroupPtrs
                )
                {
                    const objectRegistry& db = popBal.fluid().mesh();

                    if (!db.foundObject<groups>(name(popBal.name())))
                    {
                        FatalErrorInFunction
                            << "No velocity groups exist for population "
                            << "balance \"" << popBal.name() << "\""
                            << exit(FatalError);
                    }

                    groups& ps =
                        db.lookupObjectRef<groups>(name(popBal.name()));

                    velGroupPtrs.transfer(ps.velocityGroupPtrs_);
                    szGroupPtrs.transfer(ps.sizeGroups_);

                    ps.checkOut();
                }

                //- Dummy write for regIOobject
                virtual bool writeData(Ostream&) const
                {
                    return true;
                }
        };


private:

    // Private Data

        //- Reference to the phaseSystem
        const phaseSystem& fluid_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Name of the populationBalance
        const word name_;

        //- Continuous phase
        const phaseModel& continuousPhase_;

        //- Velocity groups belonging to this populationBalance
        HashTable<const velocityGroup*> velocityGroupPtrs_;

        //- Size groups belonging to this populationBalance
        UPtrList<sizeGroup> sizeGroups_;

        //- Size group boundaries
        PtrList<dimensionedScalar> v_;

        //- Section width required for binary breakup formulation
        PtrList<PtrList<dimensionedScalar>> delta_;

        //- Explicit coalescence and breakup source terms
        PtrList<volScalarField::Internal> Su_;

        //- Implicit coalescence and breakup source terms
        PtrList<volScalarField::Internal> Sp_;

        //- Coalescence and breakup mass transfer rates
        dmdtfTable dmdtfs_;

        //- Expansion mass transfer rates
        dmdtfTable expansionDmdtfs_;

        //- Model source mass transfer rates
        dmdtfTable modelSourceDmdtfs_;

        //- Rates of change of volume per unit volume
        PtrList<volScalarField::Internal> expansionRates_;

        //- Dilatation errors
        PtrList<volScalarField::Internal> dilatationErrors_;

        //- Coalescence models
        PtrList<coalescenceModel> coalescenceModels_;

        //- Coalescence rate
        autoPtr<volScalarField::Internal> coalescenceRate_;

        //- Coalescence relevant size group pairs
        List<labelPair> coalescencePairs_;

        //- Breakup models
        PtrList<breakupModel> breakupModels_;

        //- Breakup rate
        autoPtr<volScalarField::Internal> breakupRate_;

        //- Binary breakup models
        PtrList<binaryBreakupModel> binaryBreakupModels_;

        //- Binary breakup rate
        autoPtr<volScalarField::Internal> binaryBreakupRate_;

        //- Binary breakup relevant size group pairs
        List<labelPair> binaryBreakupPairs_;

        //- Total void fraction
        autoPtr<volScalarField> alphas_;

        //- Mean Sauter diameter
        autoPtr<volScalarField> dsm_;

        //- Average velocity
        autoPtr<volVectorField> U_;

        //- Counter for interval between source term updates
        label sourceUpdateCounter_;


    // Private Member Functions

        const dictionary& coeffDict() const;

        void precomputeCoalescenceAndBreakup();

        void birthByCoalescence(const label j, const label k);

        void deathByCoalescence(const label i, const label j);

        void birthByBreakup(const label k, const label model);

        void deathByBreakup(const label i);

        void birthByBinaryBreakup(const label i, const label j);

        void deathByBinaryBreakup(const label j, const label i);

        void computeCoalescenceAndBreakup();

        void precomputeExpansion();

        Pair<tmp<volScalarField::Internal>> expansionSus
        (
            const label i,
            const UPtrList<const volScalarField>& flds =
                UPtrList<const volScalarField>()
        ) const;

        void computeExpansion();

        void precomputeModelSources();

        Pair<tmp<volScalarField::Internal>> modelSourceRhoSus
        (
            const label i,
            const UPtrList<const volScalarField>& flds =
                UPtrList<const volScalarField>()
        ) const;

        void computeModelSources();

        void computeDilatationErrors();

        //- Return whether the sources should be updated on this iteration
        bool updateSources();

        //- Return the interval at which the sources are updated
        inline label sourceUpdateInterval() const;

        //- Return the coefficients of the lower half of the number allocation
        //  coefficient function. This spans the range from group i-1 to group
        //  i. The first coefficient is a constant, and the second is
        //  multiplied by v.
        Pair<dimensionedScalar> etaCoeffs0(const label i) const;

        //- Return the coefficients of the lower half of the number allocation
        //  coefficient function. This spans the range from group i to group
        //  i+1. The first coefficient is a constant, and the second is
        //  multiplied by v.
        Pair<dimensionedScalar> etaCoeffs1(const label i) const;

        //- Return the coefficients of the lower half of the volume allocation
        //  coefficient function. This spans the range from group i-1 to group
        //  i. The first coefficient is a constant, and the second is
        //  multiplied by 1/v.
        Pair<dimensionedScalar> etaVCoeffs0(const label i) const;

        //- Return the coefficients of the lower half of the volume allocation
        //  coefficient function. This spans the range from group i to group
        //  i+1. The first coefficient is a constant, and the second is
        //  multiplied by 1/v.
        Pair<dimensionedScalar> etaVCoeffs1(const label i) const;


public:

    //- Runtime type information
    TypeName("populationBalanceModel");


    // Constructors

        //- Construct for a fluid
        populationBalanceModel(const phaseSystem& fluid, const word& name);

        //- Return clone
        autoPtr<populationBalanceModel> clone() const;

        //- Return a pointer to a new populationBalanceModel object created on
        //  freestore from Istream
        class iNew
        {
            const phaseSystem& fluid_;

        public:

            iNew(const phaseSystem& fluid)
            :
                fluid_(fluid)
            {}

            autoPtr<populationBalanceModel> operator()(Istream& is) const
            {
                const word name(is);

                Info << "Setting up population balance: " << name << endl;

                return autoPtr<populationBalanceModel>
                (
                    new populationBalanceModel(fluid_, name)
                );
            }
        };


    //- Destructor
    virtual ~populationBalanceModel();


    // Member Functions

        //- Dummy write for regIOobject
        bool writeData(Ostream&) const;

        //- Return reference to the phaseSystem
        inline const phaseSystem& fluid() const;

        //- Return reference to the coalescence and breakup interfacial mass
        //  transfer rates
        inline const dmdtfTable& dmdtfs() const;

        //- Return reference to the expansion interfacial mass transfer rates
        inline const dmdtfTable& expansionDmdtfs() const;

        //- Return reference to the model source interfacial mass transfer rates
        inline const dmdtfTable& modelSourceDmdtfs() const;

        //- Return reference to the mesh
        inline const fvMesh& mesh() const;

        //- Return solution settings dictionary for this populationBalance
        inline const dictionary& solverDict() const;

        //- Solve on final pimple iteration only
        inline bool solveOnFinalIterOnly() const;

        //- Return continuous phase
        inline const phaseModel& continuousPhase() const;

        //- Return the size groups belonging to this populationBalance
        inline const UPtrList<sizeGroup>& sizeGroups() const;

        //- Access the size groups belonging to this populationBalance
        inline UPtrList<sizeGroup>& sizeGroups();

        //- Return coalescence relevant size group pairs
        inline const List<labelPair>& coalescencePairs() const;

        //- Return binary breakup relevant size group pairs
        inline const List<labelPair>& binaryBreakupPairs() const;

        //- Return total void of phases belonging to this populationBalance
        inline const volScalarField& alphas() const;

        //- Return average velocity
        inline const volVectorField& U() const;

        //- Return the number allocation coefficient for a single volume
        dimensionedScalar eta(const label i, const dimensionedScalar& v) const;

        //- Return the number allocation coefficient for a field of volumes
        tmp<volScalarField::Internal> eta
        (
            const label i,
            const volScalarField::Internal& v
        ) const;

        //- Return the volume allocation coefficient for a single volume
        dimensionedScalar etaV(const label i, const dimensionedScalar& v) const;

        //- Return the volume allocation coefficient for a field of volumes
        tmp<volScalarField::Internal> etaV
        (
            const label i,
            const volScalarField::Internal& v
        ) const;

        //- Return the volume allocation coefficient for a given distribution
        //  of diameters over a range of size-groups. The diameter distribution
        //  should be volumetrically sampled (i.e., sampleQ should equal 3).
        dimensionedScalar etaV(const labelPair is, const distribution& d) const;

        //- Return the volume allocation coefficient for a given distribution
        //  of diameters. The diameter distribution should be volumetrically
        //  sampled (i.e., sampleQ should equal 3). The result is normalised
        //  with respect to the volume allocation coefficient for the phase.
        dimensionedScalar etaV(const label i, const distribution& d) const;

        //- Return the surface tension coefficient between a given dispersed
        //  and the continuous phase
        const tmp<volScalarField> sigmaWithContinuousPhase
        (
            const phaseModel& dispersedPhase
        ) const;

        //- Return reference to momentumTransport model of the continuous phase
        const phaseCompressibleMomentumTransportModel&
            continuousTurbulence() const;

        //- Return the implicit coalescence and breakup source term
        tmp<volScalarField::Internal> Sp(const label i) const;

        //- Return the explicit expansion source term
        tmp<volScalarField::Internal> expansionSu
        (
            const label i,
            const UPtrList<const volScalarField>& flds =
                UPtrList<const volScalarField>()
        ) const;

        //- Return the implicit expansion source term
        tmp<volScalarField::Internal> expansionSp(const label i) const;

        //- Return the explicit model source source term
        tmp<volScalarField::Internal> modelSourceSu
        (
            const label i,
            const UPtrList<const volScalarField>& flds =
                UPtrList<const volScalarField>()
        ) const;

        //- Solve the population balance equation
        void solve();

        //- Correct derived quantities
        void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "populationBalanceModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
