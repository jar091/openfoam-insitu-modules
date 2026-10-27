/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenFOAM-Mesio integration project
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

#include "mesioVDB.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "clockTime.H"
#include "OStringStream.H"
#include "openfoam_mpi.H"
#include "addToRunTimeSelectionTable.H"

#include "mesio_insitu.h"

#if !defined(MESIO_INSITU_VERSION) || MESIO_INSITU_VERSION < 2
    #error "mesioVDB requires the Mesio in-situ API version 2 (error codes, several contexts): build against install/mesio-next or newer"
#endif

#include <iomanip>
#include <type_traits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mesioVDB, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        mesioVDB,
        dictionary
    );
}
}


const Foam::Enum
<
    Foam::functionObjects::mesioVDB::missingAction
>
Foam::functionObjects::mesioVDB::missingActionNames_
({
    { missingAction::ERROR, "error" },
    { missingAction::WARN, "warn" },
    { missingAction::SKIP, "skip" },
});


const Foam::Enum
<
    Foam::functionObjects::mesioVDB::staleAction
>
Foam::functionObjects::mesioVDB::staleActionNames_
({
    { staleAction::ERROR, "error" },
    { staleAction::WARN, "warn" },
    { staleAction::IGNORE, "ignore" },
});


const Foam::Enum
<
    Foam::functionObjects::mesioVDB::policyType
>
Foam::functionObjects::mesioVDB::policyTypeNames_
({
    { policyType::STANDARD, "standard" },
    { policyType::VDB_ONLY, "vdbOnly" },
});


const Foam::Enum
<
    Foam::functionObjects::mesioVDB::errorAction
>
Foam::functionObjects::mesioVDB::errorActionNames_
({
    { errorAction::ERROR, "error" },
    { errorAction::WARN, "warn" },
});


// The adapter passes OpenFOAM storage to Mesio without conversion
static_assert
(
    sizeof(Foam::label) == sizeof(MESIOInt),
    "OpenFOAM label and MESIOInt differ: build Mesio with --intwidth matching"
    " WM_LABEL_SIZE"
);

static_assert
(
    std::is_same_v<Foam::scalar, MESIOReal>,
    "mesioVDB requires a double precision (DP) OpenFOAM build"
);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Zero-padded integer
std::string zeroPadded(const Foam::label val, const Foam::label width)
{
    std::ostringstream os;
    os << std::setw(width) << std::setfill('0') << val;
    return os.str();
}


// Max over all ranks
Foam::scalar maxOverRanks(const Foam::scalar val)
{
    return Foam::returnReduce(val, Foam::maxOp<Foam::scalar>());
}


// Value of one rank (collective)
Foam::scalar valueOfRank(const Foam::scalar val, const Foam::label rank)
{
    return Foam::returnReduce
    (
        Foam::UPstream::myProcNo() == rank ? val : 0,
        Foam::sumOp<Foam::scalar>()
    );
}


// Exponent as text: integer if possible
std::string exponentName(const Foam::scalar e)
{
    if (e == Foam::label(e))
    {
        return std::to_string(Foam::label(e));
    }
    std::ostringstream os;
    os << e;
    return os.str();
}

} // End anonymous namespace


// SI units of a dimension set, e.g. "kg m^-1 s^-2"; "1" for dimensionless
Foam::string Foam::functionObjects::mesioVDB::unitsName
(
    const dimensionSet& dims
)
{
    static const char* names[dimensionSet::nDimensions] =
        {"kg", "m", "s", "K", "mol", "A", "cd"};

    std::string units;
    for (int d = 0; d < dimensionSet::nDimensions; ++d)
    {
        const scalar e = dims[d];
        if (Foam::mag(e) > dimensionSet::smallExponent)
        {
            units += (units.empty() ? "" : " ") + std::string(names[d]);
            if (Foam::mag(e - 1) > dimensionSet::smallExponent)
            {
                units += "^" + exponentName(e);
            }
        }
    }
    return units.empty() ? "1" : units;
}


// OpenFOAM dimension set, e.g. "[0 1 -1 0 0 0 0]"
Foam::string Foam::functionObjects::mesioVDB::dimensionsName
(
    const dimensionSet& dims
)
{
    std::string text("[");
    for (int d = 0; d < dimensionSet::nDimensions; ++d)
    {
        text += (d ? " " : "") + exponentName(dims[d]);
    }
    return text + "]";
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::string Foam::functionObjects::mesioVDB::settingsSignature() const
{
    OStringStream os;
    os  << flatOutput(fieldNames_) << ' '
        << flatOutput(vectorComponents_) << ' '
        << flatOutput(tensorComponents_) << ' '
        << voxelSize_ << ' ' << density_ << ' '
        << useBox_ << ' ' << box_ << ' '
        << movingMesh_ << ' ' << writerRank_ << ' '
        << writer_ << ' ' << compression_;

    return os.str();
}


void Foam::functionObjects::mesioVDB::checkOutputPolicy(const bool fatal) const
{
    DynamicList<string> issues;

    const dictionary& cd = time_.controlDict();

    // Solver (Time) write schedule inside [startTime, endTime]
    {
        const word wc = cd.getOrDefault<word>("writeControl", "timeStep");
        const scalar wi = cd.get<scalar>("writeInterval");
        const scalar t0 = time_.startTime().value();
        const scalar t1 = time_.endTime().value();
        const scalar dt = time_.deltaTValue();

        bool writes = true;

        if (wc == "timeStep")
        {
            const label i0 = time_.startTimeIndex();
            const label nSteps = label((t1 - t0)/dt + 0.5);
            const label interval = max(label(1), label(wi + 0.5));
            const label firstWrite = (i0/interval + 1)*interval;
            writes = (firstWrite <= i0 + nSteps);
        }
        else if
        (
            wc == "runTime"
         || wc == "adjustableRunTime"
         || wc == "adjustable"
        )
        {
            writes =
            (
                label(t1/wi + SMALL) > label(t0/wi + SMALL)
            );
        }

        if (writes)
        {
            issues.push_back
            (
                "controlDict: writeControl " + wc + ", writeInterval "
              + Foam::name(wi)
              + " writes OpenFOAM time directories before endTime"
                " (or cannot be verified); use writeControl timeStep with"
                " writeInterval larger than the number of steps"
            );
        }
    }

    // Other function objects must not write
    if (const dictionary* fdPtr = cd.findDict("functions"))
    {
        for (const entry& e : *fdPtr)
        {
            if (!e.isDict() || e.keyword() == name())
            {
                continue;
            }

            const dictionary& d = e.dict();

            // Other mesioVDB instances write .vdb files only
            if
            (
                !d.getOrDefault<bool>("enabled", true)
             || d.getOrDefault<word>("type", word::null) == typeName
            )
            {
                continue;
            }

            const word wc = d.getOrDefaultCompat<word>
            (
                "writeControl",
                {{"outputControl", 1606}},
                "timeStep"
            );

            if (wc != "none")
            {
                issues.push_back
                (
                    "function object " + e.keyword()
                  + " (type " + d.getOrDefault<word>("type", "?")
                  + ") uses writeControl " + wc
                  + "; set writeControl none"
                );
            }
        }
    }

    if (issues.empty())
    {
        Info<< type() << ' ' << name()
            << ": outputPolicy vdbOnly - configuration audit passed" << nl
            << "    (no solver writes and no function-object writes"
            << " before endTime)" << endl;
        return;
    }

    OStringStream msg;
    msg << "Function object " << name()
        << ": outputPolicy vdbOnly violated:" << nl;
    for (const string& s : issues)
    {
        msg << "    - " << s.c_str() << nl;
    }

    if (fatal)
    {
        FatalErrorInFunction
            << msg.str().c_str()
            << "Use 'policyAction warn' to continue anyway."
            << exit(FatalError);
    }
    else
    {
        WarningInFunction
            << msg.str().c_str() << endl;
    }
}


bool Foam::functionObjects::mesioVDB::resolveGrids()
{
    grids_.clear();

    wordHashSet seen;

    for (const word& fieldName : fieldNames_)
    {
        if (!seen.insert(fieldName))
        {
            continue;
        }

        const bool found =
        (
            addFieldGrids<scalar>(fieldName)
         || addFieldGrids<vector>(fieldName)
         || addFieldGrids<sphericalTensor>(fieldName)
         || addFieldGrids<symmTensor>(fieldName)
         || addFieldGrids<tensor>(fieldName)
        );

        const bool foundAny = returnReduceOr(found);
        const bool foundAll = returnReduceAnd(found);

        if (foundAny && !foundAll)
        {
            FatalErrorInFunction
                << "Function object " << name() << ": field " << fieldName
                << " exists only on some ranks"
                << exit(FatalError);
        }

        if (found)
        {
            continue;
        }

        // Registered, but of an unsupported type?
        const regIOobject* objPtr = obr().cfindObject<regIOobject>(fieldName);
        const bool unsupported = returnReduceOr(objPtr != nullptr);

        if (unsupported)
        {
            FatalErrorInFunction
                << "Function object " << name() << ": object " << fieldName
                << " has unsupported type "
                << (objPtr ? objPtr->type() : word("unknown")) << nl
                << "    Supported: volScalarField, volVectorField,"
                << " volSphericalTensorField, volSymmTensorField,"
                << " volTensorField and their internal fields" << nl
                << "    (surface and point fields are not supported)"
                << exit(FatalError);
        }

        switch (missingAction_)
        {
            case missingAction::ERROR:
            {
                FatalErrorInFunction
                    << "Function object " << name() << ": field "
                    << fieldName << " not found in registry "
                    << obr().name() << " at time " << time_.timeName() << nl
                    << "    Derived fields must be computed by function"
                    << " objects listed before " << name()
                    << " and executed at this time step." << nl
                    << "    Available volScalarFields: "
                    << flatOutput(obr().sortedNames<volScalarField>()) << nl
                    << "    Available volVectorFields: "
                    << flatOutput(obr().sortedNames<volVectorField>())
                    << exit(FatalError);
                break;
            }
            case missingAction::WARN:
            {
                WarningInFunction
                    << "Function object " << name() << ": field "
                    << fieldName << " not found at time "
                    << time_.timeName()
                    << "; it is not exported in this run" << endl;
                break;
            }
            case missingAction::SKIP:
            {
                break;
            }
        }
    }

    return !grids_.empty();
}


void Foam::functionObjects::mesioVDB::checkFreshness
(
    const label refEvent
) const
{
    if (staleAction_ == staleAction::IGNORE)
    {
        return;
    }

    wordHashSet checked;

    for (const gridEntry& g : grids_)
    {
        if (!checked.insert(g.field))
        {
            continue;
        }

        const regIOobject* objPtr = obr().cfindObject<regIOobject>(g.field);

        if (!objPtr || objPtr->eventNo() > refEvent)
        {
            continue;
        }

        if (staleAction_ == staleAction::ERROR)
        {
            FatalErrorInFunction
                << "Function object " << name() << ": field " << g.field
                << " was not updated in the current time step "
                << time_.timeName() << " (stale)." << nl
                << "    Align the executeControl/executeInterval of its"
                << " producer with the writeControl/writeInterval of "
                << name() << " and list the producer before " << name()
                << exit(FatalError);
        }
        else
        {
            WarningInFunction
                << "Function object " << name() << ": field " << g.field
                << " was not updated in time step " << time_.timeName()
                << " (stale data exported)" << endl;
        }
    }
}


bool Foam::functionObjects::mesioVDB::mesioOk
(
    const int code,
    const char* call
) const
{
    if (code == MESIO_INSITU_OK)
    {
        return true;
    }

    const char* msg = MESIOInSituLastError(ctx_);

    if (errorAction_ == errorAction::ERROR)
    {
        FatalErrorInFunction
            << "Function object " << name() << ": " << call
            << " failed (Mesio error " << code << "):" << nl
            << "    " << msg << nl
            << "    Use 'onError warn' to skip failed exports instead."
            << exit(FatalError);
    }

    WarningInFunction
        << "Function object " << name() << ": " << call
        << " failed (Mesio error " << code << "):" << nl
        << "    " << msg << nl
        << "    Export at time " << time_.timeName() << " skipped"
        << " (onError warn)" << endl;

    return false;
}


void Foam::functionObjects::mesioVDB::mesioLocal
(
    const int code,
    const char* call
) const
{
    if (code != MESIO_INSITU_OK)
    {
        // Reported again (on all ranks) by the next collective Mesio call
        Pout<< "--> FOAM Warning : function object " << name() << ": "
            << call << " failed (Mesio error " << code << "): "
            << MESIOInSituLastError(ctx_) << endl;
    }
}


bool Foam::functionObjects::mesioVDB::applyOptions()
{
    int code = MESIO_INSITU_OK;

    if (!writer_.empty() || !compression_.empty())
    {
        const int writer =
        (
            writer_ == "single"
          ? MESIO_INSITU_WRITER_SINGLE
          : writer_ == "parallel"
          ? MESIO_INSITU_WRITER_PARALLEL
          : MESIO_INSITU_WRITER_DISTRIBUTED
        );
        const int compression =
        (
            compression_ == "none"
          ? MESIO_INSITU_COMPRESSION_NONE
          : MESIO_INSITU_COMPRESSION_BLOSC
        );
        code = MESIOInSituSetWriterOptions(ctx_, writer, compression);
        mesioLocal(code, "MESIOInSituSetWriterOptions");
    }

    if (code == MESIO_INSITU_OK && !sampling_.empty())
    {
        code = MESIOInSituSetSampling
        (
            ctx_,
            sampling_ == "cell"
          ? MESIO_INSITU_SAMPLING_CELL
          : MESIO_INSITU_SAMPLING_CELLPOINT
        );
        mesioLocal(code, "MESIOInSituSetSampling");
    }

    if (code == MESIO_INSITU_OK)
    {
        code = MESIOInSituSetMemoryTracking(ctx_, memoryTracking_ ? 1 : 0);
        mesioLocal(code, "MESIOInSituSetMemoryTracking");
    }

    return code == MESIO_INSITU_OK;
}


void Foam::functionObjects::mesioVDB::setMesioMesh()
{
    // Local polyMesh in CSR face format
    const pointField& points = mesh_.points();
    const faceList& faces = mesh_.faces();

    labelList offsets(faces.size() + 1);
    offsets[0] = 0;
    forAll(faces, facei)
    {
        offsets[facei + 1] = offsets[facei] + faces[facei].size();
    }

    labelList faceNodes(offsets.back());
    forAll(faces, facei)
    {
        SubList<label>(faceNodes, faces[facei].size(), offsets[facei]) =
            faces[facei];
    }

    // Errors (invalid mesh) are reported by the following Commit
    mesioLocal
    (
        MESIOInSituSetFVMesh
        (
            ctx_,
            points.size(),
            reinterpret_cast<const MESIOReal*>(points.cdata()),
            faces.size(),
            offsets.cdata(),
            faceNodes.cdata(),
            mesh_.faceOwner().cdata(),
            mesh_.nInternalFaces(),
            mesh_.faceNeighbour().cdata(),
            mesh_.nCells()
        ),
        "MESIOInSituSetFVMesh"
    );
}


bool Foam::functionObjects::mesioVDB::createContext()
{
    clockTime timer;

    if (!UPstream::parRun())
    {
        // Mesio needs MPI; OpenFOAM finalizes it at exit (ourMpi)
        UPstream::initNull();
    }

    if (!resolveGrids())
    {
        WarningInFunction
            << "Function object " << name()
            << ": none of the fields " << flatOutput(fieldNames_)
            << " can be exported at time " << time_.timeName() << endl;
        return false;
    }

    MESIOInSituSettings settings;
    MESIOInSituDefaultSettings(&settings);
    settings.density = density_;
    settings.voxelSize = voxelSize_;
    settings.useBox = useBox_ ? 1 : 0;
    if (useBox_)
    {
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            settings.boxMin[d] = box_.min()[d];
            settings.boxMax[d] = box_.max()[d];
        }
    }
    ctxMoving_ = (movingMesh_ || pointsMoved_ || mesh_.moving());
    settings.movingMesh = ctxMoving_ ? 1 : 0;
    settings.writerRank = writerRank_;
    settings.verbosity = verbosity_;

    MPI_Comm comm = MPI_COMM_WORLD;
    if (UPstream::parRun())
    {
        comm = PstreamUtils::Cast::to_mpi
        (
            UPstream::Communicator::lookup(mesh_.comm())
        );
    }

    MESIOInSitu ctx = nullptr;
    const int code = MESIOInSituCreate(&ctx, comm, &settings);
    ctx_ = ctx;
    if (!mesioOk(code, "MESIOInSituCreate"))
    {
        ++failedExports_;
        return false;
    }

    applyOptions();
    setMesioMesh();

    for (const gridEntry& g : grids_)
    {
        mesioLocal
        (
            MESIOInSituAddCellField(ctx_, g.grid.c_str(), g.dimension),
            "MESIOInSituAddCellField"
        );

        const std::pair<const char*, const string*> meta[] =
        {
            {"units", &g.units},
            {"dimensions", &g.dimensions},
            {"field", &g.field},
            {"component", &g.componentName}
        };
        for (const auto& m : meta)
        {
            mesioLocal
            (
                MESIOInSituSetFieldMetadata
                (
                    ctx_, g.grid.c_str(), m.first, m.second->c_str()
                ),
                "MESIOInSituSetFieldMetadata"
            );
        }
    }

    if (!mesioOk(MESIOInSituCommit(ctx_), "MESIOInSituCommit"))
    {
        ++failedExports_;
        destroyContext();
        return false;
    }

    ctxSignature_ = settingsSignature();
    pointsMoved_ = false;
    topoChanged_ = false;
    settingsChanged_ = false;

    setupTime_ = maxOverRanks(timer.elapsedTime());

    int dims[3] = {0, 0, 0};
    double origin[3] = {0, 0, 0};
    double h = 0;
    MESIOInSituGridInfo(ctx_, dims, origin, &h);

    wordList gridNames(grids_.size());
    forAll(grids_, i)
    {
        gridNames[i] = grids_[i].grid;
    }

    Log << type() << ' ' << name() << ": Mesio context created" << nl
        << "    grids      : " << flatOutput(gridNames) << nl
        << "    voxels     : " << dims[0] << " x " << dims[1] << " x "
        << dims[2] << ", voxel size " << h << nl
        << "    origin     : (" << origin[0] << ' ' << origin[1] << ' '
        << origin[2] << ") (centre of voxel 0)" << nl
        << "    moving mesh: " << Switch(ctxMoving_) << nl
        << "    writer     : "
        << (writer_.empty() ? word("Mesio default") : writer_) << ", compression "
        << (compression_.empty() ? word("Mesio default") : compression_)
        << ", sampling "
        << (sampling_.empty() ? word("Mesio default") : sampling_) << nl
        << "    on error   : " << errorActionNames_[errorAction_] << nl
        << "    setup time : " << setupTime_ << " s" << nl
        << "    output     : " << outputDir_ << endl;

    return true;
}


bool Foam::functionObjects::mesioVDB::resetContextMesh()
{
    clockTime timer;

    // Keeps the communicators, settings, options and registered grids
    if (!mesioOk(MESIOInSituResetMesh(ctx_), "MESIOInSituResetMesh"))
    {
        destroyContext();
        return false;
    }

    setMesioMesh();

    if (!mesioOk(MESIOInSituCommit(ctx_), "MESIOInSituCommit"))
    {
        destroyContext();
        return false;
    }

    pointsMoved_ = false;
    topoChanged_ = false;

    const scalar resetTime = maxOverRanks(timer.elapsedTime());
    setupTime_ += resetTime;

    Log << type() << ' ' << name() << ": mesh changed, Mesio mesh rebuilt ("
        << returnReduce(mesh_.nCells(), sumOp<label>()) << " cells) in "
        << resetTime << " s" << endl;

    return true;
}


void Foam::functionObjects::mesioVDB::destroyContext()
{
    if (ctx_)
    {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized)
        {
            MESIOInSituDestroy(ctx_);
        }
        ctx_ = nullptr;
    }
    grids_.clear();
}


Foam::fileName Foam::functionObjects::mesioVDB::outputFile() const
{
    string file(filePattern_);
    file.replaceAll("%NAME%", name());
    file.replaceAll
    (
        "%TIMEINDEX%",
        zeroPadded(time_.timeIndex(), indexWidth_)
    );
    file.replaceAll("%TIME%", time_.timeName());
    file.replaceAll("%COUNT%", zeroPadded(exportCount_, indexWidth_));

    return outputDir_/fileName(file);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mesioVDB::mesioVDB
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldNames_(),
    vectorComponents_(one{}, "mag"),
    tensorComponents_(one{}, "mag"),
    voxelSize_(0),
    density_(100),
    box_(),
    useBox_(false),
    movingMesh_(false),
    writerRank_(0),
    writer_(),
    compression_(),
    sampling_(),
    errorAction_(errorAction::ERROR),
    memoryTracking_(false),
    verbosity_(0),
    outputDir_(),
    filePattern_("%NAME%_%TIMEINDEX%.vdb"),
    indexWidth_(6),
    missingAction_(missingAction::ERROR),
    staleAction_(staleAction::WARN),
    outputPolicy_(policyType::STANDARD),
    policyFatal_(true),
    ctx_(nullptr),
    ctxMoving_(false),
    grids_(),
    ctxSignature_(),
    pointsMoved_(false),
    topoChanged_(false),
    settingsChanged_(false),
    failedExports_(0),
    settingsRead_(false),
    curMarker_(mesh_.getEvent()),
    refMarker_(curMarker_),
    markerTimeIndex_(runTime.timeIndex()),
    lastExportIndex_(-1),
    exportCount_(0),
    buffer_(),
    totalExportTime_(0),
    setupTime_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::mesioVDB::~mesioVDB()
{
    destroyContext();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mesioVDB::readSettings(const dictionary& dict)
{
    // All entries are parsed and validated into local values first; the
    // members change only if everything is valid.

    const wordList fieldNames(dict.get<wordList>("fields"));
    const wordList vectorComponents
    (
        dict.getOrDefault<wordList>("vectorComponents", wordList(one{}, "mag"))
    );
    const wordList tensorComponents
    (
        dict.getOrDefault<wordList>("tensorComponents", wordList(one{}, "mag"))
    );

    const scalar voxelSize = dict.getOrDefault<scalar>("voxelSize", 0);
    const label density = dict.getOrDefault<label>("density", 100);
    if (voxelSize <= 0 && density <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Either voxelSize > 0 or density > 0 is required"
            << exit(FatalIOError);
    }

    // The box is read as two points: the boundBox stream operator reads raw
    // binary data from binary streams, and the dictionary that the master
    // broadcasts at a run-time re-read is such a stream on the other ranks
    boundBox box;
    const bool useBox = dict.found("box");
    if (useBox)
    {
        ITstream& is = dict.lookup("box");
        point bmin, bmax;
        is >> bmin >> bmax;
        dict.checkITstream(is, "box");
        box = boundBox(bmin, bmax);
        if (!box.good())
        {
            FatalIOErrorInFunction(dict)
                << "Invalid box " << bmin << ' ' << bmax
                << " (min must not exceed max)" << exit(FatalIOError);
        }
    }

    const bool movingMesh = dict.getOrDefault<bool>("movingMesh", false);
    const label writerRank = dict.getOrDefault<label>("writerRank", 0);
    const label verbosity = dict.getOrDefault<label>("verbosity", 0);

    const word writer(dict.getOrDefault<word>("writer", word::null));
    if
    (
        !writer.empty()
     && writer != "distributed" && writer != "parallel" && writer != "single"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unknown writer " << writer
            << " (distributed, parallel or single)" << exit(FatalIOError);
    }
    const word compression(dict.getOrDefault<word>("compression", word::null));
    if (!compression.empty() && compression != "blosc" && compression != "none")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown compression " << compression
            << " (blosc or none)" << exit(FatalIOError);
    }
    const word sampling(dict.getOrDefault<word>("sampling", word::null));
    if (!sampling.empty() && sampling != "cellPoint" && sampling != "cell")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown sampling " << sampling
            << " (cellPoint or cell)" << exit(FatalIOError);
    }
    const errorAction onError = errorActionNames_.getOrDefault
    (
        "onError",
        dict,
        errorAction::ERROR
    );
    const bool memoryTracking = dict.getOrDefault<bool>("memoryTracking", false);

    if (writerRank < 0 || writerRank >= UPstream::nProcs())
    {
        FatalIOErrorInFunction(dict)
            << "writerRank " << writerRank << " outside [0,"
            << UPstream::nProcs() << ')' << exit(FatalIOError);
    }

    fileName outputDir(dict.getOrDefault<fileName>("outputDir", "VDB"));
    outputDir.expand();
    if (!outputDir.isAbsolute())
    {
        outputDir = time_.globalPath()/outputDir;
    }
    outputDir.clean();

    const string filePattern
    (
        dict.getOrDefault<string>("filePattern", "%NAME%_%TIMEINDEX%.vdb")
    );
    const label indexWidth = dict.getOrDefault<label>("indexWidth", 6);

    const missingAction missing = missingActionNames_.getOrDefault
    (
        "missingField",
        dict,
        missingAction::ERROR
    );
    const staleAction stale = staleActionNames_.getOrDefault
    (
        "staleField",
        dict,
        staleAction::WARN
    );
    const policyType policy = policyTypeNames_.getOrDefault
    (
        "outputPolicy",
        dict,
        policyType::STANDARD
    );
    const bool policyFatal =
        (dict.getOrDefault<word>("policyAction", "error") != "warn");

    if (policy == policyType::VDB_ONLY)
    {
        checkOutputPolicy(policyFatal);
    }

    // Valid: take over
    const word oldSampling(sampling_);
    const bool oldMemoryTracking = memoryTracking_;
    const errorAction oldErrorAction = errorAction_;

    fieldNames_ = fieldNames;
    vectorComponents_ = vectorComponents;
    tensorComponents_ = tensorComponents;
    voxelSize_ = voxelSize;
    density_ = density;
    box_ = box;
    useBox_ = useBox;
    movingMesh_ = movingMesh;
    writerRank_ = writerRank;
    verbosity_ = verbosity;
    writer_ = writer;
    compression_ = compression;
    sampling_ = sampling;
    errorAction_ = onError;
    memoryTracking_ = memoryTracking;
    outputDir_ = outputDir;
    filePattern_ = filePattern;
    indexWidth_ = indexWidth;
    missingAction_ = missing;
    staleAction_ = stale;
    outputPolicy_ = policy;
    policyFatal_ = policyFatal;

    if (ctx_)
    {
        if (ctxSignature_ != settingsSignature())
        {
            // Mesh/grid settings changed: the context is re-created at the
            // next export (collectively, all ranks read the same dictionary)
            Info<< type() << ' ' << name()
                << ": settings changed, Mesio context is re-created at the"
                << " next export" << endl;
            settingsChanged_ = true;
        }
        else if
        (
            sampling_ != oldSampling
         || memoryTracking_ != oldMemoryTracking
         || errorAction_ != oldErrorAction
        )
        {
            // Applied to the existing context
            applyOptions();
            Info<< type() << ' ' << name() << ": settings re-read, sampling "
                << (sampling_.empty() ? word("default") : sampling_)
                << " applied, memoryTracking " << Switch(memoryTracking_)
                << ", onError " << errorActionNames_[errorAction_]
                << " (no re-creation)" << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::mesioVDB::read(const dictionary& dict)
{
    if (!settingsRead_)
    {
        // Construction: invalid settings are fatal
        fvMeshFunctionObject::read(dict);
        readSettings(dict);
        settingsRead_ = true;
        return true;
    }

    // Run-time re-read (runTimeModifiable): OpenFOAM does not trap errors of
    // read(), so an invalid entry would stop the run. Invalid settings are
    // reported and the previous settings are kept. All ranks read the same
    // dictionary (broadcast by the master), so the outcome is the same on
    // all ranks.
    const bool oldThrowingError = FatalError.throwing(true);
    const bool oldThrowingIOError = FatalIOError.throwing(true);

    try
    {
        fvMeshFunctionObject::read(dict);
        readSettings(dict);
    }
    catch (const Foam::error& err)
    {
        WarningInFunction
            << "Function object " << name()
            << ": invalid settings after the run-time re-read of the"
            << " dictionary; settings not changed:" << nl;
        err.write(Warning, false);
        Warning<< endl;
    }

    FatalError.throwing(oldThrowingError);
    FatalIOError.throwing(oldThrowingIOError);

    return true;
}


bool Foam::functionObjects::mesioVDB::execute()
{
    const label timeIndex = time_.timeIndex();

    if (timeIndex != markerTimeIndex_)
    {
        refMarker_ = curMarker_;
        curMarker_ = mesh_.getEvent();
        markerTimeIndex_ = timeIndex;
    }

    return true;
}


bool Foam::functionObjects::mesioVDB::write()
{
    const label timeIndex = time_.timeIndex();

    if (timeIndex == lastExportIndex_)
    {
        return true;
    }

    // Reference event: fields must have been modified after it
    const label refEvent =
    (
        markerTimeIndex_ == timeIndex ? refMarker_ : curMarker_
    );

    clockTime timer;

    if (ctx_ && (settingsChanged_ || (pointsMoved_ && !ctxMoving_)))
    {
        Log << type() << ' ' << name()
            << ": settings or mesh motion changed, re-creating the Mesio"
            << " context" << endl;
        destroyContext();
    }

    if (ctx_ && topoChanged_)
    {
        if (!resetContextMesh())
        {
            ++failedExports_;
            return true;
        }
        timer.resetTime();
    }
    else if (!ctx_)
    {
        if (!createContext())
        {
            return true;
        }
        timer.resetTime();
    }
    else if (pointsMoved_)
    {
        mesioLocal
        (
            MESIOInSituUpdatePoints
            (
                ctx_,
                reinterpret_cast<const MESIOReal*>(mesh_.points().cdata())
            ),
            "MESIOInSituUpdatePoints"
        );
        pointsMoved_ = false;
    }

    checkFreshness(refEvent);

    for (const gridEntry& g : grids_)
    {
        const bool ok =
        (
            setGridValues<scalar>(g)
         || setGridValues<vector>(g)
         || setGridValues<sphericalTensor>(g)
         || setGridValues<symmTensor>(g)
         || setGridValues<tensor>(g)
        );

        if (!ok)
        {
            if (missingAction_ == missingAction::ERROR)
            {
                FatalErrorInFunction
                    << "Function object " << name() << ": field " << g.field
                    << " disappeared from the registry at time "
                    << time_.timeName() << exit(FatalError);
            }

            WarningInFunction
                << "Function object " << name() << ": field " << g.field
                << " not found at time " << time_.timeName()
                << "; grid " << g.grid << " is written as zero" << endl;

            buffer_.resize_nocopy(g.dimension*mesh_.nCells());
            buffer_ = Zero;
            mesioLocal
            (
                MESIOInSituSetCellField(ctx_, g.grid.c_str(), buffer_.cdata()),
                "MESIOInSituSetCellField"
            );
        }
    }

    const double fieldTime = timer.elapsedTime();

    const fileName file(outputFile());

    const int code = MESIOInSituWriteVDB(ctx_, time_.value(), file.c_str());
    if (!mesioOk(code, "MESIOInSituWriteVDB"))
    {
        ++failedExports_;
        if (code != MESIO_INSITU_ERROR_IO && code != MESIO_INSITU_ERROR_ARGUMENT)
        {
            // Mesio failed internally: rebuild the mesh at the next export
            topoChanged_ = true;
        }
        return true;
    }

    const scalar totalTime = maxOverRanks(timer.elapsedTime());

    lastExportIndex_ = timeIndex;
    ++exportCount_;
    totalExportTime_ += totalTime;

    if (log)
    {
        double seconds[4] = {0, 0, 0, 0};
        long long activeVoxels = 0;
        long long fileBytes = 0;
        MESIOInSituLastWriteStats(ctx_, seconds, &activeVoxels, &fileBytes);

        scalar maxSec[4];
        for (label i = 0; i < 4; ++i)
        {
            maxSec[i] = maxOverRanks(seconds[i]);
        }
        const scalar maxField = maxOverRanks(fieldTime);

        Info<< type() << ' ' << name() << " write: " << file << nl
            << "    voxels " << int64_t(activeVoxels) << ", size "
            << scalar(fileBytes)/(1024.0*1024.0) << " MB" << nl
            << "    time [s] " << totalTime
            << " (fields " << maxField
            << ", convert " << maxSec[0]
            << ", sample " << maxSec[1]
            << ", gather " << maxSec[2]
            << ", write " << maxSec[3] << ')' << nl;

        double timings[MESIO_INSITU_TIMINGS];
        MESIOInSituLastWriteTimings(ctx_, timings);
        for (label i = 0; i < MESIO_INSITU_TIMINGS; ++i)
        {
            timings[i] = maxOverRanks(timings[i]);
        }
        Info<< "    phases [s] (max over ranks): plan "
            << timings[MESIO_INSITU_TIME_PLAN]
            << ", sample " << timings[MESIO_INSITU_TIME_SAMPLE]
            << ", exchange " << timings[MESIO_INSITU_TIME_EXCHANGE]
            << ", merge+compress " << timings[MESIO_INSITU_TIME_MERGE]
            << ", gather " << timings[MESIO_INSITU_TIME_GATHER]
            << ", file " << timings[MESIO_INSITU_TIME_WRITE] << nl;

        long long memory[MESIO_INSITU_MEMORY];
        MESIOInSituLastWriteMemory(ctx_, memory);
        const scalar MB = 1024.0*1024.0;
        const scalar rss = scalar(memory[MESIO_INSITU_MEMORY_RSS])/MB;
        const scalar peak = scalar(memory[MESIO_INSITU_MEMORY_PEAK])/MB;
        const scalar buffers = scalar(memory[MESIO_INSITU_MEMORY_BUFFERS])/MB;
        Info<< "    memory [MB]: writer rank RSS "
            << valueOfRank(rss, writerRank_)
            << ", peak " << valueOfRank(peak, writerRank_)
            << (memoryTracking_ ? " (this export)" : " (since start)")
            << ", Mesio buffers " << valueOfRank(buffers, writerRank_)
            << "; max over ranks: peak " << maxOverRanks(peak)
            << ", Mesio buffers " << maxOverRanks(buffers) << nl;
        Info<< endl;
    }

    return true;
}


bool Foam::functionObjects::mesioVDB::end()
{
    Log << type() << ' ' << name() << " end: " << exportCount_
        << " export(s), total export time " << totalExportTime_
        << " s, setup " << setupTime_ << " s";
    if (failedExports_)
    {
        Log << ", " << failedExports_ << " export(s) skipped after Mesio errors";
    }
    Log << endl;

    return true;
}


void Foam::functionObjects::mesioVDB::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == static_cast<const polyMesh*>(&mesh_))
    {
        topoChanged_ = true;
    }
}


void Foam::functionObjects::mesioVDB::movePoints(const polyMesh& mesh)
{
    if (&mesh == static_cast<const polyMesh*>(&mesh_))
    {
        pointsMoved_ = true;
    }
}


// ************************************************************************* //
