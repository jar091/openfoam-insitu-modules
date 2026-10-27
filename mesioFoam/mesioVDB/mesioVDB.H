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

Class
    Foam::functionObjects::mesioVDB

Group
    grpUtilitiesFunctionObjects

Description
    Exports cell fields directly from memory to OpenVDB (.vdb) files using
    the Mesio in-situ API (ESPRESO/Mesio, \c mesio_insitu.h).

    All MPI ranks pass their local polyMesh (points, faces, owner, neighbour;
    processor faces as ordinary boundary faces) and their local cell values
    to Mesio. Mesio redistributes the data, averages cell values to points,
    samples them at voxel centres and writes exactly \em one .vdb file per
    export step (written to a temporary name and renamed when complete). No
    intermediate per-rank field files are written.

    Every VDB grid is an \c openvdb::FloatGrid (fog volume, background 0)
    named after the exported quantity:
    \table
      OpenFOAM object                  | Grids
      volScalarField                   | \<name\>
      volVectorField                   | \<name\>_mag, \<name\>_x, ...
      volSphericalTensorField          | \<name\>
      volSymmTensorField, volTensorField | \<name\>_mag, \<name\>_xx, ...
    \endtable
    The internal fields (DimensionedField) of these types are accepted as
    well. Surface fields, point fields and other registry objects are
    rejected with a diagnostic. Every grid carries the string metadata
    \c units (SI units from the field dimensions, e.g. "m s^-1"),
    \c dimensions (OpenFOAM dimension set "[kg m s K mol A cd]"), \c field
    (source field) and \c component, besides the \c time metadata.

    Derived fields (e.g. \c Q, \c vorticity, \c mag(U)) are computed by the
    standard OpenFOAM function objects. They must be listed \em before this
    function object in the \c functions dictionary, so that they are executed
    earlier within the same time step, and they should use
    \c writeControl \c none to avoid writing their result to disk.
    Each exported field is checked for freshness: a field that was not
    modified since the previous step (or since the previous export when this
    function object does not execute every step) is reported as stale.

    The export cadence is controlled with the standard \c writeControl and
    \c writeInterval entries. Mesh motion (same topology) is supported via
    \c MESIOInSituUpdatePoints; after topology changes the Mesio mesh is
    rebuilt in the existing context (\c MESIOInSituResetMesh).

    Several mesioVDB instances can be used in one run (e.g. two boxes or
    resolutions); each has its own Mesio context.

    At a run-time re-read of the dictionary (runTimeModifiable) invalid
    entries are reported and the previous settings are kept; changed grid
    settings (fields, box, voxelSize, density, ...) re-create the Mesio
    context at the next export, other settings apply immediately.

    Mesio errors are reported with the Mesio message: \c onError \c error
    stops the run (FatalError), \c onError \c warn reports a warning, skips
    the export and retries at the next export (the context is rebuilt if
    Mesio failed internally).

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    Q1
    {
        type            Q;
        libs            (fieldFunctionObjects);
        field           U;
        result          Q;
        executeControl  timeStep;
        executeInterval 1;
        writeControl    none;
    }

    mesioVDB1
    {
        // Mandatory entries
        type            mesioVDB;
        libs            (mesioFoam);
        fields          (p U Q);

        // Optional entries
        vectorComponents (mag);
        tensorComponents (mag);
        voxelSize       0.001;
        density         256;
        box             (0 0 0) (0.1 0.1 0.1);
        movingMesh      false;
        writerRank      0;
        writer          distributed;
        compression     blosc;
        sampling        cellPoint;
        onError         error;
        memoryTracking  false;
        verbosity       0;
        outputDir       "<case>/VDB";
        filePattern     "%NAME%_%TIMEINDEX%.vdb";
        indexWidth      6;
        missingField    error;
        staleField      warn;
        outputPolicy    standard;
        policyAction    error;

        // Inherited entries
        region          region0;
        writeControl    timeStep;
        writeInterval   10;
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                              | Type | Reqd | Deflt
      type         | Type name: mesioVDB                      | word | yes  | -
      libs         | Library name: mesioFoam                  | word | yes  | -
      fields       | Names of the fields to export            | wordList | yes | -
      vectorComponents | Grids per vector field: mag x y z    | wordList | no | (mag)
      tensorComponents | Grids per (symm)tensor: mag xx xy ... | wordList | no | (mag)
      voxelSize    | Voxel size [m]; overrides density        | scalar | no | -
      density      | Voxels along the longest box edge        | label | no  | 100
      box          | Voxelized sub-box (min max)              | boundBox | no | mesh bounds
      movingMesh   | Prepare for mesh motion from the start   | bool | no   | false
      writerRank   | MPI rank writing the .vdb file           | label | no  | 0
      writer       | distributed, parallel or single (see below) | word | no | Mesio default
      compression  | blosc or none                            | word | no   | Mesio default
      sampling     | cellPoint (interpolated) or cell (cell value) | word | no | Mesio default
      onError      | Mesio error: error (stop) or warn (skip export) | word | no | error
      memoryTracking | Peak RSS per export (resets the process VmHWM) | bool | no | false
      verbosity    | Mesio progress output (0 or 1)           | label | no  | 0
      outputDir    | Output directory (relative: to case)     | fileName | no | VDB
      filePattern  | File name pattern (see below)            | string | no | %NAME%_%TIMEINDEX%.vdb
      indexWidth   | Zero padding of %TIMEINDEX% and %COUNT%  | label | no  | 6
      missingField | Missing field: error, warn or skip       | word | no   | error
      staleField   | Not updated field: error, warn or ignore | word | no   | warn
      outputPolicy | standard or vdbOnly (audit, see below)   | word | no   | standard
      policyAction | vdbOnly violation: error or warn         | word | no   | error
    \endtable

    Writers (Mesio): \c distributed (default): OpenVDB leaf nodes are built, merged
    and compressed on all ranks, the writer rank writes the file; \c parallel: as
    distributed, but one rank per node writes its part of the leaf data into the same
    file (faster for large files); \c single: legacy single-writer path.
    The default can also be set by the environment variable MESIO_VDB_WRITER.

    Sampling (Mesio): \c cellPoint (default): cell values are averaged to the
    mesh points and interpolated inside the cells (smooth, but extreme values
    in boundary cells are smoothed); \c cell: every voxel gets the value of the
    cell containing its centre (exact cell values, piecewise constant).
    The default can also be set by the environment variable MESIO_VDB_SAMPLING.

    File pattern tokens: \c %NAME% (function object name), \c %TIME%
    (time name), \c %TIMEINDEX% (time step index, zero padded),
    \c %COUNT% (export counter of this run, zero padded).

    With \c outputPolicy \c vdbOnly the configuration is audited at start-up:
    the solver write schedule (\c writeControl / \c writeInterval of the
    controlDict) must not trigger a write before \c endTime, and every other
    enabled function object (except other mesioVDB entries) must use
    \c writeControl \c none. Violations are
    reported as errors (or warnings with \c policyAction \c warn). Logs,
    pre-existing case files, processor directories and explicit file writes
    of third-party code are not controlled.

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link timeControlFunctionObject.H \endlink

Note
    MPI is required by Mesio. In a serial run the library initializes MPI
    (UPstream::initNull) and OpenFOAM finalizes it at exit.

SourceFiles
    mesioVDB.C
    mesioVDBTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_mesioVDB_H
#define Foam_functionObjects_mesioVDB_H

#include "fvMeshFunctionObject.H"
#include "boundBox.H"
#include "dimensionSet.H"
#include "Enum.H"
#include "DynamicList.H"
#include "scalarField.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Opaque Mesio context (struct MESIOInSituData* in mesio_insitu.h)
struct MESIOInSituData;

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class mesioVDB Declaration
\*---------------------------------------------------------------------------*/

class mesioVDB
:
    public fvMeshFunctionObject
{
public:

    // Public Data Types

        //- Handling of missing fields
        enum class missingAction : char { ERROR, WARN, SKIP };

        //- Handling of fields not updated since the reference event
        enum class staleAction : char { ERROR, WARN, IGNORE };

        //- Disk output policy
        enum class policyType : char { STANDARD, VDB_ONLY };

        //- Handling of Mesio errors
        enum class errorAction : char { ERROR, WARN };

        //- Names for missingAction
        static const Enum<missingAction> missingActionNames_;

        //- Names for staleAction
        static const Enum<staleAction> staleActionNames_;

        //- Names for policyType
        static const Enum<policyType> policyTypeNames_;

        //- Names for errorAction
        static const Enum<errorAction> errorActionNames_;


private:

    // Private Data Types

        //- One exported VDB grid
        struct gridEntry
        {
            //- VDB grid name
            word grid;

            //- Source field name in the registry
            word field;

            //- Component index (>= 0), FULL or MAG
            label component;

            //- Dimension passed to Mesio (1, or 3 = magnitude in Mesio)
            int dimension;

            //- Component name written as metadata (value, mag, x, xx, ...)
            word componentName;

            //- SI units of the field (e.g. "m s^-1")
            string units;

            //- OpenFOAM dimension set (e.g. "[0 1 -1 0 0 0 0]")
            string dimensions;
        };

        //- Special component selectors
        static constexpr label FULL = -1;
        static constexpr label MAG = -2;


    // Private Data

        //- Requested field names
        wordList fieldNames_;

        //- Grids per vector field
        wordList vectorComponents_;

        //- Grids per tensor/symmTensor field
        wordList tensorComponents_;

        //- Voxel size (<= 0: use density)
        scalar voxelSize_;

        //- Voxels along the longest edge
        label density_;

        //- Optional voxelization box
        boundBox box_;

        //- Use the box
        bool useBox_;

        //- Mesh motion requested by the user
        bool movingMesh_;

        //- Writer rank
        label writerRank_;

        //- Writer: distributed or single (empty: Mesio default/environment)
        word writer_;

        //- Compression: blosc or none (empty: Mesio default/environment)
        word compression_;

        //- Sampling: cellPoint or cell (empty: Mesio default/environment)
        word sampling_;

        //- Handling of Mesio errors
        errorAction errorAction_;

        //- Report the peak RSS of every export (resets VmHWM)
        bool memoryTracking_;

        //- Mesio verbosity
        label verbosity_;

        //- Output directory (absolute after read)
        fileName outputDir_;

        //- File name pattern
        string filePattern_;

        //- Zero padding width
        label indexWidth_;

        //- Missing field action
        missingAction missingAction_;

        //- Stale field action
        staleAction staleAction_;

        //- Output policy
        policyType outputPolicy_;

        //- Policy violation is fatal
        bool policyFatal_;

        //- Mesio context (nullptr until the first export)
        MESIOInSituData* ctx_;

        //- Context was created with movingMesh = 1
        bool ctxMoving_;

        //- Exported grids (fixed at context creation)
        DynamicList<gridEntry> grids_;

        //- Signature of the settings used for the current context
        string ctxSignature_;

        //- Points moved since the last export
        bool pointsMoved_;

        //- Topology changed: rebuild the Mesio mesh at the next export
        bool topoChanged_;

        //- Settings changed: re-create the context at the next export
        bool settingsChanged_;

        //- Number of exports skipped after Mesio errors (onError warn)
        label failedExports_;

        //- Settings were read (construction done): later reads are
        //- run-time re-reads
        bool settingsRead_;

        //- Event marker of the most recent execute() call
        label curMarker_;

        //- Event marker of the execute() before that
        label refMarker_;

        //- Time index of curMarker_
        label markerTimeIndex_;

        //- Time index of the last export (avoid double export)
        label lastExportIndex_;

        //- Number of exports in this run
        label exportCount_;

        //- Scratch buffer for component / magnitude extraction
        scalarField buffer_;

        //- Accumulated export wall time (max over ranks) [s]
        scalar totalExportTime_;

        //- Setup (create/mesh/commit) wall time [s]
        scalar setupTime_;



    // Private Member Functions

        //- SI units of a dimension set (e.g. "kg m^-1 s^-2", "1")
        static string unitsName(const dimensionSet& dims);

        //- OpenFOAM dimension set as text (e.g. "[0 1 -1 0 0 0 0]")
        static string dimensionsName(const dimensionSet& dims);

        //- Signature of the settings that require a new context
        string settingsSignature() const;

        //- Audit the configuration for the vdbOnly policy (violations are
        //- fatal or warnings)
        void checkOutputPolicy(const bool fatal) const;

        //- Parse and validate all settings; the members change only if
        //- all entries are valid (FatalIOError otherwise)
        bool readSettings(const dictionary& dict);

        //- Resolve the grid list from the registry (collective).
        //  Returns false if nothing can be exported.
        bool resolveGrids();

        //- Add the grids of one field of a given type, if it is one
        template<class Type>
        bool addFieldGrids(const word& fieldName);

        //- Find cell values of a field (volume or internal field)
        template<class Type>
        const Field<Type>* findCellValues
        (
            const word& fieldName,
            label& eventNo
        ) const;

        //- Dimensions of a field (volume or internal field), nullptr if
        //- not found
        template<class Type>
        const dimensionSet* findDimensions(const word& fieldName) const;

        //- Pass the values of one grid to Mesio. False if the field type
        //- does not match Type.
        template<class Type>
        bool setGridValues(const gridEntry& g);

        //- Check that the exported fields were updated after refEvent
        void checkFreshness(const label refEvent) const;

        //- Check the return code of a collective Mesio call (the same on
        //- all ranks): FatalError (onError error) or a warning (false)
        bool mesioOk(const int code, const char* call) const;

        //- Check the return code of a rank-local Mesio call: a failure is
        //- reported on this rank and by the next collective call
        void mesioLocal(const int code, const char* call) const;

        //- Pass the local polyMesh to Mesio (rank-local)
        void setMesioMesh();

        //- Create the Mesio context: mesh, grids, commit (collective)
        bool createContext();

        //- Topology change: rebuild the Mesio mesh in the context (collective)
        bool resetContextMesh();

        //- Apply the run-time options (writer, compression, sampling,
        //- memory tracking) to the context
        bool applyOptions();

        //- Destroy the Mesio context (collective)
        void destroyContext();

        //- Output file name for the current time
        fileName outputFile() const;

        //- No copy construct
        mesioVDB(const mesioVDB&) = delete;

        //- No copy assignment
        void operator=(const mesioVDB&) = delete;


public:

    //- Runtime type information
    TypeName("mesioVDB");


    // Constructors

        //- Construct from name, Time and dictionary
        mesioVDB
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor. Releases the Mesio context (collective).
    virtual ~mesioVDB();


    // Member Functions

        //- Read the function-object dictionary. At a run-time re-read,
        //- invalid settings are reported and the previous ones kept;
        //- changed grid settings re-create the context at the next export.
        virtual bool read(const dictionary& dict);

        //- Record the event marker used for the freshness check
        virtual bool execute();

        //- Export one .vdb file (collective)
        virtual bool write();

        //- Report the accumulated export time
        virtual bool end();

        //- Topology change: re-create the Mesio context at the next export
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Mesh motion: update the Mesio point coordinates at the next export
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mesioVDBTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
