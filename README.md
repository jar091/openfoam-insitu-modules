# openfoam-insitu-modules

OpenFOAM v2606 modules for in-situ and in-transit visualisation:

| Module | Contents |
|---|---|
| `mesioFoam` | Function object `mesioVDB`: exports OpenFOAM volume fields (native or derived by other function objects) as OpenVDB fog volumes through the Mesio in-situ API, one `.vdb` file per export step for all MPI ranks, without intermediate field files. |
| `anariVis` | Function object `anariVis` and application `anariVisRender`: in-situ and in-transit rendering of OpenFOAM data with ANARI devices (`barney_cpu`, `cycles` from cyclesphi-anari, `helide`, `ptc`) and optional deep compositing. |

The OpenFOAM core is not modified; both modules are built with `wmake` against an existing
OpenFOAM installation. Build and run scripts, environment setup, examples and detailed notes
live in the companion repository
[openfoam-insitu-workflow](https://github.com/jar091/openfoam-insitu-workflow).

## Requirements

- OpenFOAM v2606, `WM_LABEL_SIZE=32` (tested with foss/2025b: GCC 14.3, OpenMPI 5.0.8)
- `mesioFoam`:
  - Mesio from the espreso fork, branch `mjar/devel`, built with `--intwidth=32`;
    in-situ API version 2 (`include/mesio_insitu.h`; older headers are rejected at compile time)
  - OpenVDB 12
- `anariVis`:
  - ANARI-SDK (jar091 fork)
  - deepCompositing CPU build (jar091 fork)
  - an ANARI device: barney (CPU backend), cyclesphi-anari or helide

## Build

Source the OpenFOAM environment first, then set the variables used by `Make/options`.

```bash
# mesioFoam -> $FOAM_USER_LIBBIN/libmesioFoam.so
export MESIO_ROOT=<mesio install prefix>      # include/mesio_insitu.h, lib64/libmesioapi.so
export OPENVDB_ROOT=<openvdb install prefix>
(cd mesioFoam && ./Allwmake)

# anariVis -> $ANARIVIS_LIBBIN/libanariVis.so, $ANARIVIS_BINDIR/anariVisRender
export ANARI_ROOT=<ANARI-SDK prefix>
export DEEPCOMP_ROOT=<deepCompositing prefix>
export ANARIVIS_LIBBIN=<library output dir>
export ANARIVIS_BINDIR=<application output dir>
(cd anariVis && ./Allwmake)
```

The workflow repository provides `scripts/env_mesio_openfoam.sh`,
`scripts/env_anari_openfoam.sh`, `scripts/build_mesio_openfoam.sh` and
`scripts/build_anariVis.sh`, which set all of these for the Karolina setup.

## Usage

`system/controlDict` excerpt for `mesioVDB`. Derived fields are computed in memory by
standard function objects listed before the exporter (`writeControl none`):

```
functions
{
    Q1
    {
        type            Q;
        libs            (fieldFunctionObjects);
        field           U;
        executeControl  timeStep;
        executeInterval 50;
        writeControl    none;
    }

    mesioVDB1
    {
        type            mesioVDB;
        libs            (mesioFoam);
        writeControl    timeStep;
        writeInterval   50;
        fields          (U p Q);          // vectors are exported as U_mag (vectorComponents)
        density         256;              // or voxelSize
        // box          (0 0 0) (0.1 0.1 0.05);
        outputDir       "VDB";
        writer          distributed;      // distributed | parallel | single
        outputPolicy    vdbOnly;
    }
}
```

Full examples, including all `mesioVDB` options (`sampling`, `onError`, `compression`,
several instances) and the `anariVis` in-situ / in-transit setups, are in the workflow
repository under `data/openfoam/` and documented in
`docs/notes/mesio_openfoam_functionobject.md` and
`docs/notes/anari_openfoam_insitu_intransit.md`.

## Development

AI agent harnesses were used for selected development, testing, and optimization tasks

## License

GNU General Public License v3 or later (see `COPYING`), as required for code linked
against OpenFOAM.

Third-party code: `anariVis/render/stb_image_write.h` by Sean Barrett (public domain / MIT).
