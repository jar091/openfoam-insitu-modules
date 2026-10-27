/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 IT4Innovations, VSB-TU Ostrava
-------------------------------------------------------------------------------
License
    This file is part of the anariVis module for OpenFOAM.

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
    Foam::functionObjects::anariVis

Group
    grpUtilitiesFunctionObjects

Description
    In-situ and in-transit visualisation through ANARI (Barney CPU, helide,
    cyclesphi-anari, deepCompositing) without intermediate files.

    At every execution the function object
    - updates its sampled surfaces (any sampledSurface type: isoSurface,
      patch, cuttingPlane, ...), triangulates them and interpolates the
      colour field onto the surface points,
    - resamples volume fields of the local subdomain onto the points of a
      global structured grid that is split into one brick per render rank,
    - packs the local data into frame packets and sends them (MPI, no
      files) to the render ranks.

    Render ranks build an ANARI world from the data they receive and render
    collectively. The image is written as PNG by render rank 0.

    Modes:
    - \c inSitu    : the render ranks are a subset of the solver ranks
                     (\c renderRanksPerNode per node, -1 for all ranks).
                     Rendering is synchronous; with \c synchronise the
                     other ranks wait in a sleeping barrier and leave their
                     cores to the render threads.
    - \c inTransit : the render ranks belong to a separate OpenFOAM world
                     running the anariVisRender application (MPMD launch,
                     option -world). The solver only waits for the
                     non-blocking sends of the previous frame; if the
                     renderer is still busy the frame is dropped
                     (\c dropPolicy drop) or the solver waits
                     (\c dropPolicy block).

    Rendering libraries (\c library): barney_cpu (Barney CPU/Embree,
    data-parallel "mpi" device), barney, helide, cycles, ptc (deep
    compositing wrapper of \c wrappedLibrary). Compositing (\c compositing):
    auto, none, barney, deep, ptc.

    The Barney "mpi" device expects MPI_THREAD_MULTIPLE: start the solver
    (in-situ) or the renderer (in-transit) with option \c -mpi-threads.

Usage
    Example of function object specification:
    \verbatim
    anariVis1
    {
        type            anariVis;
        libs            (anariVis);

        executeControl  timeStep;
        executeInterval 10;
        writeControl    none;

        mode            inSitu;        // inSitu | inTransit
        renderRanksPerNode 1;          // inSitu: render ranks per node
        synchronise     true;          // inSitu: sleeping wait
        dropPolicy      drop;          // inTransit: drop | block

        library         barney_cpu;    // barney_cpu | helide | cycles | ptc
        compositing     auto;          // auto | none | barney | deep | ptc
        resolution      (1280 720);
        pixelSamples    8;
        passes          1;
        background      (0.1 0.1 0.12);
        outputDir       "<case>/postProcessing/anariVis1";

        camera
        {
            viewDirection (-1 -0.6 -1.2);
            up            (0 1 0);
            fov           35;
            zoom          1.1;
        }

        surfaces
        {
            walls
            {
                type        patch;
                patches     (fixedWalls);
                colour      (0.7 0.7 0.7);
                opacity     1;
            }
            qIso
            {
                type        isoSurface;
                isoField    Q;
                isoValue    0.5;
                isoMethod   topo;
                colourField U;
                colourMap   viridis;
            }
        }

        volumes
        {
            pressure
            {
                field       p;
                resolution  128;
                colourMap   coolToWarm;
                opacity     ((0 0) (0.5 0.05) (1 0.6));
            }
        }
    }
    \endverbatim

    Main entries:
    \table
        Property     | Description                          | Required | Default
        type         | Type name: anariVis                  | yes |
        mode         | inSitu or inTransit                  | no  | inSitu
        renderRanksPerNode | inSitu render ranks per node (-1: all) | no | 1
        synchronise  | inSitu: sleeping barrier after render | no | true
        renderAffinityAll | inSitu: render ranks may use all node cores | no | true
        dropPolicy   | inTransit: drop or block             | no  | drop
        interpolationScheme | point interpolation of colour fields | no | cellPoint
        surfaceRouting | node (to a render rank on the same node / block mapping) or spatial (triangles to the convex brick of a render rank, exact sort-last ordering) | no | auto (spatial for deep/ptc)
        surfaces     | sampledSurface dictionaries + style  | no  |
        volumes      | volume dictionaries                  | no  |
    \endtable

    Render settings: see Foam::anariVisIO.

Note
    The mesh is assumed static (volume sampling weights and routing are
    computed once).

SourceFiles
    anariVis.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_anariVis_H
#define Foam_functionObjects_anariVis_H

#include "fvMeshFunctionObject.H"
#include "sampledSurface.H"
#include "PtrList.H"
#include "anariVisVolumeSampler.H"
#include "anariVisTypes.H"

#include <memory>

// Forward declarations
namespace anariVis
{
    class Transport;
    class Renderer;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class anariVis Declaration
\*---------------------------------------------------------------------------*/

class anariVis
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Copy of the dictionary (sent to the in-transit renderer)
        dictionary dict_;

        //- Render settings
        ::anariVis::RenderSettings settings_;

        //- inSitu | inTransit
        word mode_;

        //- inSitu: render ranks per node (-1: all)
        label renderRanksPerNode_;

        //- inSitu: sleeping barrier after rendering
        bool synchronise_;

        //- inSitu: let render ranks use all cores of the node
        bool widenAffinity_;

        //- inTransit: drop | block
        word dropPolicy_;

        //- Interpolation scheme for the surface colour fields
        word interpolationScheme_;

        //- Surface routing: auto | node | spatial
        word surfaceRouting_;

        //- Surfaces are routed by triangle centroid to convex bricks
        bool spatialRouting_;

        //- Point grid and brick decomposition of the spatial routing
        ::anariVis::VolumeGrid sceneGrid_;
        ::anariVis::BrickDecomposition sceneDecomp_;

        //- Sampled surfaces
        PtrList<sampledSurface> surfaces_;

        //- Colour field per surface (empty: uniform)
        wordList colourFields_;

        //- Volume samplers
        PtrList<anariVisVolumeSampler> volumes_;

        //- Transport between producers and render ranks
        std::unique_ptr<::anariVis::Transport> transport_;

        //- Renderer (inSitu render ranks only)
        std::unique_ptr<::anariVis::Renderer> renderer_;

        //- Consumers this producer sends to every frame
        std::vector<int> destinations_;

        //- Global mesh bounds
        boundBox bounds_;

        //- Global data ranges of the last frame (surfaces, then volumes)
        std::vector<::anariVis::vec2f> dataRanges_;

        bool initialised_;

        bool ended_;

        //- MPI was initialised by this object (serial runs)
        bool ownMpi_;

        // Statistics
        label nFrames_;
        label nDropped_;
        scalar sumWait_;
        scalar sumExtract_;
        scalar sumSend_;
        scalar sumRender_;
        scalar sumTotal_;


    // Private Member Functions

        //- Set up transport, renderer and volume sampling. Collective.
        void initialise();

        //- Extract the local surface pieces. Collective.
        void extractSurfaces(std::vector<::anariVis::SurfacePiece>& pieces);

        //- Send the stop packets (inTransit)
        void sendStop();

        //- No copy construct
        anariVis(const anariVis&) = delete;

        //- No copy assignment
        void operator=(const anariVis&) = delete;


public:

    //- Runtime type information
    TypeName("anariVis");


    // Constructors

        //- Construct from Time and dictionary
        anariVis
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~anariVis();


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary& dict);

        //- Extract, send and (inSitu) render
        virtual bool execute();

        //- No-op
        virtual bool write();

        //- Send the stop packets to the in-transit renderer
        virtual bool end();

        //- Mesh changes are not supported (warning)
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Mesh motion is not supported (warning)
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
