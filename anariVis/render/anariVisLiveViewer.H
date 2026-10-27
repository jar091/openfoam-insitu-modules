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
    anariVis::LiveViewer

Description
    Interactive viewing of the frames with the Nexaris client (nexaris_live,
    see the Nexaris docs/INSITU.md). Used by the in-transit renderer
    (anariVisRender -viewer) and by the function object in-situ (viewer
    dictionary).

    Render rank 0 hosts the Nexaris server (session, token, JPEG encoder,
    video stream). Every step(), on all render ranks:

    \verbatim
    1. rank 0: client commands (layers, follow / show frame, pause) and the
       next view the Nexaris server asks for -> MPI_Bcast of a fixed command
       block (+ layer settings when they changed)
    2. all ranks: apply, render one accumulation pass (collective);
       rank 0 hands colour + ray depth to the server
    3. rank 0: publish the state to the client
    \endverbatim

    Frames are added with addFrame() (collective): a bounded in-memory
    history, no files. Only the calling (main) thread uses MPI and ANARI; the
    Nexaris threads talk to it through nexaris::live::Service.

    In-situ the solver waits while runInSituFrame() renders the new frame
    for the client: passesPerFrame passes or maxSecondsPerFrame, nothing
    without a client, and as long as the client pauses the simulation.

    Compiled only with ANARIVIS_NEXARIS (NEXARIS_LIVE_ROOT at build time).

SourceFiles
    anariVisLiveViewer.C

\*---------------------------------------------------------------------------*/

#ifndef anariVis_LiveViewer_H
#define anariVis_LiveViewer_H

#ifdef ANARIVIS_NEXARIS

#include "anariVisRenderer.H"

#include <chrono>
#include <deque>
#include <memory>
#include <string>
#include <vector>

// Forward declarations
namespace nexaris
{
namespace live
{
    class Service;
    struct Layer;
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace anariVis
{

struct ViewerOptions
{
    //- TCP port on 127.0.0.1 of render rank 0 (client via SSH tunnel)
    int port = 47800;

    //- Frames kept in memory for browsing back in time
    int historySize = 16;

    //- Memory budget of the history per render rank [MB]
    double historyMB = 16384;

    //- After the solver finished: keep serving until no client was
    //  connected for this long [s]
    double keepAlive = 600;

    std::string title;

    int jpegQuality = 85;

    float maxFps = 30;

    //- Default samples per pixel per pass (the client may override)
    int pixelSamples = 1;

    //- Follow mode (in-transit): passes of the shown frame before a newer
    //  one replaces it
    int minPasses = 8;

    //- Up axis of the case ('z' or 'y'; the client's world is +Y up)
    char upAxis = 'z';

    //- In-situ: the solver waits for the viewer (runInSituFrame)
    bool inSitu = false;

    //- In-situ: passes rendered for the client per new frame
    int passesPerFrame = 2;

    //- In-situ: upper bound of the time the solver waits per frame [s]
    double maxSecondsPerFrame = 2;
};


/*---------------------------------------------------------------------------*\
                         Class LiveViewer Declaration
\*---------------------------------------------------------------------------*/

class LiveViewer
{
public:

    //- What a step did (identical on all render ranks)
    struct StepResult
    {
        bool quit = false;
        bool rendered = false;
        bool paused = false;
        //- In-situ: the frame loop continues
        bool more = false;
    };


private:

    //- One frame (this rank's part) with the global field ranges
    struct Entry
    {
        LocalScene scene;
        std::vector<vec2f> ranges;   // surfaces, then volumes
        double bytes = 0;
    };

    typedef std::chrono::steady_clock clock_type;

    Renderer& renderer_;
    RenderSettings settings_;
    ViewerOptions opt_;
    MPI_Comm comm_;
    int rank_;
    int size_;
    size_t nSurfaces_;
    size_t nVolumes_;
    bool ok_;

    //- Rank 0 only
    std::unique_ptr<nexaris::live::Service> svc_;
    std::vector<nexaris::live::Layer> layers_;
    std::vector<nexaris::live::Layer> defaults_;
    std::vector<std::string> mapNames_;

    std::deque<std::unique_ptr<Entry>> history_;
    int shown_;
    bool follow_;
    bool paused_;
    bool solverFinished_;
    bool quit_;
    long received_;

    //- Passes of the shown frame / of the current in-situ frame
    long passes_;
    long framePasses_;
    bool inFrame_;
    clock_type::time_point frameStart_;
    clock_type::time_point lastRender_;
    clock_type::time_point lastClient_;

    std::vector<uint32_t> rgba_;
    std::vector<float> depth_;
    std::vector<float> stateBuf_;
    FrameTimings last_;
    double lastAssemble_;
    double lastFrameTotal_;
    double historyMaxMB_;


    // Private Member Functions

        void showEntry(int i);

        //- Rank 0: publish the state and the telemetry
        void publish();


public:

    // Constructors

        //- Start the viewer (Nexaris server on render rank 0). Collective
        //  over renderComm; ok() tells whether it started on rank 0.
        LiveViewer
        (
            Renderer& renderer,
            const RenderSettings& settings,
            MPI_Comm renderComm,
            const ViewerOptions& options
        );

        LiveViewer(const LiveViewer&) = delete;
        void operator=(const LiveViewer&) = delete;


    //- Destructor
    ~LiveViewer();


    // Member Functions

        bool ok() const { return ok_; }

        //- The viewer stopped (server stopped, or keepAlive expired)
        bool quit() const { return quit_; }

        long received() const { return received_; }

        //- Add a frame (this rank's part of it). Collective.
        void addFrame(LocalScene&& scene, double assembleSeconds = 0);

        //- The solver ended (no more frames). Collective.
        void setSolverFinished();

        //- One iteration: commands, the next view (rank 0 waits up to
        //  waitMs for it), one render pass. Collective.
        StepResult step(int waitMs);

        //- In-situ: render the newest frame for the client, while the
        //  solver waits. Collective.
        void runInSituFrame();

        //- Serve the history until quit (after the solver ended).
        //  Collective.
        void serveUntilQuit();
};


} // End namespace anariVis

#endif // ANARIVIS_NEXARIS

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
