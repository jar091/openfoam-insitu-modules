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

Application
    anariVisRender

Group
    grpPostProcessingUtilities

Description
    In-transit ANARI renderer for the anariVis function object.

    Runs as a separate OpenFOAM world in an MPMD launch together with the
    solver (option -world on both sides). The render ranks receive the
    render settings and, every frame, the extracted surfaces and volume
    samples from the solver ranks (MPI, no files), render them collectively
    (Barney "mpi" data-parallel, or per rank + deepCompositing) and render
    rank 0 writes the image.

Usage
    Slurm heterogeneous step (one MPI_COMM_WORLD):
    \verbatim
    srun --mpi=pmix -N1 -n128 -w cn627 icoFoam -parallel -world sim \
       : -N2 -n2 -w cn651,cn652 anariVisRender -parallel -world render \
         -mpi-threads
    \endverbatim

    Options:
      - \par -threads \<N\>
        Render threads per render rank (default: cores / render ranks on
        the node)
      - \par -viewer
        Interactive viewing with the Nexaris client instead of PNG files
        (built with NEXARIS_LIVE_ROOT): render rank 0 serves the frames on
        127.0.0.1:\<port\> (client through an SSH tunnel); the frames
        stay in memory, nothing is written. Further options -viewerPort,
        -history, -historyMB, -keepAlive, -upAxis, -viewerSpp

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IStringStream.H"
#include "anariVisSettingsIO.H"
#include "anariVisMPI.H"
#include "anariVisTransport.H"
#include "anariVisRenderer.H"
#include "anariVisScene.H"
#ifdef ANARIVIS_NEXARIS
#include "anariVisRenderViewer.H"
#endif

#include <chrono>
#include <cstdlib>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "In-transit ANARI renderer for the anariVis function object."
        " Run as a separate world (-world) of an MPMD launch."
    );
    argList::noCheckProcessorDirectories();
    argList::noFunctionObjects();
    argList::addOption
    (
        "threads",
        "N",
        "Render threads per render rank"
    );

    #ifdef ANARIVIS_NEXARIS
    argList::addBoolOption
    (
        "viewer",
        "Interactive viewing with the Nexaris client (no image files)"
    );
    argList::addOption("viewerPort", "port", "Viewer TCP port (47800)");
    argList::addOption("history", "N", "Frames kept in memory (16)");
    argList::addOption("historyMB", "MB", "History memory per rank (16384)");
    argList::addOption
    (
        "keepAlive", "s",
        "Viewer lifetime after the solver ended without a client (600)"
    );
    argList::addOption("upAxis", "y|z", "Up axis of the case (z)");
    argList::addOption("viewerSpp", "N", "Default samples per pixel (1)");
    argList::addOption
    (
        "viewerMinPasses", "N",
        "Passes of a frame before a newer one replaces it when following (8)"
    );
    argList::addOption("viewerTitle", "text", "Title shown by the client");
    #endif

    argList args(argc, argv);

    if (!UPstream::parRun() || UPstream::allWorlds().size() < 2)
    {
        FatalErrorInFunction
            << "anariVisRender must run with -parallel in its own world"
            << " (-world <name>) next to the solver world" << nl
            << exit(FatalError);
    }

    typedef std::chrono::steady_clock clock_type;
    auto seconds = [](const clock_type::time_point& t0)
    {
        return std::chrono::duration<double>(clock_type::now() - t0).count();
    };

    const MPI_Comm global = ::anariVis::toMPIComm
    (
        UPstream::Communicator::lookup(UPstream::commGlobal())
    );

    ::anariVis::Transport transport(global, false, true);

    // Settings from the simulation (the function object dictionary)
    const std::string text = transport.receiveConfig();
    IStringStream is(text);
    const dictionary dict(is);

    ::anariVis::RenderSettings settings =
        anariVisIO::readSettings(dict, ".", "anariVis");

    // all cores of the node, unless the launcher placed the rank next to its
    // GPU (ANARIVIS_KEEP_AFFINITY)
    if (!std::getenv("ANARIVIS_KEEP_AFFINITY"))
    {
        ::anariVis::widenAffinity();
    }

    if (args.found("threads"))
    {
        settings.threads = args.get<int>("threads");
    }
    if (settings.threads <= 0)
    {
        settings.threads = std::max
        (
            1,
            ::anariVis::availableCores()/transport.nConsumersOnNode()
        );
    }

    transport.setDestinations(std::vector<int>());

    #ifdef ANARIVIS_NEXARIS
    const bool viewer = args.found("viewer");
    if (viewer)
    {
        // in-transit viewing: frames stay in memory
        settings.writeImages = false;
    }
    #endif

    ::anariVis::Renderer renderer(settings, transport.renderComm());

    {
        int provided = 0;
        MPI_Query_thread(&provided);
        Info<< "anariVisRender: render ranks " << transport.nConsumers()
            << ", producers " << transport.nProducers()
            << ", renderer " << renderer.libraryName().c_str()
            << '/' << renderer.mode().c_str()
            << ", threads/rank " << settings.threads
            << ", MPI thread level " << provided
            << ", image " << settings.width << 'x' << settings.height
            << ", output " << settings.outputDir.c_str() << nl << endl;
    }

    #ifdef ANARIVIS_NEXARIS
    if (viewer)
    {
        ::anariVis::ViewerOptions vo;
        vo.port = args.getOrDefault<label>("viewerPort", 47800);
        vo.historySize = args.getOrDefault<label>("history", 16);
        vo.historyMB = args.getOrDefault<scalar>("historyMB", 16384);
        vo.keepAlive = args.getOrDefault<scalar>("keepAlive", 600);
        vo.pixelSamples = args.getOrDefault<label>("viewerSpp", 1);
        vo.minPasses = args.getOrDefault<label>("viewerMinPasses", 8);
        vo.upAxis = args.getOrDefault<word>("upAxis", "z") == "y" ? 'y' : 'z';
        std::string title = args.globalCaseName();
        title +=
            " (in-transit, " + std::to_string(transport.nProducers())
          + " solver ranks)";
        vo.title =
            args.found("viewerTitle")
          ? std::string(args.get<string>("viewerTitle")) : title;
        const long n = ::anariVis::runViewer(transport, renderer, settings, vo);
        Info<< "anariVisRender: viewer finished, " << label(n)
            << " frames received" << nl << "End" << nl << endl;
        return n < 0 ? 1 : 0;
    }
    #endif

    std::vector<std::vector<char>> packets;
    int frameIndex = 0;
    label nFrames = 0;
    double sumWait = 0, sumRender = 0, sumTotal = 0;
    auto tWait = clock_type::now();

    while (transport.receiveFrame(packets, frameIndex))
    {
        const double wait = seconds(tWait);
        const auto t0 = clock_type::now();

        size_t nBytes = 0;
        for (const auto& p : packets)
        {
            nBytes += p.size();
        }

        ::anariVis::LocalScene scene;
        ::anariVis::assembleScene
        (
            packets,
            transport.consumerIndex(),
            transport.nConsumers(),
            scene
        );
        const double tAssemble = seconds(t0);

        const ::anariVis::FrameTimings ft = renderer.render(scene);
        const double total = seconds(t0);

        // Load balance information
        double maxBytes = double(nBytes);
        MPI_Allreduce
        (
            MPI_IN_PLACE, &maxBytes, 1, MPI_DOUBLE, MPI_MAX,
            transport.renderComm()
        );

        ++nFrames;
        sumWait += wait;
        sumRender += ft.render;
        sumTotal += total;

        Info<< "anariVisRender: frame " << scene.frameIndex
            << " t=" << scene.time
            << " packets " << label(packets.size())
            << " max bytes/rank " << maxBytes
            << " wait " << wait << " s, assemble " << tAssemble
            << " s, build " << ft.build << " s, render " << ft.render
            << " s, composite " << ft.composite << " s, write " << ft.write
            << " s, total " << total << " s" << endl;

        tWait = clock_type::now();
    }

    if (nFrames)
    {
        Info<< nl << "anariVisRender: " << nFrames << " frames, mean wait "
            << sumWait/nFrames << " s, render " << sumRender/nFrames
            << " s, total " << sumTotal/nFrames << " s" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
