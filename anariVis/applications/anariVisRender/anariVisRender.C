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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IStringStream.H"
#include "anariVisSettingsIO.H"
#include "anariVisMPI.H"
#include "anariVisTransport.H"
#include "anariVisRenderer.H"
#include "anariVisScene.H"

#include <chrono>

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

    ::anariVis::widenAffinity();

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
