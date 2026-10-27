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

\*---------------------------------------------------------------------------*/

#include "anariVisRenderViewer.H"
#include "anariVisScene.H"

#include <chrono>
#include <climits>
#include <cstdio>

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

long anariVis::runViewer
(
    Transport& transport,
    Renderer& renderer,
    const RenderSettings& settings,
    const ViewerOptions& opt
)
{
    typedef std::chrono::steady_clock clock_type;

    const MPI_Comm comm = transport.renderComm();
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    LiveViewer viewer(renderer, settings, comm, opt);
    if (!viewer.ok())
    {
        return -1;
    }

    std::vector<std::vector<char>> packets;
    bool producersDone = false;
    while (!viewer.quit())
    {
        // New solver frames (non-blocking, agreed by all render ranks)
        if (!producersDone)
        {
            int idx = 0;
            const bool ready = transport.pollFrame(idx);
            // all ranks ready (a rank that is not gives INT_MIN) with the
            // same frame (min index == max index)
            int v[2] = {ready ? idx : INT_MAX, ready ? -idx : INT_MIN};
            MPI_Allreduce(MPI_IN_PLACE, v, 2, MPI_INT, MPI_MIN, comm);
            if (v[1] != INT_MIN && v[0] == -v[1])
            {
                const auto t0 = clock_type::now();
                if (!transport.takeFrame(v[0], packets))
                {
                    // cannot happen: every rank reported this frame complete
                    std::fprintf
                    (
                        stderr, "anariVisRender: render rank %d lost frame %d\n",
                        rank, v[0]
                    );
                }
                LocalScene scene;
                assembleScene
                (
                    packets, transport.consumerIndex(),
                    transport.nConsumers(), scene
                );
                packets.clear();
                viewer.addFrame
                (
                    std::move(scene),
                    std::chrono::duration<double>(clock_type::now() - t0).count()
                );
            }
            int fin = transport.finished() ? 1 : 0;
            MPI_Allreduce(MPI_IN_PLACE, &fin, 1, MPI_INT, MPI_MIN, comm);
            if (fin)
            {
                producersDone = true;
                viewer.setSolverFinished();
            }
        }

        // Client commands, the next view, one render pass
        viewer.step(producersDone ? 50 : 10);
    }

    if (rank == 0)
    {
        std::printf
        (
            "anariVisRender: viewer stopped (%ld frames received)\n",
            viewer.received()
        );
        std::fflush(stdout);
    }
    return viewer.received();
}


// ************************************************************************* //
