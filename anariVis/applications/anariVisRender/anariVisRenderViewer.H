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

Description
    Interactive viewing of the in-transit frames with the Nexaris client
    (nexaris_live, see the Nexaris docs/INSITU.md).

    Render rank 0 hosts the Nexaris server (session, token, JPEG encoder,
    video stream). Every loop iteration, on all render ranks:

    \verbatim
    1. receive solver packets without blocking (Transport::pollFrame);
       a frame that is complete on every rank (allreduce) is assembled
       and kept in a bounded in-memory history (no files)
    2. rank 0: client commands (layers, follow / show frame) and the next
       view the Nexaris server asks for -> MPI_Bcast of a fixed command
       block (+ layer settings when they changed)
    3. all ranks: apply, render one accumulation pass (collective);
       rank 0 hands colour + ray depth to the server
    \endverbatim

    Only the main thread calls MPI and ANARI; the Nexaris threads talk to
    it through nexaris::live::Service.

SourceFiles
    anariVisRenderViewer.C

\*---------------------------------------------------------------------------*/

#ifndef anariVis_RenderViewer_H
#define anariVis_RenderViewer_H

#include "anariVisLiveViewer.H"
#include "anariVisTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace anariVis
{

//- Run the interactive loop until the solver finished and no client is
//  connected (keepAlive) or the server stopped. Collective over the render
//  ranks. Returns the number of frames received.
long runViewer
(
    Transport& transport,
    Renderer& renderer,
    const RenderSettings& settings,
    const ViewerOptions& options
);

} // End namespace anariVis

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
