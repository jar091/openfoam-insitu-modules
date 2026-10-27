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
    Assembly of the local scene of a render (consumer) rank from the frame
    packets it received: surface fragments of the same surface id are
    merged into one triangle mesh, volume samples are scattered into this
    consumer's brick of the global grid.

SourceFiles
    anariVisScene.C

\*---------------------------------------------------------------------------*/

#ifndef anariVis_Scene_H
#define anariVis_Scene_H

#include "anariVisRenderer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace anariVis
{

//- Assemble the scene of consumer consumerIdx (of nConsumers)
void assembleScene
(
    const std::vector<std::vector<char>>& packets,
    int consumerIdx,
    int nConsumers,
    LocalScene& scene
);

} // End namespace anariVis

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
