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
    Assembly of the local scene of a render rank from received packets.

\*---------------------------------------------------------------------------*/

#include "anariVisScene.H"

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void anariVis::assembleScene
(
    const std::vector<std::vector<char>>& packets,
    int consumerIdx,
    int nConsumers,
    LocalScene& scene
)
{
    scene.surfaces.clear();
    scene.bricks.clear();

    std::map<int, SurfacePiece> surfaces;
    std::map<int, VolumeBrick> bricks;
    std::map<int, BrickDecomposition> decomp;

    bool first = true;
    FramePacket pkt;

    for (const auto& buf : packets)
    {
        if (!pkt.unpack(buf.data(), buf.size()))
        {
            throw std::runtime_error("anariVis: corrupt frame packet");
        }
        if (first)
        {
            scene.frameIndex = pkt.header.frameIndex;
            scene.time = pkt.header.time;
            scene.bounds = pkt.header.bounds;
            scene.hasBrick = false;
            first = false;
        }
        if (pkt.header.hasBrick)
        {
            scene.hasBrick = true;
            scene.brick = pkt.header.brick;
        }

        for (SurfacePiece& s : pkt.surfaces)
        {
            if (s.triangles.empty())
            {
                // Still record the surface (range) but without geometry
                surfaces[s.id].id = s.id;
                surfaces[s.id].range = s.range;
                continue;
            }
            SurfacePiece& m = surfaces[s.id];
            m.id = s.id;
            m.range = s.range;

            const uint32_t offset = uint32_t(m.nPoints());
            const bool scalars =
                (s.scalars.size() == s.nPoints())
             && (m.points.empty() || m.scalars.size() == m.nPoints());

            m.points.insert(m.points.end(), s.points.begin(), s.points.end());
            if (scalars)
            {
                m.scalars.insert
                (
                    m.scalars.end(), s.scalars.begin(), s.scalars.end()
                );
            }
            else
            {
                m.scalars.clear();
            }
            m.triangles.reserve(m.triangles.size() + s.triangles.size());
            for (const uint32_t i : s.triangles)
            {
                m.triangles.push_back(i + offset);
            }
        }

        for (const VolumePiece& v : pkt.volumes)
        {
            auto it = bricks.find(v.id);
            if (it == bricks.end())
            {
                VolumeBrick vb;
                vb.id = v.id;
                vb.grid = v.grid;
                vb.range = v.range;

                auto d = decomp.emplace
                (
                    v.id, BrickDecomposition(v.grid.dims, nConsumers)
                ).first;
                d->second.range(consumerIdx, vb.lo, vb.hi);

                const auto n = vb.size();
                vb.data.assign
                (
                    size_t(n[0])*n[1]*n[2],
                    std::numeric_limits<float>::quiet_NaN()
                );
                it = bricks.emplace(v.id, std::move(vb)).first;
            }

            VolumeBrick& vb = it->second;
            const auto n = vb.size();
            const int64_t nx = v.grid.dims[0];
            const int64_t ny = v.grid.dims[1];

            for (size_t s = 0; s < v.index.size(); ++s)
            {
                const int64_t g = v.index[s];
                const int64_t i = g % nx - vb.lo[0];
                const int64_t j = (g/nx) % ny - vb.lo[1];
                const int64_t k = g/(nx*ny) - vb.lo[2];
                if
                (
                    i < 0 || j < 0 || k < 0
                 || i >= n[0] || j >= n[1] || k >= n[2]
                )
                {
                    continue;
                }
                vb.data[size_t(i + n[0]*(j + n[1]*k))] = v.value[s];
            }
        }
    }

    for (auto& s : surfaces)
    {
        scene.surfaces.push_back(std::move(s.second));
    }

    for (auto& b : bricks)
    {
        // Points not covered by any cell: lowest value of the range
        // (transparent with the usual opacity maps)
        const float fill = b.second.range[0];
        for (float& x : b.second.data)
        {
            if (std::isnan(x))
            {
                x = fill;
            }
        }
        scene.bricks.push_back(std::move(b.second));
    }
}


// ************************************************************************* //
