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
    Serialisation of frame packets, brick decomposition and colour maps.

    Wire format (native endianness, all ranks run on the same platform):

    \verbatim
    FrameHeader (fixed size, see anariVisTypes.H)
    int32 nSurfaces
    int32 nVolumes
    per surface:
        int32 id, float range[2],
        uint64 nPoints, uint64 nScalars, uint64 nTriangles,
        float points[3*nPoints], float scalars[nScalars],
        uint32 triangles[3*nTriangles]
    per volume:
        int32 id, int32 dims[3], float origin[3], float spacing[3],
        float range[2], uint64 nSamples,
        uint32 index[nSamples], float value[nSamples]
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "anariVisTypes.H"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

template<class T>
void put(std::vector<char>& buf, const T& v)
{
    const size_t n = buf.size();
    buf.resize(n + sizeof(T));
    std::memcpy(buf.data() + n, &v, sizeof(T));
}

template<class T>
void putArray(std::vector<char>& buf, const std::vector<T>& v)
{
    const size_t n = buf.size();
    const size_t nb = v.size()*sizeof(T);
    buf.resize(n + nb);
    if (nb)
    {
        std::memcpy(buf.data() + n, v.data(), nb);
    }
}


class Reader
{
    const char* data_;
    size_t size_;
    size_t pos_;
    bool ok_;

public:

    Reader(const char* data, size_t size)
    :
        data_(data), size_(size), pos_(0), ok_(true)
    {}

    bool ok() const { return ok_; }

    template<class T>
    T get()
    {
        T v{};
        if (pos_ + sizeof(T) > size_)
        {
            ok_ = false;
            return v;
        }
        std::memcpy(&v, data_ + pos_, sizeof(T));
        pos_ += sizeof(T);
        return v;
    }

    template<class T>
    void getArray(std::vector<T>& v, uint64_t n)
    {
        const size_t nb = n*sizeof(T);
        if (!ok_ || pos_ + nb > size_)
        {
            ok_ = false;
            v.clear();
            return;
        }
        v.resize(n);
        if (nb)
        {
            std::memcpy(v.data(), data_ + pos_, nb);
        }
        pos_ += nb;
    }
};


//- Prime factors in descending order
std::vector<int> primeFactors(int n)
{
    std::vector<int> f;
    for (int p = 2; p*p <= n; ++p)
    {
        while (n % p == 0)
        {
            f.push_back(p);
            n /= p;
        }
    }
    if (n > 1)
    {
        f.push_back(n);
    }
    std::sort(f.rbegin(), f.rend());
    return f;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * FramePacket  * * * * * * * * * * * * * * //

void anariVis::FramePacket::pack(std::vector<char>& buf) const
{
    buf.clear();

    size_t nBytes = sizeof(FrameHeader) + 64;
    for (const auto& s : surfaces)
    {
        nBytes += 64 + 4*(s.points.size() + s.scalars.size()
            + s.triangles.size());
    }
    for (const auto& v : volumes)
    {
        nBytes += 64 + 8*v.index.size();
    }
    buf.reserve(nBytes);

    put(buf, header);
    put(buf, int32_t(surfaces.size()));
    put(buf, int32_t(volumes.size()));

    for (const auto& s : surfaces)
    {
        put(buf, s.id);
        put(buf, s.range);
        put(buf, uint64_t(s.nPoints()));
        put(buf, uint64_t(s.scalars.size()));
        put(buf, uint64_t(s.nTriangles()));
        putArray(buf, s.points);
        putArray(buf, s.scalars);
        putArray(buf, s.triangles);
    }

    for (const auto& v : volumes)
    {
        put(buf, v.id);
        put(buf, v.grid.dims);
        put(buf, v.grid.origin);
        put(buf, v.grid.spacing);
        put(buf, v.range);
        put(buf, uint64_t(v.index.size()));
        putArray(buf, v.index);
        putArray(buf, v.value);
    }
}


bool anariVis::FramePacket::peekHeader
(
    const char* data,
    size_t nBytes,
    FrameHeader& h
)
{
    Reader r(data, nBytes);
    h = r.get<FrameHeader>();
    return
        r.ok()
     && h.magic == FrameHeader::magicValue
     && h.version == FrameHeader::versionValue;
}


bool anariVis::FramePacket::unpack(const char* data, size_t nBytes)
{
    surfaces.clear();
    volumes.clear();

    Reader r(data, nBytes);
    header = r.get<FrameHeader>();
    if
    (
        !r.ok()
     || header.magic != FrameHeader::magicValue
     || header.version != FrameHeader::versionValue
    )
    {
        return false;
    }

    const int32_t nSurf = r.get<int32_t>();
    const int32_t nVol = r.get<int32_t>();
    if (!r.ok() || nSurf < 0 || nVol < 0)
    {
        return false;
    }

    surfaces.resize(nSurf);
    for (auto& s : surfaces)
    {
        s.id = r.get<int32_t>();
        s.range = r.get<vec2f>();
        const uint64_t nPts = r.get<uint64_t>();
        const uint64_t nSca = r.get<uint64_t>();
        const uint64_t nTri = r.get<uint64_t>();
        r.getArray(s.points, 3*nPts);
        r.getArray(s.scalars, nSca);
        r.getArray(s.triangles, 3*nTri);
    }

    volumes.resize(nVol);
    for (auto& v : volumes)
    {
        v.id = r.get<int32_t>();
        v.grid.dims = r.get<std::array<int32_t, 3>>();
        v.grid.origin = r.get<vec3f>();
        v.grid.spacing = r.get<vec3f>();
        v.range = r.get<vec2f>();
        const uint64_t n = r.get<uint64_t>();
        r.getArray(v.index, n);
        r.getArray(v.value, n);
    }

    return r.ok();
}


// * * * * * * * * * * * * * * BrickDecomposition  * * * * * * * * * * * * * //

anariVis::BrickDecomposition::BrickDecomposition
(
    const std::array<int32_t, 3>& dims,
    int nBricks
)
:
    dims_(dims),
    parts_{{1, 1, 1}}
{
    // Greedy: give each prime factor to the axis with the largest
    // number of cells per part
    for (const int p : primeFactors(std::max(nBricks, 1)))
    {
        int best = 0;
        double bestLen = -1;
        for (int a = 0; a < 3; ++a)
        {
            const double len = double(dims_[a] - 1)/parts_[a];
            if (len > bestLen)
            {
                bestLen = len;
                best = a;
            }
        }
        parts_[best] *= p;
    }

    for (int a = 0; a < 3; ++a)
    {
        const int np = parts_[a];
        const int nCells = std::max(dims_[a] - 1, 1);
        start_[a].resize(np);
        end_[a].resize(np);
        for (int b = 0; b < np; ++b)
        {
            start_[a][b] = int32_t((int64_t(b)*nCells)/np);
            end_[a][b] = int32_t((int64_t(b + 1)*nCells)/np);
        }
    }
}


void anariVis::BrickDecomposition::range
(
    int b,
    std::array<int32_t, 3>& lo,
    std::array<int32_t, 3>& hi
) const
{
    const int bi = b % parts_[0];
    const int bj = (b/parts_[0]) % parts_[1];
    const int bk = b/(parts_[0]*parts_[1]);
    const int idx[3] = {bi, bj, bk};
    for (int a = 0; a < 3; ++a)
    {
        lo[a] = start_[a][idx[a]];
        hi[a] = end_[a][idx[a]];
    }
}


void anariVis::BrickDecomposition::bricksOf
(
    int32_t i,
    int32_t j,
    int32_t k,
    std::vector<int>& bricks
) const
{
    bricks.clear();

    const int32_t ijk[3] = {i, j, k};
    int cand[3][2];
    int nCand[3];

    for (int a = 0; a < 3; ++a)
    {
        nCand[a] = 0;
        const auto& s = start_[a];
        const auto& e = end_[a];

        // First part whose end >= index
        auto it = std::lower_bound(e.begin(), e.end(), ijk[a]);
        int b = int(it - e.begin());
        if (b >= parts_[a])
        {
            b = parts_[a] - 1;
        }
        if (s[b] <= ijk[a] && ijk[a] <= e[b])
        {
            cand[a][nCand[a]++] = b;
        }
        if (b + 1 < parts_[a] && s[b + 1] <= ijk[a] && ijk[a] <= e[b + 1])
        {
            cand[a][nCand[a]++] = b + 1;
        }
    }

    for (int x = 0; x < nCand[0]; ++x)
    {
        for (int y = 0; y < nCand[1]; ++y)
        {
            for (int z = 0; z < nCand[2]; ++z)
            {
                bricks.push_back
                (
                    cand[0][x]
                  + parts_[0]*(cand[1][y] + parts_[1]*cand[2][z])
                );
            }
        }
    }
}


int anariVis::BrickDecomposition::brickOfCell
(
    int32_t i,
    int32_t j,
    int32_t k
) const
{
    const int32_t ijk[3] = {i, j, k};
    int b[3];
    for (int a = 0; a < 3; ++a)
    {
        const int32_t c =
            std::min(std::max(ijk[a], int32_t(0)), std::max(dims_[a] - 2, 0));
        // First part whose end > c (cell c spans points c and c+1)
        const auto& e = end_[a];
        auto it = std::upper_bound(e.begin(), e.end(), c);
        b[a] = std::min(int(it - e.begin()), parts_[a] - 1);
    }
    return b[0] + parts_[0]*(b[1] + parts_[1]*b[2]);
}


std::array<float, 6> anariVis::BrickDecomposition::box
(
    int b,
    const VolumeGrid& grid
) const
{
    std::array<int32_t, 3> lo, hi;
    range(b, lo, hi);
    std::array<float, 6> bb;
    for (int a = 0; a < 3; ++a)
    {
        bb[a] = grid.origin[a] + lo[a]*grid.spacing[a];
        bb[3 + a] = grid.origin[a] + hi[a]*grid.spacing[a];
    }
    return bb;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

std::string anariVis::compositingMode
(
    const RenderSettings& s,
    int nRenderRanks
)
{
    const bool isBarney = (s.library.rfind("barney", 0) == 0);
    std::string comp = s.compositing;

    if (s.library == "ptc")
    {
        return "ptc";
    }
    if (comp == "auto" || comp == "none")
    {
        if (nRenderRanks == 1)
        {
            return "none";
        }
        return isBarney ? "barney" : "deep";
    }
    if (comp == "barney" && !isBarney)
    {
        return "deep";
    }
    return comp;
}


void anariVis::splitByBrick
(
    const SurfacePiece& piece,
    const BrickDecomposition& decomp,
    const VolumeGrid& grid,
    std::vector<std::pair<int, SurfacePiece>>& out
)
{
    out.clear();
    std::map<int, size_t> slot;
    std::map<int, std::vector<int64_t>> remap;
    const bool scalars = (piece.scalars.size() == piece.nPoints());

    for (size_t t = 0; t < piece.nTriangles(); ++t)
    {
        const uint32_t* tri = &piece.triangles[3*t];
        int32_t c[3];
        for (int a = 0; a < 3; ++a)
        {
            const float x =
            (
                piece.points[3*tri[0] + a]
              + piece.points[3*tri[1] + a]
              + piece.points[3*tri[2] + a]
            )/3.f;
            c[a] = int32_t(std::floor((x - grid.origin[a])/grid.spacing[a]));
        }
        const int b = decomp.brickOfCell(c[0], c[1], c[2]);

        auto it = slot.find(b);
        if (it == slot.end())
        {
            it = slot.emplace(b, out.size()).first;
            out.emplace_back(b, SurfacePiece());
            out.back().second.id = piece.id;
            out.back().second.range = piece.range;
            remap[b].assign(piece.nPoints(), -1);
        }
        SurfacePiece& dst = out[it->second].second;
        std::vector<int64_t>& map = remap[b];

        for (int v = 0; v < 3; ++v)
        {
            const uint32_t p = tri[v];
            if (map[p] < 0)
            {
                map[p] = int64_t(dst.nPoints());
                dst.points.push_back(piece.points[3*p]);
                dst.points.push_back(piece.points[3*p + 1]);
                dst.points.push_back(piece.points[3*p + 2]);
                if (scalars)
                {
                    dst.scalars.push_back(piece.scalars[p]);
                }
            }
            dst.triangles.push_back(uint32_t(map[p]));
        }
    }
}


// * * * * * * * * * * * * * * * * Colour maps * * * * * * * * * * * * * * * //

std::vector<anariVis::ColourPoint>
anariVis::namedColourMap(const std::string& name)
{
    if (name == "viridis")
    {
        return
        {
            {0.00f, {{0.267f, 0.005f, 0.329f}}},
            {0.25f, {{0.229f, 0.322f, 0.546f}}},
            {0.50f, {{0.128f, 0.567f, 0.551f}}},
            {0.75f, {{0.369f, 0.789f, 0.383f}}},
            {1.00f, {{0.993f, 0.906f, 0.144f}}}
        };
    }
    else if (name == "inferno")
    {
        return
        {
            {0.00f, {{0.001f, 0.000f, 0.014f}}},
            {0.25f, {{0.341f, 0.062f, 0.429f}}},
            {0.50f, {{0.735f, 0.216f, 0.330f}}},
            {0.75f, {{0.978f, 0.557f, 0.035f}}},
            {1.00f, {{0.988f, 0.998f, 0.645f}}}
        };
    }
    else if (name == "jet")
    {
        return
        {
            {0.000f, {{0.0f, 0.0f, 0.5f}}},
            {0.125f, {{0.0f, 0.0f, 1.0f}}},
            {0.375f, {{0.0f, 1.0f, 1.0f}}},
            {0.625f, {{1.0f, 1.0f, 0.0f}}},
            {0.875f, {{1.0f, 0.0f, 0.0f}}},
            {1.000f, {{0.5f, 0.0f, 0.0f}}}
        };
    }
    else if (name == "grey" || name == "gray")
    {
        return {{0.f, {{0.f, 0.f, 0.f}}}, {1.f, {{1.f, 1.f, 1.f}}}};
    }
    else if (name == "blueWhiteRed")
    {
        return
        {
            {0.0f, {{0.0f, 0.0f, 1.0f}}},
            {0.5f, {{1.0f, 1.0f, 1.0f}}},
            {1.0f, {{1.0f, 0.0f, 0.0f}}}
        };
    }

    // Default: ParaView "Cool to Warm" (Moreland), sampled
    return
    {
        {0.00f, {{0.230f, 0.299f, 0.754f}}},
        {0.25f, {{0.552f, 0.690f, 0.996f}}},
        {0.50f, {{0.865f, 0.865f, 0.865f}}},
        {0.75f, {{0.958f, 0.603f, 0.482f}}},
        {1.00f, {{0.706f, 0.016f, 0.150f}}}
    };
}


const std::vector<std::string>& anariVis::colourMapNames()
{
    static const std::vector<std::string> names
    {
        "coolToWarm", "viridis", "inferno", "jet", "grey", "blueWhiteRed"
    };
    return names;
}


anariVis::vec3f anariVis::evalColourMap
(
    const std::vector<ColourPoint>& map,
    float x
)
{
    if (map.empty())
    {
        return {{x, x, x}};
    }
    if (!(x > map.front().x))
    {
        return map.front().rgb;
    }
    if (!(x < map.back().x))
    {
        return map.back().rgb;
    }
    for (size_t i = 1; i < map.size(); ++i)
    {
        if (x <= map[i].x)
        {
            const float d = map[i].x - map[i-1].x;
            const float t = d > 0 ? (x - map[i-1].x)/d : 0.f;
            vec3f c;
            for (int k = 0; k < 3; ++k)
            {
                c[k] = (1 - t)*map[i-1].rgb[k] + t*map[i].rgb[k];
            }
            return c;
        }
    }
    return map.back().rgb;
}


float anariVis::evalOpacityMap
(
    const std::vector<OpacityPoint>& map,
    float x
)
{
    if (map.empty())
    {
        return x;
    }
    if (!(x > map.front().x))
    {
        return map.front().alpha;
    }
    if (!(x < map.back().x))
    {
        return map.back().alpha;
    }
    for (size_t i = 1; i < map.size(); ++i)
    {
        if (x <= map[i].x)
        {
            const float d = map[i].x - map[i-1].x;
            const float t = d > 0 ? (x - map[i-1].x)/d : 0.f;
            return (1 - t)*map[i-1].alpha + t*map[i].alpha;
        }
    }
    return map.back().alpha;
}


// ************************************************************************* //
