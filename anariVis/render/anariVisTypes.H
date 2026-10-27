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

Namespace
    anariVis

Description
    Renderer-agnostic data types of the anariVis in-situ / in-transit
    visualisation layer: render settings, surface and volume packets
    exchanged between OpenFOAM producer ranks and render (consumer) ranks,
    and the brick decomposition of the structured volume grid.

    These types deliberately do not depend on OpenFOAM headers so that the
    renderer side (ANARI, Barney, deepCompositing) is isolated from the
    solver side, as proposed in the design document.

SourceFiles
    anariVisPacket.C

\*---------------------------------------------------------------------------*/

#ifndef anariVis_Types_H
#define anariVis_Types_H

#include <array>
#include <cstdint>
#include <string>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace anariVis
{

typedef std::array<float, 2> vec2f;
typedef std::array<float, 3> vec3f;
typedef std::array<float, 4> vec4f;

//- Colour map control point: x in [0,1] and an RGB colour
struct ColourPoint
{
    float x;
    vec3f rgb;
};

//- Opacity control point: x in [0,1] and opacity
struct OpacityPoint
{
    float x;
    float alpha;
};


/*---------------------------------------------------------------------------*\
                          Render settings (POD)
\*---------------------------------------------------------------------------*/

//- Style of one surface (iso-surface, patch, plane ...)
struct SurfaceStyle
{
    std::string name;

    //- Name of the colour field (informational, "" = uniform colour)
    std::string field;

    //- Colour by a scalar field (else uniform colour)
    bool useField = false;

    //- Uniform colour
    vec3f colour{{0.8f, 0.8f, 0.8f}};

    //- Colour map control points (if useField)
    std::vector<ColourPoint> colourMap;

    //- Name of the colour map ("" = explicit control points)
    std::string colourMapName;

    //- Use the global min/max of the colour field as range
    bool autoRange = true;

    //- Fixed range (if !autoRange)
    vec2f range{{0.f, 1.f}};

    //- Surface opacity
    float opacity = 1.f;

    //- Smooth (per-vertex) normals
    bool smoothNormals = true;
};


//- Style of one volume
struct VolumeStyle
{
    std::string name;

    //- Sampled field (informational)
    std::string field;

    std::vector<ColourPoint> colourMap;

    //- Name of the colour map ("" = explicit control points)
    std::string colourMapName;

    std::vector<OpacityPoint> opacityMap;

    bool autoRange = true;

    vec2f range{{0.f, 1.f}};

    //- ANARI transferFunction1D unitDistance (<= 0: automatic)
    float unitDistance = 0.f;
};


//- Camera
struct CameraSettings
{
    //- Place the camera automatically from the scene bounds
    bool autoPosition = true;

    //- Use the centre of the scene bounds as look-at point
    bool autoLookAt = true;

    vec3f position{{0.f, 0.f, 1.f}};
    vec3f lookAt{{0.f, 0.f, 0.f}};
    vec3f up{{0.f, 1.f, 0.f}};

    //- Direction from the camera towards the look-at point (autoPosition)
    vec3f viewDirection{{-1.f, -0.6f, -1.2f}};

    //- Vertical field of view [deg]
    float fovy = 35.f;

    //- Zoom factor for the automatic placement (> 1: closer)
    float zoom = 1.f;
};


//- Light
struct LightSettings
{
    //- directional | point
    std::string type = "directional";

    //- Direction (directional light)
    vec3f direction{{0.f, 0.f, -1.f}};

    //- Position (point light)
    vec3f position{{0.f, 0.f, 0.f}};

    //- Interpret direction as (right, up, forward) in camera space
    bool cameraRelative = false;

    vec3f colour{{1.f, 1.f, 1.f}};

    float intensity = 3.f;
};


//- Complete set of rendering settings
struct RenderSettings
{
    //- ANARI library: barney_cpu | barney | helide | cycles | ptc
    std::string library = "barney_cpu";

    //- Wrapped library for library 'ptc'
    std::string wrappedLibrary = "helide";

    //- Compositing: auto | none | barney | deep | ptc
    std::string compositing = "auto";

    //- ANARI renderer subtype
    std::string rendererType = "default";

    int width = 1024;
    int height = 768;

    //- Samples per pixel per pass
    int pixelSamples = 4;

    //- Accumulation passes per image
    int passes = 1;

    vec4f background{{0.12f, 0.12f, 0.15f, 1.f}};

    float ambientRadiance = 0.3f;

    vec3f ambientColour{{1.f, 1.f, 1.f}};

    //- Renderer denoising (if supported by the device)
    bool denoise = false;

    //- Render threads per render rank (0: automatic)
    int threads = 0;

    //- Let the render calls (and the threads they create) use all cores of
    //  the node; the caller's CPU affinity is restored afterwards
    bool widenAffinity = false;

    //- Fragment depth in deep compositing: auto | device | brick
    std::string depthMode = "auto";

    CameraSettings camera;

    std::vector<LightSettings> lights;

    std::vector<SurfaceStyle> surfaces;

    std::vector<VolumeStyle> volumes;

    //- Directory for the images
    std::string outputDir = ".";

    //- Image file prefix
    std::string prefix = "anariVis";

    bool writeImages = true;

    //- Additional log output
    bool verbose = false;
};


/*---------------------------------------------------------------------------*\
                       Data packets (producer -> consumer)
\*---------------------------------------------------------------------------*/

//- Definition of a global structured volume grid (point based)
struct VolumeGrid
{
    std::array<int32_t, 3> dims{{0, 0, 0}};
    vec3f origin{{0.f, 0.f, 0.f}};
    vec3f spacing{{1.f, 1.f, 1.f}};

    uint64_t nPoints() const
    {
        return uint64_t(dims[0])*uint64_t(dims[1])*uint64_t(dims[2]);
    }
};


//- One triangulated surface fragment of one producer
struct SurfacePiece
{
    int32_t id = 0;

    //- Global colour range
    vec2f range{{0.f, 1.f}};

    //- xyz triplets
    std::vector<float> points;

    //- One scalar per point (empty: uniform colour)
    std::vector<float> scalars;

    //- Three point indices per triangle
    std::vector<uint32_t> triangles;

    size_t nPoints() const { return points.size()/3; }
    size_t nTriangles() const { return triangles.size()/3; }
};


//- Samples of a volume grid sent to one brick owner
struct VolumePiece
{
    int32_t id = 0;

    VolumeGrid grid;

    //- Global value range
    vec2f range{{0.f, 1.f}};

    //- Global linear point index (i + nx*(j + ny*k))
    std::vector<uint32_t> index;

    std::vector<float> value;
};


//- Header of a frame packet
struct FrameHeader
{
    static constexpr int32_t magicValue = 0x53495641;   // "AVIS"
    static constexpr int32_t versionValue = 2;

    //- Flag: no more frames from this producer
    static constexpr int32_t stopFlag = 1;

    int32_t magic = magicValue;
    int32_t version = versionValue;

    //- Frame index (the OpenFOAM time index)
    int32_t frameIndex = 0;

    int32_t flags = 0;

    //- Producer rank (in the visualisation communicator)
    int32_t producer = 0;

    //- Simulation time
    double time = 0;

    //- Global scene bounds (min xyz, max xyz)
    std::array<float, 6> bounds{{0.f, 0.f, 0.f, 1.f, 1.f, 1.f}};

    //- Data of this packet was routed spatially: it lies in the convex
    //  brick 'brick' (min xyz, max xyz) of the destination render rank
    int32_t hasBrick = 0;

    std::array<float, 6> brick{{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};
};


//- A frame packet: header + surfaces + volume samples
struct FramePacket
{
    FrameHeader header;
    std::vector<SurfacePiece> surfaces;
    std::vector<VolumePiece> volumes;

    //- Serialise into a contiguous byte buffer
    void pack(std::vector<char>& buf) const;

    //- Deserialise from a byte buffer. Returns false on a format error
    bool unpack(const char* data, size_t nBytes);

    //- Peek at the header only
    static bool peekHeader(const char* data, size_t nBytes, FrameHeader& h);
};


/*---------------------------------------------------------------------------*\
                      Brick decomposition of a volume grid
\*---------------------------------------------------------------------------*/

//- Split the point grid 'dims' into nBricks convex bricks.
//  Neighbouring bricks share their boundary point plane, so that
//  trilinear interpolation is continuous across bricks.
class BrickDecomposition
{
    std::array<int32_t, 3> dims_{{0, 0, 0}};
    std::array<int32_t, 3> parts_{{1, 1, 1}};

    //- Point range [start, end] (inclusive) per axis and part
    std::array<std::vector<int32_t>, 3> start_, end_;

public:

    BrickDecomposition() = default;

    BrickDecomposition(const std::array<int32_t, 3>& dims, int nBricks);

    int nBricks() const { return parts_[0]*parts_[1]*parts_[2]; }

    const std::array<int32_t, 3>& parts() const { return parts_; }

    //- Point range of brick b: lo (inclusive), hi (inclusive)
    void range
    (
        int b,
        std::array<int32_t, 3>& lo,
        std::array<int32_t, 3>& hi
    ) const;

    //- All bricks that contain the grid point (i,j,k) (1 to 8 bricks)
    void bricksOf(int32_t i, int32_t j, int32_t k, std::vector<int>& b) const;

    //- The brick that owns the grid cell (i,j,k) (cell i between points
    //  i and i+1; indices are clamped)
    int brickOfCell(int32_t i, int32_t j, int32_t k) const;

    //- Physical box (min xyz, max xyz) of brick b of grid
    std::array<float, 6> box(int b, const VolumeGrid& grid) const;
};


//- The compositing mode that a renderer with these settings uses on
//  nRenderRanks ranks: none | barney | deep | ptc
std::string compositingMode(const RenderSettings& s, int nRenderRanks);

//- Split a surface piece by triangle centroid into the bricks of decomp
//  (grid = the point grid the decomposition refers to).
//  Returns brick index -> piece (with compact point lists).
void splitByBrick
(
    const SurfacePiece& piece,
    const BrickDecomposition& decomp,
    const VolumeGrid& grid,
    std::vector<std::pair<int, SurfacePiece>>& out
);

//- Colour map by name (coolToWarm, viridis, inferno, jet, grey, blueWhiteRed)
std::vector<ColourPoint> namedColourMap(const std::string& name);

//- The names namedColourMap() knows, in a fixed order
const std::vector<std::string>& colourMapNames();

//- Evaluate a colour map at x in [0,1]
vec3f evalColourMap(const std::vector<ColourPoint>& map, float x);

//- Evaluate an opacity map at x in [0,1]
float evalOpacityMap(const std::vector<OpacityPoint>& map, float x);


} // End namespace anariVis

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
