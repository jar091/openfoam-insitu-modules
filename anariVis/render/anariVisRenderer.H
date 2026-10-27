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
    anariVis::Renderer

Description
    ANARI scene manager and renderer of the anariVis layer.

    Every render rank owns one Renderer. It converts its local scene
    (triangulated surfaces and structured volume bricks) into ANARI objects,
    renders them and produces the final image on rank 0 of the render
    communicator. Compositing modes:

    - \c none   : a single render rank, any ANARI library.
    - \c barney : Barney device subtype "mpi" (data-parallel, Barney traces
                  rays across ranks and composites the image itself).
    - \c deep   : every rank renders its local data with any ANARI library
                  (helide, cycles, barney_cpu) and the colour + depth images
                  are sort-last composited with deepCompositing
                  (dc::Compositor).
    - \c ptc    : the deepCompositing ANARI wrapper device "ptc" wraps
                  another library and composites inside anariRenderFrame.

    All calls of render() are collective over the render communicator.

SourceFiles
    anariVisRenderer.C

\*---------------------------------------------------------------------------*/

#ifndef anariVis_Renderer_H
#define anariVis_Renderer_H

#include "anariVisTypes.H"

#include <anari/anari.h>
#include <mpi.h>

#include <memory>

// Forward declarations
namespace dc
{
    struct Compositor;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace anariVis
{

//- A structured volume brick owned by one render rank
struct VolumeBrick
{
    int32_t id = 0;

    //- The global grid
    VolumeGrid grid;

    //- Point range of this brick (inclusive)
    std::array<int32_t, 3> lo{{0, 0, 0}};
    std::array<int32_t, 3> hi{{0, 0, 0}};

    vec2f range{{0.f, 1.f}};

    //- Point values (x fastest), size = product of (hi - lo + 1)
    std::vector<float> data;

    std::array<int32_t, 3> size() const
    {
        return {{hi[0] - lo[0] + 1, hi[1] - lo[1] + 1, hi[2] - lo[2] + 1}};
    }
};


//- The scene of one render rank for one frame
struct LocalScene
{
    int32_t frameIndex = 0;

    double time = 0;

    std::array<float, 6> bounds{{0.f, 0.f, 0.f, 1.f, 1.f, 1.f}};

    //- All data of this rank lies in a convex brick (spatial routing)
    bool hasBrick = false;

    //- The brick (min xyz, max xyz)
    std::array<float, 6> brick{{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};

    //- Surfaces, merged per surface id
    std::vector<SurfacePiece> surfaces;

    std::vector<VolumeBrick> bricks;
};


//- Interactive viewing: how the user shows one layer (surface or volume)
struct LayerView
{
    bool visible = true;
    //- Surface: opacity [0,1]; volume: density scale (1 = as configured)
    float opacity = 1.f;
    //- Index into colourMapNames(); -1 = the configured map; -2 = uniform
    //  colour (surfaces)
    int colourMap = -1;
    bool autoRange = true;
    vec2f range{{0.f, 1.f}};
};


//- Interactive viewing: the user's settings of all layers
struct ViewState
{
    std::vector<LayerView> surfaces;
    std::vector<LayerView> volumes;
};


//- Interactive viewing: one camera (simulation coordinates)
struct ViewCamera
{
    vec3f eye{{0.f, 0.f, 1.f}};
    vec3f dir{{0.f, 0.f, -1.f}};
    vec3f up{{0.f, 1.f, 0.f}};
    //- Look-at point (defaultCamera only)
    vec3f target{{0.f, 0.f, 0.f}};
    //- Vertical field of view [rad]
    float fovy = 0.6f;
    float aspect = 1.f;
    int width = 640;
    int height = 480;
    int pixelSamples = 1;
};


//- Timings of one frame [s]
struct FrameTimings
{
    double build = 0;
    double render = 0;
    double composite = 0;
    double write = 0;
    size_t nTriangles = 0;
};


/*---------------------------------------------------------------------------*\
                          Class Renderer Declaration
\*---------------------------------------------------------------------------*/

class Renderer
{
    // Private Data

        RenderSettings settings_;

        //- Render communicator (duplicate, owned)
        MPI_Comm comm_;

        int rank_;
        int size_;

        //- Compositing mode: none | barney | deep | ptc
        std::string mode_;

        //- Name of the loaded ANARI library
        std::string libName_;

        ANARILibrary library_;
        ANARIDevice device_;
        ANARIRenderer renderer_;
        ANARICamera camera_;
        ANARIFrame frame_;

        std::unique_ptr<dc::Compositor> compositor_;

        //- Camera frame used for the last render (for brick depth)
        vec3f eye_, dir_, up_;
        float fovy_;
        float aspect_;

        // Interactive viewing (setScene / setViewState / renderView)

            //- Scene of the interactive world (not owned)
            const LocalScene* iScene_;

            ViewState iState_;
            bool haveState_;

            //- World must be rebuilt before the next renderView
            bool iWorldDirty_;

            ANARIWorld iWorld_;

            //- Lights of iWorld_ that follow the camera, with their
            //  camera-space direction
            std::vector<std::pair<ANARILight, vec3f>> iCamLights_;

            //- Last committed camera / frame size / samples
            ViewCamera iCam_;
            bool iCamValid_;
            int iSpp_;

            //- channel.depth was added to the frame
            bool iDepthRequested_;


    // Private Member Functions

        void setupCamera(const LocalScene& scene);

        //- Build the ANARI world; with a view state its layer settings
        //  replace the configured style, with camLights the lights that
        //  follow the camera are kept (referenced) in it
        ANARIWorld buildWorld
        (
            const LocalScene& scene,
            size_t& nTris,
            const ViewState* view = nullptr,
            std::vector<std::pair<ANARILight, vec3f>>* camLights = nullptr
        );

        ANARIArray1D newArray1D
        (
            ANARIDataType type,
            const void* data,
            size_t n,
            size_t elemBytes
        );

        void compositeDeep(const LocalScene& scene, std::vector<uint32_t>&);

        void writeImage
        (
            const std::vector<uint32_t>& rgba,
            const LocalScene& scene
        ) const;


public:

    // Constructors

        //- Construct from settings. Collective over renderComm
        Renderer(const RenderSettings& settings, MPI_Comm renderComm);

        //- No copy
        Renderer(const Renderer&) = delete;
        void operator=(const Renderer&) = delete;


    //- Destructor
    ~Renderer();


    // Member Functions

        //- The compositing mode in use
        const std::string& mode() const { return mode_; }

        //- The ANARI library in use
        const std::string& libraryName() const { return libName_; }

        int rank() const { return rank_; }
        int size() const { return size_; }

        //- Render one frame. Collective over the render communicator.
        //  The image is written by render rank 0.
        FrameTimings render(const LocalScene& scene);


        // Interactive viewing (anariVisRender with the Nexaris viewer).
        // All calls are collective over the render communicator and must
        // get identical arguments on every render rank (except the scene,
        // which is the rank's own part).

            //- The scene to show; it must stay alive until the next
            //  setScene (the world is rebuilt at the next renderView)
            void setScene(const LocalScene& scene);

            //- The user's layer settings
            void setViewState(const ViewState& state);

            //- Render one accumulation pass of the current scene. Unchanged
            //  camera, scene and settings refine the previous image.
            //  Render rank 0 gets the image (RGBA8 sRGB, row 0 = bottom)
            //  and, where the device delivers it, the ray distance per
            //  pixel (empty otherwise).
            FrameTimings renderView
            (
                const ViewCamera& cam,
                std::vector<uint32_t>& rgba,
                std::vector<float>& depth
            );

            //- The camera the configured settings give for these bounds
            ViewCamera defaultCamera(const std::array<float, 6>& bounds) const;
};


//- Write an RGBA8 image (row 0 = bottom) as PNG
bool writePNG
(
    const std::string& fileName,
    int width,
    int height,
    const uint32_t* rgba
);


} // End namespace anariVis

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
