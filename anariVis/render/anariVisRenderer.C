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

#include "anariVisRenderer.H"

#include <deepCompositing/deepCompositing.h>

#include <fenv.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <sched.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <limits>
#include <stdexcept>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

typedef std::chrono::steady_clock clock_type;

double seconds(const clock_type::time_point& t0)
{
    return std::chrono::duration<double>(clock_type::now() - t0).count();
}

anariVis::vec3f sub(const anariVis::vec3f& a, const anariVis::vec3f& b)
{
    return {{a[0] - b[0], a[1] - b[1], a[2] - b[2]}};
}

anariVis::vec3f cross(const anariVis::vec3f& a, const anariVis::vec3f& b)
{
    return
    {{
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    }};
}

float dot(const anariVis::vec3f& a, const anariVis::vec3f& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

anariVis::vec3f normalised(const anariVis::vec3f& a)
{
    const float l = std::sqrt(dot(a, a));
    if (l > 0)
    {
        return {{a[0]/l, a[1]/l, a[2]/l}};
    }
    return {{0.f, 0.f, 1.f}};
}


int g_statusRank = 0;


//- Masks floating-point exception traps (OpenFOAM's FOAM_SIGFPE) while the
//  ANARI devices run: Embree/SIMD code raises benign FP exceptions. Worker
//  threads created inside the guard inherit the masked environment.
class FpeGuard
{
    int saved_;

public:

    FpeGuard()
    :
        saved_(fegetexcept())
    {
        fedisableexcept(FE_ALL_EXCEPT);
    }

    ~FpeGuard()
    {
        feclearexcept(FE_ALL_EXCEPT);
        if (saved_ > 0)
        {
            feenableexcept(saved_);
        }
    }
};

void statusFunc
(
    const void*,
    ANARIDevice,
    ANARIObject,
    ANARIDataType,
    ANARIStatusSeverity severity,
    ANARIStatusCode,
    const char* message
)
{
    const char* tag = nullptr;
    switch (severity)
    {
        case ANARI_SEVERITY_FATAL_ERROR: tag = "FATAL"; break;
        case ANARI_SEVERITY_ERROR: tag = "ERROR"; break;
        case ANARI_SEVERITY_WARNING: tag = "WARN "; break;
        default: return;
    }
    std::fprintf
    (
        stderr, "[anariVis render rank %d][%s] %s\n",
        g_statusRank, tag, message
    );
}


//- Distance at which a ray enters box [lo,hi] (+inf: miss)
float boxEntry
(
    const anariVis::vec3f& org,
    const anariVis::vec3f& dir,
    const anariVis::vec3f& lo,
    const anariVis::vec3f& hi
)
{
    float t0 = 0.f;
    float t1 = std::numeric_limits<float>::infinity();
    for (int k = 0; k < 3; ++k)
    {
        const float inv = 1.f/dir[k];
        float a = (lo[k] - org[k])*inv;
        float b = (hi[k] - org[k])*inv;
        if (a > b) std::swap(a, b);
        t0 = std::max(t0, a);
        t1 = std::min(t1, b);
    }
    return t0 <= t1 ? t0 : std::numeric_limits<float>::infinity();
}


uint8_t toSrgb8(float l)
{
    l = std::clamp(l, 0.f, 1.f);
    const float s =
        l <= 0.0031308f ? 12.92f*l : 1.055f*std::pow(l, 1.f/2.4f) - 0.055f;
    return uint8_t(std::clamp(s, 0.f, 1.f)*255.f + 0.5f);
}


//- Temporarily allow the calling thread to run on all cores of the node
//  (in-situ render ranks are usually bound to one core by the launcher).
//  Threads created inside the scope (Barney's CPU launch system, Embree,
//  helide/cycles workers) keep the wide mask; the caller's own mask is
//  restored on exit so that the solver part of the rank stays on its core.
class AffinityScope
{
    cpu_set_t saved_;
    bool active_;

public:

    explicit AffinityScope(bool on)
    :
        active_(false)
    {
        if (!on)
        {
            return;
        }
        CPU_ZERO(&saved_);
        if (sched_getaffinity(0, sizeof(saved_), &saved_) != 0)
        {
            return;
        }
        cpu_set_t all;
        CPU_ZERO(&all);
        const long n = sysconf(_SC_NPROCESSORS_CONF);
        for (long i = 0; i < n && i < CPU_SETSIZE; ++i)
        {
            CPU_SET(int(i), &all);
        }
        active_ = (sched_setaffinity(0, sizeof(all), &all) == 0);
    }

    ~AffinityScope()
    {
        if (active_)
        {
            sched_setaffinity(0, sizeof(saved_), &saved_);
        }
    }
};

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

anariVis::Renderer::Renderer
(
    const RenderSettings& settings,
    MPI_Comm renderComm
)
:
    settings_(settings),
    comm_(MPI_COMM_NULL),
    rank_(0),
    size_(1),
    mode_("none"),
    libName_(settings.library),
    library_(nullptr),
    device_(nullptr),
    renderer_(nullptr),
    camera_(nullptr),
    frame_(nullptr),
    compositor_(nullptr),
    eye_{{0.f, 0.f, 1.f}},
    dir_{{0.f, 0.f, -1.f}},
    up_{{0.f, 1.f, 0.f}},
    fovy_(0.6f),
    aspect_(1.f),
    iScene_(nullptr),
    haveState_(false),
    iWorldDirty_(true),
    iWorld_(nullptr),
    iCamValid_(false),
    iSpp_(-1),
    iDepthRequested_(false)
{
    FpeGuard fpeGuard;
    AffinityScope affinity(settings.widenAffinity);

    MPI_Comm_dup(renderComm, &comm_);
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
    g_statusRank = rank_;

    // Select the compositing mode
    mode_ = compositingMode(settings_, size_);
    if
    (
        rank_ == 0
     && settings_.compositing == "barney"
     && mode_ != "barney"
    )
    {
        std::fprintf
        (
            stderr, "anariVis: compositing 'barney' requires a barney"
            " library, using '%s' for '%s'\n", mode_.c_str(), libName_.c_str()
        );
    }

    // Render threads: Barney CPU reads BARNEY_MAX_CPU_THREADS (Embree
    // threads) at device creation; OpenMP is used by the deep compositing
    if (settings_.threads > 0)
    {
        const std::string n = std::to_string(settings_.threads);
        setenv("BARNEY_MAX_CPU_THREADS", n.c_str(), 1);
        #ifdef _OPENMP
        omp_set_num_threads(settings_.threads);
        #endif
    }

    // Load library and create device
    std::string loadName = libName_;
    std::string wrapped = settings_.wrappedLibrary;
    if (mode_ == "ptc" && libName_ != "ptc")
    {
        wrapped = libName_;
        loadName = "ptc";
    }

    library_ = anariLoadLibrary(loadName.c_str(), statusFunc, nullptr);
    if (!library_)
    {
        throw std::runtime_error
        (
            "anariVis: could not load ANARI library '" + loadName + "'"
        );
    }

    device_ = anariNewDevice
    (
        library_,
        mode_ == "barney" ? "mpi" : "default"
    );
    if (!device_)
    {
        throw std::runtime_error("anariVis: could not create ANARI device");
    }

    if (mode_ == "barney")
    {
        // Barney dereferences the pointer at commit; comm_ outlives device
        const uint64_t p = uint64_t(reinterpret_cast<uintptr_t>(&comm_));
        anariSetParameter
        (
            device_, device_, "pointer_to_mpi_communicator",
            ANARI_UINT64, &p
        );
    }
    else if (mode_ == "ptc")
    {
        anariSetParameter
        (
            device_, device_, "wrappedLibrary", ANARI_STRING, wrapped.c_str()
        );
        // VOID_POINTER: 'mem' is the pointer value itself
        anariSetParameter
        (
            device_, device_, "mpiCommunicator", ANARI_VOID_POINTER, &comm_
        );
        libName_ = "ptc(" + wrapped + ")";
    }
    anariCommitParameters(device_, device_);

    // Renderer
    renderer_ = anariNewRenderer(device_, settings_.rendererType.c_str());
    {
        const bool transparent = (mode_ == "deep" || mode_ == "ptc");
        const vec4f bg =
            transparent ? vec4f{{0.f, 0.f, 0.f, 0.f}} : settings_.background;
        anariSetParameter
        (
            device_, renderer_, "background", ANARI_FLOAT32_VEC4, bg.data()
        );
        anariSetParameter
        (
            device_, renderer_, "pixelSamples", ANARI_INT32,
            &settings_.pixelSamples
        );
        anariSetParameter
        (
            device_, renderer_, "ambientRadiance", ANARI_FLOAT32,
            &settings_.ambientRadiance
        );
        anariSetParameter
        (
            device_, renderer_, "ambientColor", ANARI_FLOAT32_VEC3,
            settings_.ambientColour.data()
        );
        if (settings_.denoise)
        {
            const bool on = true;
            anariSetParameter(device_, renderer_, "denoise", ANARI_BOOL, &on);
        }
        if (settings_.library == "cycles" || wrapped == "cycles")
        {
            // Offline frames: full resolution and full sample count
            const bool off = false;
            anariSetParameter
            (
                device_, renderer_, "interactiveScaling", ANARI_BOOL, &off
            );
            anariSetParameter(device_, renderer_, "runAsync", ANARI_BOOL, &off);
        }
        anariCommitParameters(device_, renderer_);
    }

    camera_ = anariNewCamera(device_, "perspective");

    frame_ = anariNewFrame(device_);
    {
        const uint32_t sz[2] =
        {
            uint32_t(settings_.width), uint32_t(settings_.height)
        };
        anariSetParameter(device_, frame_, "size", ANARI_UINT32_VEC2, sz);

        const ANARIDataType colourType =
            (mode_ == "deep") ? ANARI_FLOAT32_VEC4 : ANARI_UFIXED8_RGBA_SRGB;
        anariSetParameter
        (
            device_, frame_, "channel.color", ANARI_DATA_TYPE, &colourType
        );
        if (mode_ == "deep")
        {
            const ANARIDataType depthType = ANARI_FLOAT32;
            anariSetParameter
            (
                device_, frame_, "channel.depth", ANARI_DATA_TYPE, &depthType
            );
        }
        if (mode_ == "ptc")
        {
            anariSetParameter
            (
                device_, frame_, "ptc.background", ANARI_FLOAT32_VEC4,
                settings_.background.data()
            );
        }
        const bool accum = true;
        anariSetParameter(device_, frame_, "accumulation", ANARI_BOOL, &accum);
        anariSetParameter(device_, frame_, "renderer", ANARI_RENDERER, &renderer_);
        anariSetParameter(device_, frame_, "camera", ANARI_CAMERA, &camera_);
        anariCommitParameters(device_, frame_);
    }

    if (mode_ == "deep")
    {
        compositor_.reset(new dc::Compositor(1, comm_));
        compositor_->resize(settings_.width, settings_.height);
    }

    if (settings_.writeImages && !settings_.outputDir.empty() && rank_ == 0)
    {
        std::error_code ec;
        std::filesystem::create_directories(settings_.outputDir, ec);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

anariVis::Renderer::~Renderer()
{
    FpeGuard fpeGuard;
    AffinityScope affinity(settings_.widenAffinity);

    compositor_.reset(nullptr);

    if (device_)
    {
        for (auto& l : iCamLights_) anariRelease(device_, l.first);
        if (iWorld_) anariRelease(device_, iWorld_);
        if (frame_) anariRelease(device_, frame_);
        if (camera_) anariRelease(device_, camera_);
        if (renderer_) anariRelease(device_, renderer_);
        anariRelease(device_, device_);
    }
    if (library_)
    {
        anariUnloadLibrary(library_);
    }

    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized && comm_ != MPI_COMM_NULL)
    {
        MPI_Comm_free(&comm_);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

ANARIArray1D anariVis::Renderer::newArray1D
(
    ANARIDataType type,
    const void* data,
    size_t n,
    size_t elemBytes
)
{
    // Managed array: the device owns a copy, no lifetime constraints
    ANARIArray1D arr =
        anariNewArray1D(device_, nullptr, nullptr, nullptr, type, n);
    void* dst = anariMapArray(device_, arr);
    if (dst && n)
    {
        std::memcpy(dst, data, n*elemBytes);
    }
    anariUnmapArray(device_, arr);
    return arr;
}


void anariVis::Renderer::setupCamera(const LocalScene& scene)
{
    const CameraSettings& c = settings_.camera;
    const auto& b = scene.bounds;

    const vec3f lo{{b[0], b[1], b[2]}};
    const vec3f hi{{b[3], b[4], b[5]}};
    const vec3f centre
    {{
        0.5f*(lo[0] + hi[0]), 0.5f*(lo[1] + hi[1]), 0.5f*(lo[2] + hi[2])
    }};
    const vec3f diag = sub(hi, lo);
    const float radius = std::max(0.5f*std::sqrt(dot(diag, diag)), 1e-12f);

    aspect_ = float(settings_.width)/float(settings_.height);
    fovy_ = c.fovy*float(M_PI)/180.f;

    const vec3f target = c.autoLookAt ? centre : c.lookAt;

    if (c.autoPosition)
    {
        dir_ = normalised(c.viewDirection);
        const float fovx = 2.f*std::atan(std::tan(0.5f*fovy_)*aspect_);
        const float halfFov = 0.5f*std::min(fovy_, fovx);
        const float dist = radius/std::sin(halfFov)/std::max(c.zoom, 1e-3f);
        eye_ =
        {{
            target[0] - dir_[0]*dist,
            target[1] - dir_[1]*dist,
            target[2] - dir_[2]*dist
        }};
    }
    else
    {
        eye_ = c.position;
        dir_ = normalised(sub(target, eye_));
    }

    up_ = normalised(c.up);
    if (std::abs(dot(up_, dir_)) > 0.999f)
    {
        up_ = std::abs(dir_[2]) < 0.9f
            ? vec3f{{0.f, 0.f, 1.f}} : vec3f{{0.f, 1.f, 0.f}};
    }

    anariSetParameter(device_, camera_, "position", ANARI_FLOAT32_VEC3, eye_.data());
    anariSetParameter(device_, camera_, "direction", ANARI_FLOAT32_VEC3, dir_.data());
    anariSetParameter(device_, camera_, "up", ANARI_FLOAT32_VEC3, up_.data());
    anariSetParameter(device_, camera_, "fovy", ANARI_FLOAT32, &fovy_);
    anariSetParameter(device_, camera_, "aspect", ANARI_FLOAT32, &aspect_);
    anariCommitParameters(device_, camera_);
}


ANARIWorld anariVis::Renderer::buildWorld
(
    const LocalScene& scene,
    size_t& nTris,
    const ViewState* view,
    std::vector<std::pair<ANARILight, vec3f>>* camLights
)
{
    nTris = 0;
    ANARIWorld world = anariNewWorld(device_);

    std::vector<ANARISurface> surfaces;
    std::vector<ANARIVolume> volumes;
    std::vector<ANARILight> lights;

    // Surfaces
    for (const SurfacePiece& s : scene.surfaces)
    {
        if (s.triangles.empty() || s.id < 0)
        {
            continue;
        }
        const SurfaceStyle* style =
            (size_t(s.id) < settings_.surfaces.size())
          ? &settings_.surfaces[s.id] : nullptr;
        const LayerView* lv =
            (view && size_t(s.id) < view->surfaces.size())
          ? &view->surfaces[s.id] : nullptr;
        if (lv && !lv->visible)
        {
            continue;
        }

        const size_t nPts = s.nPoints();
        const size_t nTri = s.nTriangles();
        nTris += nTri;

        ANARIGeometry geom = anariNewGeometry(device_, "triangle");

        ANARIArray1D arr = newArray1D
        (
            ANARI_FLOAT32_VEC3, s.points.data(), nPts, 3*sizeof(float)
        );
        anariSetParameter(device_, geom, "vertex.position", ANARI_ARRAY1D, &arr);
        anariRelease(device_, arr);

        arr = newArray1D
        (
            ANARI_UINT32_VEC3, s.triangles.data(), nTri, 3*sizeof(uint32_t)
        );
        anariSetParameter(device_, geom, "primitive.index", ANARI_ARRAY1D, &arr);
        anariRelease(device_, arr);

        // Area-weighted vertex normals that agree with the triangle winding
        if (!style || style->smoothNormals)
        {
            std::vector<float> nrm(3*nPts, 0.f);
            for (size_t t = 0; t < nTri; ++t)
            {
                const uint32_t* tri = &s.triangles[3*t];
                const float* p0 = &s.points[3*tri[0]];
                const float* p1 = &s.points[3*tri[1]];
                const float* p2 = &s.points[3*tri[2]];
                const vec3f e1{{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]}};
                const vec3f e2{{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]}};
                const vec3f n = cross(e1, e2);
                for (int v = 0; v < 3; ++v)
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        nrm[3*tri[v] + k] += n[k];
                    }
                }
            }
            for (size_t i = 0; i < nPts; ++i)
            {
                const vec3f n = normalised
                (
                    vec3f{{nrm[3*i], nrm[3*i + 1], nrm[3*i + 2]}}
                );
                nrm[3*i] = n[0];
                nrm[3*i + 1] = n[1];
                nrm[3*i + 2] = n[2];
            }
            arr = newArray1D
            (
                ANARI_FLOAT32_VEC3, nrm.data(), nPts, 3*sizeof(float)
            );
            anariSetParameter(device_, geom, "vertex.normal", ANARI_ARRAY1D, &arr);
            anariRelease(device_, arr);
        }

        const bool byField =
            style && style->useField && s.scalars.size() == nPts
         && !(lv && lv->colourMap == -2);
        const float opacity =
            lv ? std::clamp(lv->opacity, 0.f, 1.f)
               : (style ? style->opacity : 1.f);

        if (byField)
        {
            const vec2f range =
                (lv && !lv->autoRange) ? lv->range : s.range;
            const float lo = range[0];
            const float d = (range[1] > lo) ? (range[1] - lo) : 1.f;
            const std::vector<ColourPoint>& cmap =
                (lv && lv->colourMap >= 0)
              ? namedColourMap(colourMapNames()[lv->colourMap])
              : style->colourMap;
            std::vector<float> col(4*nPts);
            for (size_t i = 0; i < nPts; ++i)
            {
                const vec3f c =
                    evalColourMap(cmap, (s.scalars[i] - lo)/d);
                col[4*i] = c[0];
                col[4*i + 1] = c[1];
                col[4*i + 2] = c[2];
                col[4*i + 3] = opacity;
            }
            arr = newArray1D
            (
                ANARI_FLOAT32_VEC4, col.data(), nPts, 4*sizeof(float)
            );
            anariSetParameter(device_, geom, "vertex.color", ANARI_ARRAY1D, &arr);
            anariRelease(device_, arr);
        }
        anariCommitParameters(device_, geom);

        ANARIMaterial mat = anariNewMaterial(device_, "matte");
        if (byField)
        {
            anariSetParameter(device_, mat, "color", ANARI_STRING, "color");
        }
        else
        {
            const vec3f c = style ? style->colour : vec3f{{0.8f, 0.8f, 0.8f}};
            anariSetParameter(device_, mat, "color", ANARI_FLOAT32_VEC3, c.data());
        }
        if (opacity < 1.f)
        {
            anariSetParameter(device_, mat, "opacity", ANARI_FLOAT32, &opacity);
            anariSetParameter(device_, mat, "alphaMode", ANARI_STRING, "blend");
        }
        anariCommitParameters(device_, mat);

        ANARISurface surf = anariNewSurface(device_);
        anariSetParameter(device_, surf, "geometry", ANARI_GEOMETRY, &geom);
        anariSetParameter(device_, surf, "material", ANARI_MATERIAL, &mat);
        anariCommitParameters(device_, surf);
        anariRelease(device_, geom);
        anariRelease(device_, mat);
        surfaces.push_back(surf);
    }

    // Volumes
    for (const VolumeBrick& vb : scene.bricks)
    {
        if (vb.data.empty() || vb.id < 0)
        {
            continue;
        }
        const VolumeStyle* style =
            (size_t(vb.id) < settings_.volumes.size())
          ? &settings_.volumes[vb.id] : nullptr;
        const LayerView* lv =
            (view && size_t(vb.id) < view->volumes.size())
          ? &view->volumes[vb.id] : nullptr;
        if (lv && !lv->visible)
        {
            continue;
        }

        const auto n = vb.size();
        ANARISpatialField field =
            anariNewSpatialField(device_, "structuredRegular");

        ANARIArray3D arr = anariNewArray3D
        (
            device_, nullptr, nullptr, nullptr, ANARI_FLOAT32,
            uint64_t(n[0]), uint64_t(n[1]), uint64_t(n[2])
        );
        void* dst = anariMapArray(device_, arr);
        std::memcpy(dst, vb.data.data(), vb.data.size()*sizeof(float));
        anariUnmapArray(device_, arr);
        anariSetParameter(device_, field, "data", ANARI_ARRAY3D, &arr);
        anariRelease(device_, arr);

        const vec3f origin
        {{
            vb.grid.origin[0] + vb.lo[0]*vb.grid.spacing[0],
            vb.grid.origin[1] + vb.lo[1]*vb.grid.spacing[1],
            vb.grid.origin[2] + vb.lo[2]*vb.grid.spacing[2]
        }};
        anariSetParameter(device_, field, "origin", ANARI_FLOAT32_VEC3, origin.data());
        anariSetParameter
        (
            device_, field, "spacing", ANARI_FLOAT32_VEC3, vb.grid.spacing.data()
        );
        anariSetParameter(device_, field, "filter", ANARI_STRING, "linear");
        anariCommitParameters(device_, field);

        // Transfer function
        const int nTF = 256;
        std::vector<float> tfColour(3*nTF);
        std::vector<float> tfOpacity(nTF);
        const std::vector<ColourPoint> cmap =
            (lv && lv->colourMap >= 0)
          ? namedColourMap(colourMapNames()[lv->colourMap])
          : (style ? style->colourMap : namedColourMap("coolToWarm"));
        const std::vector<OpacityPoint> omap =
            style ? style->opacityMap : std::vector<OpacityPoint>();
        for (int i = 0; i < nTF; ++i)
        {
            const float x = float(i)/(nTF - 1);
            const vec3f c = evalColourMap(cmap, x);
            tfColour[3*i] = c[0];
            tfColour[3*i + 1] = c[1];
            tfColour[3*i + 2] = c[2];
            tfOpacity[i] = evalOpacityMap(omap, x);
        }

        ANARIVolume vol = anariNewVolume(device_, "transferFunction1D");
        anariSetParameter(device_, vol, "value", ANARI_SPATIAL_FIELD, &field);
        anariRelease(device_, field);

        ANARIArray1D carr =
            newArray1D(ANARI_FLOAT32_VEC3, tfColour.data(), nTF, 3*sizeof(float));
        anariSetParameter(device_, vol, "color", ANARI_ARRAY1D, &carr);
        anariRelease(device_, carr);
        ANARIArray1D oarr =
            newArray1D(ANARI_FLOAT32, tfOpacity.data(), nTF, sizeof(float));
        anariSetParameter(device_, vol, "opacity", ANARI_ARRAY1D, &oarr);
        anariRelease(device_, oarr);

        float vr[2] = {vb.range[0], vb.range[1]};
        if (lv && !lv->autoRange)
        {
            vr[0] = lv->range[0];
            vr[1] = lv->range[1];
        }
        if (!(vr[1] > vr[0]))
        {
            vr[1] = vr[0] + 1.f;
        }
        anariSetParameter(device_, vol, "valueRange", ANARI_FLOAT32_BOX1, vr);

        float unitDistance = style ? style->unitDistance : 0.f;
        if (!(unitDistance > 0))
        {
            const auto& b = scene.bounds;
            const float dx = b[3] - b[0], dy = b[4] - b[1], dz = b[5] - b[2];
            unitDistance = 0.02f*std::sqrt(dx*dx + dy*dy + dz*dz);
        }
        if (lv)
        {
            // density scale: a shorter unit distance is more opaque
            unitDistance /= std::max(lv->opacity, 1e-3f);
        }
        anariSetParameter(device_, vol, "unitDistance", ANARI_FLOAT32, &unitDistance);
        anariCommitParameters(device_, vol);
        volumes.push_back(vol);
    }

    // Lights
    {
        const vec3f right = normalised(cross(dir_, up_));
        const vec3f camUp = cross(right, dir_);

        std::vector<LightSettings> ls = settings_.lights;
        if (ls.empty())
        {
            LightSettings key;
            key.cameraRelative = true;
            key.direction = {{0.3f, -0.5f, 1.f}};
            key.intensity = 3.f;
            ls.push_back(key);
        }
        for (const LightSettings& l : ls)
        {
            vec3f d = l.direction;
            if (l.cameraRelative)
            {
                d =
                {{
                    d[0]*right[0] + d[1]*camUp[0] + d[2]*dir_[0],
                    d[0]*right[1] + d[1]*camUp[1] + d[2]*dir_[1],
                    d[0]*right[2] + d[1]*camUp[2] + d[2]*dir_[2]
                }};
            }
            d = normalised(d);

            ANARILight light = anariNewLight(device_, l.type.c_str());
            anariSetParameter(device_, light, "color", ANARI_FLOAT32_VEC3, l.colour.data());
            if (l.type == "point")
            {
                anariSetParameter
                (
                    device_, light, "position", ANARI_FLOAT32_VEC3,
                    l.position.data()
                );
                anariSetParameter
                (
                    device_, light, "intensity", ANARI_FLOAT32, &l.intensity
                );
            }
            else
            {
                anariSetParameter(device_, light, "direction", ANARI_FLOAT32_VEC3, d.data());
                anariSetParameter
                (
                    device_, light, "irradiance", ANARI_FLOAT32, &l.intensity
                );
            }
            anariCommitParameters(device_, light);
            lights.push_back(light);
            if (camLights && l.cameraRelative && l.type != "point")
            {
                anariRetain(device_, light);
                camLights->push_back({light, l.direction});
            }
        }
    }

    if (!surfaces.empty())
    {
        ANARIArray1D a = newArray1D
        (
            ANARI_SURFACE, surfaces.data(), surfaces.size(), sizeof(ANARISurface)
        );
        anariSetParameter(device_, world, "surface", ANARI_ARRAY1D, &a);
        anariRelease(device_, a);
    }
    if (!volumes.empty())
    {
        ANARIArray1D a = newArray1D
        (
            ANARI_VOLUME, volumes.data(), volumes.size(), sizeof(ANARIVolume)
        );
        anariSetParameter(device_, world, "volume", ANARI_ARRAY1D, &a);
        anariRelease(device_, a);
    }
    if (!lights.empty())
    {
        ANARIArray1D a = newArray1D
        (
            ANARI_LIGHT, lights.data(), lights.size(), sizeof(ANARILight)
        );
        anariSetParameter(device_, world, "light", ANARI_ARRAY1D, &a);
        anariRelease(device_, a);
    }
    anariCommitParameters(device_, world);

    for (auto o : surfaces) anariRelease(device_, o);
    for (auto o : volumes) anariRelease(device_, o);
    for (auto o : lights) anariRelease(device_, o);

    return world;
}


void anariVis::Renderer::compositeDeep
(
    const LocalScene& scene,
    std::vector<uint32_t>& image
)
{
    const int w = settings_.width;
    const int h = settings_.height;
    const size_t nPix = size_t(w)*h;

    dc::DeviceInterface di = compositor_->prepare();

    uint32_t mw = 0, mh = 0;
    ANARIDataType t = ANARI_UNKNOWN;
    const float* colour = static_cast<const float*>
    (
        anariMapFrame(device_, frame_, "channel.color", &mw, &mh, &t)
    );
    const float* depth = static_cast<const float*>
    (
        anariMapFrame(device_, frame_, "channel.depth", &mw, &mh, &t)
    );

    if (colour)
    {
        std::vector<float> z(nPix, std::numeric_limits<float>::infinity());
        if (depth)
        {
            std::copy(depth, depth + nPix, z.begin());
        }

        // Fragment depth: the device depth is exact for opaque surfaces.
        // For volumes use the distance at which the pixel ray enters the
        // (convex) brick of this rank, which orders disjoint bricks
        // correctly; some devices (cycles) write depth 0 for volumes.
        const bool haveBox = scene.hasBrick || !scene.bricks.empty();
        const bool brickDepth =
            haveBox
         && (
                settings_.depthMode == "brick"
             || settings_.depthMode == "auto"
             || !depth
            );

        if (brickDepth)
        {
            const vec3f right = normalised(cross(dir_, up_));
            const vec3f camUp = cross(right, dir_);
            const float tanHalf = std::tan(0.5f*fovy_);

            #pragma omp parallel for schedule(static)
            for (int64_t i = 0; i < int64_t(nPix); ++i)
            {
                const int x = int(i % w);
                const int y = int(i / w);
                const float sx = (2.f*(x + 0.5f)/w - 1.f)*tanHalf*aspect_;
                const float sy = (2.f*(y + 0.5f)/h - 1.f)*tanHalf;
                const vec3f d = normalised
                (
                    vec3f
                    {{
                        dir_[0] + sx*right[0] + sy*camUp[0],
                        dir_[1] + sx*right[1] + sy*camUp[1],
                        dir_[2] + sx*right[2] + sy*camUp[2]
                    }}
                );
                float zb = std::numeric_limits<float>::infinity();
                if (scene.hasBrick)
                {
                    const auto& bb = scene.brick;
                    zb = boxEntry
                    (
                        eye_, d,
                        vec3f{{bb[0], bb[1], bb[2]}},
                        vec3f{{bb[3], bb[4], bb[5]}}
                    );
                }
                else
                for (const VolumeBrick& vb : scene.bricks)
                {
                    const vec3f lo
                    {{
                        vb.grid.origin[0] + vb.lo[0]*vb.grid.spacing[0],
                        vb.grid.origin[1] + vb.lo[1]*vb.grid.spacing[1],
                        vb.grid.origin[2] + vb.lo[2]*vb.grid.spacing[2]
                    }};
                    const vec3f hi
                    {{
                        vb.grid.origin[0] + vb.hi[0]*vb.grid.spacing[0],
                        vb.grid.origin[1] + vb.hi[1]*vb.grid.spacing[1],
                        vb.grid.origin[2] + vb.hi[2]*vb.grid.spacing[2]
                    }};
                    zb = std::min(zb, boxEntry(eye_, d, lo, hi));
                }
                if (std::isfinite(zb))
                {
                    z[i] = zb;
                }
            }
        }

        dc::writeImageFragments(di, colour, z.data());
    }
    else if (rank_ == 0)
    {
        std::fprintf(stderr, "anariVis: device returned no colour channel\n");
    }

    if (colour) anariUnmapFrame(device_, frame_, "channel.color");
    if (depth) anariUnmapFrame(device_, frame_, "channel.depth");

    std::vector<uint32_t> linear(rank_ == 0 ? nPix : 0);
    compositor_->finish(rank_ == 0 ? linear.data() : nullptr);

    if (rank_ == 0)
    {
        // Blend background under the premultiplied linear result, sRGB encode
        const vec4f& bg = settings_.background;
        image.resize(nPix);
        #pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < int64_t(nPix); ++i)
        {
            const uint32_t p = linear[i];
            const float a = ((p >> 24) & 0xff)/255.f;
            uint32_t out = 0xff000000u;
            for (int k = 0; k < 3; ++k)
            {
                const float c = ((p >> (8*k)) & 0xff)/255.f + (1.f - a)*bg[k];
                out |= uint32_t(toSrgb8(c)) << (8*k);
            }
            image[i] = out;
        }
    }
}


void anariVis::Renderer::writeImage
(
    const std::vector<uint32_t>& rgba,
    const LocalScene& scene
) const
{
    char name[64];
    std::snprintf(name, sizeof(name), "_%06d.png", int(scene.frameIndex));
    std::string file = settings_.prefix + name;
    if (!settings_.outputDir.empty())
    {
        file = settings_.outputDir + "/" + file;
    }
    if (!writePNG(file, settings_.width, settings_.height, rgba.data()))
    {
        std::fprintf(stderr, "anariVis: could not write %s\n", file.c_str());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

anariVis::FrameTimings anariVis::Renderer::render(const LocalScene& scene)
{
    FpeGuard fpeGuard;
    AffinityScope affinity(settings_.widenAffinity);

    FrameTimings t;

    auto t0 = clock_type::now();

    setupCamera(scene);
    ANARIWorld world = buildWorld(scene, t.nTriangles);
    anariSetParameter(device_, frame_, "world", ANARI_WORLD, &world);
    anariCommitParameters(device_, frame_);
    anariRelease(device_, world);

    t.build = seconds(t0);

    auto t1 = clock_type::now();
    for (int pass = 0; pass < std::max(settings_.passes, 1); ++pass)
    {
        anariRenderFrame(device_, frame_);
        anariFrameReady(device_, frame_, ANARI_WAIT);
    }
    t.render = seconds(t1);

    auto t2 = clock_type::now();
    std::vector<uint32_t> image;
    const size_t nPix = size_t(settings_.width)*settings_.height;

    if (mode_ == "deep")
    {
        compositeDeep(scene, image);
    }
    else if (rank_ == 0)
    {
        uint32_t w = 0, h = 0;
        ANARIDataType type = ANARI_UNKNOWN;
        const uint32_t* px = static_cast<const uint32_t*>
        (
            anariMapFrame(device_, frame_, "channel.color", &w, &h, &type)
        );
        if (px && size_t(w)*h == nPix)
        {
            image.assign(px, px + nPix);
        }
        else
        {
            std::fprintf(stderr, "anariVis: could not map colour channel\n");
        }
        anariUnmapFrame(device_, frame_, "channel.color");

        // The final image is opaque
        for (auto& p : image)
        {
            p |= 0xff000000u;
        }
    }
    t.composite = seconds(t2);

    auto t3 = clock_type::now();
    if (rank_ == 0 && settings_.writeImages && !image.empty())
    {
        writeImage(image, scene);
    }
    t.write = seconds(t3);

    return t;
}


// * * * * * * * * * * * * * * * Interactive viewing * * * * * * * * * * * * //

void anariVis::Renderer::setScene(const LocalScene& scene)
{
    iScene_ = &scene;
    iWorldDirty_ = true;
}


void anariVis::Renderer::setViewState(const ViewState& state)
{
    iState_ = state;
    haveState_ = true;
    iWorldDirty_ = true;
}


anariVis::ViewCamera anariVis::Renderer::defaultCamera
(
    const std::array<float, 6>& b
) const
{
    const CameraSettings& c = settings_.camera;
    const vec3f lo{{b[0], b[1], b[2]}};
    const vec3f hi{{b[3], b[4], b[5]}};
    const vec3f centre
    {{
        0.5f*(lo[0] + hi[0]), 0.5f*(lo[1] + hi[1]), 0.5f*(lo[2] + hi[2])
    }};
    const vec3f diag = sub(hi, lo);
    const float radius = std::max(0.5f*std::sqrt(dot(diag, diag)), 1e-12f);

    ViewCamera v;
    v.width = settings_.width;
    v.height = settings_.height;
    v.aspect = float(settings_.width)/float(settings_.height);
    v.fovy = c.fovy*float(M_PI)/180.f;
    const vec3f target = c.autoLookAt ? centre : c.lookAt;
    v.target = target;
    if (c.autoPosition)
    {
        v.dir = normalised(c.viewDirection);
        const float fovx = 2.f*std::atan(std::tan(0.5f*v.fovy)*v.aspect);
        const float halfFov = 0.5f*std::min(v.fovy, fovx);
        const float dist = radius/std::sin(halfFov)/std::max(c.zoom, 1e-3f);
        v.eye =
        {{
            target[0] - v.dir[0]*dist,
            target[1] - v.dir[1]*dist,
            target[2] - v.dir[2]*dist
        }};
    }
    else
    {
        v.eye = c.position;
        v.dir = normalised(sub(target, v.eye));
    }
    v.up = normalised(c.up);
    if (std::abs(dot(v.up, v.dir)) > 0.999f)
    {
        v.up = std::abs(v.dir[2]) < 0.9f
            ? vec3f{{0.f, 0.f, 1.f}} : vec3f{{0.f, 1.f, 0.f}};
    }
    return v;
}


anariVis::FrameTimings anariVis::Renderer::renderView
(
    const ViewCamera& cam,
    std::vector<uint32_t>& rgba,
    std::vector<float>& depth
)
{
    FpeGuard fpeGuard;
    AffinityScope affinity(settings_.widenAffinity);

    FrameTimings t;
    rgba.clear();
    depth.clear();
    if (!iScene_)
    {
        return t;
    }

    auto t0 = clock_type::now();

    bool camChanged = !iCamValid_;
    if (iWorldDirty_)
    {
        for (auto& l : iCamLights_) anariRelease(device_, l.first);
        iCamLights_.clear();
        if (iWorld_) anariRelease(device_, iWorld_);
        iWorld_ = buildWorld
        (
            *iScene_, t.nTriangles, haveState_ ? &iState_ : nullptr,
            &iCamLights_
        );
        anariSetParameter(device_, frame_, "world", ANARI_WORLD, &iWorld_);
        anariCommitParameters(device_, frame_);
        iWorldDirty_ = false;
        camChanged = true;
    }

    const int w = std::max(cam.width, 1);
    const int h = std::max(cam.height, 1);
    const bool sizeChanged =
        !iCamValid_ || w != iCam_.width || h != iCam_.height
     || !iDepthRequested_;
    if (sizeChanged)
    {
        const uint32_t sz[2] = {uint32_t(w), uint32_t(h)};
        anariSetParameter(device_, frame_, "size", ANARI_UINT32_VEC2, sz);
        if (!iDepthRequested_ && mode_ != "deep")
        {
            // ray distance per pixel for the client's reprojection
            const ANARIDataType depthType = ANARI_FLOAT32;
            anariSetParameter
            (
                device_, frame_, "channel.depth", ANARI_DATA_TYPE, &depthType
            );
        }
        iDepthRequested_ = true;
        anariCommitParameters(device_, frame_);
        settings_.width = w;
        settings_.height = h;
        if (compositor_)
        {
            compositor_->resize(w, h);
        }
    }

    if
    (
        !camChanged
     && (cam.eye != iCam_.eye || cam.dir != iCam_.dir || cam.up != iCam_.up
      || cam.fovy != iCam_.fovy || cam.aspect != iCam_.aspect)
    )
    {
        camChanged = true;
    }
    if (camChanged)
    {
        eye_ = cam.eye;
        dir_ = normalised(cam.dir);
        const vec3f right = normalised(cross(dir_, cam.up));
        up_ = cross(right, dir_);
        fovy_ = cam.fovy;
        aspect_ = cam.aspect;
        anariSetParameter(device_, camera_, "position", ANARI_FLOAT32_VEC3, eye_.data());
        anariSetParameter(device_, camera_, "direction", ANARI_FLOAT32_VEC3, dir_.data());
        anariSetParameter(device_, camera_, "up", ANARI_FLOAT32_VEC3, up_.data());
        anariSetParameter(device_, camera_, "fovy", ANARI_FLOAT32, &fovy_);
        anariSetParameter(device_, camera_, "aspect", ANARI_FLOAT32, &aspect_);
        anariCommitParameters(device_, camera_);

        // Lights that follow the camera (the default key light)
        for (const auto& l : iCamLights_)
        {
            const vec3f& c = l.second;
            const vec3f d = normalised
            (
                vec3f
                {{
                    c[0]*right[0] + c[1]*up_[0] + c[2]*dir_[0],
                    c[0]*right[1] + c[1]*up_[1] + c[2]*dir_[1],
                    c[0]*right[2] + c[1]*up_[2] + c[2]*dir_[2]
                }}
            );
            anariSetParameter(device_, l.first, "direction", ANARI_FLOAT32_VEC3, d.data());
            anariCommitParameters(device_, l.first);
        }
        if (!iCamLights_.empty() && std::getenv("ANARIVIS_LIGHT_WORLD_COMMIT"))
        {
            anariCommitParameters(device_, iWorld_);
        }
    }

    const int spp = std::max(cam.pixelSamples, 1);
    if (spp != iSpp_)
    {
        anariSetParameter(device_, renderer_, "pixelSamples", ANARI_INT32, &spp);
        anariCommitParameters(device_, renderer_);
        iSpp_ = spp;
    }

    iCam_ = cam;
    iCam_.width = w;
    iCam_.height = h;
    iCamValid_ = true;
    t.build = seconds(t0);

    auto t1 = clock_type::now();
    anariRenderFrame(device_, frame_);
    anariFrameReady(device_, frame_, ANARI_WAIT);
    t.render = seconds(t1);

    auto t2 = clock_type::now();
    const size_t nPix = size_t(w)*h;
    if (mode_ == "deep")
    {
        compositeDeep(*iScene_, rgba);
    }
    else if (rank_ == 0)
    {
        uint32_t mw = 0, mh = 0;
        ANARIDataType type = ANARI_UNKNOWN;
        const uint32_t* px = static_cast<const uint32_t*>
        (
            anariMapFrame(device_, frame_, "channel.color", &mw, &mh, &type)
        );
        if (px && size_t(mw)*mh == nPix)
        {
            rgba.assign(px, px + nPix);
        }
        if (px) anariUnmapFrame(device_, frame_, "channel.color");

        const float* z = static_cast<const float*>
        (
            anariMapFrame(device_, frame_, "channel.depth", &mw, &mh, &type)
        );
        if (z && type == ANARI_FLOAT32 && size_t(mw)*mh == nPix)
        {
            depth.assign(z, z + nPix);
        }
        if (z) anariUnmapFrame(device_, frame_, "channel.depth");
    }
    t.composite = seconds(t2);
    return t;
}


// ************************************************************************* //
