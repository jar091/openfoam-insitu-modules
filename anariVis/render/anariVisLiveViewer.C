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

#include "anariVisLiveViewer.H"

#ifdef ANARIVIS_NEXARIS

#include <nexaris/live/live.h>

#include <sched.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <thread>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

typedef std::chrono::steady_clock clock_type;

double seconds(const clock_type::time_point& t0)
{
    return std::chrono::duration<double>(clock_type::now() - t0).count();
}


double sceneBytes(const anariVis::LocalScene& s)
{
    double b = 0;
    for (const auto& p : s.surfaces)
    {
        b += 4.0*(p.points.size() + p.scalars.size() + p.triangles.size());
    }
    for (const auto& v : s.bricks)
    {
        b += 4.0*v.data.size();
    }
    return b;
}


//- The fixed command block render rank 0 broadcasts every step
struct Command
{
    int32_t quit = 0;
    int32_t follow = 1;
    int32_t paused = 0;
    //- In-situ: the frame loop continues
    int32_t more = 0;
    //- History entry to show (-1: unchanged)
    int32_t show = -1;
    int32_t stateChanged = 0;
    int32_t render = 0;
    float eye[3] = {0, 0, 1};
    float dir[3] = {0, 0, -1};
    float up[3] = {0, 1, 0};
    float fovy = 0.6f;
    float aspect = 1.f;
    int32_t width = 16;
    int32_t height = 16;
    int32_t pixelSamples = 1;
};


// Layer settings <-> the client's layers. Colour maps: the named maps,
// then "configured" (the style's own map, also explicit control points).

int configuredMapIndex()
{
    return int(anariVis::colourMapNames().size());
}


int namedIndex(const std::string& name)
{
    const auto& n = anariVis::colourMapNames();
    for (size_t i = 0; i < n.size(); ++i)
    {
        if (n[i] == name)
        {
            return int(i);
        }
    }
    return configuredMapIndex();
}


std::vector<nexaris::live::Layer> defaultLayers
(
    const anariVis::RenderSettings& settings
)
{
    std::vector<nexaris::live::Layer> layers;
    for (const auto& st : settings.surfaces)
    {
        nexaris::live::Layer l;
        l.name = st.name;
        l.field = st.field;
        l.opacity = st.opacity;
        l.colourMap = st.useField ? namedIndex(st.colourMapName) : -1;
        l.autoRange = st.autoRange;
        l.range[0] = st.range[0];
        l.range[1] = st.range[1];
        layers.push_back(l);
    }
    for (const auto& st : settings.volumes)
    {
        nexaris::live::Layer l;
        l.volume = true;
        l.name = st.name;
        l.field = st.field;
        l.opacity = 1.f;
        l.colourMap = namedIndex(st.colourMapName);
        l.autoRange = st.autoRange;
        l.range[0] = st.range[0];
        l.range[1] = st.range[1];
        layers.push_back(l);
    }
    return layers;
}


anariVis::LayerView toLayerView(const nexaris::live::Layer& l)
{
    anariVis::LayerView v;
    v.visible = l.visible;
    v.opacity = l.opacity;
    if (l.colourMap < 0)
    {
        v.colourMap = l.volume ? -1 : -2;
    }
    else if (l.colourMap >= configuredMapIndex())
    {
        v.colourMap = -1;
    }
    else
    {
        v.colourMap = l.colourMap;
    }
    v.autoRange = l.autoRange;
    v.range = {{l.range[0], l.range[1]}};
    return v;
}


// 6 floats per layer: visible, opacity, colourMap, autoRange, range
void packState
(
    const std::vector<nexaris::live::Layer>& layers,
    std::vector<float>& buf
)
{
    buf.clear();
    for (const auto& l : layers)
    {
        const anariVis::LayerView v = toLayerView(l);
        buf.insert
        (
            buf.end(),
            {
                v.visible ? 1.f : 0.f, v.opacity, float(v.colourMap),
                v.autoRange ? 1.f : 0.f, v.range[0], v.range[1]
            }
        );
    }
}


anariVis::ViewState unpackState
(
    const std::vector<float>& buf,
    size_t nSurfaces,
    size_t nVolumes
)
{
    anariVis::ViewState s;
    for (size_t i = 0; i < nSurfaces + nVolumes; ++i)
    {
        const float* f = &buf[6*i];
        anariVis::LayerView v;
        v.visible = f[0] > 0.5f;
        v.opacity = f[1];
        v.colourMap = int(std::lround(f[2]));
        v.autoRange = f[3] > 0.5f;
        v.range = {{f[4], f[5]}};
        (i < nSurfaces ? s.surfaces : s.volumes).push_back(v);
    }
    return s;
}


//- The Nexaris threads of an in-situ render rank must not share the
//  solver rank's core: they are created with the node's full CPU set
class WideAffinity
{
    cpu_set_t saved_;
    bool active_;

public:

    WideAffinity()
    :
        active_(false)
    {
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

    ~WideAffinity()
    {
        if (active_)
        {
            sched_setaffinity(0, sizeof(saved_), &saved_);
        }
    }
};

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

anariVis::LiveViewer::LiveViewer
(
    Renderer& renderer,
    const RenderSettings& settings,
    MPI_Comm renderComm,
    const ViewerOptions& options
)
:
    renderer_(renderer),
    settings_(settings),
    opt_(options),
    comm_(renderComm),
    rank_(0),
    size_(1),
    nSurfaces_(settings.surfaces.size()),
    nVolumes_(settings.volumes.size()),
    ok_(false),
    shown_(-1),
    follow_(true),
    paused_(false),
    solverFinished_(false),
    quit_(false),
    received_(0),
    passes_(0),
    framePasses_(0),
    inFrame_(false),
    frameStart_(clock_type::now()),
    lastRender_(clock_type::now()),
    lastClient_(clock_type::now()),
    lastAssemble_(0),
    lastFrameTotal_(0),
    historyMaxMB_(0)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    int ok = 1;
    if (rank_ == 0)
    {
        layers_ = defaultLayers(settings_);
        defaults_ = layers_;
        mapNames_ = colourMapNames();
        mapNames_.push_back("configured");

        nexaris::live::Config cfg;
        cfg.port = uint16_t(opt_.port);
        cfg.title = opt_.title;
        cfg.rendererName =
            renderer_.libraryName() + " " + renderer_.mode()
          + " x" + std::to_string(size_);
        cfg.upAxis = opt_.upAxis;
        cfg.jpegQuality = opt_.jpegQuality;
        cfg.maxFps = opt_.maxFps;
        cfg.pixelSamples = opt_.pixelSamples;
        std::string err;
        {
            WideAffinity wide;
            svc_ = nexaris::live::Service::start(cfg, err);
        }
        if (!svc_)
        {
            std::fprintf
            (
                stderr, "anariVis: Nexaris viewer failed to start: %s\n",
                err.c_str()
            );
            ok = 0;
        }
        else
        {
            std::printf
            (
                "anariVis: Nexaris viewer (%s) on %s:%u (loopback),"
                " run dir %s (token file there)\n"
                "  connect: scripts\\win\\connect_live.ps1 -Node %s -Port %u\n",
                opt_.inSitu ? "in-situ" : "in-transit",
                svc_->host().c_str(), unsigned(svc_->port()),
                svc_->runDir().c_str(), svc_->host().c_str(),
                unsigned(svc_->port())
            );
            std::fflush(stdout);
        }
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, comm_);
    ok_ = (ok != 0);
    if (ok_ && rank_ == 0)
    {
        publish();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

anariVis::LiveViewer::~LiveViewer()
{
    svc_.reset();
}


// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * //

void anariVis::LiveViewer::showEntry(int i)
{
    shown_ = i;
    passes_ = 0;
    renderer_.setScene(history_[i]->scene);
}


void anariVis::LiveViewer::publish()
{
    nexaris::live::SceneInfo info;
    info.follow = follow_;
    info.framesReceived = uint64_t(received_);
    info.simulationRunning = !solverFinished_;
    info.canPause = opt_.inSitu;
    info.paused = paused_;
    info.colourMaps = mapNames_;
    for (const auto& h : history_)
    {
        info.historyTimes.push_back(h->scene.time);
    }
    const size_t nLayers = nSurfaces_ + nVolumes_;
    if (shown_ >= 0 && shown_ < int(history_.size()))
    {
        const Entry& e = *history_[shown_];
        const auto& b = e.scene.bounds;
        for (int k = 0; k < 3; ++k)
        {
            info.boundsMin[k] = b[k];
            info.boundsMax[k] = b[k + 3];
        }
        const ViewCamera home = renderer_.defaultCamera(b);
        for (int k = 0; k < 3; ++k)
        {
            info.homeEye[k] = home.eye[k];
            info.homeTarget[k] = home.target[k];
            // linear -> display (sRGB)
            const float c = std::clamp(settings_.background[k], 0.f, 1.f);
            info.background[k] =
                c <= 0.0031308f
              ? 12.92f*c : 1.055f*std::pow(c, 1.f/2.4f) - 0.055f;
        }
        info.homeFovy = home.fovy;
        info.time = e.scene.time;
        info.frameIndex = e.scene.frameIndex;
        info.historyIndex = uint32_t(shown_);
        for (size_t k = 0; k < nLayers && k < layers_.size(); ++k)
        {
            layers_[k].dataRange[0] = e.ranges[k][0];
            layers_[k].dataRange[1] = e.ranges[k][1];
            if (layers_[k].autoRange)
            {
                layers_[k].range[0] = e.ranges[k][0];
                layers_[k].range[1] = e.ranges[k][1];
            }
        }
    }
    info.layers = layers_;
    svc_->publish(info);
    svc_->setTelemetry
    ({
        {"live_render_ranks", double(size_)},
        {"live_frames_received", double(received_)},
        {"live_history_frames", double(history_.size())},
        {"live_history_mb_max", historyMaxMB_},
        {"live_assemble_ms", 1e3*lastAssemble_},
        {"live_build_ms", 1e3*last_.build},
        {"live_trace_ms", 1e3*last_.render},
        {"live_frame_ms", 1e3*lastFrameTotal_},
        {"live_solver_running", solverFinished_ ? 0.0 : 1.0},
        {"live_in_situ", opt_.inSitu ? 1.0 : 0.0}
    });
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void anariVis::LiveViewer::addFrame(LocalScene&& scene, double assembleSeconds)
{
    const size_t nLayers = nSurfaces_ + nVolumes_;
    auto e = std::make_unique<Entry>();
    e->scene = std::move(scene);

    // Global field range per layer (min of lo, min of -hi): a rank may hold
    // no piece of a surface
    const float inf = std::numeric_limits<float>::infinity();
    std::vector<float> r(2*nLayers, inf);
    for (const auto& p : e->scene.surfaces)
    {
        if (p.id >= 0 && size_t(p.id) < nSurfaces_ && !p.triangles.empty())
        {
            r[2*p.id] = std::min(r[2*p.id], p.range[0]);
            r[2*p.id + 1] = std::min(r[2*p.id + 1], -p.range[1]);
        }
    }
    for (const auto& b : e->scene.bricks)
    {
        if (b.id >= 0 && size_t(b.id) < nVolumes_)
        {
            const size_t k = nSurfaces_ + b.id;
            r[2*k] = std::min(r[2*k], b.range[0]);
            r[2*k + 1] = std::min(r[2*k + 1], -b.range[1]);
        }
    }
    if (!r.empty())
    {
        MPI_Allreduce
        (
            MPI_IN_PLACE, r.data(), int(r.size()), MPI_FLOAT, MPI_MIN, comm_
        );
    }
    e->ranges.resize(nLayers);
    for (size_t k = 0; k < nLayers; ++k)
    {
        e->ranges[k] = {{r[2*k], -r[2*k + 1]}};
    }
    e->bytes = sceneBytes(e->scene);
    history_.push_back(std::move(e));
    ++received_;

    // Bounded history: count and the largest rank's memory
    while (history_.size() > 1)
    {
        double mb = 0;
        for (const auto& h : history_)
        {
            mb += h->bytes/1048576.0;
        }
        MPI_Allreduce(MPI_IN_PLACE, &mb, 1, MPI_DOUBLE, MPI_MAX, comm_);
        historyMaxMB_ = mb;
        if
        (
            int(history_.size()) <= std::max(opt_.historySize, 1)
         && mb <= opt_.historyMB
        )
        {
            break;
        }
        // the oldest frame goes, except the shown one (a held frame stays
        // until the user leaves it)
        const int victim = (shown_ == 0) ? 1 : 0;
        history_.erase(history_.begin() + victim);
        if (shown_ > victim)
        {
            --shown_;
        }
    }

    // In-situ the solver waits for the viewer: follow shows the new frame at
    // once (in-transit a newer frame waits for minPasses, see step)
    if (shown_ < 0 || (opt_.inSitu && follow_))
    {
        showEntry(int(history_.size()) - 1);
    }
    lastAssemble_ = assembleSeconds;
    if (rank_ == 0)
    {
        const LocalScene& s = history_.back()->scene;
        std::printf
        (
            "anariVis viewer: frame %d t=%g received (assemble %.3f s),"
            " history %d frames, %.0f MB/rank max\n",
            s.frameIndex, s.time, assembleSeconds,
            int(history_.size()), historyMaxMB_
        );
        std::fflush(stdout);
    }
}


void anariVis::LiveViewer::setSolverFinished()
{
    if (solverFinished_)
    {
        return;
    }
    solverFinished_ = true;
    paused_ = false;
    lastClient_ = clock_type::now();
    if (rank_ == 0)
    {
        std::printf
        (
            "anariVis viewer: the solver finished after %ld frames;"
            " the viewer stays until no client for %.0f s\n",
            received_, opt_.keepAlive
        );
        std::fflush(stdout);
        publish();
    }
}


anariVis::LiveViewer::StepResult anariVis::LiveViewer::step(int waitMs)
{
    const size_t nLayers = nSurfaces_ + nVolumes_;

    // 1. Render rank 0: client commands and the next view
    Command cmd;
    cmd.follow = follow_ ? 1 : 0;
    cmd.paused = paused_ ? 1 : 0;
    if (rank_ == 0)
    {
        for (const auto& c : svc_->takeCommands())
        {
            if (nexaris::live::applyLayerCommand(layers_, c, mapNames_.size()))
            {
                cmd.stateChanged = 1;
                continue;
            }
            switch (c.op)
            {
                case nexaris::live::Command::Follow:
                    cmd.follow = c.flag ? 1 : 0;
                    if (c.flag && !history_.empty())
                    {
                        cmd.show = int(history_.size()) - 1;
                    }
                    break;
                case nexaris::live::Command::ShowFrame:
                    if (c.index < history_.size())
                    {
                        cmd.follow = 0;
                        cmd.show = int(c.index);
                    }
                    break;
                case nexaris::live::Command::Reset:
                    layers_ = defaults_;
                    cmd.stateChanged = 1;
                    cmd.follow = 1;
                    if (!history_.empty())
                    {
                        cmd.show = int(history_.size()) - 1;
                    }
                    break;
                case nexaris::live::Command::PauseSimulation:
                    cmd.paused = (opt_.inSitu && c.flag && !solverFinished_) ? 1 : 0;
                    break;
                default:
                    break;
            }
        }

        // in-transit follow: a newer frame after minPasses passes (or when
        // nobody renders)
        if
        (
            cmd.follow && cmd.show < 0 && !history_.empty()
         && shown_ != int(history_.size()) - 1
         && (passes_ >= opt_.minPasses || seconds(lastRender_) > 0.5)
        )
        {
            cmd.show = int(history_.size()) - 1;
        }

        const bool client = svc_->clientConnected();
        if (client)
        {
            lastClient_ = clock_type::now();
        }
        if
        (
            svc_->stopped()
         || (solverFinished_ && seconds(lastClient_) > opt_.keepAlive)
        )
        {
            cmd.quit = 1;
            cmd.paused = 0;
        }
        else if (!history_.empty())
        {
            nexaris::live::View view;
            if (svc_->nextRender(view, waitMs))
            {
                cmd.render = 1;
                std::copy(view.eye, view.eye + 3, cmd.eye);
                std::copy(view.dir, view.dir + 3, cmd.dir);
                std::copy(view.up, view.up + 3, cmd.up);
                cmd.fovy = view.fovy;
                cmd.aspect = view.aspect;
                cmd.width = int32_t(view.width);
                cmd.height = int32_t(view.height);
                cmd.pixelSamples = view.pixelSamples;
            }
        }
        else if (waitMs > 0)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(waitMs));
        }

        // in-situ: does the solver keep waiting? While paused; else until
        // the frame got its passes (only with a client) or the time is up
        if (inFrame_ && !cmd.quit)
        {
            const long done = framePasses_ + (cmd.render ? 1 : 0);
            cmd.more =
                cmd.paused
             || (
                    (client || cmd.render)
                 && done < opt_.passesPerFrame
                 && seconds(frameStart_) < opt_.maxSecondsPerFrame
                );
        }
        if (cmd.stateChanged)
        {
            packState(layers_, stateBuf_);
        }
    }
    MPI_Bcast(&cmd, int(sizeof(Command)), MPI_BYTE, 0, comm_);

    // 2. All render ranks: apply and render
    follow_ = cmd.follow != 0;
    paused_ = cmd.paused != 0;
    if (cmd.stateChanged && nLayers)
    {
        stateBuf_.resize(6*nLayers);
        MPI_Bcast(stateBuf_.data(), int(stateBuf_.size()), MPI_FLOAT, 0, comm_);
        renderer_.setViewState(unpackState(stateBuf_, nSurfaces_, nVolumes_));
    }
    if (cmd.show >= 0 && cmd.show < int(history_.size()))
    {
        showEntry(cmd.show);
    }
    if (cmd.render)
    {
        const auto t0 = clock_type::now();
        ViewCamera cam;
        cam.eye = {{cmd.eye[0], cmd.eye[1], cmd.eye[2]}};
        cam.dir = {{cmd.dir[0], cmd.dir[1], cmd.dir[2]}};
        cam.up = {{cmd.up[0], cmd.up[1], cmd.up[2]}};
        cam.fovy = cmd.fovy;
        cam.aspect = cmd.aspect;
        cam.width = cmd.width;
        cam.height = cmd.height;
        cam.pixelSamples = cmd.pixelSamples;
        last_ = renderer_.renderView(cam, rgba_, depth_);
        lastFrameTotal_ = seconds(t0);
        ++passes_;
        ++framePasses_;
        lastRender_ = clock_type::now();
        if (rank_ == 0)
        {
            if (rgba_.size() == size_t(cam.width)*cam.height)
            {
                nexaris::live::Image img;
                img.width = uint32_t(cam.width);
                img.height = uint32_t(cam.height);
                img.rgba.resize(4*rgba_.size());
                std::memcpy(img.rgba.data(), rgba_.data(), img.rgba.size());
                img.depth = std::move(depth_);
                svc_->finishRender(std::move(img));
            }
            else
            {
                svc_->failRender("the renderer returned no image");
            }
        }
    }

    // 3. Rank 0: the state for the client
    if (rank_ == 0)
    {
        publish();
    }

    quit_ = cmd.quit != 0;
    StepResult r;
    r.quit = quit_;
    r.rendered = cmd.render != 0;
    r.paused = paused_;
    r.more = cmd.more != 0;
    return r;
}


void anariVis::LiveViewer::runInSituFrame()
{
    if (quit_)
    {
        return;
    }
    inFrame_ = true;
    framePasses_ = 0;
    frameStart_ = clock_type::now();
    // the first step only takes a view that is already waiting (a client
    // that kept the Nexaris render thread busy while the solver ran)
    int waitMs = 2;
    while (true)
    {
        const StepResult r = step(waitMs);
        if (r.quit || !r.more)
        {
            break;
        }
        waitMs = r.paused ? 50 : 20;
    }
    inFrame_ = false;
}


void anariVis::LiveViewer::serveUntilQuit()
{
    while (!quit_)
    {
        step(50);
    }
    if (rank_ == 0)
    {
        std::printf
        (
            "anariVis viewer: stopped (%ld frames received)\n", received_
        );
        std::fflush(stdout);
    }
}


#endif // ANARIVIS_NEXARIS

// ************************************************************************* //
