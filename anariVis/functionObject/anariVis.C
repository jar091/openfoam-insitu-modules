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

#include "anariVis.H"
#include "anariVisSettingsIO.H"
#include "anariVisMPI.H"
#include "anariVisTransport.H"
#include "anariVisRenderer.H"
#include "anariVisScene.H"
#include "interpolation.H"
#include "OStringStream.H"
#include "addToRunTimeSelectionTable.H"

#include <algorithm>
#include <chrono>
#include <map>
#include <set>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(anariVis, 0);
    addToRunTimeSelectionTable(functionObject, anariVis, dictionary);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

typedef std::chrono::steady_clock clock_type;

double seconds(const clock_type::time_point& t0)
{
    return std::chrono::duration<double>(clock_type::now() - t0).count();
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::anariVis::initialise()
{
    if (initialised_)
    {
        return;
    }

    const auto tInit = clock_type::now();

    // MPI communicator of the visualisation layer
    MPI_Comm parent = MPI_COMM_NULL;

    if (!UPstream::parRun())
    {
        int flag = 0;
        MPI_Initialized(&flag);
        if (!flag)
        {
            int provided = 0;
            MPI_Init_thread
            (
                nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided
            );
            ownMpi_ = true;
        }
        parent = MPI_COMM_WORLD;
    }
    else if (mode_ == "inTransit")
    {
        if (UPstream::allWorlds().size() < 2)
        {
            FatalErrorInFunction
                << "mode inTransit requires a separate render world:" << nl
                << "    launch the solver with -world <name> and"
                << " anariVisRender with -world <renderName>" << nl
                << exit(FatalError);
        }
        parent = ::anariVis::toMPIComm
        (
            UPstream::Communicator::lookup(UPstream::commGlobal())
        );
    }
    else
    {
        parent = ::anariVis::toMPIComm
        (
            UPstream::Communicator::lookup(UPstream::worldComm)
        );
    }

    if (mode_ == "inTransit" && !UPstream::parRun())
    {
        FatalErrorInFunction
            << "mode inTransit requires a parallel run" << nl
            << exit(FatalError);
    }

    bool consumer = false;
    if (mode_ == "inSitu")
    {
        consumer =
            renderRanksPerNode_ < 0
         || ::anariVis::Transport::nodeLocalRank(parent) < renderRanksPerNode_;
    }

    transport_.reset(new ::anariVis::Transport(parent, true, consumer));

    const int nConsumers = transport_->nConsumers();

    bounds_ = boundBox(mesh_.points(), true);

    // In-transit: hand the settings to the render world
    if (mode_ == "inTransit")
    {
        dictionary cfg(dict_);
        cfg.set("outputDir", fileName(settings_.outputDir));
        cfg.set("prefix", word(settings_.prefix));
        OStringStream os;
        cfg.write(os, false);
        transport_->sendConfig(os.str());
    }

    // Volume sampling and brick routing
    forAll(volumes_, i)
    {
        volumes_[i].setup(nConsumers);
    }

    // Surface routing
    const std::string compMode =
        ::anariVis::compositingMode(settings_, nConsumers);
    spatialRouting_ =
        (surfaceRouting_ == "spatial")
     || (
            surfaceRouting_ == "auto"
         && nConsumers > 1
         && (compMode == "deep" || compMode == "ptc")
        );

    std::set<int> spatialDest;
    if (spatialRouting_)
    {
        if (volumes_.size())
        {
            // Same bricks as the (first) volume
            sceneGrid_ = volumes_[0].grid();
        }
        else
        {
            const vector span = bounds_.span();
            const scalar h = cmptMax(span)/1023;
            for (direction a = 0; a < 3; ++a)
            {
                sceneGrid_.dims[a] = int32_t
                (
                    max(label(2), label(std::ceil(span[a]/h - 1e-6)) + 1)
                );
                sceneGrid_.origin[a] = float(bounds_.min()[a]);
                sceneGrid_.spacing[a] = float(h);
            }
        }
        sceneDecomp_ =
            ::anariVis::BrickDecomposition(sceneGrid_.dims, nConsumers);

        // Static candidate bricks: those overlapping the local mesh bounds
        boundBox lbb(mesh_.points(), false);
        if (lbb.good())
        {
            lbb.inflate(1e-3);
            for (int b = 0; b < nConsumers; ++b)
            {
                const auto bb = sceneDecomp_.box(b, sceneGrid_);
                const boundBox bbox
                (
                    point(bb[0], bb[1], bb[2]),
                    point(bb[3], bb[4], bb[5])
                );
                if (bbox.overlaps(lbb))
                {
                    spatialDest.insert(b);
                }
            }
        }
    }

    // Destinations of every frame
    {
        std::set<int> dest(spatialDest);
        if (!spatialRouting_)
        {
            dest.insert(transport_->surfaceConsumer());
        }
        forAll(volumes_, i)
        {
            for (const int b : volumes_[i].destinations())
            {
                dest.insert(b);
            }
        }
        if (transport_->isLeader())
        {
            for (int c = 0; c < nConsumers; ++c)
            {
                dest.insert(c);
            }
        }
        destinations_.assign(dest.begin(), dest.end());
    }
    transport_->setDestinations(destinations_);

    // In-situ renderer on the render ranks
    if (transport_->isConsumer())
    {
        settings_.widenAffinity = widenAffinity_;
        if (settings_.threads <= 0)
        {
            const int cores =
                widenAffinity_
              ? int(sysconf(_SC_NPROCESSORS_ONLN))
              : ::anariVis::availableCores();
            settings_.threads =
                std::max(1, cores/transport_->nConsumersOnNode());
        }
        renderer_.reset
        (
            new ::anariVis::Renderer(settings_, transport_->renderComm())
        );

        #ifdef ANARIVIS_NEXARIS
        if (viewerEnabled_)
        {
            if (viewerOptions_.title.empty())
            {
                viewerOptions_.title =
                    std::string(time_.globalCaseName()) + " (in-situ, "
                  + std::to_string(UPstream::nProcs()) + " ranks, "
                  + std::to_string(nConsumers) + " render ranks)";
            }
            viewer_.reset
            (
                new ::anariVis::LiveViewer
                (
                    *renderer_, settings_, transport_->renderComm(),
                    viewerOptions_
                )
            );
            if (!viewer_->ok())
            {
                viewer_.reset(nullptr);
            }
        }
        #endif

        if (renderer_->rank() == 0)
        {
            int provided = 0;
            MPI_Query_thread(&provided);
            if
            (
                renderer_->mode() == "barney"
             && provided < MPI_THREAD_MULTIPLE
            )
            {
                WarningInFunction
                    << "Barney 'mpi' device with MPI thread level "
                    << provided << " < MPI_THREAD_MULTIPLE;"
                    << " consider running with -mpi-threads" << endl;
            }
        }
    }

    // Summary
    {
        int haveRenderer = renderer_ ? 1 : 0;
        int r0 = (renderer_ && renderer_->rank() == 0) ? 1 : 0;
        std::string desc;
        if (r0)
        {
            desc = renderer_->libraryName() + "/" + renderer_->mode();
        }
        reduce(haveRenderer, sumOp<int>());
        Info<< type() << ' ' << name() << ": mode " << mode_
            << ", producers " << transport_->nProducers()
            << ", render ranks " << nConsumers;
        if (mode_ == "inSitu")
        {
            Info<< ", renderer " << desc.c_str()
                << ", threads/render rank " << settings_.threads;
        }
        Info<< ", surface routing "
            << (spatialRouting_ ? "spatial" : "node")
            << ", surfaces " << surfaces_.size()
            << ", volumes " << volumes_.size()
            << ", images " << settings_.outputDir.c_str()
            << ", setup " << seconds(tInit) << " s" << endl;
    }

    initialised_ = true;
}


void Foam::functionObjects::anariVis::extractSurfaces
(
    std::vector<::anariVis::SurfacePiece>& pieces
)
{
    pieces.clear();
    dataRanges_.clear();

    forAll(surfaces_, i)
    {
        sampledSurface& surf = surfaces_[i];
        surf.update();

        ::anariVis::SurfacePiece piece;
        piece.id = int32_t(i);

        const pointField& pts = surf.points();
        const faceList& faces = surf.faces();

        piece.points.resize(3*pts.size());
        forAll(pts, pi)
        {
            piece.points[3*pi] = float(pts[pi].x());
            piece.points[3*pi + 1] = float(pts[pi].y());
            piece.points[3*pi + 2] = float(pts[pi].z());
        }

        // Fan triangulation (sampled surface faces are convex polygons);
        // keeps the face orientation
        label nTri = 0;
        for (const face& f : faces)
        {
            nTri += max(f.size() - 2, label(0));
        }
        piece.triangles.reserve(3*nTri);
        for (const face& f : faces)
        {
            for (label k = 1; k + 1 < f.size(); ++k)
            {
                piece.triangles.push_back(uint32_t(f[0]));
                piece.triangles.push_back(uint32_t(f[k]));
                piece.triangles.push_back(uint32_t(f[k + 1]));
            }
        }

        // Colour field interpolated to the surface points
        const ::anariVis::SurfaceStyle& style = settings_.surfaces[i];
        scalar minV = GREAT;
        scalar maxV = -GREAT;

        if (!colourFields_[i].empty())
        {
            tmp<volScalarField> tfld =
                anariVisVolumeSampler::lookupScalar(mesh_, colourFields_[i]);

            if (tfld.valid())
            {
                autoPtr<interpolation<scalar>> interp =
                    interpolation<scalar>::New(interpolationScheme_, tfld());

                tmp<scalarField> tvals = surf.interpolate(*interp);
                const scalarField& vals = tvals();

                if (vals.size() == pts.size())
                {
                    piece.scalars.resize(vals.size());
                    forAll(vals, pi)
                    {
                        piece.scalars[pi] = float(vals[pi]);
                        minV = min(minV, vals[pi]);
                        maxV = max(maxV, vals[pi]);
                    }
                }
            }
            else if (nFrames_ <= 1)
            {
                WarningInFunction
                    << "Colour field " << colourFields_[i]
                    << " not found for surface " << surf.name() << endl;
            }
        }

        reduce(minV, minOp<scalar>());
        reduce(maxV, maxOp<scalar>());

        if (maxV >= minV)
        {
            dataRanges_.push_back({{float(minV), float(maxV)}});
        }
        else
        {
            dataRanges_.push_back({{0.f, 0.f}});
        }

        if (style.autoRange)
        {
            if (maxV < minV)
            {
                minV = 0;
                maxV = 1;
            }
            piece.range = {{float(minV), float(maxV)}};
        }
        else
        {
            piece.range = style.range;
        }

        pieces.push_back(std::move(piece));
    }
}


void Foam::functionObjects::anariVis::sendStop()
{
    if (ended_ || !transport_ || mode_ != "inTransit")
    {
        ended_ = true;
        return;
    }

    transport_->waitSends();

    ::anariVis::FramePacket pkt;
    pkt.header.frameIndex = int32_t(time_.timeIndex());
    pkt.header.time = time_.value();
    pkt.header.flags = ::anariVis::FrameHeader::stopFlag;
    pkt.header.producer = transport_->visRank();

    for (const int c : destinations_)
    {
        std::vector<char> buf;
        pkt.pack(buf);
        transport_->post(c, std::move(buf));
    }
    transport_->waitSends();

    ended_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::anariVis::anariVis
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dict_(),
    settings_(),
    mode_("inSitu"),
    renderRanksPerNode_(1),
    synchronise_(true),
    widenAffinity_(true),
    dropPolicy_("drop"),
    interpolationScheme_("cellPoint"),
    surfaceRouting_("auto"),
    spatialRouting_(false),
    sceneGrid_(),
    sceneDecomp_(),
    surfaces_(),
    colourFields_(),
    volumes_(),
    transport_(nullptr),
    renderer_(nullptr),
    viewerEnabled_(false),
    destinations_(),
    bounds_(),
    initialised_(false),
    ended_(false),
    ownMpi_(false),
    nFrames_(0),
    nDropped_(0),
    sumWait_(0),
    sumExtract_(0),
    sumSend_(0),
    sumRender_(0),
    sumTotal_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::anariVis::~anariVis()
{
    int finalized = 0;
    MPI_Finalized(&finalized);

    if (!finalized)
    {
        sendStop();
    }

    #ifdef ANARIVIS_NEXARIS
    viewer_.reset(nullptr);
    #endif
    renderer_.reset(nullptr);
    transport_.reset(nullptr);

    if (ownMpi_ && !finalized)
    {
        MPI_Finalize();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::anariVis::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (initialised_)
    {
        WarningInFunction
            << "Changes of the " << type() << " settings are ignored after"
            << " the first frame" << endl;
        return true;
    }

    dict_ = dict;

    mode_ = dict.getOrDefault<word>("mode", "inSitu");
    if (mode_ != "inSitu" && mode_ != "inTransit")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown mode " << mode_ << ", use inSitu or inTransit" << nl
            << exit(FatalIOError);
    }
    renderRanksPerNode_ = dict.getOrDefault<label>("renderRanksPerNode", 1);
    synchronise_ = dict.getOrDefault<bool>("synchronise", true);
    widenAffinity_ = dict.getOrDefault<bool>("renderAffinityAll", true);
    dropPolicy_ = dict.getOrDefault<word>("dropPolicy", "drop");
    interpolationScheme_ =
        dict.getOrDefault<word>("interpolationScheme", "cellPoint");
    surfaceRouting_ = dict.getOrDefault<word>("surfaceRouting", "auto");

    settings_ = anariVisIO::readSettings
    (
        dict,
        time_.globalPath()/functionObject::outputPrefix/name(),
        name()
    );

    // Interactive viewing in-situ (in-transit: anariVisRender -viewer)
    viewerEnabled_ = false;
    if (const dictionary* vd = dict.findDict("viewer"))
    {
        viewerEnabled_ = vd->getOrDefault<bool>("enabled", true);
        if (viewerEnabled_ && mode_ != "inSitu")
        {
            Info<< type() << ' ' << name() << ": the viewer dictionary"
                << " applies in-situ; in-transit use anariVisRender -viewer"
                << endl;
            viewerEnabled_ = false;
        }
        #ifdef ANARIVIS_NEXARIS
        if (viewerEnabled_)
        {
            ::anariVis::ViewerOptions& o = viewerOptions_;
            o.inSitu = true;
            o.port = vd->getOrDefault<label>("port", 47800);
            o.historySize = vd->getOrDefault<label>("history", 4);
            o.historyMB = vd->getOrDefault<scalar>("historyMB", 4096);
            o.passesPerFrame = vd->getOrDefault<label>("passesPerFrame", 2);
            o.maxSecondsPerFrame =
                vd->getOrDefault<scalar>("maxSecondsPerFrame", 2);
            o.keepAlive = vd->getOrDefault<scalar>("keepAlive", 120);
            o.pixelSamples = vd->getOrDefault<label>("pixelSamples", 1);
            o.upAxis =
                vd->getOrDefault<word>("upAxis", "z") == "y" ? 'y' : 'z';
            o.title = vd->getOrDefault<string>("title", string::null);
            // frames stay in memory
            settings_.writeImages = false;
        }
        #else
        if (viewerEnabled_)
        {
            WarningInFunction
                << "anariVis was built without Nexaris (NEXARIS_LIVE_ROOT):"
                << " the viewer is ignored" << endl;
            viewerEnabled_ = false;
        }
        #endif
    }

    // Surfaces (sampledSurface selection), same order as the styles
    surfaces_.clear();
    colourFields_.clear();
    if (const dictionary* sdicts = dict.findDict("surfaces"))
    {
        for (const entry& e : *sdicts)
        {
            if (!e.isDict())
            {
                continue;
            }
            surfaces_.push_back
            (
                sampledSurface::New(e.keyword(), mesh_, e.dict())
            );
            colourFields_.push_back
            (
                e.dict().getOrDefault<word>("colourField", word::null)
            );
        }
    }

    // Volumes
    volumes_.clear();
    if (const dictionary* vdicts = dict.findDict("volumes"))
    {
        for (const entry& e : *vdicts)
        {
            if (!e.isDict())
            {
                continue;
            }
            volumes_.push_back
            (
                new anariVisVolumeSampler(mesh_, e.keyword(), e.dict())
            );
        }
    }

    return true;
}


bool Foam::functionObjects::anariVis::execute()
{
    if (ended_)
    {
        return true;
    }

    initialise();

    const auto tStart = clock_type::now();

    // Bounded in-flight data: at most one frame per producer
    if (mode_ == "inTransit")
    {
        bool done = transport_->sendsComplete();
        if (dropPolicy_ == "block")
        {
            if (!done)
            {
                transport_->waitSends();
            }
        }
        else if (!returnReduceAnd(done))
        {
            ++nDropped_;
            Info<< type() << ' ' << name() << ": renderer busy, frame "
                << time_.timeIndex() << " dropped (" << nDropped_
                << " dropped)" << endl;
            return true;
        }
    }
    else
    {
        transport_->waitSends();
    }
    const double tWait = seconds(tStart);
    const auto tExtractStart = clock_type::now();

    ++nFrames_;

    // Extraction ------------------------------------------------------------
    std::vector<::anariVis::SurfacePiece> surfPieces;
    extractSurfaces(surfPieces);

    List<scalarField> volValues(volumes_.size());
    List<bool> volValid(volumes_.size(), false);
    std::vector<::anariVis::vec2f> volRange(volumes_.size());

    forAll(volumes_, i)
    {
        volValid[i] = volumes_[i].sample(volValues[i]);

        scalar minV = GREAT;
        scalar maxV = -GREAT;
        for (const scalar v : volValues[i])
        {
            minV = min(minV, v);
            maxV = max(maxV, v);
        }
        reduce(minV, minOp<scalar>());
        reduce(maxV, maxOp<scalar>());
        dataRanges_.push_back
        (
            maxV >= minV
          ? ::anariVis::vec2f{{float(minV), float(maxV)}}
          : ::anariVis::vec2f{{0.f, 0.f}}
        );
        if (maxV < minV)
        {
            minV = 0;
            maxV = 1;
        }
        const ::anariVis::VolumeStyle& style = settings_.volumes[i];
        volRange[i] =
            style.autoRange
          ? ::anariVis::vec2f{{float(minV), float(maxV)}}
          : style.range;

        if (!volValid[i] && nFrames_ <= 1)
        {
            WarningInFunction
                << "Field " << volumes_[i].fieldName()
                << " not found for volume " << volumes_[i].name() << endl;
        }
    }

    const double tExtract = seconds(tExtractStart);

    if (settings_.verbose || nFrames_ == 1)
    {
        // Data ranges of the colour fields (surfaces) and volume fields
        Info<< type() << ' ' << name() << ": data ranges";
        forAll(surfaces_, i)
        {
            Info<< ' ' << surfaces_[i].name() << " ["
                << dataRanges_[i][0] << ' ' << dataRanges_[i][1] << ']';
        }
        forAll(volumes_, i)
        {
            const auto& r = dataRanges_[surfaces_.size() + i];
            Info<< ' ' << volumes_[i].name() << " ["
                << r[0] << ' ' << r[1] << ']';
        }
        Info<< endl;
    }

    // Packing and sending ---------------------------------------------------
    const auto tSendStart = clock_type::now();
    {
        ::anariVis::FrameHeader header;
        header.frameIndex = int32_t(time_.timeIndex());
        header.time = time_.value();
        header.producer = transport_->visRank();
        for (direction a = 0; a < 3; ++a)
        {
            header.bounds[a] = float(bounds_.min()[a]);
            header.bounds[3 + a] = float(bounds_.max()[a]);
        }

        const int surfConsumer = transport_->surfaceConsumer();

        // Spatial routing: split the local surfaces by brick
        std::map<int, std::vector<::anariVis::SurfacePiece>> routed;
        if (spatialRouting_)
        {
            std::vector<std::pair<int, ::anariVis::SurfacePiece>> parts;
            label nLost = 0;
            for (const auto& sp : surfPieces)
            {
                ::anariVis::splitByBrick(sp, sceneDecomp_, sceneGrid_, parts);
                for (auto& bp : parts)
                {
                    if
                    (
                        std::binary_search
                        (
                            destinations_.begin(), destinations_.end(),
                            bp.first
                        )
                    )
                    {
                        routed[bp.first].push_back(std::move(bp.second));
                    }
                    else
                    {
                        nLost += label(bp.second.nTriangles());
                    }
                }
            }
            if (nLost)
            {
                WarningInFunction
                    << nLost << " triangles outside the candidate bricks"
                    << endl;
            }
        }

        for (const int c : destinations_)
        {
            ::anariVis::FramePacket pkt;
            pkt.header = header;

            if (spatialRouting_)
            {
                pkt.header.hasBrick = 1;
                pkt.header.brick = sceneDecomp_.box(c, sceneGrid_);
                auto it = routed.find(c);
                if (it != routed.end())
                {
                    pkt.surfaces = std::move(it->second);
                }
            }
            else if (c == surfConsumer)
            {
                pkt.surfaces = std::move(surfPieces);
            }

            forAll(volumes_, i)
            {
                if (!volValid[i])
                {
                    continue;
                }
                const auto& d = volumes_[i].destinations();
                if (std::find(d.begin(), d.end(), c) == d.end())
                {
                    continue;
                }
                ::anariVis::VolumePiece vp;
                vp.id = int32_t(i);
                vp.range = volRange[i];
                volumes_[i].fillPiece(c, volValues[i], vp);
                pkt.volumes.push_back(std::move(vp));
            }

            std::vector<char> buf;
            pkt.pack(buf);
            transport_->post(c, std::move(buf));
        }
    }
    const double tSend = seconds(tSendStart);

    // In-situ rendering -----------------------------------------------------
    ::anariVis::FrameTimings ft;
    double tRender = 0;
    if (renderer_)
    {
        const auto tr = clock_type::now();
        std::vector<std::vector<char>> packets;
        int frameIndex = 0;
        if (transport_->receiveFrame(packets, frameIndex))
        {
            ::anariVis::LocalScene scene;
            ::anariVis::assembleScene
            (
                packets,
                transport_->consumerIndex(),
                transport_->nConsumers(),
                scene
            );
            #ifdef ANARIVIS_NEXARIS
            if (viewer_)
            {
                // the solver waits while the client sees the new frame
                viewer_->addFrame(std::move(scene), seconds(tr));
                viewer_->runInSituFrame();
                if (viewer_->quit())
                {
                    viewer_.reset(nullptr);
                }
            }
            else
            #endif
            {
                ft = renderer_->render(scene);
            }
        }
        tRender = seconds(tr);
    }

    if (mode_ == "inSitu" && synchronise_)
    {
        transport_->sleepBarrier();
    }

    const double tTotal = seconds(tStart);

    // Statistics ------------------------------------------------------------
    scalar maxExtract = tExtract;
    scalar maxSend = tSend;
    scalar maxWait = tWait;
    reduce(maxExtract, maxOp<scalar>());
    reduce(maxSend, maxOp<scalar>());
    reduce(maxWait, maxOp<scalar>());

    sumWait_ += maxWait;
    sumExtract_ += maxExtract;
    sumSend_ += maxSend;
    sumRender_ += tRender;
    sumTotal_ += tTotal;

    Info<< type() << ' ' << name() << ": frame " << time_.timeIndex()
        << " t=" << time_.value()
        << " wait " << maxWait
        << " s, extract " << maxExtract << " s, send " << maxSend << " s";
    if (renderer_ && renderer_->rank() == 0)
    {
        Info<< ", render " << tRender << " s (build " << ft.build
            << " render " << ft.render << " composite " << ft.composite
            << " write " << ft.write << ", local triangles "
            << label(ft.nTriangles) << ")";
    }
    Info<< ", total " << tTotal << " s" << endl;

    return true;
}


bool Foam::functionObjects::anariVis::write()
{
    return true;
}


bool Foam::functionObjects::anariVis::end()
{
    sendStop();

    #ifdef ANARIVIS_NEXARIS
    if (viewer_)
    {
        // the run is over: the client can still browse the frames in memory
        viewer_->setSolverFinished();
        viewer_->serveUntilQuit();
        viewer_.reset(nullptr);
    }
    if (initialised_ && mode_ == "inSitu" && viewerEnabled_)
    {
        // the other ranks wait without spinning
        transport_->sleepBarrier();
    }
    #endif

    if (nFrames_)
    {
        Info<< type() << ' ' << name() << ": " << nFrames_ << " frames, "
            << nDropped_ << " dropped; mean per frame: wait "
            << sumWait_/nFrames_ << " s, extract "
            << sumExtract_/nFrames_ << " s, send " << sumSend_/nFrames_
            << " s, render (master) " << sumRender_/nFrames_
            << " s, total (master) " << sumTotal_/nFrames_ << " s" << endl;
    }
    return true;
}


void Foam::functionObjects::anariVis::updateMesh(const mapPolyMesh&)
{
    WarningInFunction
        << "Topology changes are not supported by " << type() << endl;
}


void Foam::functionObjects::anariVis::movePoints(const polyMesh&)
{
    if (initialised_ && volumes_.size())
    {
        WarningInFunction
            << "Mesh motion is not supported by the volume sampling of "
            << type() << endl;
    }
}


// ************************************************************************* //
