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

#include "anariVisSettingsIO.H"
#include "vector.H"
#include "FixedList.H"
#include "ITstream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

::anariVis::vec3f toVec3f(const Foam::vector& v)
{
    return {{float(v.x()), float(v.y()), float(v.z())}};
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

std::vector<::anariVis::ColourPoint> Foam::anariVisIO::readColourMap
(
    const dictionary& dict,
    const word& key,
    const word& defaultName
)
{
    if (!dict.found(key))
    {
        return ::anariVis::namedColourMap(defaultName);
    }

    ITstream& is = dict.lookup(key);
    if (is.peek().isWord())
    {
        const word name(is);
        return ::anariVis::namedColourMap(name);
    }

    const List<FixedList<scalar, 4>> pts(is);
    std::vector<::anariVis::ColourPoint> map;
    for (const auto& p : pts)
    {
        map.push_back
        (
            {float(p[0]), {{float(p[1]), float(p[2]), float(p[3])}}}
        );
    }
    return map;
}


::anariVis::RenderSettings Foam::anariVisIO::readSettings
(
    const dictionary& dict,
    const fileName& defaultOutputDir,
    const word& defaultPrefix
)
{
    ::anariVis::RenderSettings s;

    s.library = dict.getOrDefault<word>("library", "barney_cpu");
    s.wrappedLibrary = dict.getOrDefault<word>("wrappedLibrary", "helide");
    s.compositing = dict.getOrDefault<word>("compositing", "auto");
    s.rendererType = dict.getOrDefault<word>("renderer", "default");

    if (dict.found("resolution"))
    {
        const FixedList<label, 2> res(dict.get<FixedList<label, 2>>("resolution"));
        s.width = int(res[0]);
        s.height = int(res[1]);
    }
    s.pixelSamples = int(dict.getOrDefault<label>("pixelSamples", 4));
    s.passes = int(dict.getOrDefault<label>("passes", 1));

    if (dict.found("background"))
    {
        const vector bg(dict.get<vector>("background"));
        s.background = {{float(bg.x()), float(bg.y()), float(bg.z()), 1.f}};
    }
    s.ambientRadiance = float(dict.getOrDefault<scalar>("ambientRadiance", 0.3));
    s.denoise = dict.getOrDefault<bool>("denoise", false);
    s.threads = int(dict.getOrDefault<label>("threads", 0));
    s.depthMode = dict.getOrDefault<word>("depthMode", "auto");

    s.outputDir = dict.getOrDefault<fileName>("outputDir", defaultOutputDir);
    s.outputDir = fileName(s.outputDir).expand();
    s.prefix = dict.getOrDefault<word>("prefix", defaultPrefix);
    s.writeImages = dict.getOrDefault<bool>("writeImages", true);
    s.verbose = dict.getOrDefault<bool>("verbose", false);

    // Camera
    if (const dictionary* cdict = dict.findDict("camera"))
    {
        ::anariVis::CameraSettings& c = s.camera;
        vector v;
        if (cdict->readIfPresent("position", v))
        {
            c.position = toVec3f(v);
            c.autoPosition = false;
        }
        if (cdict->readIfPresent("lookAt", v))
        {
            c.lookAt = toVec3f(v);
            c.autoLookAt = false;
        }
        if (cdict->readIfPresent("up", v))
        {
            c.up = toVec3f(v);
        }
        if (cdict->readIfPresent("viewDirection", v))
        {
            c.viewDirection = toVec3f(v);
        }
        c.fovy = float(cdict->getOrDefault<scalar>("fov", c.fovy));
        c.zoom = float(cdict->getOrDefault<scalar>("zoom", c.zoom));
    }

    // Lights
    if (const dictionary* ldicts = dict.findDict("lights"))
    {
        for (const entry& e : *ldicts)
        {
            if (!e.isDict())
            {
                continue;
            }
            const dictionary& ld = e.dict();
            ::anariVis::LightSettings l;
            l.type = ld.getOrDefault<word>("type", "directional");
            vector v;
            if (ld.readIfPresent("direction", v)) l.direction = toVec3f(v);
            if (ld.readIfPresent("position", v)) l.position = toVec3f(v);
            if (ld.readIfPresent("colour", v)) l.colour = toVec3f(v);
            l.intensity = float(ld.getOrDefault<scalar>("intensity", 3));
            l.cameraRelative = ld.getOrDefault<bool>("cameraRelative", false);
            s.lights.push_back(l);
        }
    }

    // Surface styles, in dictionary order
    if (const dictionary* sdicts = dict.findDict("surfaces"))
    {
        for (const entry& e : *sdicts)
        {
            if (!e.isDict())
            {
                continue;
            }
            const dictionary& sd = e.dict();
            ::anariVis::SurfaceStyle st;
            st.name = e.keyword();
            st.useField = sd.found("colourField");
            vector v;
            if (sd.readIfPresent("colour", v)) st.colour = toVec3f(v);
            st.colourMap = readColourMap(sd, "colourMap", "coolToWarm");
            if (sd.found("range"))
            {
                const FixedList<scalar, 2> r(sd.get<FixedList<scalar, 2>>("range"));
                st.range = {{float(r[0]), float(r[1])}};
                st.autoRange = false;
            }
            st.opacity = float(sd.getOrDefault<scalar>("opacity", 1));
            st.smoothNormals = sd.getOrDefault<bool>("smoothNormals", true);
            s.surfaces.push_back(st);
        }
    }

    // Volume styles, in dictionary order
    if (const dictionary* vdicts = dict.findDict("volumes"))
    {
        for (const entry& e : *vdicts)
        {
            if (!e.isDict())
            {
                continue;
            }
            const dictionary& vd = e.dict();
            ::anariVis::VolumeStyle st;
            st.name = e.keyword();
            st.colourMap = readColourMap(vd, "colourMap", "coolToWarm");
            if (vd.found("opacity"))
            {
                const List<FixedList<scalar, 2>> op
                (
                    vd.get<List<FixedList<scalar, 2>>>("opacity")
                );
                for (const auto& p : op)
                {
                    st.opacityMap.push_back({float(p[0]), float(p[1])});
                }
            }
            else
            {
                st.opacityMap = {{0.f, 0.f}, {1.f, 0.8f}};
            }
            if (vd.found("range"))
            {
                const FixedList<scalar, 2> r(vd.get<FixedList<scalar, 2>>("range"));
                st.range = {{float(r[0]), float(r[1])}};
                st.autoRange = false;
            }
            st.unitDistance = float(vd.getOrDefault<scalar>("unitDistance", 0));
            s.volumes.push_back(st);
        }
    }

    return s;
}


// ************************************************************************* //
