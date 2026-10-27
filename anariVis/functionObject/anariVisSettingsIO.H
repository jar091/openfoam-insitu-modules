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
    Foam::anariVisIO

Description
    Conversion of the anariVis function object dictionary into the
    renderer-agnostic anariVis::RenderSettings. Used both by the function
    object (in-situ) and by the anariVisRender application (in-transit),
    which receives the dictionary text from the simulation.

    Dictionary entries (all optional):
    \table
        Property     | Description                          | Default
        library      | barney_cpu, barney, helide, cycles, ptc | barney_cpu
        wrappedLibrary | library wrapped by ptc             | helide
        compositing  | auto, none, barney, deep, ptc        | auto
        renderer     | ANARI renderer subtype               | default
        resolution   | image size (width height)            | (1024 768)
        pixelSamples | samples per pixel per pass           | 4
        passes       | accumulated passes per image         | 1
        background   | background colour (r g b)            | (0.12 0.12 0.15)
        ambientRadiance | ambient light                     | 0.3
        denoise      | device denoiser                      | false
        threads      | render threads per render rank (0: auto) | 0
        depthMode    | deep compositing depth: auto, device, brick | auto
        outputDir    | image directory                      | \<case\>/postProcessing/\<name\>
        prefix       | image file prefix                    | \<name\>
        writeImages  | write PNG images                     | true
        camera       | dictionary, see below                |
        lights       | dictionary of light dictionaries     | camera key light
        surfaces     | dictionary of surface dictionaries   |
        volumes      | dictionary of volume dictionaries    |
    \endtable

    camera: position, lookAt, up, viewDirection, fov [deg], zoom.
    If position is missing the camera is placed along viewDirection so
    that the scene bounds fit into the image.

    lights: type (directional|point), direction, position, colour,
    intensity, cameraRelative.

    Per surface (besides the sampledSurface entries): colourField,
    colour, colourMap (name or list of (x r g b)), range (min max),
    opacity, smoothNormals.

    Per volume: field, resolution, colourMap, opacity (list of (x alpha)),
    range (min max), unitDistance.

SourceFiles
    anariVisSettingsIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_anariVisSettingsIO_H
#define Foam_anariVisSettingsIO_H

#include "dictionary.H"
#include "fileName.H"
#include "anariVisTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace anariVisIO
{

//- Read the render settings from the function object dictionary
::anariVis::RenderSettings readSettings
(
    const dictionary& dict,
    const fileName& defaultOutputDir,
    const word& defaultPrefix
);

//- Name of a colour map entry ("" for a list of control points)
std::string colourMapName
(
    const dictionary& dict,
    const word& key,
    const word& defaultName
);

//- Read a colour map entry (name or list of (x r g b))
std::vector<::anariVis::ColourPoint> readColourMap
(
    const dictionary& dict,
    const word& key,
    const word& defaultName
);

} // End namespace anariVisIO
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
