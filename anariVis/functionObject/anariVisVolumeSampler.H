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
    Foam::anariVisVolumeSampler

Description
    Resamples a cell field of the local (sub)domain onto the points of a
    global structured grid that covers the global mesh bounds.

    Every rank samples only the grid points that lie inside its own cells
    (no global gather). The global grid is split into convex bricks, one per
    render rank (anariVis::BrickDecomposition, neighbouring bricks share one
    point plane); each sample is routed to the brick(s) that contain it.
    The sample locations, containing cells and interpolation weights are
    computed once (static mesh) and reused every frame.

    Dictionary entries:
    \table
        Property    | Description                           | Required | Default
        field       | volScalarField or volVectorField (magnitude) | yes |
        resolution  | grid points along the longest axis    | no  | 128
        bounds      | sub-box (min max)                     | no  | mesh bounds
        interpolationScheme | cellPoint or cell             | no  | cellPoint
    \endtable

SourceFiles
    anariVisVolumeSampler.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_anariVisVolumeSampler_H
#define Foam_anariVisVolumeSampler_H

#include "fvMesh.H"
#include "volFields.H"
#include "boundBox.H"
#include "cellPointWeight.H"
#include "anariVisTypes.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class anariVisVolumeSampler Declaration
\*---------------------------------------------------------------------------*/

class anariVisVolumeSampler
{
    // Private Data

        const fvMesh& mesh_;

        word name_;

        word fieldName_;

        word scheme_;

        label resolution_;

        bool hasBounds_;

        boundBox bounds_;

        ::anariVis::VolumeGrid grid_;

        //- Interpolation weights per local sample (cellPoint)
        std::vector<cellPointWeight> weights_;

        //- Containing cell per local sample
        labelList sampleCells_;

        //- Global grid point index per local sample
        std::vector<uint32_t> sampleIndex_;

        //- Local sample ids per consumer brick
        std::vector<std::vector<label>> destSamples_;

        //- Consumers that receive samples from this rank
        std::vector<int> destinations_;


    // Private Member Functions

        //- Inclusive point-in-cell test with a relative tolerance
        bool insideCell(const point& p, const label celli) const;


public:

    // Constructors

        //- Construct from dictionary
        anariVisVolumeSampler
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& dict
        );


    // Static Functions

        //- Look up a volScalarField, or the magnitude of a volVectorField.
        //  Returns an invalid tmp if neither exists.
        static tmp<volScalarField> lookupScalar
        (
            const fvMesh& mesh,
            const word& fieldName
        );


    // Member Functions

        const word& name() const noexcept { return name_; }

        const word& fieldName() const noexcept { return fieldName_; }

        const ::anariVis::VolumeGrid& grid() const noexcept { return grid_; }

        label nSamples() const noexcept { return sampleCells_.size(); }

        //- Build grid, local samples and brick routing. Collective.
        void setup(const label nConsumers);

        //- Consumers that receive samples from this rank
        const std::vector<int>& destinations() const noexcept
        {
            return destinations_;
        }

        //- Interpolate the field at the local samples.
        //  Returns false if the field does not exist. Collective.
        bool sample(scalarField& values) const;

        //- Fill the packet for consumer brick b
        void fillPiece
        (
            const int b,
            const scalarField& values,
            ::anariVis::VolumePiece& piece
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
