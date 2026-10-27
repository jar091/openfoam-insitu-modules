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

#include "anariVisVolumeSampler.H"
#include "interpolationCellPoint.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anariVisVolumeSampler::anariVisVolumeSampler
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& dict
)
:
    mesh_(mesh),
    name_(name),
    fieldName_(dict.get<word>("field")),
    scheme_(dict.getOrDefault<word>("interpolationScheme", "cellPoint")),
    resolution_(max(dict.getOrDefault<label>("resolution", 128), 2)),
    hasBounds_(dict.found("bounds")),
    bounds_()
{
    if (hasBounds_)
    {
        bounds_ = dict.get<boundBox>("bounds");
    }
    if (scheme_ != "cellPoint" && scheme_ != "cell")
    {
        FatalIOErrorInFunction(dict)
            << "interpolationScheme " << scheme_
            << " not supported, use cellPoint or cell" << nl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::anariVisVolumeSampler::lookupScalar
(
    const fvMesh& mesh,
    const word& fieldName
)
{
    if (const auto* sf = mesh.cfindObject<volScalarField>(fieldName))
    {
        return tmp<volScalarField>(*sf);
    }
    if (const auto* vf = mesh.cfindObject<volVectorField>(fieldName))
    {
        return mag(*vf);
    }
    return tmp<volScalarField>();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::anariVisVolumeSampler::insideCell
(
    const point& p,
    const label celli
) const
{
    const cell& cFaces = mesh_.cells()[celli];
    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Cf = mesh_.faceCentres();
    const labelUList& own = mesh_.faceOwner();

    // Tolerance relative to the cell size
    const scalar tol = 1e-6*Foam::cbrt(mag(mesh_.cellVolumes()[celli]));

    for (const label facei : cFaces)
    {
        vector n = Sf[facei];
        if (own[facei] != celli)
        {
            n = -n;
        }
        const scalar magN = mag(n);
        if (((p - Cf[facei]) & n) > tol*magN)
        {
            return false;
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::anariVisVolumeSampler::setup(const label nConsumers)
{
    // Global grid
    boundBox bb = hasBounds_ ? bounds_ : boundBox(mesh_.points(), true);
    const vector span = bb.span();
    const scalar maxL = cmptMax(span);
    const scalar h = maxL/scalar(resolution_ - 1);

    for (direction a = 0; a < 3; ++a)
    {
        grid_.dims[a] =
            int32_t(max(label(2), label(std::ceil(span[a]/h - 1e-6)) + 1));
        grid_.origin[a] = float(bb.min()[a]);
        grid_.spacing[a] = float(h);
    }

    if (grid_.nPoints() > uint64_t(std::numeric_limits<uint32_t>::max()))
    {
        FatalErrorInFunction
            << "Volume grid of " << name_ << " too large" << exit(FatalError);
    }

    // Local samples: loop over cells and test the grid points inside the
    // cell bounding box. First cell to claim a point wins.
    const boundBox localBb(mesh_.points(), false);

    // cellPointWeight uses the lazily built polyMesh::tetBasePtIs(), whose
    // construction is collective (syncTools face sync). Build it on all
    // ranks here: ranks without any sample point would otherwise never
    // enter it and the other ranks would deadlock.
    if (scheme_ == "cellPoint")
    {
        (void)mesh_.tetBasePtIs();
    }

    weights_.clear();
    DynamicList<label> cells;
    sampleIndex_.clear();

    if (mesh_.nCells() && localBb.good())
    {
        label lo[3], n[3];
        for (direction a = 0; a < 3; ++a)
        {
            lo[a] = max
            (
                label(0),
                label(std::floor((localBb.min()[a] - grid_.origin[a])/h)) - 1
            );
            const label hi = min
            (
                label(grid_.dims[a] - 1),
                label(std::ceil((localBb.max()[a] - grid_.origin[a])/h)) + 1
            );
            n[a] = max(label(0), hi - lo[a] + 1);
        }

        bitSet claimed(n[0]*n[1]*n[2]);

        const pointField& pts = mesh_.points();
        const labelListList& cellPts = mesh_.cellPoints();

        forAll(cellPts, celli)
        {
            boundBox cb(pts, cellPts[celli], false);
            cb.inflate(1e-6);

            label c0[3], c1[3];
            for (direction a = 0; a < 3; ++a)
            {
                c0[a] = max
                (
                    lo[a],
                    label(std::ceil((cb.min()[a] - grid_.origin[a])/h))
                );
                c1[a] = min
                (
                    lo[a] + n[a] - 1,
                    label(std::floor((cb.max()[a] - grid_.origin[a])/h))
                );
            }

            for (label k = c0[2]; k <= c1[2]; ++k)
            {
                for (label j = c0[1]; j <= c1[1]; ++j)
                {
                    for (label i = c0[0]; i <= c1[0]; ++i)
                    {
                        const label l =
                            (i - lo[0])
                          + n[0]*((j - lo[1]) + n[1]*(k - lo[2]));
                        if (claimed.test(l))
                        {
                            continue;
                        }
                        const point p
                        (
                            grid_.origin[0] + i*h,
                            grid_.origin[1] + j*h,
                            grid_.origin[2] + k*h
                        );
                        if (!insideCell(p, celli))
                        {
                            continue;
                        }
                        claimed.set(l);
                        cells.push_back(celli);
                        sampleIndex_.push_back
                        (
                            uint32_t
                            (
                                i + grid_.dims[0]*(j + label(grid_.dims[1])*k)
                            )
                        );
                        if (scheme_ == "cellPoint")
                        {
                            weights_.emplace_back(mesh_, p, celli);
                        }
                    }
                }
            }
        }
    }
    sampleCells_.transfer(cells);

    // Brick routing
    const ::anariVis::BrickDecomposition decomp(grid_.dims, int(nConsumers));
    destSamples_.assign(nConsumers, std::vector<label>());
    std::vector<int> bricks;
    const int64_t nx = grid_.dims[0];
    const int64_t ny = grid_.dims[1];

    forAll(sampleCells_, s)
    {
        const int64_t g = sampleIndex_[s];
        decomp.bricksOf
        (
            int32_t(g % nx), int32_t((g/nx) % ny), int32_t(g/(nx*ny)), bricks
        );
        for (const int b : bricks)
        {
            destSamples_[b].push_back(s);
        }
    }

    destinations_.clear();
    for (label b = 0; b < nConsumers; ++b)
    {
        if (!destSamples_[b].empty())
        {
            destinations_.push_back(int(b));
        }
    }

    Info<< "anariVis: volume " << name_ << " field " << fieldName_
        << " grid " << grid_.dims[0] << 'x' << grid_.dims[1] << 'x'
        << grid_.dims[2] << " spacing " << h
        << " bricks " << decomp.parts()[0] << 'x' << decomp.parts()[1]
        << 'x' << decomp.parts()[2]
        << " samples " << returnReduce(sampleCells_.size(), sumOp<label>())
        << endl;
}


bool Foam::anariVisVolumeSampler::sample(scalarField& values) const
{
    tmp<volScalarField> tfld = lookupScalar(mesh_, fieldName_);
    if (!tfld.valid())
    {
        return false;
    }
    const volScalarField& fld = tfld();

    values.resize(sampleCells_.size());

    if (scheme_ == "cellPoint")
    {
        const interpolationCellPoint<scalar> interp(fld);
        forAll(values, s)
        {
            values[s] = interp.interpolate(weights_[s]);
        }
    }
    else
    {
        forAll(values, s)
        {
            values[s] = fld[sampleCells_[s]];
        }
    }
    return true;
}


void Foam::anariVisVolumeSampler::fillPiece
(
    const int b,
    const scalarField& values,
    ::anariVis::VolumePiece& piece
) const
{
    piece.grid = grid_;
    piece.index.clear();
    piece.value.clear();
    if (b < 0 || size_t(b) >= destSamples_.size())
    {
        return;
    }
    const std::vector<label>& ids = destSamples_[b];
    piece.index.reserve(ids.size());
    piece.value.reserve(ids.size());
    for (const label s : ids)
    {
        piece.index.push_back(sampleIndex_[s]);
        piece.value.push_back(float(values[s]));
    }
}


// ************************************************************************* //
