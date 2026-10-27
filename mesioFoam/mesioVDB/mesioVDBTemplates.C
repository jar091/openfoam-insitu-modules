/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenFOAM-Mesio integration project
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "volFields.H"
#include "mesio_insitu.h"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
const Foam::Field<Type>* Foam::functionObjects::mesioVDB::findCellValues
(
    const word& fieldName,
    label& eventNo
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef DimensionedField<Type, volMesh> IntFieldType;

    if (const auto* fldPtr = obr().cfindObject<VolFieldType>(fieldName))
    {
        eventNo = fldPtr->eventNo();
        return &(fldPtr->primitiveField());
    }

    if (const auto* fldPtr = obr().cfindObject<IntFieldType>(fieldName))
    {
        eventNo = fldPtr->eventNo();
        return &(fldPtr->field());
    }

    return nullptr;
}


template<class Type>
const Foam::dimensionSet* Foam::functionObjects::mesioVDB::findDimensions
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef DimensionedField<Type, volMesh> IntFieldType;

    if (const auto* fldPtr = obr().cfindObject<VolFieldType>(fieldName))
    {
        return &(fldPtr->dimensions());
    }

    if (const auto* fldPtr = obr().cfindObject<IntFieldType>(fieldName))
    {
        return &(fldPtr->dimensions());
    }

    return nullptr;
}


template<class Type>
bool Foam::functionObjects::mesioVDB::addFieldGrids(const word& fieldName)
{
    label eventNo = 0;
    const dimensionSet* dims = findDimensions<Type>(fieldName);
    if (!findCellValues<Type>(fieldName, eventNo) || !dims)
    {
        return false;
    }

    // Grid metadata: the units of all grids of a field are the field units
    const string units(unitsName(*dims));
    const string dimensions(dimensionsName(*dims));

    constexpr direction nCmpt = pTraits<Type>::nComponents;

    if constexpr (nCmpt == 1)
    {
        // Scalar or spherical tensor: one grid, the value itself
        grids_.append
        (
            gridEntry{fieldName, fieldName, FULL, 1, word::null, units, dimensions}
        );
        return true;
    }
    else
    {
        const wordList& selection =
        (
            nCmpt == 3 ? vectorComponents_ : tensorComponents_
        );

        for (const word& cmptName : selection)
        {
            if (cmptName == "mag")
            {
                if constexpr (nCmpt == 3)
                {
                    // Interleaved xyz, magnitude formed by Mesio
                    grids_.append
                    (
                        gridEntry
                        {
                            fieldName + "_mag", fieldName, FULL, 3,
                            "mag", units, dimensions
                        }
                    );
                }
                else
                {
                    grids_.append
                    (
                        gridEntry
                        {
                            fieldName + "_mag", fieldName, MAG, 1,
                            "mag", units, dimensions
                        }
                    );
                }
                continue;
            }

            label cmpt = -1;
            for (direction d = 0; d < nCmpt; ++d)
            {
                if (cmptName == pTraits<Type>::componentNames[d])
                {
                    cmpt = d;
                    break;
                }
            }

            if (cmpt < 0)
            {
                wordList valid(nCmpt + 1);
                valid[0] = "mag";
                for (direction d = 0; d < nCmpt; ++d)
                {
                    valid[d + 1] = pTraits<Type>::componentNames[d];
                }

                FatalErrorInFunction
                    << "Function object " << name()
                    << ": unknown component '" << cmptName
                    << "' for field " << fieldName
                    << " of type " << pTraits<Type>::typeName << nl
                    << "    Valid: " << flatOutput(valid)
                    << exit(FatalError);
            }

            grids_.append
            (
                gridEntry
                {
                    fieldName + '_' + cmptName, fieldName, cmpt, 1,
                    cmptName, units, dimensions
                }
            );
        }
        return true;
    }
}


template<class Type>
bool Foam::functionObjects::mesioVDB::setGridValues
(
    const gridEntry& g
)
{
    label eventNo = 0;
    const Field<Type>* fldPtr = findCellValues<Type>(g.field, eventNo);

    if (!fldPtr)
    {
        return false;
    }

    const Field<Type>& fld = *fldPtr;

    const MESIOReal* values = nullptr;

    if (g.component == FULL)
    {
        // Contiguous scalar (dimension 1) or vector (dimension 3) storage
        values = reinterpret_cast<const MESIOReal*>(fld.cdata());
    }
    else if (g.component == MAG)
    {
        buffer_.resize_nocopy(fld.size());
        forAll(fld, i)
        {
            buffer_[i] = Foam::mag(fld[i]);
        }
        values = buffer_.cdata();
    }
    else
    {
        buffer_.resize_nocopy(fld.size());
        const direction d = direction(g.component);
        forAll(fld, i)
        {
            buffer_[i] = component(fld[i], d);
        }
        values = buffer_.cdata();
    }

    // A failure is reported by the following MESIOInSituWriteVDB on all ranks
    mesioLocal
    (
        MESIOInSituSetCellField(ctx_, g.grid.c_str(), values),
        "MESIOInSituSetCellField"
    );

    return true;
}


// ************************************************************************* //
