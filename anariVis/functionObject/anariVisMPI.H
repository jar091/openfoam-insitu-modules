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

Description
    Helpers between OpenFOAM's Pstream layer and raw MPI:
    the MPI_Comm behind an OpenFOAM communicator and the number of cores
    available to this process.

\*---------------------------------------------------------------------------*/

#ifndef anariVis_MPI_H
#define anariVis_MPI_H

#include "UPstream.H"

#include <mpi.h>
#include <sched.h>
#include <unistd.h>

#include <cstdint>
#include <thread>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace anariVis
{

//- Convert the wrapped value of a UPstream::Communicator into MPI_Comm
//  (MPI_Comm is a pointer in Open MPI and an integer in MPICH)
template<class CommType = MPI_Comm>
inline CommType toMPIComm(const Foam::UPstream::Communicator& c)
{
    if constexpr (std::is_pointer_v<CommType>)
    {
        return reinterpret_cast<CommType>(c.value());
    }
    else
    {
        return static_cast<CommType>(c.value());
    }
}


//- Number of cores in the affinity mask of this process
inline int availableCores()
{
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0)
    {
        const int n = CPU_COUNT(&set);
        if (n > 0)
        {
            return n;
        }
    }
    const unsigned n = std::thread::hardware_concurrency();
    return n ? int(n) : 1;
}


//- Allow the calling thread (and the threads it creates later) to run on
//  all cores of the node that the cgroup/cpuset permits. Render ranks call
//  this before creating the ANARI device, so that the render threads are not
//  confined to the core(s) the launcher bound the solver rank to.
inline bool widenAffinity()
{
    const long n = sysconf(_SC_NPROCESSORS_CONF);
    cpu_set_t set;
    CPU_ZERO(&set);
    for (long i = 0; i < n && i < CPU_SETSIZE; ++i)
    {
        CPU_SET(int(i), &set);
    }
    // The kernel intersects the mask with the permitted cpuset
    return sched_setaffinity(0, sizeof(set), &set) == 0;
}

} // End namespace anariVis

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
