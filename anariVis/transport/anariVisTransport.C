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

#include "anariVisTransport.H"
#include "anariVisTypes.H"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <time.h>

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

int anariVis::Transport::nodeLocalRank(MPI_Comm comm)
{
    MPI_Comm node;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node);
    int r = 0;
    MPI_Comm_rank(node, &r);
    MPI_Comm_free(&node);
    return r;
}


int anariVis::Transport::nodeLocalSize(MPI_Comm comm)
{
    MPI_Comm node;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node);
    int n = 0;
    MPI_Comm_size(node, &n);
    MPI_Comm_free(&node);
    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

anariVis::Transport::Transport
(
    MPI_Comm parentComm,
    bool producer,
    bool consumer
)
:
    visComm_(MPI_COMM_NULL),
    renderComm_(MPI_COMM_NULL),
    visRank_(0),
    visSize_(1),
    producer_(producer),
    consumer_(consumer),
    nExpected_(0),
    nStopped_(0)
{
    MPI_Comm_dup(parentComm, &visComm_);
    MPI_Comm_rank(visComm_, &visRank_);
    MPI_Comm_size(visComm_, &visSize_);

    // Node id = vis rank of the lowest rank on the same node
    int nodeId = visRank_;
    {
        MPI_Comm node;
        MPI_Comm_split_type
        (
            visComm_, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node
        );
        MPI_Allreduce(MPI_IN_PLACE, &nodeId, 1, MPI_INT, MPI_MIN, node);
        MPI_Comm_free(&node);
    }

    const int mine[2] = {(producer_ ? 1 : 0) | (consumer_ ? 2 : 0), nodeId};
    std::vector<int> all(2*visSize_);
    MPI_Allgather(mine, 2, MPI_INT, all.data(), 2, MPI_INT, visComm_);

    nodeOf_.resize(visSize_);
    for (int r = 0; r < visSize_; ++r)
    {
        if (all[2*r] & 1) producers_.push_back(r);
        if (all[2*r] & 2) consumers_.push_back(r);
        nodeOf_[r] = all[2*r + 1];
    }

    if (consumers_.empty() || producers_.empty())
    {
        throw std::runtime_error
        (
            "anariVis::Transport: need at least one producer and one consumer"
        );
    }

    MPI_Comm_split
    (
        visComm_, consumer_ ? 0 : MPI_UNDEFINED, visRank_, &renderComm_
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

anariVis::Transport::~Transport()
{
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized)
    {
        return;
    }

    if (!requests_.empty())
    {
        MPI_Waitall(int(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE);
    }
    if (renderComm_ != MPI_COMM_NULL)
    {
        MPI_Comm_free(&renderComm_);
    }
    if (visComm_ != MPI_COMM_NULL)
    {
        MPI_Comm_free(&visComm_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int anariVis::Transport::producerIndex() const
{
    auto it = std::lower_bound(producers_.begin(), producers_.end(), visRank_);
    return (it != producers_.end() && *it == visRank_)
        ? int(it - producers_.begin()) : -1;
}


int anariVis::Transport::consumerIndex() const
{
    auto it = std::lower_bound(consumers_.begin(), consumers_.end(), visRank_);
    return (it != consumers_.end() && *it == visRank_)
        ? int(it - consumers_.begin()) : -1;
}


bool anariVis::Transport::isLeader() const
{
    return producer_ && visRank_ == producers_.front();
}


int anariVis::Transport::nConsumersOnNode() const
{
    int n = 0;
    for (const int c : consumers_)
    {
        if (nodeOf_[c] == nodeOf_[visRank_]) ++n;
    }
    return n;
}


int anariVis::Transport::surfaceConsumer() const
{
    if (!producer_)
    {
        return -1;
    }
    const int node = nodeOf_[visRank_];

    std::vector<int> localConsumers;
    for (int c = 0; c < int(consumers_.size()); ++c)
    {
        if (nodeOf_[consumers_[c]] == node)
        {
            localConsumers.push_back(c);
        }
    }

    if (!localConsumers.empty())
    {
        // Index of this producer among the producers of this node
        int k = 0;
        for (const int p : producers_)
        {
            if (p == visRank_) break;
            if (nodeOf_[p] == node) ++k;
        }
        return localConsumers[k % localConsumers.size()];
    }

    // Block mapping of all producers onto all consumers
    const long long p = producerIndex();
    return int((p*consumers_.size())/producers_.size());
}


void anariVis::Transport::sendConfig(const std::string& text)
{
    if (!isLeader())
    {
        return;
    }
    for (const int c : consumers_)
    {
        if (c == visRank_)
        {
            continue;
        }
        MPI_Send
        (
            text.data(), int(text.size()), MPI_CHAR, c, tagConfig, visComm_
        );
    }
}


std::string anariVis::Transport::receiveConfig()
{
    MPI_Status status;
    MPI_Probe(producers_.front(), tagConfig, visComm_, &status);
    int n = 0;
    MPI_Get_count(&status, MPI_CHAR, &n);
    std::string text(size_t(n), '\0');
    MPI_Recv
    (
        &text[0], n, MPI_CHAR, producers_.front(), tagConfig, visComm_,
        MPI_STATUS_IGNORE
    );
    return text;
}


void anariVis::Transport::setDestinations(const std::vector<int>& dest)
{
    std::vector<int> counts(visSize_, 0);
    if (producer_)
    {
        for (const int c : dest)
        {
            counts[consumers_[c]] = 1;
        }
    }
    int n = 0;
    MPI_Reduce_scatter_block
    (
        counts.data(), &n, 1, MPI_INT, MPI_SUM, visComm_
    );
    nExpected_ = consumer_ ? n : 0;
}


bool anariVis::Transport::sendsComplete()
{
    if (requests_.empty())
    {
        return true;
    }
    int flag = 0;
    MPI_Testall
    (
        int(requests_.size()), requests_.data(), &flag, MPI_STATUSES_IGNORE
    );
    if (flag)
    {
        requests_.clear();
        sendBuffers_.clear();
    }
    return flag != 0;
}


void anariVis::Transport::waitSends()
{
    if (!requests_.empty())
    {
        MPI_Waitall
        (
            int(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE
        );
    }
    requests_.clear();
    sendBuffers_.clear();
}


void anariVis::Transport::post(int consumerIdx, std::vector<char>&& buffer)
{
    sendBuffers_.push_back(std::move(buffer));
    const std::vector<char>& buf = sendBuffers_.back();

    if (buf.size() > size_t(std::numeric_limits<int>::max()))
    {
        throw std::runtime_error("anariVis::Transport: packet exceeds 2 GB");
    }

    MPI_Request req;
    MPI_Issend
    (
        buf.data(), int(buf.size()), MPI_BYTE, consumers_[consumerIdx],
        tagFrame, visComm_, &req
    );
    requests_.push_back(req);
}


bool anariVis::Transport::receiveFrame
(
    std::vector<std::vector<char>>& packets,
    int& frameIndex
)
{
    packets.clear();

    while (true)
    {
        // A complete frame (the lowest index first)?
        if (!pending_.empty())
        {
            auto it = pending_.begin();
            if (int(it->second.size()) >= nExpected_)
            {
                frameIndex = it->first;
                packets = std::move(it->second);
                pending_.erase(it);
                return true;
            }
        }

        if (nStopped_ >= nExpected_)
        {
            return false;
        }

        MPI_Message msg;
        MPI_Status status;
        MPI_Mprobe(MPI_ANY_SOURCE, tagFrame, visComm_, &msg, &status);
        int n = 0;
        MPI_Get_count(&status, MPI_BYTE, &n);
        std::vector<char> buf(static_cast<size_t>(n));
        MPI_Mrecv(buf.data(), n, MPI_BYTE, &msg, MPI_STATUS_IGNORE);

        FrameHeader h;
        if (!FramePacket::peekHeader(buf.data(), buf.size(), h))
        {
            throw std::runtime_error("anariVis::Transport: corrupt packet");
        }
        if (h.flags & FrameHeader::stopFlag)
        {
            ++nStopped_;
        }
        else
        {
            pending_[h.frameIndex].push_back(std::move(buf));
        }
    }
}


bool anariVis::Transport::pollFrame(int& frameIndex)
{
    while (nStopped_ < nExpected_)
    {
        int flag = 0;
        MPI_Message msg;
        MPI_Status status;
        MPI_Improbe(MPI_ANY_SOURCE, tagFrame, visComm_, &flag, &msg, &status);
        if (!flag)
        {
            break;
        }
        int n = 0;
        MPI_Get_count(&status, MPI_BYTE, &n);
        std::vector<char> buf(static_cast<size_t>(n));
        MPI_Mrecv(buf.data(), n, MPI_BYTE, &msg, MPI_STATUS_IGNORE);

        FrameHeader h;
        if (!FramePacket::peekHeader(buf.data(), buf.size(), h))
        {
            throw std::runtime_error("anariVis::Transport: corrupt packet");
        }
        if (h.flags & FrameHeader::stopFlag)
        {
            ++nStopped_;
        }
        else
        {
            pending_[h.frameIndex].push_back(std::move(buf));
        }
    }

    if (!pending_.empty())
    {
        auto it = pending_.begin();
        if (int(it->second.size()) >= nExpected_)
        {
            frameIndex = it->first;
            return true;
        }
    }
    return false;
}


bool anariVis::Transport::takeFrame
(
    int frameIndex,
    std::vector<std::vector<char>>& packets
)
{
    auto it = pending_.find(frameIndex);
    if (it == pending_.end() || int(it->second.size()) < nExpected_)
    {
        return false;
    }
    packets = std::move(it->second);
    pending_.erase(it);
    return true;
}


void anariVis::Transport::sleepBarrier(int sleepMicroSeconds)
{
    MPI_Request req;
    MPI_Ibarrier(visComm_, &req);
    int flag = 0;
    struct timespec ts;
    ts.tv_sec = 0;
    ts.tv_nsec = long(sleepMicroSeconds)*1000L;
    while (true)
    {
        MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
        if (flag)
        {
            break;
        }
        nanosleep(&ts, nullptr);
    }
}


// ************************************************************************* //
