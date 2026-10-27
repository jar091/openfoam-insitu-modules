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
    anariVis::Transport

Description
    MPI transport between OpenFOAM producer ranks and render (consumer)
    ranks of the anariVis layer. Renderer-agnostic, no OpenFOAM types.

    Communicators:
    - visComm    : duplicate of the parent communicator. In-situ the parent
                   is the OpenFOAM world, in-transit it is the global
                   MPI_COMM_WORLD that spans the simulation world and the
                   render world (OpenFOAM option -world).
    - renderComm : the consumer ranks only (Barney / deepCompositing).

    Protocol (per frame f, all messages with tag tagFrame on visComm):
    \verbatim
    setup (collective over visComm, once):
        MPI_Comm_dup, MPI_Allgather(role, node), MPI_Comm_split(consumer)
        in-transit: leader producer -> every consumer: settings (tagConfig)
        MPI_Reduce_scatter_block: every consumer learns how many producers
        send to it per frame (nExpected)
    producer, every frame f:
        if previous sends not complete (MPI_Testall, agreed over all
        producers by the caller): drop frame f (latest-frame policy)
        else: pack FramePacket per destination consumer, MPI_Issend
        (synchronous mode: completion implies the consumer matched the
        message, which bounds the in-flight data to one frame)
    consumer:
        MPI_Mprobe/MPI_Mrecv(ANY_SOURCE) until all nExpected packets of the
        lowest pending frame index have arrived, then assemble and render
    end:
        every producer sends a packet with FrameHeader::stopFlag to all of
        its destinations; a consumer stops after nExpected stop packets
    \endverbatim

    The producer leader always sends to every consumer, so that every
    consumer receives at least one packet (header, time, bounds) per frame.

SourceFiles
    anariVisTransport.C

\*---------------------------------------------------------------------------*/

#ifndef anariVis_Transport_H
#define anariVis_Transport_H

#include <mpi.h>

#include <map>
#include <string>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace anariVis
{

/*---------------------------------------------------------------------------*\
                          Class Transport Declaration
\*---------------------------------------------------------------------------*/

class Transport
{
public:

    static constexpr int tagConfig = 7301;
    static constexpr int tagFrame = 7302;


private:

    // Private Data

        MPI_Comm visComm_;
        MPI_Comm renderComm_;

        int visRank_;
        int visSize_;

        bool producer_;
        bool consumer_;

        //- Vis ranks of all producers / consumers (ascending)
        std::vector<int> producers_;
        std::vector<int> consumers_;

        //- Node id (vis rank of the node leader) of every vis rank
        std::vector<int> nodeOf_;

        //- Number of packets a consumer receives per frame
        int nExpected_;

        //- Outstanding sends (producer)
        std::vector<MPI_Request> requests_;
        std::vector<std::vector<char>> sendBuffers_;

        //- Received packets of not yet complete frames (consumer)
        std::map<int, std::vector<std::vector<char>>> pending_;

        int nStopped_;


public:

    // Constructors

        //- Construct, collective over parentComm
        Transport(MPI_Comm parentComm, bool producer, bool consumer);

        Transport(const Transport&) = delete;
        void operator=(const Transport&) = delete;


    //- Destructor
    ~Transport();


    // Static Functions

        //- Rank within the node (shared memory domain) of comm
        static int nodeLocalRank(MPI_Comm comm);

        //- Number of ranks of comm on this node
        static int nodeLocalSize(MPI_Comm comm);


    // Access

        MPI_Comm visComm() const { return visComm_; }

        //- The consumer communicator (MPI_COMM_NULL on pure producers)
        MPI_Comm renderComm() const { return renderComm_; }

        int visRank() const { return visRank_; }
        int visSize() const { return visSize_; }

        bool isProducer() const { return producer_; }
        bool isConsumer() const { return consumer_; }

        int nProducers() const { return int(producers_.size()); }
        int nConsumers() const { return int(consumers_.size()); }

        //- Index of this rank among the producers (-1 if none)
        int producerIndex() const;

        //- Index of this rank among the consumers (-1 if none)
        int consumerIndex() const;

        //- Vis rank of consumer idx
        int consumerRank(int idx) const { return consumers_[idx]; }

        //- True on the producer with the lowest vis rank
        bool isLeader() const;

        //- Consumer index that receives this producer's surfaces:
        //  a consumer on the same node if there is one, otherwise a block
        //  mapping of producers onto consumers
        int surfaceConsumer() const;

        //- Packets per frame expected by this consumer
        int nExpected() const { return nExpected_; }

        //- Number of consumers on the node of this rank
        int nConsumersOnNode() const;


    // Setup

        //- Leader: send the settings text to all non-producer consumers
        void sendConfig(const std::string& text);

        //- Consumer (non-producer): receive the settings text
        std::string receiveConfig();

        //- Collective over visComm: declare the consumers this producer
        //  sends a packet to every frame (consumers pass an empty list)
        void setDestinations(const std::vector<int>& consumerIdx);


    // Producer

        //- Non-blocking test whether all previous sends completed
        bool sendsComplete();

        //- Wait for all previous sends
        void waitSends();

        //- Post a packet to consumer idx (MPI_Issend, takes the buffer)
        void post(int consumerIdx, std::vector<char>&& buffer);


    // Consumer

        //- Receive all packets of the next frame. Returns false once all
        //  producers have sent their stop packet.
        bool receiveFrame
        (
            std::vector<std::vector<char>>& packets,
            int& frameIndex
        );


    // Synchronisation

        //- Barrier over visComm that sleeps while waiting, so that waiting
        //  producers leave their cores to the render threads
        void sleepBarrier(int sleepMicroSeconds = 200);
};


} // End namespace anariVis

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
