#pragma once

#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Common/Timer.h" 
#include "psi/Defines.h"

namespace CmpFuzzyPSI{

    class PeqtSender : public oc::TimerAdapter
    {
    public:
        u64 mDataSize = 0;
        u64 mEqLength = 0;

        BitVector rShare;
        std::vector<u64> tShare;
        std::vector<std::vector<u8>> U, V;
        std::vector<u32> vose_val;
        BitVector vose_table;

        Proto setUp(u64 datasize, u64 eqlength, PRNG& prng, Socket& chl);
        Proto run(BitVector& data, BitVector& output, Socket& chl, u64 mNumThreads = 1);

    };

    class PeqtReceiver : public oc::TimerAdapter
    {
    public:
        u64 mDataSize = 0;
        u64 mEqLength = 0;

        BitVector rShare;
        std::vector<u64> tShare;
        std::vector<u32> eps;
        std::vector<std::vector<u8>> W;
        BitVector vose_table;

        Proto setUp(u64 datasize, u64 eqlength, PRNG& prng, Socket& chl);
        Proto run(BitVector& data, BitVector& output, Socket& chl, u64 mNumThreads = 1);

    };
}
