#pragma once

// #include "volePSI/Defines.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Common/Timer.h" 
#include "psi/Defines.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h"
#include "andpair/triple.h"

using namespace oc;

namespace CmpFuzzyPSI{

    class mIMTSender : public oc::TimerAdapter
    {
    public:
        u64 mTableSize = 0;
        u64 mCmpsize = 0;
        u64 mDim = 0;
        u64 mDelta = 0;
        u64 mMetric = 0;
        u64 mCmp_len = 0;
        u64 mShareSize = 0;
        u64 mOutputSize = 0;

        BitVector e;
        BitVector d;
        BitVector e_and;
        BitVector d_and;
        BitVector mImt_e_Share;
        BitVector mImt_d_Share;

        Proto setUp(u64 tablesize, u64 dim, u64 delta, u64 metric, u64 Cmp_len, PRNG& prng, Socket& chl, u64 mNumThreads = 0);
        // For L_1, the first output is x \in [y_0, y_2], second output is x \in [y_1, y_2]
        Proto run(span<u64> Val, BitVector& v_s, BitVector& output, Socket& chl, u64 mNumThreads = 1); // Receiver sends the first message using OPPRF
    };

    class mIMTReceiver : public oc::TimerAdapter
    {
    public:
        u64 mTableSize = 0;
        u64 mCmpsize = 0;
        u64 mDim = 0;
        u64 mDelta = 0;
        u64 mMetric = 0;
        u64 mCmp_len = 0;
        u64 mShareSize = 0;
        u64 mOutputSize = 0;

        BitVector e;
        BitVector d;
        BitVector e_and;
        BitVector d_and;
        BitVector mImt_e_Share;
        BitVector mImt_d_Share;

        Proto setUp(u64 tablesize, u64 dim, u64 delta, u64 metric, u64 Cmp_len, PRNG& prng, Socket& chl, u64 mNumThreads = 0);
        Proto run(BitVector& output, Socket& chl, u64 mNumThreads = 1);
    };
}