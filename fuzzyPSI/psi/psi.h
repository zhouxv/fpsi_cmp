#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <future>
#include <cstdint>
#include <stdexcept>
#include "volePSI/RsOpprf.h"
#include "sparsehash/dense_hash_map"
#include "cryptoTools/Common/Timer.h" 
#include "cryptoTools/Common/CuckooIndex.h"
#include "volePSI/SimpleIndex.h"
#include "fmap.h"
#include "mImt.h"
#include "Defines.h"
// #include "peqt/neweq.h"
#include "mPeqt/peqt.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h"

using namespace volePSI;
using namespace oc;

namespace CmpFuzzyPSI
{
    namespace details
    {
        struct FuzzyPsiBase
        {

            u64 mSenderSize = 0;
            u64 mRecverSize = 0;
            u64 mSsp = 0;
            u64 mDim = 0;
            u64 mMetric = 0;
            u64 mDelta = 0;
            u64 mLorH = 0;

            PRNG mPrng;
            bool mCompress = true;
            u64 mNumThreads = 0;
            u64 mMaskSize = 0;
            bool mUseReducedRounds = false;
            bool mDebug = false;

            void init(u64 senderSize, u64 recverSize, u64 statSecParam, u64 dim, u64 metric, u64 delta, u64 LorH, block seed, u64 numThreads, bool useReducedRounds = false);

        };
    }

    class FuzzyPsiSender : public details::FuzzyPsiBase, public oc::TimerAdapter
    {
    public:

        Proto run(span<block> inputs, Socket& chl);
        Proto runLinfty(span<block> inputs, Socket& chl);
        Proto runL1(span<block> inputs, Socket& chl);
        Proto runL2(span<block> inputs, Socket& chl);
    };


    class FuzzyPsiReceiver : public details::FuzzyPsiBase, public oc::TimerAdapter
    {
    public:

        std::vector<u64> mIntersection;

        Proto run(span<block> inputs, Socket& chl);
        Proto runLinfty(span<block> inputs, Socket& chl);
        Proto runL1(span<block> inputs, Socket& chl);
        Proto runL2(span<block> inputs, Socket& chl);
    };

    struct ModOps {
        u64 N;

        ModOps(u64 n) : N(n) {
            if (n == 0) throw std::invalid_argument("N must be > 0");
        }

        u64 add(u64 a, u64 b) const {
            a %= N;
            b %= N;
            u64 s = a + b;
            return (s >= N) ? s - N : s;
        }

        u64 sub(u64 a, u64 b) const {
            a %= N;
            b %= N;
            return (a >= b) ? a - b : N - (b - a);
        }

        u64 mul(u64 a, u64 b) const {
            a %= N;
            b %= N;
    #ifdef __SIZEOF_INT128__
            return (u64)(((__uint128_t)a * b) % N);
    #else
            // fallback to binary multiplication
            a %= N; b %= N;
            u64 res = 0;
            while (b) {
                if (b & 1) res = add(res, a);
                a = add(a, a);
                b >>= 1;
            }
            return res;
    #endif
        }

        u64 pow(u64 a, u64 e) const { // 快速幂
            u64 res = 1 % N;
            a %= N;
            while (e) {
                if (e & 1) res = mul(res, a);
                a = mul(a, a);
                e >>= 1;
            }
            return res;
        }

        u64 div2(u64 a) const {
            a %= N;
            if (N%2 == 1){
                if (a % 2 == 0) {
                    return a / 2;
                } else {
                    return (a / 2) + (N / 2) + 1;
                }
            }
            else
                throw std::invalid_argument("N must be odd");
        }
    };
}