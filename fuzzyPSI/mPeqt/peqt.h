#pragma once

#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Common/Timer.h" 
#include "psi/Defines.h"
#include "n1not.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h"

namespace CmpFuzzyPSI{

    class PeqtSender : public oc::TimerAdapter
    {
    public:
        u64 mDataSize = 0;
        u64 mEqLength = 0;

        BitVector convert_bit; // convert_bit_r[i] ^ convert_bit_s[i] = convert_val_r[i] + convert_val_s[i]
        std::vector<u8> convert_val;
        std::vector<u8> vose_val; // vose_table[vose_val_r[i] + vose_val_s[i]]=1, else 0
        BitVector vose_table;

        Proto setUp(u64 datasize, u64 eqlength, PRNG& prng, Socket& chl);
        Proto vose(BitVector &output, int nums, int size, PRNG& prng, Socket& chl, u64 mNumThreads = 1);
        Proto convert(u64 p, int size, PRNG& prng, Socket& chl);
        Proto run(BitVector& data, BitVector& output, Socket& chl, u64 mNumThreads = 1);

    };

    class PeqtReceiver : public oc::TimerAdapter
    {
    public:
        u64 mDataSize = 0;
        u64 mEqLength = 0;

        BitVector convert_bit;
        std::vector<u8> convert_val;
        std::vector<u8> vose_val;
        BitVector vose_table;

        Proto setUp(u64 datasize, u64 eqlength, PRNG& prng, Socket& chl);
        Proto vose(BitVector &output, int nums, int size, PRNG& prng, Socket& chl, u64 mNumThreads = 1);
        Proto convert(u64 p, int size, PRNG& prng, Socket& chl);
        Proto run(BitVector& data, BitVector& output, Socket& chl, u64 mNumThreads = 1);

    };
}