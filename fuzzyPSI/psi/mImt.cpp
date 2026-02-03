#include "mImt.h"
#include <array>
#include <future>

namespace CmpFuzzyPSI
{
    Proto mIMTSender::setUp(u64 tablesize, u64 dim, u64 delta, u64 metric, u64 Cmp_len, PRNG& prng, Socket& chl, u64 mNumThreads)
    {
        // DEBUG_LOG("mIMTSender::setUp called with tablesize=" << tablesize << ", dim=" << dim << ", metric=" << metric << ", Cmp_len=" << Cmp_len);
        mTableSize = tablesize;
        mDim = dim;
        mDelta = delta;
        mMetric = metric;
        mCmp_len = Cmp_len;

        if (mMetric==0 || mMetric==2){
            mCmpsize = mTableSize*mDim*2;
            mShareSize = mTableSize*mCmp_len*mDim*2;
            mOutputSize = mTableSize*mDim;
        } else if (mMetric==1){
            mCmpsize = mTableSize*mDim*3;
            mShareSize = mTableSize*mCmp_len*mDim*3;
            mOutputSize = mTableSize*mDim*2;
        } else{
            throw std::runtime_error("mIMTSender::setUp unsupported metric");
        }
        // e_r * e_s = d_r + d_s
        e.resize(mShareSize); 
        d.resize(mShareSize);

        // AND pairs, sender order is eq_2k, cmp_2k+1, eq_2k+1
        e_and.resize(3*mShareSize);
        d_and.resize(3*mShareSize);

        mImt_e_Share.resize(mOutputSize*2);
        mImt_d_Share.resize(mOutputSize*2);

        u64 num_triples = mShareSize + 3*mShareSize + mOutputSize*2;
        bool silent = true;
        Triples triples_0(num_triples, silent);

        macoro::sync_wait(triples_0.gen0(chl));
        macoro::sync_wait(trans_andpair0(chl, triples_0));

        e.copy(triples_0.b, 0, mShareSize);
        d.copy(triples_0.c, 0, mShareSize);
        e_and.copy(triples_0.b, mShareSize, 3*mShareSize);
        d_and.copy(triples_0.c, mShareSize, 3*mShareSize);
        mImt_e_Share.copy(triples_0.b, 4*mShareSize, mOutputSize*2);
        mImt_d_Share.copy(triples_0.c, 4*mShareSize, mOutputSize*2);

        if (isImplement){
            // e,d e_and,d_and mImt_e_Share, mImt_d_Share are random AND shares

            // generate correlated e, d
            u64 masksize = mTableSize*mDim*mCmp_len;
            if (mMetric==0 || mMetric==2){
                auto recvmask0 = BitVector(masksize);
                // DEBUG_LOG("IMT [Sender] setUp step 1");
                sync_wait(chl.recv(recvmask0));
                for (u64 i=0; i < masksize; i++){
                    d[2*i+1] = (recvmask0[i] & e[2*i+1]) ^ d[2*i+1];
                }
            } else if (mMetric==1){
                auto recvmask0 = BitVector(2*masksize);
                sync_wait(chl.recv(recvmask0));
                for (u64 i=0; i < masksize; i++){
                    d[3*i+1] = (recvmask0[2*i] & e[3*i+1]) ^ d[3*i+1];
                    d[3*i+2] = (recvmask0[2*i+1] & e[3*i+2]) ^ d[3*i+2];
                }
            }


            // generate correlated e_and, d_and
            auto mask = BitVector(mShareSize);
            auto recvmask = BitVector(mShareSize);
            for (u64 i=0; i<mShareSize; i++){
                mask[i] = e_and[3*i] ^ e_and[3*i+2];
            }
            sync_wait(chl.send(std::move(mask)));
            sync_wait(chl.recv(recvmask));
            for (u64 i=0; i<mShareSize; i++){
                d_and[3*i+2] = (recvmask[i] & e_and[3*i+2]);
            }

        }

        co_return;
    }

    Proto mIMTReceiver::setUp(u64 tablesize, u64 dim, u64 delta, u64 metric, u64 Cmp_len, PRNG& prng, Socket& chl, u64 mNumThreads)
    {
        mTableSize = tablesize;
        mDim = dim;
        mDelta = delta;
        mMetric = metric;
        mCmp_len = Cmp_len;

        if (mMetric==0 || mMetric==2){
            mCmpsize = mTableSize*mDim*2;
            mShareSize = mTableSize*mCmp_len*mDim*2;
            mOutputSize = mTableSize*mDim;
        } else if (mMetric==1){
            mCmpsize = mTableSize*mDim*3;
            mShareSize = mTableSize*mCmp_len*mDim*3;
            mOutputSize = mTableSize*mDim*2;
        } else{
            throw std::runtime_error("mIMTSender::setUp unsupported metric");
        }
        // e_r * e_s = d_r + d_s, e[2*mTableSize*mDim]=e[2*mTableSize*mDim+1]
        e.resize(mShareSize);
        d.resize(mShareSize);

        // AND pairs, receiver order is cmp_2k+1, eq_2k, eq_2k+1, d_and_3k+2 = {e_and_3k+1}_r * {e_and_3k}_s + {e_and_3k+1}_r * {e_and_3k+2}_s
        e_and.resize(3*mShareSize);
        d_and.resize(3*mShareSize);

        mImt_e_Share.resize(mOutputSize*2);
        mImt_d_Share.resize(mOutputSize*2);

        u64 num_triples = mShareSize + 3*mShareSize + mOutputSize*2;
        bool silent = true;
        Triples triples_1(num_triples, silent);

        macoro::sync_wait(triples_1.gen1(chl));
        macoro::sync_wait(trans_andpair1(chl, triples_1));

        e.copy(triples_1.a, 0, mShareSize);
        d.copy(triples_1.c, 0, mShareSize);
        e_and.copy(triples_1.a, mShareSize, 3*mShareSize);
        d_and.copy(triples_1.c, mShareSize, 3*mShareSize);
        mImt_e_Share.copy(triples_1.a, 4*mShareSize, mOutputSize*2);
        mImt_d_Share.copy(triples_1.c, 4*mShareSize, mOutputSize*2);

        if (isImplement){
            // e,d e_and,d_and mImt_e_Share, mImt_d_Share are AND shares, (num: mShareSize+ 3*mShareSize+ mOutputSize*2)

            // generate correlated e, d
            u64 masksize = mTableSize*mDim*mCmp_len;
            if (mMetric==0 || mMetric==2){
                auto mask0 = BitVector(masksize);
                for (u64 i=0; i<masksize; i++){
                    mask0[i] = e[2*i] ^ e[2*i+1];
                    e[2*i] = e[2*i+1];
                }
                sync_wait(chl.send(std::move(mask0)));
            } else if (mMetric==1){
                auto mask0 = BitVector(2*masksize);
                for (u64 i=0; i < masksize; i++){
                    mask0[2*i] = e[3*i] ^ e[3*i+1];
                    mask0[2*i+1] = e[3*i] ^ e[3*i+2];
                    e[3*i+1] = e[3*i];
                    e[3*i+2] = e[3*i];
                }
                sync_wait(chl.send(std::move(mask0)));
            }

            // generate correlated e_and, d_and
            auto mask = BitVector(mShareSize);
            auto recvmask = BitVector(mShareSize);
            for (u64 i=0; i<mShareSize; i++){
                mask[i] = e_and[3*i+1] ^ e_and[3*i+2];
            }
            sync_wait(chl.recv(recvmask));
            sync_wait(chl.send(std::move(mask)));
            for (u64 i=0; i<mShareSize; i++){
                d_and[3*i+2] = (recvmask[i] & e_and[3*i+2]);
            }
        }

        co_return;
    }

    Proto mIMTSender::run(span<u64> Val, BitVector& v_s, BitVector& output, Socket& chl, u64 mNumThreads)
    {
        // run bCMP
        auto v = BitVector(mShareSize);
        auto cmp = BitVector(mShareSize);
        auto eq = BitVector(mShareSize);
        auto v_tmp = BitVector(mShareSize);
        auto cmp_output = BitVector(mOutputSize*2);
        output.resize(mOutputSize);

        // Generating phase
        for (u64 i=0; i < mCmpsize; i++){
            for (u64 j = 0; j < mCmp_len; j++){
                bool tmp = (Val[i] >> j) & 1;
                if ((mMetric==0 || mMetric==2) && i %2 ==1){ // > using !b_i
                    tmp = !tmp;
                } else if (mMetric==1 && (i % 3 == 2)){ // > using !b_i
                    tmp = !tmp;
                } else {
                    v_s[i*mCmp_len + j] = !v_s[i*mCmp_len + j];
                }
                v_tmp[i*mCmp_len + j] = tmp;
            }
        }
        v = v_tmp ^ e;
        cmp = (v_tmp & v_s) ^ d;
        eq = v_tmp ^ v_s;

        sync_wait(chl.send(std::move(v)));
        // Iteration phase

        for (u64 i=0; i < mShareSize/2; i++){
            bool tmp = cmp[i];
            cmp[i] = cmp[mShareSize-i-1];
            cmp[mShareSize-i-1] = tmp;

            tmp = eq[i];
            eq[i] = eq[mShareSize-i-1];
            eq[mShareSize-i-1] = tmp;
        }

        u64 IterSize = osuCrypto::log2ceil(mCmp_len);
        u64 eAndIndex = 0;
        u64 dAndIndex = 0;
        u64 PreStep = mCmp_len;
        for (u64 j=0; j < IterSize; j++){
            u64 step = PreStep/2;
            auto tem_cmp = BitVector(3 * mCmpsize * step);
            auto tem_cmpshare = BitVector(3 * mCmpsize * step);
            for (u64 i=0; i < mCmpsize; i++){
                for (u64 k=0; k < step; k++){
                    tem_cmp[3*(i*step + k)] = eq[(i*mCmp_len) + 2*k] ^ e_and[eAndIndex];
                    tem_cmp[3*(i*step + k)+1] = cmp[(i*mCmp_len) + 2*k+1] ^ e_and[eAndIndex + 1];
                    tem_cmp[3*(i*step + k)+2] = eq[(i*mCmp_len) + 2*k+1] ^ e_and[eAndIndex + 2];
                    eAndIndex+=3;
                }
            }
            sync_wait(chl.recv(tem_cmpshare));
            sync_wait(chl.send(std::move(tem_cmp)));
            for (u64 i=0; i < mCmpsize; i++){
                for (u64 k=0; k < step; k++){
                    cmp[(i*mCmp_len) + k] = 
                        (tem_cmpshare[3*(i*step + k)] & e_and[dAndIndex]) ^ d_and[dAndIndex] ^
                        (tem_cmpshare[3*(i*step + k)+1] & e_and[dAndIndex + 1]) ^ d_and[dAndIndex + 1] ^
                        (cmp[(i*mCmp_len) + (2*k)+1] & eq[(i*mCmp_len) + (2*k)]) ^ cmp[(i*mCmp_len) + (2*k)];
                    eq[(i*mCmp_len) + k] = 
                        (tem_cmpshare[3*(i*step + k)+2] & e_and[dAndIndex]) ^
                        (tem_cmpshare[3*(i*step + k)+1] & e_and[dAndIndex + 2]) ^ d_and[dAndIndex+2] ^
                        (eq[(i*mCmp_len) + 2*k+1] & eq[(i*mCmp_len) + 2*k]);
                    dAndIndex+=3;
                }
            }
        
            if (PreStep %2 == 1) {
                for (u64 i=0; i < mCmpsize; i++){
                    cmp[(i*mCmp_len) + step] = cmp[(i*mCmp_len) + 2 *step]; ;
                    eq[(i*mCmp_len) + step] = eq[(i*mCmp_len) + 2 *step];
                }
            }
            PreStep = (PreStep +1)/2;
        }

        if (mMetric==0 || mMetric==2){
            for (u64 i=0; i < mOutputSize; i++){
                if (Val[2*(mOutputSize - i -1)] > Val[2*(mOutputSize - i -1)+1]){ // reverse the output
                    // exchange the position of 2 CMP results in an IMT
                    cmp_output[2*i+1] = cmp[2*i*mCmp_len];
                    cmp_output[2*i] = cmp[(2*i+1)*mCmp_len];
                } else{ // reverse the input
                    cmp_output[2*i+1] = !cmp[2*i*mCmp_len];
                    cmp_output[2*i] = !cmp[(2*i+1)*mCmp_len];
                }
            }
        } else {
            for (u64 i=0; i < mOutputSize/2; i++){

                if (Val[3*(mOutputSize/2 - i -1)+1] > Val[3*(mOutputSize/2 - i -1)+2]){ // reverse the output
                    cmp_output[4*i+1] = cmp[3*i*mCmp_len];
                    cmp_output[4*i] = cmp[(3*i+1)*mCmp_len];

                } else{ // reverse the input
                    cmp_output[4*i+1] = !cmp[3*i*mCmp_len];
                    cmp_output[4*i] = !cmp[(3*i+1)*mCmp_len];
                }

                if (Val[3*(mOutputSize/2 - i -1)] > Val[3*(mOutputSize/2 - i -1)+2]){ // reverse the output
                    cmp_output[4*i+3] = cmp[(3*i)*mCmp_len];
                    cmp_output[4*i+2] = cmp[(3*i+2)*mCmp_len];
                } else{ // reverse the input
                    cmp_output[4*i+3] = !cmp[(3*i)*mCmp_len];
                    cmp_output[4*i+2] = !cmp[(3*i+2)*mCmp_len];
                }
            }
        }

        for (u64 i=0; i < mOutputSize; i++){
            bool tmp = cmp_output[i];
            cmp_output[i] = cmp_output[mOutputSize*2-i-1];
            cmp_output[mOutputSize*2-i-1] = tmp;
        }

        // output phase
        auto tem_Imtshare = BitVector(mOutputSize*2);
        sync_wait(chl.recv(tem_Imtshare));
        sync_wait(chl.send(std::move(cmp_output ^ mImt_e_Share)));

        if (mMetric==0 || mMetric==2){
            for (u64 i=0; i < mOutputSize; i++){
                output[i] = (tem_Imtshare[2*i] & cmp_output[2*i]) ^ mImt_d_Share[2*i] ^
                        (tem_Imtshare[2*i+1] & cmp_output[2*i+1]) ^ mImt_d_Share[2*i+1]^
                        (cmp_output[2*i] & cmp_output[2*i+1]);
                if (Val[2*i] > Val[2*i+1]){ // reverse the output
                    output[i] = !output[i];
                }
            }
        } else {
            for (u64 i=0; i < mOutputSize/2; i++){
                output[2*i] = (tem_Imtshare[4*i] & cmp_output[4*i]) ^ mImt_d_Share[4*i] ^
                        (tem_Imtshare[4*i+1] & cmp_output[4*i+1]) ^ mImt_d_Share[4*i+1] ^
                        (cmp_output[4*i] & cmp_output[4*i+1]);
                output[2*i+1] = (tem_Imtshare[4*i+2] & cmp_output[4*i+2]) ^ mImt_d_Share[4*i+2] ^
                        (tem_Imtshare[4*i+3] & cmp_output[4*i+3]) ^ mImt_d_Share[4*i+3] ^
                        (cmp_output[4*i+2] & cmp_output[4*i+3]);
                if (Val[3*i] > Val[3*i+2]){ // reverse the output
                    output[2*i] = !output[2*i];
                }
                if (Val[3*i+1] > Val[3*i+2]){ // reverse the output
                    output[2*i+1] = !output[2*i+1];
                }
            }
        }
        
        co_return;
    }

    Proto mIMTReceiver::run(BitVector& output, Socket& chl, u64 mNumThreads)
    {
        // run bCMP
        auto v = BitVector(mShareSize);
        auto cmp = BitVector(mShareSize);
        auto eq = BitVector(mShareSize);
        auto cmp_output = BitVector(mOutputSize*2);
        output.resize(mOutputSize);

        // generating phase
        sync_wait(chl.recv(v));
        cmp = (e & v) ^ d;
        eq = e;
        // Iteration phase

        for (u64 i=0; i < mShareSize/2; i++){
            bool tmp = cmp[i];
            cmp[i] = cmp[mShareSize-i-1];
            cmp[mShareSize-i-1] = tmp;

            tmp = eq[i];
            eq[i] = eq[mShareSize-i-1];
            eq[mShareSize-i-1] = tmp;
        }

        u64 IterSize = osuCrypto::log2ceil(mCmp_len);
        u64 eAndIndex = 0;
        u64 dAndIndex = 0;
        u64 PreStep = mCmp_len;
        for (u64 j=0; j < IterSize; j++){
            u64 step = PreStep/2;
            auto tem_cmp = BitVector(3 * mCmpsize * step);
            auto tem_cmpshare = BitVector(3 * mCmpsize * step);
            for (u64 i=0; i < mCmpsize; i++){
                for (u64 k=0; k < step; k++){
                    tem_cmp[3*(i*step + k)] = cmp[(i*mCmp_len) + 2*k+1] ^ e_and[eAndIndex];
                    tem_cmp[3*(i*step + k)+1] = eq[(i*mCmp_len) + 2*k] ^ e_and[eAndIndex + 1];
                    tem_cmp[3*(i*step + k)+2] = eq[(i*mCmp_len) + 2*k+1] ^ e_and[eAndIndex + 2];
                    eAndIndex+=3;
                }
            }
            sync_wait(chl.send(std::move(tem_cmp)));
            sync_wait(chl.recv(tem_cmpshare));
            for (u64 i=0; i < mCmpsize; i++){
                for (u64 k=0; k < step; k++){
                    cmp[(i*mCmp_len) + k] = 
                        (tem_cmpshare[3*(i*step + k)] & cmp[(i*mCmp_len) + 2*k+1]) ^ d_and[dAndIndex] ^
                        (tem_cmpshare[3*(i*step + k)+1] & eq[(i*mCmp_len) + 2*k]) ^ d_and[dAndIndex + 1] ^
                        (cmp[(i*mCmp_len) + (2*k)+1] & eq[(i*mCmp_len) + (2*k)]) ^ cmp[(i*mCmp_len) + (2*k)];

                    eq[(i*mCmp_len) + k] = 
                        (tem_cmpshare[3*(i*step + k)+2] & eq[(i*mCmp_len) + 2*k+1]) ^
                        (tem_cmpshare[3*(i*step + k)] & eq[(i*mCmp_len) + 2*k+1]) ^ d_and[dAndIndex+2] ^
                        (eq[(i*mCmp_len) + 2*k+1] & eq[(i*mCmp_len) + 2*k]);
                    dAndIndex+=3;
                }
            }
            if (PreStep %2 == 1) {
                for (u64 i=0; i < mCmpsize; i++){
                    cmp[(i*mCmp_len) + step] = cmp[(i*mCmp_len) + 2 *step]; ;
                    eq[(i*mCmp_len) + step] = eq[(i*mCmp_len) + 2 *step];
                }
            }
            PreStep = (PreStep +1)/2;
        }

        if (mMetric==0 || mMetric==2){
            for (u64 i=0; i < mCmpsize; i++){
                cmp_output[i] = cmp[i*mCmp_len];
            }
        } else {
            for (u64 i=0; i < mOutputSize/2; i++){
                cmp_output[4*i] = cmp[3*i*mCmp_len];
                cmp_output[4*i+1] = cmp[(3*i+1)*mCmp_len];
                cmp_output[4*i+2] = cmp[(3*i)*mCmp_len];
                cmp_output[4*i+3] = cmp[(3*i+2)*mCmp_len];
            }
        }

        for (u64 i=0; i < mOutputSize; i++){
            bool tmp = cmp_output[i];
            cmp_output[i] = cmp_output[mOutputSize*2-i-1];
            cmp_output[mOutputSize*2-i-1] = tmp;
        }

        // output phase
        auto tem_Imtshare = BitVector(mOutputSize*2);
        sync_wait(chl.send(std::move(cmp_output ^ mImt_e_Share)));
        sync_wait(chl.recv(tem_Imtshare));

        for (u64 i=0; i < mOutputSize; i++){
            output[i] = (tem_Imtshare[2*i] & mImt_e_Share[2*i]) ^ mImt_d_Share[2*i] ^
                (tem_Imtshare[2*i+1] & mImt_e_Share[2*i+1]) ^ mImt_d_Share[2*i+1]^
                (cmp_output[2*i] & cmp_output[2*i+1]);
        }

        co_return;
    }
}