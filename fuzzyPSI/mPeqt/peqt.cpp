#include "peqt.h"
#include <array>
#include <future>


namespace CmpFuzzyPSI
{
    Proto PeqtSender::setUp(u64 datasize, u64 eqlength, PRNG& prng, Socket& chl){
        mDataSize = datasize;
        mEqLength = eqlength;
        bool doFakeSetup = true;
        u64 modlength = oc::log2ceil(eqlength);
        u64 mod = 1 << modlength;

        convert_bit.resize(mDataSize*mEqLength);
        convert_val.resize(mDataSize*mEqLength);
        vose_val.resize(mDataSize);
        vose_table.resize(mDataSize * mod, 0);

        if (doFakeSetup){
            prng.get(convert_bit.data(), convert_bit.sizeBytes());
            prng.get(convert_val.data(), convert_val.size());
            BitVector convert_bit_recv(mDataSize*mEqLength);
            std::vector<u8> convert_val_send(mDataSize*mEqLength);
            macoro::sync_wait(chl.recv(convert_bit_recv));
            for (u64 i=0; i<mDataSize*mEqLength; i++){
                convert_val_send[i] = (convert_bit_recv[i]^convert_bit[i]) - convert_val[i];
            }
            macoro::sync_wait(chl.send(std::move(convert_val_send)));

            std::vector<u8> vose_val_recv(mDataSize);
            prng.get(vose_val.data(), vose_val.size());
            macoro::sync_wait(chl.recv(vose_val_recv));
            macoro::sync_wait(chl.send(vose_val));
            for (u64 i=0; i<mDataSize; i++){
                u8 tmp = (vose_val_recv[i] + vose_val[i]) % mod;
                vose_table[i*mod+tmp] = 1;
            }
            BitVector vose_table_recv(mDataSize*mod);
            macoro::sync_wait(chl.recv(vose_table_recv));
            vose_table = vose_table ^ vose_table_recv;
        } else {
            convert(mod, mDataSize*mEqLength, prng, chl);
            vose(vose_table, mDataSize, mEqLength, prng, chl, 1);
        }


        co_return;
    }

    Proto PeqtReceiver::setUp(u64 datasize, u64 eqlength, PRNG& prng, Socket& chl){
        mDataSize = datasize;
        mEqLength = eqlength;
        bool doFakeSetup = true;
        u64 modlength = oc::log2ceil(eqlength);
        u64 mod = 1 << modlength;

        convert_bit.resize(mDataSize*mEqLength);
        convert_val.resize(mDataSize*mEqLength);
        vose_val.resize(mDataSize);
        vose_table.resize(mDataSize * mod, 0);

        if (doFakeSetup){
            prng.get(convert_bit.data(), convert_bit.sizeBytes());
            macoro::sync_wait(chl.send(convert_bit));
            macoro::sync_wait(chl.recv(convert_val));

            std::vector<u8> vose_val_recv(mDataSize);
            prng.get(vose_val.data(), vose_val.size());
            macoro::sync_wait(chl.send(vose_val));
            macoro::sync_wait(chl.recv(vose_val_recv));
            prng.get(vose_table.data(), vose_table.sizeBytes());
            macoro::sync_wait(chl.send(vose_table));
        } else {
            convert(mod, mDataSize*mEqLength, prng, chl);
            vose(vose_table, mDataSize, mEqLength, prng, chl, 1);
        }

        // for (u64 i=0; i < 5; i++){
        //     std::cout << "data " << i << " : ";
        //     for (u64 j=0; j < 5; j++){
        //         std::cout << (u64)(convert_val[i*mEqLength+j]) << " " << convert_bit[i*mEqLength+j] << "; "; 
        //     }
        //     std::cout << std::endl;
        // }

        co_return;
    }

    Proto PeqtSender::convert(u64 p, int size, PRNG& prng, Socket& chl) {
        prng.get(convert_bit.data(), convert_bit.sizeBytes());
        IknpOtExtSender sender;
        IknpOtExtReceiver recver;

        u64 *pad = new u64[2*size];
        coproto::span<u64> pat(pad, 2*size);
        prng.get(osuCrypto::span<u8>(convert_val));

        std::vector<std::array<block, 2>> data(size);
        sync_wait(sender.send(data, prng, chl));

        //#pragma omp parallel for
        for (int64_t i = 0; i < size; ++i) {
            pat[2*i] = *(u64*)(data[i][0].data()) ^ ((convert_val[i] - convert_bit[i]) % p);
            pat[2*i+1] = *(u64*)(data[i][1].data()) ^ ((convert_val[i] - 1 + convert_bit[i]) % p);
        }

        sync_wait(chl.send(pat));
        delete pad;

        co_return;
    }

    Proto PeqtReceiver::convert(u64 p, int size, PRNG& prng, Socket& chl) {
        prng.get(convert_bit.data(), convert_bit.sizeBytes());
        IknpOtExtSender sender;
        IknpOtExtReceiver recver;

        vector<block> data(size);
        sync_wait(recver.receive(convert_bit, data, prng, chl));
    
        u64 *pad = new u64[2*size];
        coproto::span<u64> pat(pad, 2 * size);
        sync_wait(chl.recv(pat));

        //#pragma omp parallel for
        for (int i = 0; i < size; ++i) {
            convert_val[i] = *(u64*)data[i].data() ^ pat[2*i+convert_bit[i]];
        }

        co_return;
    }

    Proto PeqtSender::vose(BitVector &output, int nums, int size, PRNG& prng, Socket& chl, u64 mNumThreads) {
        int simd_size = 16;
        uint64_t** seed = new uint64_t*[nums];
        for(int i = 0; i < nums; ++i) {
            seed[i] = new uint64_t[size]();
        }
        int simd_num = (size + 15)/simd_size;

        int ell2 = 1 << size;
        int sumsize = nums * ell2;
        auto bits = BitVector(sumsize);
        prng.get(osuCrypto::span<u8>(vose_val));
        // //#pragma omp parallel for
        for (int i = 0; i < nums; ++i) {
            vose_val[i] %= size;
            bits[i*size + vose_val[i]] = 1;
        }
        
        // TODO: 
        IknpOtExtSender *sender;
        IknpOtExtReceiver *recver;
        NCN1OT<uint64_t> ot(Role::Sender, nums, size, sender, recver, chl);

        BitVector u;
        u.reset(nums*size);

        block* v_simd = new block[nums*simd_num];
        block* u_simd = new block[nums*simd_num];

        int sizebyte = size / 8;

        ot.send(seed);

        //#pragma omp parallel for
        for (int i = 0; i < nums; ++i) {
            int tmp_i = i * size;
            int tmp, tmp_k0, tmp_jk;

            BitVector temp0(size), temp1;
            // bool* temp0 = new bool [size]();
            // bool* temp1 = new bool [size]();
            for (int j = 0; j < size; ++j) {
                block seed_one = toBlock(0, seed[i][j]);
                // printf("[%d][%d], seed = %lX\n",i,j,seed[i][j]);
                PRNG prng;
                prng.SetSeed(seed_one);
                prng.get<uint8_t>(temp0.data(), temp0.sizeBytes());
                for(int k = 0; k < simd_num; k++){
                    v_simd[i*simd_num + k] ^= *(temp0.blocks() + k);
                }
                
                
                temp1.copy(temp0, 0, size-j);
                temp1.append(temp0, j, size-j);
                // // printf("i=%d, &tmp=%d\n",i,&temp);
                // print(temp1,simd_size);
                for(int k = 0; k < simd_num; k++){
                    u_simd[i*simd_num + k] ^= *(temp1.blocks() + k);
                }
                
            }
            memcpy(output.data()+i*sizebyte, v_simd+i*simd_num, sizebyte);
            memcpy(u.data()+i*sizebyte, u_simd+i*simd_num, sizebyte);
        }

        for (int i = 0; i < bits.sizeBlocks(); ++i) {
            *(u.blocks() + i) ^= *(bits.blocks() + i);
        }

        // cout << u.sizeBytes() << endl;
        sync_wait(chl.send(u));
        sync_wait(chl.flush());

        free(u_simd);
        free(v_simd);

        co_return;
    }

    Proto PeqtReceiver::vose(BitVector &output, int nums, int size, PRNG& prng, Socket& chl, u64 mNumThreads) {
        int simd_size = 16;
        PRNG prg;
        uint64_t** seed = new uint64_t*[nums];
        for(int i = 0; i < nums; ++i) {
            seed[i] = new uint64_t[size]();
        }
        int simd_num = (size + 15)/simd_size;
        
        // TODO: 
        IknpOtExtSender *sender;
        IknpOtExtReceiver *recver;
        NCN1OT<uint64_t> ot(Role::Receiver, nums, size, sender, recver, chl);

        BitVector u;
        u.reset(nums*size);

        block* v_simd = new block[nums*simd_num];
        block* u_simd = new block[nums*simd_num];

        int sizebyte = size / 8;

        ot.recv(seed, vose_val);

        for (int i = 0; i < nums; ++i) {
            int tmp_i = i * size;
            int tmp, tmp_k0, tmp_jk;
            
            BitVector temp0(size), temp1;
            for (int j = 0; j < size; ++j) {
                block seed_one = toBlock(0, seed[i][j]);
                // printf("[%d][%d], seed = %lX\n",i,j,seed[i][j]);
                prg.SetSeed(seed_one);
                prg.get<uint8_t>(temp0.data(), temp0.sizeBytes());
                for(int k = 0; k < simd_num; k++){
                    v_simd[i*simd_num + k] ^= *(temp0.blocks() + k);
                }

                int tmp_jk = (j+size-vose_val[i])%size; // shold be +
                
                temp1.copy(temp0, size-tmp_jk, tmp_jk);
                temp1.append(temp0, size-tmp_jk);
                // memcpy(temp1 + ((tmp_jk)), temp0, ((size-tmp_jk)) * sizeof(bool));
                // memcpy(temp1, temp0 + ((size-tmp_jk)), (tmp_jk) * sizeof(bool));
                // // printf("i=%d, &tmp=%d\n",i,&temp);
                // print(temp1,simd_size);
                for(int k = 0; k < simd_num; k++){
                    v_simd[i*simd_num + k] ^= *(temp1.blocks() + k);
                }
            }
            memcpy(output.data()+i*sizebyte, v_simd+i*simd_num, sizebyte);
        }

        u.reset(nums*size);
        // cout << u.sizeBytes() << endl;
        sync_wait(chl.recv(u));
        sync_wait(chl.flush());

        for (int i = 0; i < nums; ++i) {
            int tmp_i = i*size;
            for (int j = 0; j < size; ++j) {
                output[tmp_i + j] = u[tmp_i + ((j + vose_val[i])%size)] - output[tmp_i + j]; //bug: XOR
            }
        }

        free(u_simd);
        free(v_simd);

        co_return;
    }

    Proto PeqtSender::run(BitVector& data, BitVector& output, Socket& chl, u64 mNumThreads){
        // Set mod p to be 2^[log(eqlength)], then need to compute mod during execution
        auto w = BitVector(mDataSize*mEqLength);
        auto w_recv = BitVector(mDataSize*mEqLength);
        std::vector<u8> modp_Share(mDataSize,0);
        std::vector<u8> modp_Share_recv(mDataSize,0);
        u64 modlength = oc::log2ceil(mEqLength);
        u64 mod = 1 << modlength;
        u64 mod_min_1 = mod-1;

        w = data ^ convert_bit;
        macoro::sync_wait(chl.send(w));
        macoro::sync_wait(chl.recv(w_recv));
        w = w ^ w_recv;

        for (u64 i=0; i < mDataSize; i++){
            for (u64 j=0; j < mEqLength; j++){
                if (w[i*mEqLength+j]==0)
                    modp_Share[i] += convert_val[i*mEqLength+j];
                else
                    modp_Share[i] += (1-convert_val[i*mEqLength+j]);
            }
            modp_Share[i] = modp_Share[i] + vose_val[i];
        }

        auto modp_Share_send = modp_Share;
        macoro::sync_wait(chl.send(std::move(modp_Share_send)));
        macoro::sync_wait(chl.recv(modp_Share_recv));
        for (u64 i=0; i < mDataSize; i++){
            u8 idx = (modp_Share_recv[i] + modp_Share[i]) & mod_min_1;
            output[i] = vose_table[i*mod + idx];
        }

        co_return;
    }

    Proto PeqtReceiver::run(BitVector& data, BitVector& output, Socket& chl, u64 mNumThreads){
        auto w = BitVector(mDataSize*mEqLength);
        auto w_recv = BitVector(mDataSize*mEqLength);
        std::vector<u8> modp_Share(mDataSize,0);
        std::vector<u8> modp_Share_recv(mDataSize,0);
        u64 modlength = oc::log2ceil(mEqLength);
        u64 mod = 1 << modlength;
        u64 mod_min_1 = mod-1;

        w = data ^ convert_bit;
        macoro::sync_wait(chl.recv(w_recv));
        macoro::sync_wait(chl.send(w));
        w = w ^ w_recv;

        for (u64 i=0; i < mDataSize; i++){
            for (u64 j=0; j < mEqLength; j++){
                if (w[i*mEqLength+j]==0)
                    modp_Share[i] += convert_val[i*mEqLength+j];
                else
                    modp_Share[i] -= convert_val[i*mEqLength+j];
            }
            modp_Share[i] = vose_val[i] + modp_Share[i];
        }

        auto modp_Share_send = modp_Share;
        macoro::sync_wait(chl.recv(modp_Share_recv));
        macoro::sync_wait(chl.send(std::move(modp_Share_send)));
        for (u64 i=0; i < mDataSize; i++){
            u8 idx = (modp_Share_recv[i] + modp_Share[i]) & mod_min_1;
            output[i] = vose_table[i*mod + idx];
        }

        co_return;
    }
    

}