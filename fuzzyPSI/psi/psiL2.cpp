#include "psi.h"
#include <array>
#include <future>

namespace CmpFuzzyPSI
{
    Proto FuzzyPsiSender::runL2(span<block> inputs, Socket& chl){
        Timer timer;
        u64 offlineComm = 0;
        u64 onlineComm = 0;
        u64 modL = mDelta * mDelta * mDim + 1;
        u64 modLLength = oc::log2ceil(modL);
        modL = 1 << modLLength;
        ModOps mod(modL);
        long long offlineTime = 0;
        long long onlineTime = 0;
        DEBUG_LOG("begin");

        auto Begin = timer.setTimePoint("FuzzyPsiSender::set-up begin");
        //fmap setup
        FmapSender mFmapSender;
        mFmapSender.setTimer(timer);
        macoro::sync_wait(mFmapSender.setUp(mSenderSize, mRecverSize, mDim, mDelta, mLorH, mPrng, chl, mNumThreads));
        //cuckoo setup
        block cuckooSeed = mPrng.get();
        sync_wait(chl.send(cuckooSeed));
        sync_wait(chl.flush());
        CuckooIndex<> cuckoo;
        cuckoo.init(mSenderSize, mSsp, 0, 3);
        u64 TableSize = cuckoo.mBins.size();
        //mIMT setup
        u64 Cmp_len = mFmapSender.orgSize;
        u64 m_min_1 = (1 << Cmp_len) - 1;
        mIMTSender mmIMTSender;
        mmIMTSender.setTimer(timer);
        macoro::sync_wait(mmIMTSender.setUp(TableSize, mDim, mDelta, mMetric, Cmp_len, mPrng, chl, mNumThreads));
        //Second mIMT setup
        mIMTSender mmIMTSender2;
        mmIMTSender2.setTimer(timer);
        macoro::sync_wait(mmIMTSender2.setUp(TableSize, 1, mDelta, 0, modLLength, mPrng, chl, mNumThreads));
        //PEQT setup
        u64 peqtLength = kappa + oc::log2ceil(TableSize);
        PeqtSender mPeqtSender;
        mPeqtSender.setTimer(timer);
        macoro::sync_wait(mPeqtSender.setUp(TableSize, peqtLength, mPrng, chl));
        //offline ole, ab=c+d mod L, toBeAdded, modL
        std::vector<u64> ole_b(TableSize * mDim, 0);
        std::vector<u64> ole_d(TableSize * mDim, 0);
        
        for (u64 i = 0; i < 5* mDim; ++i) {
            ole_b[i] = i;
            ole_d[i] = i+1;
        }

        //offline ROT output
        std::vector<std::array<block, 2>> rotOutput(TableSize * mDim);
        IknpOtExtSender senderOutput;
        sync_wait(senderOutput.send(rotOutput, mPrng, chl));
        sync_wait(chl.flush());
        auto End = timer.setTimePoint("FuzzyPsiSender::set-up end");

        offlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        offlineComm = chl.bytesSent() + chl.bytesReceived() - onlineComm;
        DEBUG_LOG("fmap setup done");

        Begin = timer.setTimePoint("FuzzyPsiSender::run-fuzzy mapping begin");
        std::vector<block> Identifiers(mSenderSize);
        std::vector<block> oringins(mSenderSize*mDim);
        macoro::sync_wait(mFmapSender.fuzzyMap(inputs, Identifiers, oringins, mPrng, chl, mNumThreads));
        End = timer.setTimePoint("FuzzyPsiSender::fuzzy mapping end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("fmap done");

        Begin = timer.setTimePoint("FuzzyPsiSender::build cuckoo table begin");

        cuckoo.insert(Identifiers, cuckooSeed);

        std::vector<block> CuckooKeys(mSenderSize); // store bin index

        for (u64 i = 0; i < TableSize; ++i) {
            if (!cuckoo.mBins[i].isEmpty()) {
                u64 b = cuckoo.mBins[i].idx();
                CuckooKeys[b] = _mm_set_epi64x(i, ((u64*)&Identifiers[b])[0]); // i||identifier
            }
        }
        End = timer.setTimePoint("FuzzyPsiSender::build cuckoo table end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("cuckoo done");

        Begin = timer.setTimePoint("FuzzyPsiSender::run-opprf begin");
        u64 valsize = peqtLength - mDim -1 + Cmp_len * mDim + modLLength * mDim;
        u64 valsize_byte = (valsize + 7) / 8;
        std::vector<u8> Opprf_val_data(mSenderSize * valsize_byte);// length d vector in each posotion, length of vector is peqtLength-d + Cmp_len*d
        MatrixView<u8> Opprf_val(Opprf_val_data.data(), mSenderSize, valsize_byte);
        auto v_s = BitVector(mmIMTSender.mShareSize); // copy the opprf value from TableSize*d*cmp_len to 3*TableSize*d*cmp_len
        auto masks = BitVector(TableSize*(peqtLength-mDim-1));
        std::vector<u64> cmp_inputs_bin(TableSize*mDim*2);//[y-delta,y+delta] for L_2
        std::vector<u64> b_prime(TableSize*mDim);
        std::vector<u64> distance_share(2*TableSize, 0);

        DEBUG_LOG("Opprf inputs defined.");
        RsOpprfReceiver mOpprfReceiver;
        mOpprfReceiver.setTimer(timer);
        macoro::sync_wait(mOpprfReceiver.receive(mRecverSize*3, CuckooKeys, Opprf_val, mPrng, mNumThreads, chl));

        for (u64 i=0; i<TableSize; i++){
            // copy the opprf value to v_s
            if (!cuckoo.mBins[i].isEmpty()) {
                u64 b = cuckoo.mBins[i].idx();
                for (u64 j=0; j < mDim; j++){
                    for (u64 k=0; k< Cmp_len; k++){
                        u64 byteIndex = (j*Cmp_len + k) / 8;
                        u64 bitIndex = (j*Cmp_len + k) % 8;
                        bool bit = (Opprf_val(b, byteIndex) >> bitIndex) & 1;
                        v_s[i*2*mDim*Cmp_len + 2*j*Cmp_len + k] = bit;
                        v_s[i*2*mDim*Cmp_len + (2*j+1)*Cmp_len + k] = bit;
                    }

                    u64 y_hat = ((u64*)&inputs[b*mDim + j])[0];
                    u64 y_org = ((u64*)&oringins[b*mDim + j])[0];
                    if ((y_hat >= y_org) && (y_hat - y_org >= mDelta))
                        cmp_inputs_bin[2*i*mDim + 2*j] = (y_hat -mDelta) & m_min_1; 
                    // & m_min_1 can be ommited, since M=2^{Cmp_len}, and (y_hat -mDelta) mod M = ((y_hat -mDelta) mod 2^64) mod M
                    else 
                        cmp_inputs_bin[2*i*mDim + 2*j] = y_org & m_min_1;

                    if ((y_hat >= y_org) && (y_hat - y_org + mDelta >= m_min_1))
                        cmp_inputs_bin[2*i*mDim + 2*j + 1] = (m_min_1 + y_org) & m_min_1;
                    else 
                        cmp_inputs_bin[2*i*mDim + 2*j + 1] = (y_hat + mDelta) & m_min_1;
                }

                u64 preidx = Cmp_len * mDim;
                for (u64 j= 0; j < mDim; j++){
                    u64 a_prime = 0;
                    for (u64 k=0; k < modLLength; k++){
                        u64 byteIndex = (j*modLLength+k+preidx) / 8;
                        u64 bitIndex = (j*modLLength+k+preidx) % 8;

                        bool bit = (Opprf_val(b, byteIndex) >> bitIndex) & 1;
                        a_prime ^= (bit << k);
                    }
                    // if (b<1)
                    //     std::cout << b << " dim " << j << ": " << (a_prime % modL) << " " << ((u64*)&inputs[b*mDim + j])[0] << " " <<
                    //     mod.sub(a_prime, ((u64*)&inputs[b*mDim + j])[0]) << std::endl;

                    b_prime[i*mDim+j] = mod.sub(a_prime, ((u64*)&inputs[b*mDim + j])[0]);
                    u64 tmp = mod.sub(b_prime[i*mDim+j]*b_prime[i*mDim+j], 2*ole_d[i*mDim+j]);
                    distance_share[2*i+1] = mod.add(distance_share[2*i+1], tmp);
                    b_prime[i*mDim+j] = mod.sub(b_prime[i*mDim+j], ole_b[i*mDim+j]);
                }

                preidx = Cmp_len * mDim + modLLength * mDim;
                for (u64 j=preidx; j < valsize; j++){
                    u64 byteIndex = j / 8;
                    u64 bitIndex = j % 8;
                    bool bit = (Opprf_val(b, byteIndex) >> bitIndex) & 1;
                    masks[i*(peqtLength - mDim - 1) + (j - preidx)] = bit;
                }
            } else {
                //randomize the empty bins, v_s[i*2*mDim*Cmp_len---(i+1)*2*mDim*Cmp_len], cmp_inputs_bin[i*mDim---(i+1)*mDim]
            }
            distance_share[2*i] = mod.sub(distance_share[2*i+1], mDelta*mDelta);
        }
        End = timer.setTimePoint("FuzzyPsiSender::run-opprf end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("Opprf done.");

        Begin = timer.setTimePoint("FuzzyPsiSender::run both mIMT begin");
        auto output = BitVector();
        macoro::sync_wait(mmIMTSender.run(cmp_inputs_bin, v_s, output, chl, mNumThreads));

        macoro::sync_wait(chl.send(std::move(b_prime)));
        BitVector v_s_2(TableSize * modLLength);
        sync_wait(chl.recv(v_s_2));
        BitVector v_s_2_input(2 * TableSize * modLLength);
        for (u64 i=1; i < TableSize; i++){
            for (u64 j=1; j < modLLength; j++){
                v_s_2_input[2*(i*modLLength) + j] = v_s_2[i*modLLength+j];
                v_s_2_input[2*(i*modLLength) + j + modLLength] = v_s_2[i*modLLength+j];
            }
        }
        auto output2 = BitVector();
        macoro::sync_wait(mmIMTSender2.run(distance_share, v_s_2_input, output2, chl, mNumThreads));
        End = timer.setTimePoint("FuzzyPsiSender::run both mIMT  end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("Both mIMT done.");

        Begin = timer.setTimePoint("FuzzyPsiSender::run PEQT begin");

        auto peqtSend = BitVector(TableSize * peqtLength);
        auto peqtoutput = BitVector(TableSize);
        for (u64 i=0; i<TableSize; i++){
            for (u64 j=0; j<peqtLength; j++){
                if (j < (peqtLength - mDim -1)){
                    peqtSend[i*peqtLength + j] = masks[i*(peqtLength - mDim-1) + j];
                } else if (j>=(peqtLength - mDim -1) && j < (peqtLength - 1)) {
                    peqtSend[i*peqtLength + j] = output[i*mDim + (j - (peqtLength - mDim - 1))];
                } else {
                    peqtSend[i*peqtLength + j] = output2[i];
                }
            }
        }
        macoro::sync_wait(mPeqtSender.run(peqtSend, peqtoutput, chl));

        End = timer.setTimePoint("FuzzyPsiSender::run PEQT end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("PEQT done.");

        //output OT begin
        Begin = timer.setTimePoint("FuzzyPsiSender:: Output OT begin");
        auto chosbitsRecv= BitVector(TableSize*mDim);
        macoro::sync_wait(chl.recv(chosbitsRecv));
        std::vector<block> maskedData(2 * TableSize * mDim);
        for (u64 i=0; i<TableSize; i++){
            if (!cuckoo.mBins[i].isEmpty()) {
                u64 b = cuckoo.mBins[i].idx();
                for (u64 j=0; j<mDim; j++){
                    if (chosbitsRecv[i*mDim + j] ^ peqtoutput[i]){
                        maskedData[2*(i*mDim+j)] = rotOutput[i*mDim+j][0] ^ inputs[b*mDim+j];
                        maskedData[2*(i*mDim+j)+1] = rotOutput[i*mDim+j][1];
                    }
                    else{
                        maskedData[2*(i*mDim+j)] = rotOutput[i*mDim+j][0];
                        maskedData[2*(i*mDim+j)+1] = rotOutput[i*mDim+j][1] ^ inputs[b*mDim+j];
                    }
                }
            } else {
                for (u64 j=0; j<mDim; j++){
                    maskedData[2*(i*mDim+j)] = rotOutput[i*mDim+j][0];
                    maskedData[2*(i*mDim+j)+1] = rotOutput[i*mDim+j][1];
                }
            }
        }
        macoro::sync_wait(chl.send(std::move(maskedData)));
        macoro::sync_wait(chl.flush());
        
        End = timer.setTimePoint("FuzzyPsiSender:: Output OT end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;

        std::cout <<"Offline time: " << offlineTime << " ms " << "Online time: " << onlineTime << " ms " << std::endl;
        std::cout << "Offline comm: " << offlineComm << " Bytes " << "Online comm: " << onlineComm << " Bytes " << std::endl;
        // DEBUG_LOG("Offline time: " << offlineTime << " ms " << "Online time: " << onlineTime << " ms ");
        // DEBUG_LOG("Offline comm: " << offlineComm << " Bytes " << "Online comm: " << onlineComm << " Bytes ");
    
        // for (u64 i=0; i<TableSize; i++){
        //     // copy the opprf value to v_s
        //     if (!cuckoo.mBins[i].isEmpty()) {
        //         u64 b = cuckoo.mBins[i].idx();
        //         if (b < 5){
        //             std::cout << "bin" << i << " item " << b << ": " ;
        //             for (u64 j=0; j<peqtLength; j++){
        //                 std::cout <<  peqtSend[i*peqtLength + j];
        //             }
        //             // std::cout << distance_share[2*i] << " " << distance_share[2*i+1];
        //             std::cout << std::endl;
        //         }
        //     }
        // }
        co_return;
    }

    Proto FuzzyPsiReceiver::runL2(span<block> inputs, Socket& chl){
        u64 offlineComm = 0;
        u64 onlineComm = 0;
        u64 modL = mDelta * mDelta * mDim + 1;
        u64 modLLength = oc::log2ceil(modL);
        modL = 1 << modLLength;
        ModOps mod(modL);
        long long offlineTime = 0;
        long long onlineTime = 0;
        Timer timer;


        DEBUG_LOG("begin");

        auto Begin = timer.setTimePoint("FuzzyPsiReceiver::set-up begin");
        //fmap setup
        FmapReceiver mFmapReceiver;
        mFmapReceiver.setTimer(timer);
        macoro::sync_wait(mFmapReceiver.setUp(mSenderSize, mRecverSize, mDim, mDelta, mLorH, mPrng, chl, mNumThreads));
        //simple setup
        block cuckooSeed;
        sync_wait(chl.recv(cuckooSeed));
        auto params = oc::CuckooIndex<>::selectParams(mRecverSize, mSsp, 0, 3);
        u64 TableSize = params.numBins();
        SimpleIndex sIdx;
        sIdx.init(TableSize, mSenderSize, mSsp, 3);
        //mIMT setup
        u64 Cmp_len = mFmapReceiver.orgSize;
        mIMTReceiver mmIMTReceiver;
        mmIMTReceiver.setTimer(timer);
        macoro::sync_wait(mmIMTReceiver.setUp(TableSize, mDim, mDelta, mMetric, Cmp_len, mPrng, chl, mNumThreads));
        //Second mIMT setup
        mIMTReceiver mmIMTReceiver2;
        mmIMTReceiver2.setTimer(timer);
        macoro::sync_wait(mmIMTReceiver2.setUp(TableSize, 1, mDelta, 0, modLLength, mPrng, chl, mNumThreads));
        //PEQT setup
        u64 peqtLength = kappa + oc::log2ceil(TableSize);
        PeqtReceiver mPeqtReceiver;
        mPeqtReceiver.setTimer(timer);
        macoro::sync_wait(mPeqtReceiver.setUp(TableSize, peqtLength, mPrng, chl));
        //offline ole, ab=c+d mod L, toBeAdded, modL
        std::vector<u64> ole_a(TableSize * mDim, 0);
        std::vector<u64> ole_c(TableSize * mDim, 0);

        for (u64 i = 0; i < 5* mDim; ++i) {
            ole_a[i] = i;
            ole_c[i] = i*i - i - 1;
        }
        //offline Output ROT
        std::vector<block> rotOutput(TableSize * mDim);
        BitVector sOutput(TableSize * mDim);
        IknpOtExtReceiver recverOutput;
        mPrng.get(sOutput.data(), sOutput.sizeBytes());
        sync_wait(recverOutput.receive(sOutput, rotOutput, mPrng, chl));
        sync_wait(chl.flush());
        auto End = timer.setTimePoint("FuzzyPsiReceiver::set-up end");

        offlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        offlineComm = chl.bytesSent() + chl.bytesReceived() - onlineComm;
        DEBUG_LOG("fmap setup done");

        Begin = timer.setTimePoint("FuzzyPsiReceiver::run-fuzzy mapping begin");
        std::vector<block> Identifiers(mRecverSize);
        // std::vector<block> oringins(mRecverSize*mDim);
        macoro::sync_wait(mFmapReceiver.fuzzyMap(inputs, Identifiers, mPrng, chl, mNumThreads));
        End = timer.setTimePoint("FuzzyPsiReceiver::fuzzy mapping end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("fmap done");

        // build simple table
        // DEBUG_LOG("Fmap Communication: " << chl.bytesSent() + chl.bytesReceived()  << " bytes");
        Begin = timer.setTimePoint("FuzzyPsiSender::build simple table begin");
        sIdx.insertItems(Identifiers, cuckooSeed);
        End = timer.setTimePoint("FuzzyPsiSender::build simple table end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("simple done");

        Begin = timer.setTimePoint("FuzzyPsiReceiver::run-opprf begin");
        u64 valsize = peqtLength - mDim -1 + Cmp_len * mDim + modLLength * mDim;
        u64 valsize_byte = (valsize + 7) / 8;
        std::vector<u8> Opprf_val_data(Identifiers.size()*valsize_byte*3, 0);// length d vector in each posotion, length of vector is peqtLength-d + Cmp_len*d
        MatrixView<u8> Opprf_val(Opprf_val_data.data(), Identifiers.size()*3, valsize_byte);
        std::vector<block> Opprf_key(Identifiers.size() * 3); // store bin index
        auto masks = BitVector(TableSize*(peqtLength-mDim-1));
        mPrng.get((u8*)masks.data(), masks.sizeBytes());

        DEBUG_LOG("OPPRF inputs defined.");
        for (u64 i = 0; i < TableSize; ++i){
            auto bin = sIdx.mBins[i];
            auto size = sIdx.mBinSizes[i];
            for (u64 p = 0; p < size; ++p){
                auto hidx = bin[p].hashIdx();
                auto vidx = bin[p].idx();
                Opprf_key[vidx * 3 + hidx] = _mm_set_epi64x(i, ((u64*)&Identifiers[vidx])[0]); // i||identifier

                for (u64 j=0; j < Cmp_len * mDim; j++){
                    u64 byteIndex = j / 8;
                    u64 bitIndex = j % 8;
                    u64 dimIndex = j / Cmp_len;
                    u64 bitPosindex = j % Cmp_len;

                    bool bit = (((u64*)&inputs[vidx*mDim + dimIndex])[0] >> bitPosindex) & 1;
                    bit ^= mmIMTReceiver.e[2*i*mDim*Cmp_len + 2*j];
                    Opprf_val(vidx * 3 + hidx, byteIndex) ^= bit << bitIndex;
                }
                u64 preidx = Cmp_len * mDim;
                for (u64 j= 0; j < mDim; j++){
                    u64 tmp = mod.add(((u64*)&inputs[vidx*mDim + j])[0], ole_a[i*mDim + j]);
                    for (u64 k=0; k < modLLength; k++){
                        u64 byteIndex = (j*modLLength+k+preidx) / 8;
                        u64 bitIndex = (j*modLLength+k+preidx) % 8;
                        bool bit = (tmp >> k) & 1;
                        Opprf_val(vidx * 3 + hidx, byteIndex) ^= bit << bitIndex;
                    }

                    // if (vidx < 1)
                    //     std::cout << vidx << " dim " << j  << ": " << ((u64*)&inputs[vidx*mDim + j])[0] << " " <<
                    //     mod.add(((u64*)&inputs[vidx*mDim + j])[0], ole_a[i*mDim + j]) << std::endl;

                }

                preidx = Cmp_len * mDim + modLLength * mDim;
                for (u64 j= preidx; j < valsize; j++){
                    u64 byteIndex = j / 8;
                    u64 bitIndex = j % 8;

                    bool bit = masks[i*(peqtLength - mDim -1) + j - preidx];
                    Opprf_val(vidx * 3 + hidx, byteIndex) ^= bit << bitIndex;
                }
            }
        }
        DEBUG_LOG("Opprf inputs prepared.");

        RsOpprfSender mOpprfSender;
        mOpprfSender.setTimer(timer);
        // DEBUG_LOG("Opprf inputs size: " << Opprf_key.size() << " value rows: " << Opprf_val.rows() << " value cols: " << Opprf_val.cols());
        macoro::sync_wait(mOpprfSender.send(mSenderSize, Opprf_key, Opprf_val, mPrng, mNumThreads, chl));
        End = timer.setTimePoint("FuzzyPsiReceiver::run-opprf end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("Opprf done.");

        Begin = timer.setTimePoint("FuzzyPsiSender::run both mIMT begin");
        auto output = BitVector();
        macoro::sync_wait(mmIMTReceiver.run(output, chl, mNumThreads));

        std::vector<u64> b_prime(TableSize*mDim);
        std::vector<u64> distance_share(TableSize, 0);
        BitVector v_s_2(TableSize * modLLength);
        macoro::sync_wait(chl.recv(b_prime));
        for (u64 i = 0; i < TableSize; ++i){
            for (u64 j= 0; j < mDim; j++){
                distance_share[i] = mod.sub(distance_share[i], ole_a[i*mDim+j]*ole_a[i*mDim+j]);
                distance_share[i] = mod.sub(distance_share[i], 2*ole_c[i*mDim+j]);
                distance_share[i] = mod.add(distance_share[i], 2*ole_a[i*mDim+j]*b_prime[i*mDim+j]);
            }
            for (u64 k=0; k< modLLength; k++){
                v_s_2[i*modLLength+k] = (distance_share[i] >> k) & 1;
            }
        }
        sync_wait(chl.send(std::move(v_s_2)));
        auto output2 = BitVector();
        macoro::sync_wait(mmIMTReceiver2.run(output2, chl, mNumThreads));

        End = timer.setTimePoint("FuzzyPsiSender::run both mIMT  end");

        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("Both mIMT done.");

        // DEBUG_LOG(output.size() << " " << output2.size() << " " << TableSize);

        Begin = timer.setTimePoint("FuzzyPsiSender::run PEQT begin");
        auto peqtSend = BitVector(TableSize * peqtLength);
        auto peqtoutput = BitVector(TableSize);
        for (u64 i=0; i<TableSize; i++){
            for (u64 j=0; j<peqtLength; j++){
                if (j < (peqtLength - mDim -1)){
                    peqtSend[i*peqtLength + j] = masks[i*(peqtLength - mDim - 1) + j];
                } else if (j>=(peqtLength - mDim -1) && j < (peqtLength - 1)) {
                    peqtSend[i*peqtLength + j] = !output[i*mDim + (j - (peqtLength - mDim - 1))];
                } else {
                    peqtSend[i*peqtLength + j] = !output2[i];
                }
            }
        }
        macoro::sync_wait(mPeqtReceiver.run(peqtSend, peqtoutput, chl));

        End = timer.setTimePoint("FuzzyPsiSender::run PEQT end");
        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
        DEBUG_LOG("PEQT done.");

        //output OT begin
        Begin = timer.setTimePoint("FuzzyPsiSender:: Output OT begin");
        auto chosbitsSend = BitVector(TableSize*mDim);
        for (u64 i=0; i<TableSize; i++){
            for (u64 j=0; j<mDim; j++){
                chosbitsSend[i*mDim+j] = peqtoutput[i] ^ sOutput[i*mDim+j];
            }
        }
        macoro::sync_wait(chl.send(std::move(chosbitsSend)));
        std::vector<block> maskedData(2 * TableSize * mDim);
        macoro::sync_wait(chl.recv(maskedData));

        std::vector<block> outputPSI;
        for (u64 i=0; i<TableSize; i++){
            bool inSet = false;
            std::vector<block> outputRecv(mDim);
            for (u64 j=0; j<mDim; j++){
                outputRecv[j] = maskedData[2*(i*mDim+j) + sOutput[i*mDim+j]] ^ rotOutput[i*mDim+j];
                if (outputRecv[j] != ZeroBlock)
                    inSet = true;
            }
            if (inSet)
                outputPSI.insert(outputPSI.end(), outputRecv.begin(), outputRecv.end());
        }

        End = timer.setTimePoint("FuzzyPsiSender:: Output OT end");
        onlineTime += std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin).count();
        onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;

        std::cout << "Intersection Size: " << (outputPSI.size()/mDim) << std::endl;

        std::cout <<"Offline time: " << offlineTime << " ms " << "Online time: " << onlineTime << " ms " << std::endl;
        std::cout << "Offline comm: " << offlineComm << " Bytes " << "Online comm: " << onlineComm << " Bytes " << std::endl;
        // DEBUG_LOG("Offline time: " << offlineTime << " ms " << "Online time: " << onlineTime << " ms ");
        // DEBUG_LOG("Offline comm: " << offlineComm << " Bytes " << "Online comm: " << onlineComm << " Bytes ");
    
        // for (u64 i = 0; i < TableSize; ++i){
        //     auto bin = sIdx.mBins[i];
        //     auto size = sIdx.mBinSizes[i];
        //     for (u64 p = 0; p < size; ++p){
        //         if (bin[p].idx() < 5){
        //             std::cout << "bin" << i << " item " << bin[p].idx() << ": " ;
        //             for (u64 j=0; j<peqtLength; j++){
        //                 std::cout <<  peqtSend[i*peqtLength + j];
        //             }
        //             // std::cout << distance_share[i];
        //             std::cout << std::endl;
        //         }
        //     }
        // }
        co_return;
    }
}