#include "psi.h"
#include <array>
#include <future>

namespace CmpFuzzyPSI {
Proto FuzzyPsiSender::runL1(span<block> inputs, Socket &chl) {
  Timer timer;
  u64 offlineComm = 0;
  [[maybe_unused]] u64 onlineComm = 0;
  u64 modL = ((mDelta * mDim + 1) % 2 == 0) ? (mDelta * mDim + 2)
                                            : (mDelta * mDim + 1);
  u64 modLLength = oc::log2ceil(modL) + 1;
  ModOps mod(modL);
  [[maybe_unused]] long long offlineTime = 0;
  [[maybe_unused]] long long onlineTime = 0;
  DEBUG_LOG("begin");

  auto offlineBegin =
      timer.setTimePoint("FuzzyPsiSender::offline begin");
  // famp setup
  FmapSender mFmapSender;
  mFmapSender.setTimer(timer);
  macoro::sync_wait(mFmapSender.setUp(mSenderSize, mRecverSize, mDim, mDelta,
                                      mPrng, chl, mNumThreads));
  // cuckoo setup
  block cuckooSeed = mPrng.get();
  sync_wait(chl.send(cuckooSeed));
  sync_wait(chl.flush());
  CuckooIndex<> cuckoo;
  cuckoo.init(mSenderSize * mFmapSender.myExpansionRate, mSsp, 0, 3);
  u64 TableSize = cuckoo.mBins.size();
  // mIMT setup
  u64 Cmp_len = mFmapSender.orgSize;
  u64 m_min_1 = (1 << Cmp_len) - 1;
  mIMTSender mmIMTSender;
  mmIMTSender.setTimer(timer);
  macoro::sync_wait(mmIMTSender.setUp(TableSize, mDim, mDelta, mMetric, Cmp_len,
                                      mPrng, chl, mNumThreads));
  // Second mIMT setup
  mIMTSender mmIMTSender2;
  mmIMTSender2.setTimer(timer);
  macoro::sync_wait(mmIMTSender2.setUp(TableSize, 1, mDelta, 0, modLLength,
                                       mPrng, chl, mNumThreads));
  // PEQT setup
  u64 peqtLength = kappa + oc::log2ceil(TableSize);
  PeqtSender mPeqtSender;
  mPeqtSender.setTimer(timer);
  macoro::sync_wait(mPeqtSender.setUp(TableSize, peqtLength, mPrng, chl));
  // offline ROT
  std::vector<std::vector<u64>> rot_prime(TableSize * mDim,
                                          std::vector<u64>(2, 0));
  std::vector<u64> rot(TableSize * mDim, 0);
  BitVector s(TableSize * mDim);
  IknpOtExtSender sender;
  IknpOtExtReceiver recver;

  mPrng.get(s.data(), s.sizeBytes());
  vector<block> data_recv(TableSize * mDim);
  sync_wait(recver.receive(s, data_recv, mPrng, chl));

  std::vector<std::array<block, 2>> data_send(TableSize * mDim);
  sync_wait(sender.send(data_send, mPrng, chl));

  for (u64 i = 0; i < TableSize * mDim; i++) {
    rot_prime[i][0] = ((u64 *)&data_send[i][0])[0];
    rot_prime[i][1] = ((u64 *)&data_send[i][1])[0];
    rot[i] = ((u64 *)&data_recv[i])[0];
  }

  std::vector<u64> rot0prime_sub_rot1prime(TableSize * mDim);
  for (u64 i = 0; i < TableSize * mDim; i++) {
    rot0prime_sub_rot1prime[i] = mod.sub(rot_prime[i][0], rot_prime[i][1]);
    // rot0prime_sub_rot1prime[i] = mod.div2(rot0prime_sub_rot1prime[i]);
  }
  std::vector<u64> rot_sum_d(TableSize, 0);
  for (u64 i = 0; i < TableSize; i++) {
    for (u64 j = 0; j < mDim; j++) {
      rot_sum_d[i] = mod.add(rot_sum_d[i], rot[i * mDim + j]);
    }
  }
  // offline ROT output
  std::vector<std::array<block, 2>> rotOutput(TableSize * mDim);
  IknpOtExtSender senderOutput;
  sync_wait(senderOutput.send(rotOutput, mPrng, chl));
  sync_wait(chl.flush());
  auto offlineEnd = timer.setTimePoint("FuzzyPsiSender::offline end");
  offlineTime = std::chrono::duration_cast<std::chrono::milliseconds>(
                    offlineEnd - offlineBegin)
                    .count();
  offlineComm = chl.bytesSent() + chl.bytesReceived();
  DEBUG_LOG("fmap setup done");

  auto onlineBegin = timer.setTimePoint("FuzzyPsiSender::online begin");

  // auto fuzzyMappingBegin =
  //     timer.setTimePoint("FuzzyPsiSender::fuzzy mapping begin");
  // auto fuzzyMappingCommBegin = chl.bytesSent() + chl.bytesReceived();

  std::vector<block> Identifiers(mSenderSize * mFmapSender.myExpansionRate); // store the identifier for each input, length is mSenderSize*expansion rate
  std::vector<block> oringins(mSenderSize * mDim * mFmapSender.myExpansionRate); // store the original vector for each input, length is mSenderSize*dim*expansion rate
  macoro::sync_wait(mFmapSender.fuzzyMap(inputs, Identifiers, oringins, mPrng,
                                         chl, mNumThreads));

  // auto fuzzyMappingEnd =
  //     timer.setTimePoint("FuzzyPsiSender::fuzzy mapping end");
  // auto fuzzyMappingTime =
  //     std::chrono::duration_cast<std::chrono::milliseconds>(
  //         fuzzyMappingEnd - fuzzyMappingBegin)
  //         .count();
  // auto fuzzyMappingComm = chl.bytesSent() + chl.bytesReceived() -
  //                         fuzzyMappingCommBegin;

  DEBUG_LOG("fmap done");

  // auto cuckooTableBuildBegin =
  //     timer.setTimePoint("FuzzyPsiSender::cuckoo table build begin");
  // auto cuckooTableBuildCommBegin = chl.bytesSent() + chl.bytesReceived();


  cuckoo.insert(Identifiers, cuckooSeed);
  std::vector<block> CuckooKeys(mSenderSize * mFmapSender.myExpansionRate); // store bin index
  for (u64 i = 0; i < TableSize; ++i) {
    if (!cuckoo.mBins[i].isEmpty()) {
      u64 b = cuckoo.mBins[i].idx();
      CuckooKeys[b] =
          _mm_set_epi64x(i, ((u64 *)&Identifiers[b])[0]); // i||identifier
    }
  }

  // auto cuckooTableBuildEnd =
  //     timer.setTimePoint("FuzzyPsiSender::cuckoo table build end");
  // auto cuckooTableBuildTime =
  //     std::chrono::duration_cast<std::chrono::milliseconds>(
  //         cuckooTableBuildEnd - cuckooTableBuildBegin)
  //         .count();
  // auto cuckooTableBuildComm = chl.bytesSent() + chl.bytesReceived() -
  //                             cuckooTableBuildCommBegin;

  DEBUG_LOG("cuckoo done");

  // auto opprfBegin = timer.setTimePoint("FuzzyPsiSender::OPPRF begin");
  // auto opprfCommBegin = chl.bytesSent() + chl.bytesReceived();

  u64 valsize = peqtLength - mDim - 1 + Cmp_len * mDim;
  u64 valsize_byte = (valsize + 7) / 8;
  std::vector<u8> Opprf_val_data(
      mSenderSize * valsize_byte * mFmapSender.myExpansionRate); // length d vector in each posotion, length
                                   // of vector is peqtLength-d + Cmp_len*d
  MatrixView<u8> Opprf_val(Opprf_val_data.data(), mSenderSize * mFmapSender.myExpansionRate, valsize_byte);
  auto v_s = BitVector(
      mmIMTSender.mShareSize); // copy the opprf value from TableSize*d*cmp_len
                               // to 3*TableSize*d*cmp_len
  auto masks = BitVector(TableSize * (peqtLength - mDim - 1));
  std::vector<u64> cmp_inputs_bin(TableSize * mDim *
                                  3); // y-delta, y ,y+delta for L1

  DEBUG_LOG("Opprf inputs defined.");

  RsOpprfReceiver mOpprfReceiver;
  mOpprfReceiver.setTimer(timer);
  macoro::sync_wait(mOpprfReceiver.receive(mRecverSize * 3 * mFmapSender.anotherExpansionRate, CuckooKeys,
                                           Opprf_val, mPrng, mNumThreads, chl));

  for (u64 i = 0; i < TableSize; i++) {
    // copy the opprf value to v_s
    if (!cuckoo.mBins[i].isEmpty()) {
      u64 b = cuckoo.mBins[i].idx();
      u64 dataindex = b / mFmapSender.myExpansionRate;
      for (u64 j = 0; j < mDim; j++) {
        for (u64 k = 0; k < Cmp_len; k++) {
          u64 byteIndex = (j * Cmp_len + k) / 8;
          u64 bitIndex = (j * Cmp_len + k) % 8;
          bool bit = (Opprf_val(b, byteIndex) >> bitIndex) & 1;
          v_s[i * 3 * mDim * Cmp_len + 3 * j * Cmp_len + k] = bit;
          v_s[i * 3 * mDim * Cmp_len + (3 * j + 1) * Cmp_len + k] = bit;
          v_s[i * 3 * mDim * Cmp_len + (3 * j + 2) * Cmp_len + k] = bit;
        }

        u64 y_hat = ((u64 *)&inputs[dataindex * mDim + j])[0];
        u64 y_org = ((u64 *)&oringins[b * mDim + j])[0];
        if (y_hat - y_org >= mDelta)
          cmp_inputs_bin[3 * i * mDim + 3 * j] = (y_hat - mDelta) & m_min_1;
        else
          cmp_inputs_bin[3 * i * mDim + 3 * j] = y_org & m_min_1;

        cmp_inputs_bin[3 * i * mDim + 3 * j + 1] = y_hat & m_min_1;

        if (y_hat - y_org + mDelta > m_min_1)
          cmp_inputs_bin[3 * i * mDim + 3 * j + 2] =
              (m_min_1 + y_org) & m_min_1;
        else
          cmp_inputs_bin[3 * i * mDim + 3 * j + 2] = (y_hat + mDelta) & m_min_1;
      }

      for (u64 j = mDim * Cmp_len; j < valsize; j++) {
        u64 byteIndex = j / 8;
        u64 bitIndex = j % 8;
        bool bit = (Opprf_val(b, byteIndex) >> bitIndex) & 1;
        masks[i * (peqtLength - mDim - 1) + (j - mDim * Cmp_len)] = bit;
      }
    } else {
      // randomize the empty bins, v_s[i*2*mDim*Cmp_len---(i+1)*2*mDim*Cmp_len],
      // cmp_inputs_bin[i*mDim---(i+1)*mDim]
    }
  }

  // auto opprfEnd = timer.setTimePoint("FuzzyPsiSender::OPPRF end");
  // auto opprfTime = std::chrono::duration_cast<std::chrono::milliseconds>(
  //                      opprfEnd - opprfBegin)
  //                      .count();
  // auto opprfComm =
  //     chl.bytesSent() + chl.bytesReceived() - opprfCommBegin;

  DEBUG_LOG("Opprf done.");

  // auto imtBegin = timer.setTimePoint("FuzzyPsiSender::mIMT begin");
  // auto imtCommBegin = chl.bytesSent() + chl.bytesReceived();

  auto output = BitVector();
  macoro::sync_wait(
      mmIMTSender.run(cmp_inputs_bin, v_s, output, chl, mNumThreads));

  // auto imtEnd = timer.setTimePoint("FuzzyPsiSender::mIMT end");
  // auto imtTime = std::chrono::duration_cast<std::chrono::milliseconds>(
  //                    imtEnd - imtBegin)
  //                    .count();
  // auto imtComm = chl.bytesSent() + chl.bytesReceived() - imtCommBegin;

  DEBUG_LOG("mIMT done.");

  // auto secondOpprfAndOtBegin =
  //     timer.setTimePoint("FuzzyPsiSender::second OPPRF and OT begin");
  // auto secondOpprfAndOtCommBegin = chl.bytesSent() + chl.bytesReceived();

  auto e_prime_s = BitVector(TableSize * mDim);
  auto t_i_prime = BitVector(TableSize * mDim);
  auto t_i = BitVector(TableSize * mDim);
  for (u64 i = 0; i < TableSize * mDim; ++i) {
    t_i[i] = output[2 * i + 1] ^ s[i];
    e_prime_s[i] = output[2 * i];
  }
  sync_wait(chl.recv(t_i_prime));
  sync_wait(chl.send(std::move(t_i)));

  u64 valsize_prime = modLLength * mDim;
  u64 valsize_prime_byte = (valsize_prime + 7) / 8;
  std::vector<u8> Opprf_val_data_prime(
      mSenderSize * mFmapSender.anotherExpansionRate *
      valsize_prime_byte); // length d vector in each posotion, length of vector
                           // is peqtLength-d + Cmp_len*d
  MatrixView<u8> Opprf_val_prime(Opprf_val_data.data(), mSenderSize * mFmapSender.anotherExpansionRate,
                                 valsize_prime_byte);
  std::vector<u64> m_hat_r(TableSize * mDim);
  std::vector<u64> m_hat_s(TableSize * mDim, 0);
  std::vector<u8> m_hat_s_message(TableSize * valsize_prime_byte);
  std::vector<u64> distance_share(2 * TableSize, 0); // input to mIMT

  RsOpprfReceiver mOpprfReceiver2;
  mOpprfReceiver2.setTimer(timer);
  macoro::sync_wait(mOpprfReceiver2.receive(
      mRecverSize * 3 * mFmapSender.anotherExpansionRate, CuckooKeys, Opprf_val_prime, mPrng, mNumThreads, chl));

  DEBUG_LOG("Second OPPRF done.");
  for (u64 i = 0; i < TableSize; i++) {
    // copy the opprf value to v_s
    if (!cuckoo.mBins[i].isEmpty()) {
      u64 b = cuckoo.mBins[i].idx();
      u64 dataindex = b / mFmapSender.myExpansionRate;
      for (u64 j = 0; j < mDim; j++) {
        for (u64 k = 0; k < modLLength; k++) {
          u64 byteIndex = (j * modLLength + k) / 8;
          u64 bitIndex = (j * modLLength + k) % 8;
          bool bit = (Opprf_val_prime(b, byteIndex) >> bitIndex) & 1;
          m_hat_r[i * mDim + j] ^= (u64)bit << k;
        }

        u64 m_i_s =
            mod.sub(((u64 *)&inputs[dataindex * mDim + j])[0], m_hat_r[i * mDim + j]);

        if (output[2 * (i * mDim + j) + 1]) {
          distance_share[2 * i + 1] = mod.add(distance_share[2 * i + 1], m_i_s);
          // distance_share[2*i+1] = mod.add(distance_share[2*i+1], m_i_s);
          m_hat_s[i * mDim + j] = mod.add(m_hat_s[i * mDim + j], 2 * m_i_s);
        } else {
          distance_share[2 * i + 1] = mod.sub(distance_share[2 * i + 1], m_i_s);
          // distance_share[2*i+1] = mod.sub(distance_share[2*i+1], m_i_s);
          m_hat_s[i * mDim + j] = mod.sub(m_hat_s[i * mDim + j], 2 * m_i_s);
        }

        if (t_i_prime[i * mDim + j]) {
          distance_share[2 * i + 1] =
              mod.sub(distance_share[2 * i + 1], rot_prime[i * mDim + j][0]);
          m_hat_s[i * mDim + j] = mod.sub(
              m_hat_s[i * mDim + j], rot0prime_sub_rot1prime[i * mDim + j]);
        } else {
          distance_share[2 * i + 1] =
              mod.sub(distance_share[2 * i + 1], rot_prime[i * mDim + j][1]);
          m_hat_s[i * mDim + j] = mod.add(
              m_hat_s[i * mDim + j], rot0prime_sub_rot1prime[i * mDim + j]);
        }

        for (u64 k = 0; k < modLLength; k++) {
          u64 byteIndex = (j * modLLength + k) / 8;
          u64 bitIndex = (j * modLLength + k) % 8;
          bool bit = (m_hat_s[i * mDim + j] >> k) & 1;
          m_hat_s_message[i * valsize_prime_byte + byteIndex] ^= (u64)bit
                                                                 << bitIndex;
        }
      }
    } else {
      distance_share[2 * i + 1] = mPrng.get();
    }
    distance_share[2 * i + 1] =
        mod.sub(distance_share[2 * i + 1], rot_sum_d[i]);
    distance_share[2 * i] = mod.sub(distance_share[2 * i + 1], mDelta);
  }
  sync_wait(chl.send(std::move(m_hat_s_message)));

  BitVector v_s_2(TableSize * modLLength);
  BitVector v_s_2_input(2 * TableSize * modLLength);
  sync_wait(chl.recv(v_s_2));
  for (u64 i = 0; i < TableSize; i++) {
    for (u64 j = 0; j < modLLength; j++) {
      v_s_2_input[2 * (i * modLLength) + j] = v_s_2[i * modLLength + j];
      v_s_2_input[2 * (i * modLLength) + j + modLLength] =
          v_s_2[i * modLLength + j];
    }
  }
  auto output2 = BitVector();
  macoro::sync_wait(
      mmIMTSender2.run(distance_share, v_s_2_input, output2, chl, mNumThreads));

  // auto secondOpprfAndOtEnd =
  //     timer.setTimePoint("FuzzyPsiSender::second OPPRF and OT end");
  // auto secondOpprfAndOtTime =
  //     std::chrono::duration_cast<std::chrono::milliseconds>(
  //         secondOpprfAndOtEnd - secondOpprfAndOtBegin)
  //         .count();
  // auto secondOpprfAndOtComm = chl.bytesSent() + chl.bytesReceived() -
  //                             secondOpprfAndOtCommBegin;

  DEBUG_LOG("Second mIMT done.");

  // auto peqtBegin = timer.setTimePoint("FuzzyPsiSender::PEQT begin");
  // auto peqtCommBegin = chl.bytesSent() + chl.bytesReceived();


  auto peqtSend = BitVector(TableSize * peqtLength);
  auto peqtoutput = BitVector(TableSize);
  for (u64 i = 0; i < TableSize; i++) {
    for (u64 j = 0; j < peqtLength; j++) {
      if (j < (peqtLength - mDim - 1)) {
        peqtSend[i * peqtLength + j] = masks[i * (peqtLength - mDim - 1) + j];
      } else if (j >= (peqtLength - mDim - 1) && j < (peqtLength - 1)) {
        peqtSend[i * peqtLength + j] =
            output[2 * (i * mDim + (j - (peqtLength - mDim - 1)))];
      } else {
        peqtSend[i * peqtLength + j] = output2[i];
      }
    }
  }
  macoro::sync_wait(mPeqtSender.run(peqtSend, peqtoutput, chl));

  // auto peqtEnd = timer.setTimePoint("FuzzyPsiSender::PEQT end");
  // auto peqtTime = std::chrono::duration_cast<std::chrono::milliseconds>(
  //                     peqtEnd - peqtBegin)
  //                     .count();
  // auto peqtComm = chl.bytesSent() + chl.bytesReceived() - peqtCommBegin;

  DEBUG_LOG("PEQT done.");

  // output OT begin

  // auto outputOtBegin = timer.setTimePoint("FuzzyPsiSender::output OT begin");
  // auto outputOtCommBegin = chl.bytesSent() + chl.bytesReceived();

  auto chosbitsRecv = BitVector(TableSize * mDim);
  macoro::sync_wait(chl.recv(chosbitsRecv));
  std::vector<block> maskedData(2 * TableSize * mDim);
  for (u64 i = 0; i < TableSize; i++) {
    if (!cuckoo.mBins[i].isEmpty()) {
      u64 b = cuckoo.mBins[i].idx();
      u64 dataindex = b / mFmapSender.myExpansionRate;
      for (u64 j = 0; j < mDim; j++) {
        if (chosbitsRecv[i * mDim + j] ^ peqtoutput[i]) {
          maskedData[2 * (i * mDim + j)] =
              rotOutput[i * mDim + j][0] ^ inputs[dataindex * mDim + j];
          maskedData[2 * (i * mDim + j) + 1] = rotOutput[i * mDim + j][1];
        } else {
          maskedData[2 * (i * mDim + j)] = rotOutput[i * mDim + j][0];
          maskedData[2 * (i * mDim + j) + 1] =
              rotOutput[i * mDim + j][1] ^ inputs[dataindex * mDim + j];
        }
      }
    } else {
      for (u64 j = 0; j < mDim; j++) {
        maskedData[2 * (i * mDim + j)] = rotOutput[i * mDim + j][0];
        maskedData[2 * (i * mDim + j) + 1] = rotOutput[i * mDim + j][1];
      }
    }
  }
  macoro::sync_wait(chl.send(std::move(maskedData)));
  macoro::sync_wait(chl.flush());

  // auto outputOtEnd = timer.setTimePoint("FuzzyPsiSender::output OT end");
  // auto outputOtTime = std::chrono::duration_cast<std::chrono::milliseconds>(
  //                         outputOtEnd - outputOtBegin)
  //                         .count();
  // auto outputOtComm =
  //     chl.bytesSent() + chl.bytesReceived() - outputOtCommBegin;

  auto onlineEnd = timer.setTimePoint("FuzzyPsiSender::online end");
  onlineTime = std::chrono::duration_cast<std::chrono::milliseconds>(
                   onlineEnd - onlineBegin)
                   .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;

  co_return;
}

Proto FuzzyPsiReceiver::runL1(span<block> inputs, Socket &chl) {
  u64 offlineComm = 0;
  u64 onlineComm = 0;
  u64 modL = ((mDelta * mDim + 1) % 2 == 0) ? (mDelta * mDim + 2)
                                            : (mDelta * mDim + 1);
  u64 modLLength = oc::log2ceil(modL) + 1;
  ModOps mod(modL);
  long long offlineTime = 0;
  long long onlineTime = 0;
  Timer timer;

  DEBUG_LOG("begin");

  auto offlineBegin =
      timer.setTimePoint("FuzzyPsiReceiver::offline begin");
  // famp setup
  FmapReceiver mFmapReceiver;
  mFmapReceiver.setTimer(timer);
  macoro::sync_wait(mFmapReceiver.setUp(mSenderSize, mRecverSize, mDim, mDelta,
                                        mPrng, chl, mNumThreads));
  // simple setup
  block cuckooSeed;
  sync_wait(chl.recv(cuckooSeed));
  auto params = oc::CuckooIndex<>::selectParams(mSenderSize * mFmapReceiver.anotherExpansionRate, mSsp, 0, 3);
  u64 TableSize = params.numBins();
  SimpleIndex sIdx;
  sIdx.init(TableSize, mRecverSize * mFmapReceiver.myExpansionRate, mSsp, 3);
  // mIMT setup
  u64 Cmp_len = mFmapReceiver.orgSize;
  mIMTReceiver mmIMTReceiver;
  mmIMTReceiver.setTimer(timer);
  macoro::sync_wait(mmIMTReceiver.setUp(TableSize, mDim, mDelta, mMetric,
                                        Cmp_len, mPrng, chl, mNumThreads));
  // Second mIMT setup
  mIMTReceiver mmIMTReceiver2;
  mmIMTReceiver2.setTimer(timer);
  macoro::sync_wait(mmIMTReceiver2.setUp(TableSize, 1, mDelta, 0, modLLength,
                                         mPrng, chl, mNumThreads));
  // PEQT setup
  u64 peqtLength = kappa + oc::log2ceil(TableSize);
  PeqtReceiver mPeqtReceiver;
  mPeqtReceiver.setTimer(timer);
  macoro::sync_wait(mPeqtReceiver.setUp(TableSize, peqtLength, mPrng, chl));
  // offline ROT
  std::vector<std::vector<u64>> rot(TableSize * mDim, std::vector<u64>(2, 0));
  std::vector<u64> rot_prime(TableSize * mDim, 0);
  BitVector s_prime(TableSize * mDim);
  IknpOtExtSender sender;
  IknpOtExtReceiver recver;

  std::vector<std::array<block, 2>> data_send(TableSize * mDim);
  sync_wait(sender.send(data_send, mPrng, chl));

  mPrng.get(s_prime.data(), s_prime.sizeBytes());
  vector<block> data_recv(TableSize * mDim);
  sync_wait(recver.receive(s_prime, data_recv, mPrng, chl));

  for (u64 i = 0; i < TableSize * mDim; i++) {
    rot[i][0] = ((u64 *)&data_send[i][0])[0];
    rot[i][1] = ((u64 *)&data_send[i][1])[0];
    rot_prime[i] = ((u64 *)&data_recv[i])[0];
  }

  std::vector<u64> rot0_sub_rot1(TableSize * mDim);
  std::vector<u64> rot0_add_rot1(TableSize * mDim);
  for (u64 i = 0; i < TableSize * mDim; i++) {
    rot0_sub_rot1[i] = mod.sub(rot[i][0], rot[i][1]);
    rot0_sub_rot1[i] = mod.div2(rot0_sub_rot1[i]);
    rot0_add_rot1[i] = mod.add(rot[i][0], rot[i][1]);
    rot0_add_rot1[i] = mod.div2(rot0_add_rot1[i]);
  }
  std::vector<u64> rot_prime_sum_d(TableSize, 0);
  for (u64 i = 0; i < TableSize; i++) {
    for (u64 j = 0; j < mDim; j++) {
      rot_prime_sum_d[i] = mod.add(rot_prime_sum_d[i], rot_prime[i * mDim + j]);
    }
  }
  // offline Output ROT
  std::vector<block> rotOutput(TableSize * mDim);
  BitVector sOutput(TableSize * mDim);
  IknpOtExtReceiver recverOutput;
  mPrng.get(sOutput.data(), sOutput.sizeBytes());
  sync_wait(recverOutput.receive(sOutput, rotOutput, mPrng, chl));
  sync_wait(chl.flush());
  auto offlineEnd = timer.setTimePoint("FuzzyPsiReceiver::offline end");
  offlineTime = std::chrono::duration_cast<std::chrono::milliseconds>(
                    offlineEnd - offlineBegin)
                    .count();
  offlineComm = chl.bytesSent() + chl.bytesReceived();
  DEBUG_LOG("fmap setup done");

  auto onlineBegin = timer.setTimePoint("FuzzyPsiReceiver::online begin");

  // auto fuzzyMappingBegin =
  //     timer.setTimePoint("FuzzyPsiReceiver::fuzzy mapping begin");
  // auto fuzzyMappingCommBegin = chl.bytesSent() + chl.bytesReceived();

  std::vector<block> Identifiers(mRecverSize * mFmapReceiver.myExpansionRate); // store the identifier for each input, length is mRecverSize*expansion rate
  // std::vector<block> oringins(mRecverSize*mDim);
  macoro::sync_wait(
      mFmapReceiver.fuzzyMap(inputs, Identifiers, mPrng, chl, mNumThreads));

  // auto fuzzyMappingEnd =
  //     timer.setTimePoint("FuzzyPsiReceiver::fuzzy mapping end");
  // auto fuzzyMappingTime =
  //     std::chrono::duration_cast<std::chrono::milliseconds>(
  //         fuzzyMappingEnd - fuzzyMappingBegin)
  //         .count();
  // auto fuzzyMappingComm = chl.bytesSent() + chl.bytesReceived() -
  //                         fuzzyMappingCommBegin;

  DEBUG_LOG("fmap done");

  // auto simpleTableBuildBegin =
  //     timer.setTimePoint("FuzzyPsiReceiver::simple table build begin");
  // auto simpleTableBuildCommBegin = chl.bytesSent() + chl.bytesReceived();

  sIdx.insertItems(Identifiers, cuckooSeed);

  // auto simpleTableBuildEnd =
  //     timer.setTimePoint("FuzzyPsiReceiver::simple table build end");
  // auto simpleTableBuildTime =
  //     std::chrono::duration_cast<std::chrono::milliseconds>(
  //         simpleTableBuildEnd - simpleTableBuildBegin)
  //         .count();
  // auto simpleTableBuildComm = chl.bytesSent() + chl.bytesReceived() -
  //                             simpleTableBuildCommBegin;

  DEBUG_LOG("simple done");

  // auto opprfBegin = timer.setTimePoint("FuzzyPsiReceiver::OPPRF begin");
  // auto opprfCommBegin = chl.bytesSent() + chl.bytesReceived();

  u64 valsize = peqtLength - mDim - 1 + Cmp_len * mDim;
  u64 valsize_byte = (valsize + 7) / 8;
  std::vector<u8> Opprf_val_data(Identifiers.size() * valsize_byte * 3,
                                 0); // length d vector in each posotion, length
                                     // of vector is peqtLength-d + Cmp_len*d
  MatrixView<u8> Opprf_val(Opprf_val_data.data(), Identifiers.size() * 3,
                           valsize_byte);
  std::vector<block> Opprf_key(Identifiers.size() * 3); // store bin index
  auto masks = BitVector(TableSize * (peqtLength - mDim - 1));
  mPrng.get((u8 *)masks.data(), masks.sizeBytes());

  DEBUG_LOG("OPPRF inputs defined.");
  for (u64 i = 0; i < TableSize; ++i) {
    auto bin = sIdx.mBins[i];
    auto size = sIdx.mBinSizes[i];
    for (u64 p = 0; p < size; ++p) {
      auto hidx = bin[p].hashIdx();
      auto vidx = bin[p].idx();
      u64 dataindex = vidx / mFmapReceiver.myExpansionRate;
      Opprf_key[vidx * 3 + hidx] =
          _mm_set_epi64x(i, ((u64 *)&Identifiers[vidx])[0]); // i||identifier

      for (u64 j = 0; j < Cmp_len * mDim; j++) {
        u64 byteIndex = j / 8;
        u64 bitIndex = j % 8;
        u64 dimIndex = j / Cmp_len;
        u64 bitPosindex = j % Cmp_len;

        bool bit =
            (((u64 *)&inputs[dataindex * mDim + dimIndex])[0] >> bitPosindex) & 1;
        bit ^= mmIMTReceiver.e[3 * i * mDim * Cmp_len + 3 * j];
        Opprf_val(vidx * 3 + hidx, byteIndex) ^= bit << bitIndex;
      }
      u64 preidx = Cmp_len * mDim;
      for (u64 j = preidx; j < valsize; j++) {
        u64 byteIndex = j / 8;
        u64 bitIndex = j % 8;

        bool bit = masks[i * (peqtLength - mDim - 1) + j - preidx];
        Opprf_val(vidx * 3 + hidx, byteIndex) ^= bit << bitIndex;
      }
    }
  }
  DEBUG_LOG("Opprf inputs prepared.");

  RsOpprfSender mOpprfSender;
  mOpprfSender.setTimer(timer);
  macoro::sync_wait(mOpprfSender.send(mSenderSize * mFmapReceiver.anotherExpansionRate, Opprf_key, Opprf_val, mPrng,
                                      mNumThreads, chl));

  // auto opprfEnd = timer.setTimePoint("FuzzyPsiReceiver::OPPRF end");
  // auto opprfTime = std::chrono::duration_cast<std::chrono::milliseconds>(
  //                      opprfEnd - opprfBegin)
  //                      .count();
  // auto opprfComm =
  //     chl.bytesSent() + chl.bytesReceived() - opprfCommBegin;

  DEBUG_LOG("Opprf done.");

  // auto imtBegin = timer.setTimePoint("FuzzyPsiReceiver::mIMT begin");
  // auto imtCommBegin = chl.bytesSent() + chl.bytesReceived();

  auto output = BitVector();
  macoro::sync_wait(mmIMTReceiver.run(output, chl, mNumThreads));

  // auto imtEnd = timer.setTimePoint("FuzzyPsiReceiver::mIMT end");
  // auto imtTime = std::chrono::duration_cast<std::chrono::milliseconds>(
  //                    imtEnd - imtBegin)
  //                    .count();
  // auto imtComm = chl.bytesSent() + chl.bytesReceived() - imtCommBegin;

  DEBUG_LOG("mIMT done.");

  // auto secondOpprfAndOtBegin =
  //     timer.setTimePoint("FuzzyPsiReceiver::second OPPRF and OT begin");
  // auto secondOpprfAndOtCommBegin = chl.bytesSent() + chl.bytesReceived();

  auto e_prime_r = BitVector(TableSize * mDim);
  auto t_i_prime = BitVector(TableSize * mDim);
  auto t_i = BitVector(TableSize * mDim);
  for (u64 i = 0; i < TableSize * mDim; ++i) {
    t_i_prime[i] = output[2 * i + 1] ^ s_prime[i];
    e_prime_r[i] = output[2 * i];
  }
  sync_wait(chl.send(std::move(t_i_prime)));
  sync_wait(chl.recv(t_i));

  u64 valsize_prime = modLLength * mDim;
  u64 valsize_prime_byte = (valsize_prime + 7) / 8;
  std::vector<u8> Opprf_val_data_prime(
      Identifiers.size() * valsize_prime_byte * 3,
      0); // length d vector in each posotion, length of vector is peqtLength-d
          // + Cmp_len*d
  MatrixView<u8> Opprf_val_prime(Opprf_val_data_prime.data(),
                                 Identifiers.size() * 3, valsize_prime_byte);
  std::vector<u64> m_i_r(TableSize * mDim);

  DEBUG_LOG("Sencond OPPRF inputs defined.");
  for (u64 i = 0; i < TableSize; ++i) {
    auto bin = sIdx.mBins[i];
    auto size = sIdx.mBinSizes[i];
    for (u64 p = 0; p < size; ++p) {
      auto hidx = bin[p].hashIdx();
      auto vidx = bin[p].idx();
      u64 dataindex = vidx / mFmapReceiver.myExpansionRate;
      for (u64 j = 0; j < mDim; j++) {
        if (output[2 * (i * mDim + j) + 1] ^ t_i[i * mDim + j]) {
          m_i_r[i * mDim + j] = mod.sub(0, rot0_sub_rot1[i * mDim + j]);
        } else {
          m_i_r[i * mDim + j] = rot0_sub_rot1[i * mDim + j];
        }
        u64 tmp =
            mod.sub(((u64 *)&inputs[dataindex * mDim + j])[0], m_i_r[i * mDim + j]);

        for (u64 k = 0; k < modLLength; k++) {
          u64 byteIndex = (j * modLLength + k) / 8;
          u64 bitIndex = (j * modLLength + k) % 8;

          bool bit = (tmp >> k) & 1;
          Opprf_val_prime(vidx * 3 + hidx, byteIndex) ^= (bit << bitIndex);
        }
      }
    }
  }
  RsOpprfSender mOpprfSender2;
  mOpprfSender2.setTimer(timer);
  macoro::sync_wait(mOpprfSender2.send(mSenderSize * mFmapReceiver.anotherExpansionRate, Opprf_key, Opprf_val_prime,
                                       mPrng, mNumThreads, chl));

  DEBUG_LOG("Second Opprf done.");

  std::vector<u64> m_hat_s(TableSize * mDim, 0);
  std::vector<u8> m_hat_s_message(TableSize * valsize_prime_byte);
  std::vector<u64> distance_share(TableSize, 0);
  BitVector v_s_2(TableSize * modLLength);

  sync_wait(chl.recv(m_hat_s_message));
  // sync_wait(chl.recv(m_hat_s));
  for (u64 i = 0; i < TableSize; ++i) {
    for (u64 j = 0; j < mDim; j++) {
      for (u64 k = 0; k < modLLength; k++) {
        u64 byteIndex = (j * modLLength + k) / 8;
        u64 bitIndex = (j * modLLength + k) % 8;
        bool bit = (m_hat_s_message[i * valsize_prime_byte + byteIndex] >>
                    (bitIndex)) &
                   1;
        m_hat_s[i * mDim + j] ^= (bit << k);
      }
      distance_share[i] =
          mod.sub(distance_share[i], rot0_add_rot1[i * mDim + j]);

      if (!output[2 * (i * mDim + j) + 1]) {
        distance_share[i] = mod.add(distance_share[i], m_hat_s[i * mDim + j]);
      }
    }

    distance_share[i] = mod.sub(distance_share[i], rot_prime_sum_d[i]);

    for (u64 k = 0; k < modLLength; k++) {
      bool bit = (distance_share[i] >> k) & 1;
      bit ^= mmIMTReceiver2.e[2 * i * modLLength + 2 * k];
      v_s_2[i * modLLength + k] = bit;
    }
  }

  sync_wait(chl.send(std::move(v_s_2)));
  auto output2 = BitVector();
  macoro::sync_wait(mmIMTReceiver2.run(output2, chl, mNumThreads));

  // auto secondOpprfAndOtEnd =
  //     timer.setTimePoint("FuzzyPsiReceiver::second OPPRF and OT end");
  // auto secondOpprfAndOtTime =
  //     std::chrono::duration_cast<std::chrono::milliseconds>(
  //         secondOpprfAndOtEnd - secondOpprfAndOtBegin)
  //         .count();
  // auto secondOpprfAndOtComm = chl.bytesSent() + chl.bytesReceived() -
  //                             secondOpprfAndOtCommBegin;

  DEBUG_LOG("Second mIMT done.");

  // auto peqtBegin = timer.setTimePoint("FuzzyPsiReceiver::PEQT begin");
  // auto peqtCommBegin = chl.bytesSent() + chl.bytesReceived();

  auto peqtSend = BitVector(TableSize * peqtLength);
  auto peqtoutput = BitVector(TableSize);
  for (u64 i = 0; i < TableSize; i++) {
    for (u64 j = 0; j < peqtLength; j++) {
      if (j < (peqtLength - mDim - 1)) {
        peqtSend[i * peqtLength + j] = masks[i * (peqtLength - mDim - 1) + j];
      } else if (j >= (peqtLength - mDim - 1) && j < (peqtLength - 1)) {
        peqtSend[i * peqtLength + j] =
            !output[2 * (i * mDim + (j - (peqtLength - mDim - 1)))];
      } else {
        peqtSend[i * peqtLength + j] = !output2[i];
      }
    }
  }
  macoro::sync_wait(mPeqtReceiver.run(peqtSend, peqtoutput, chl));

  // auto peqtEnd = timer.setTimePoint("FuzzyPsiReceiver::PEQT end");
  // auto peqtTime = std::chrono::duration_cast<std::chrono::milliseconds>(
  //                     peqtEnd - peqtBegin)
  //                     .count();
  // auto peqtComm = chl.bytesSent() + chl.bytesReceived() - peqtCommBegin;

  DEBUG_LOG("PEQT done.");

  // output OT begin

  // auto outputOtBegin =
  //     timer.setTimePoint("FuzzyPsiReceiver::output OT begin");
  // auto outputOtCommBegin = chl.bytesSent() + chl.bytesReceived();

  auto chosbitsSend = BitVector(TableSize * mDim);
  for (u64 i = 0; i < TableSize; i++) {
    for (u64 j = 0; j < mDim; j++) {
      chosbitsSend[i * mDim + j] = peqtoutput[i] ^ sOutput[i * mDim + j];
    }
  }
  macoro::sync_wait(chl.send(std::move(chosbitsSend)));
  std::vector<block> maskedData(2 * TableSize * mDim);
  macoro::sync_wait(chl.recv(maskedData));

  std::vector<block> outputPSI;
  for (u64 i = 0; i < TableSize; i++) {
    bool inSet = false;
    std::vector<block> outputRecv(mDim);
    for (u64 j = 0; j < mDim; j++) {
      outputRecv[j] = maskedData[2 * (i * mDim + j) + sOutput[i * mDim + j]] ^
                      rotOutput[i * mDim + j];
      if (outputRecv[j] != ZeroBlock)
        inSet = true;
    }
    if (inSet)
      outputPSI.insert(outputPSI.end(), outputRecv.begin(), outputRecv.end());
  }

  // auto outputOtEnd =
  //     timer.setTimePoint("FuzzyPsiReceiver::output OT end");
  // auto outputOtTime = std::chrono::duration_cast<std::chrono::milliseconds>(
  //                         outputOtEnd - outputOtBegin)
  //                         .count();
  // auto outputOtComm =
  //     chl.bytesSent() + chl.bytesReceived() - outputOtCommBegin;

  auto onlineEnd = timer.setTimePoint("FuzzyPsiReceiver::online end");
  onlineTime = std::chrono::duration_cast<std::chrono::milliseconds>(
                   onlineEnd - onlineBegin)
                   .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;

  online_time = onlineTime;
  online_commu = onlineComm;
  offline_commu = offlineComm;
  offline_time = offlineTime;

  co_return;
}
} // namespace CmpFuzzyPSI
