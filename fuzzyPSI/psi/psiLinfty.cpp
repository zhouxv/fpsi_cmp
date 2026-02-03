#include "psi.h"
#include <array>
#include <future>

namespace CmpFuzzyPSI {
namespace {
struct NoHash {
  inline size_t operator()(const block &v) const { return v.get<size_t>(0); }
};
} // namespace

void details::FuzzyPsiBase::init(u64 senderSize, u64 recverSize,
                                 u64 statSecParam, u64 dim, u64 metric,
                                 u64 delta, u64 LorH, block seed,
                                 u64 numThreads, bool useReducedRounds) {
  mSenderSize = senderSize;
  mRecverSize = recverSize;
  mSsp = statSecParam;
  mDim = dim;
  mMetric = metric;
  mDelta = delta;
  mLorH = LorH;
  u64 buffersize = 256;
  if (senderSize >= 65536)
    // buffersize = u64(16384);
    buffersize = u64(1024);
  else if (senderSize >= 1048576)
    buffersize = u64(16384);
  mPrng.SetSeed(seed, buffersize);
  mNumThreads = numThreads;
  mUseReducedRounds = useReducedRounds;
}

Proto FuzzyPsiSender::run(span<block> inputs, Socket &chl) {
  if (mMetric == 0) {
    macoro::sync_wait(runLinfty(inputs, chl));
  } else if (mMetric == 1) {
    macoro::sync_wait(runL1(inputs, chl));
  } else if (mMetric == 2) {
    macoro::sync_wait(runL2(inputs, chl));
  } else {
    throw std::runtime_error("Unknown metric type");
  }

  co_return;
}

Proto FuzzyPsiReceiver::run(span<block> inputs, Socket &chl) {
  if (mMetric == 0) {
    macoro::sync_wait(runLinfty(inputs, chl));
  } else if (mMetric == 1) {
    macoro::sync_wait(runL1(inputs, chl));
  } else if (mMetric == 2) {
    macoro::sync_wait(runL2(inputs, chl));
  } else {
    throw std::runtime_error("Unknown metric type");
  }

  co_return;
}

Proto FuzzyPsiSender::runLinfty(span<block> inputs, Socket &chl) {

  Timer timer;
  u64 offlineComm = 0;
  u64 onlineComm = 0;
  long long offlineTime = 0;
  long long onlineTime = 0;
  DEBUG_LOG("begin");

  auto Begin = timer.setTimePoint("FuzzyPsiSender::set-up begin");
  // fmap setup
  FmapSender mFmapSender;
  mFmapSender.setTimer(timer);
  macoro::sync_wait(mFmapSender.setUp(mSenderSize, mRecverSize, mDim, mDelta,
                                      mLorH, mPrng, chl, mNumThreads));
  // cuckoo setup
  block cuckooSeed = mPrng.get();
  macoro::sync_wait(chl.send(cuckooSeed));
  macoro::sync_wait(chl.flush());
  CuckooIndex<> cuckoo;
  cuckoo.init(mRecverSize, mSsp, 0, 3);
  u64 TableSize = cuckoo.mBins.size();
  // mIMT setup
  u64 Cmp_len = mFmapSender.orgSize;
  u64 m_min_1 = (1 << Cmp_len) - 1;
  mIMTSender mmIMTSender;
  mmIMTSender.setTimer(timer);
  macoro::sync_wait(mmIMTSender.setUp(TableSize, mDim, mDelta, mMetric, Cmp_len,
                                      mPrng, chl, mNumThreads));
  // peqt setup
  u64 peqtLength = kappa + oc::log2ceil(TableSize);
  PeqtSender mPeqtSender;
  mPeqtSender.setTimer(timer);
  macoro::sync_wait(mPeqtSender.setUp(TableSize, peqtLength, mPrng, chl));

  // offline ROT output
  std::vector<std::array<block, 2>> rotOutput(TableSize * mDim);
  IknpOtExtSender senderOutput;
  sync_wait(senderOutput.send(rotOutput, mPrng, chl));
  sync_wait(chl.flush());
  auto End = timer.setTimePoint("FuzzyPsiSender::set-up end");

  offlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  offlineComm = chl.bytesSent() + chl.bytesReceived() - onlineComm;
  DEBUG_LOG("fmap setup done");

  Begin = timer.setTimePoint("FuzzyPsiSender::run-fuzzy mapping begin");
  std::vector<block> Identifiers(mSenderSize);
  std::vector<block> oringins(mSenderSize * mDim);
  macoro::sync_wait(mFmapSender.fuzzyMap(inputs, Identifiers, oringins, mPrng,
                                         chl, mNumThreads));
  End = timer.setTimePoint("FuzzyPsiSender::fuzzy mapping end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("fmap done");

  // build cuckoo table
  Begin = timer.setTimePoint("FuzzyPsiSender::build cuckoo table begin");
  cuckoo.insert(Identifiers, cuckooSeed);

  std::vector<block> CuckooKeys(mRecverSize); // store bin index

  for (u64 i = 0; i < TableSize; ++i) {
    if (!cuckoo.mBins[i].isEmpty()) {
      u64 b = cuckoo.mBins[i].idx();
      CuckooKeys[b] =
          _mm_set_epi64x(i, ((u64 *)&Identifiers[b])[0]); // i||identifier
    }
  }
  End = timer.setTimePoint("FuzzyPsiSender::build cuckoo table end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("cuckoo done");

  Begin = timer.setTimePoint("FuzzyPsiSender::run-opprf begin");
  u64 valsize = peqtLength - mDim + Cmp_len * mDim;
  u64 valsize_byte = (valsize + 7) / 8;
  std::vector<u8> Opprf_val_data(
      CuckooKeys.size() *
      valsize_byte); // length d vector in each posotion, length of vector is
                     // peqtLength-d + Cmp_len*d
  MatrixView<u8> Opprf_val(Opprf_val_data.data(), CuckooKeys.size(),
                           valsize_byte);
  auto v_s = BitVector(
      mmIMTSender.mShareSize); // copy the opprf value from TableSize*d*cmp_len
                               // to 2*TableSize*d*cmp_len
  auto masks = BitVector(TableSize * (peqtLength - mDim));
  std::vector<u64> cmp_inputs_bin(TableSize * mDim *
                                  2); //[y-delta,y+delta] for L_infty

  DEBUG_LOG("Opprf inputs defined.");

  RsOpprfReceiver mOpprfReceiver;
  mOpprfReceiver.setTimer(timer);
  macoro::sync_wait(mOpprfReceiver.receive(mRecverSize * 3, CuckooKeys,
                                           Opprf_val, mPrng, mNumThreads, chl));

  for (u64 i = 0; i < TableSize; i++) {
    // copy the opprf value to v_s
    if (!cuckoo.mBins[i].isEmpty()) {
      u64 b = cuckoo.mBins[i].idx();
      for (u64 j = 0; j < mDim; j++) {
        for (u64 k = 0; k < Cmp_len; k++) {
          u64 byteIndex = (j * Cmp_len + k) / 8;
          u64 bitIndex = (j * Cmp_len + k) % 8;
          bool bit = (Opprf_val(b, byteIndex) >> bitIndex) & 1;
          v_s[i * 2 * mDim * Cmp_len + 2 * j * Cmp_len + k] = bit;
          v_s[i * 2 * mDim * Cmp_len + (2 * j + 1) * Cmp_len + k] = bit;
        }

        u64 y_hat = ((u64 *)&inputs[b * mDim + j])[0];
        u64 y_org = ((u64 *)&oringins[b * mDim + j])[0];
        if ((y_hat >= y_org) && (y_hat - y_org >= mDelta))
          cmp_inputs_bin[2 * i * mDim + 2 * j] = (y_hat - mDelta) & m_min_1;
        else
          cmp_inputs_bin[2 * i * mDim + 2 * j] = y_org & m_min_1;

        if ((y_hat >= y_org) && (y_hat - y_org + mDelta >= m_min_1))
          cmp_inputs_bin[2 * i * mDim + 2 * j + 1] =
              (m_min_1 + y_org) & m_min_1;
        else
          cmp_inputs_bin[2 * i * mDim + 2 * j + 1] = (y_hat + mDelta) & m_min_1;
      }
      for (u64 j = mDim * Cmp_len; j < valsize; j++) {
        u64 byteIndex = j / 8;
        u64 bitIndex = j % 8;
        bool bit = (Opprf_val(b, byteIndex) >> bitIndex) & 1;
        masks[i * (peqtLength - mDim) + (j - mDim * Cmp_len)] = bit;
      }
    } else {
      // randomize the empty bins, v_s[i*2*mDim*Cmp_len---(i+1)*2*mDim*Cmp_len],
      // cmp_inputs_bin[i*mDim---(i+1)*mDim]
    }
  }
  End = timer.setTimePoint("FuzzyPsiSender::run-opprf end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("Opprf done.");

  Begin = timer.setTimePoint("FuzzyPsiSender::run-mIMT begin");
  auto output = BitVector();
  macoro::sync_wait(
      mmIMTSender.run(cmp_inputs_bin, v_s, output, chl, mNumThreads));
  End = timer.setTimePoint("FuzzyPsiSender::run-mIMT end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("mIMT done.");

  Begin = timer.setTimePoint("FuzzyPsiSender::run PEQT begin");

  auto peqtSend = BitVector(TableSize * peqtLength);
  auto peqtoutput = BitVector(TableSize);
  for (u64 i = 0; i < TableSize; i++) {
    for (u64 j = 0; j < peqtLength; j++) {
      if (j < (peqtLength - mDim)) {
        peqtSend[i * peqtLength + j] = masks[i * (peqtLength - mDim) + j];
      } else {
        peqtSend[i * peqtLength + j] =
            output[i * mDim + (j - (peqtLength - mDim))];
      }
    }
  }

  macoro::sync_wait(mPeqtSender.run(peqtSend, peqtoutput, chl));

  End = timer.setTimePoint("FuzzyPsiSender::run PEQT end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("PEQT done.");

  // output OT begin
  Begin = timer.setTimePoint("FuzzyPsiSender:: Output OT begin");
  auto chosbitsRecv = BitVector(TableSize * mDim);
  macoro::sync_wait(chl.recv(chosbitsRecv));
  std::vector<block> maskedData(2 * TableSize * mDim);
  for (u64 i = 0; i < TableSize; i++) {
    if (!cuckoo.mBins[i].isEmpty()) {
      u64 b = cuckoo.mBins[i].idx();
      for (u64 j = 0; j < mDim; j++) {
        if (chosbitsRecv[i * mDim + j] ^ peqtoutput[i]) {
          maskedData[2 * (i * mDim + j)] =
              rotOutput[i * mDim + j][0] ^ inputs[b * mDim + j];
          maskedData[2 * (i * mDim + j) + 1] = rotOutput[i * mDim + j][1];
        } else {
          maskedData[2 * (i * mDim + j)] = rotOutput[i * mDim + j][0];
          maskedData[2 * (i * mDim + j) + 1] =
              rotOutput[i * mDim + j][1] ^ inputs[b * mDim + j];
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

  End = timer.setTimePoint("FuzzyPsiSender:: Output OT end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;


  co_return;
}

Proto FuzzyPsiReceiver::runLinfty(span<block> inputs, Socket &chl) {

  u64 offlineComm = 0;
  u64 onlineComm = 0;
  long long offlineTime = 0;
  long long onlineTime = 0;
  Timer timer;

  DEBUG_LOG("begin");

  auto Begin = timer.setTimePoint("FuzzyPsiReceiver::set-up begin");
  // Fmap setup
  FmapReceiver mFmapReceiver;
  mFmapReceiver.setTimer(timer);
  macoro::sync_wait(mFmapReceiver.setUp(mSenderSize, mRecverSize, mDim, mDelta,
                                        mLorH, mPrng, chl, mNumThreads));
  // simple hash setup
  block cuckooSeed;
  sync_wait(chl.recv(cuckooSeed));
  auto params = oc::CuckooIndex<>::selectParams(mRecverSize, mSsp, 0, 3);
  u64 TableSize = params.numBins();
  SimpleIndex sIdx;
  sIdx.init(TableSize, mSenderSize, mSsp, 3);
  // mIMT set up
  u64 Cmp_len = mFmapReceiver.orgSize;
  // u64 m_min_1 = (1 << Cmp_len) - 1;
  mIMTReceiver mmIMTReceiver;
  mmIMTReceiver.setTimer(timer);
  macoro::sync_wait(mmIMTReceiver.setUp(TableSize, mDim, mDelta, mMetric,
                                        Cmp_len, mPrng, chl, mNumThreads));
  // peqt setup
  u64 peqtLength = kappa + oc::log2ceil(TableSize);
  PeqtReceiver mPeqtReceiver;
  mPeqtReceiver.setTimer(timer);
  macoro::sync_wait(mPeqtReceiver.setUp(TableSize, peqtLength, mPrng, chl));
  // offline Output ROT
  std::vector<block> rotOutput(TableSize * mDim);
  BitVector sOutput(TableSize * mDim);
  IknpOtExtReceiver recverOutput;
  mPrng.get(sOutput.data(), sOutput.sizeBytes());
  sync_wait(recverOutput.receive(sOutput, rotOutput, mPrng, chl));
  sync_wait(chl.flush());
  auto End = timer.setTimePoint("FuzzyPsiReceiver::set-up end");

  offlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  offlineComm = chl.bytesSent() + chl.bytesReceived() - onlineComm;
  DEBUG_LOG("fmap setup done");

  Begin = timer.setTimePoint("FuzzyPsiReceiver::run-fuzzy mapping begin");
  std::vector<block> Identifiers(mRecverSize);
  // std::vector<block> oringins(mRecverSize*mDim);
  macoro::sync_wait(
      mFmapReceiver.fuzzyMap(inputs, Identifiers, mPrng, chl, mNumThreads));
  End = timer.setTimePoint("FuzzyPsiReceiver::fuzzy mapping end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("fmap done");

  // build simple table
  DEBUG_LOG("Fmap Communication: " << chl.bytesSent() + chl.bytesReceived()
                                   << " bytes");
  Begin = timer.setTimePoint("FuzzyPsiSender::build simple table begin");
  sIdx.insertItems(Identifiers, cuckooSeed);
  End = timer.setTimePoint("FuzzyPsiSender::build simple table end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("simple done");

  Begin = timer.setTimePoint("FuzzyPsiReceiver::run-opprf begin");
  u64 valsize = peqtLength - mDim + Cmp_len * mDim;
  u64 valsize_byte = (valsize + 7) / 8;
  std::vector<u8> Opprf_val_data(Identifiers.size() * valsize_byte * 3,
                                 0); // length d vector in each posotion, length
                                     // of vector is peqtLength-d + Cmp_len*d
  MatrixView<u8> Opprf_val(Opprf_val_data.data(), Identifiers.size() * 3,
                           valsize_byte);
  std::vector<block> Opprf_key(Identifiers.size() * 3); // store bin index
  auto masks = BitVector(TableSize * (peqtLength - mDim));
  mPrng.get((u8 *)masks.data(), masks.sizeBytes());

  DEBUG_LOG("OPPRF inputs defined.");

  for (u64 i = 0; i < TableSize; ++i) {
    auto bin = sIdx.mBins[i];
    auto size = sIdx.mBinSizes[i];
    for (u64 p = 0; p < size; ++p) {
      auto hidx = bin[p].hashIdx();
      auto vidx = bin[p].idx();
      Opprf_key[vidx * 3 + hidx] =
          _mm_set_epi64x(i, ((u64 *)&Identifiers[vidx])[0]); // i||identifier

      for (u64 j = 0; j < Cmp_len * mDim; j++) {
        u64 byteIndex = j / 8;
        u64 bitIndex = j % 8;
        u64 dimIndex = j / Cmp_len;
        u64 bitPosindex = j % Cmp_len;

        bool bit =
            (((u64 *)&inputs[vidx * mDim + dimIndex])[0] >> bitPosindex) & 1;
        bit ^= mmIMTReceiver.e[2 * i * mDim * Cmp_len + 2 * j];

        Opprf_val(vidx * 3 + hidx, byteIndex) ^= bit << bitIndex;
      }

      for (u64 j = Cmp_len * mDim; j < valsize; j++) {
        u64 byteIndex = j / 8;
        u64 bitIndex = j % 8;

        bool bit = masks[i * (peqtLength - mDim) + j - Cmp_len * mDim];
        Opprf_val(vidx * 3 + hidx, byteIndex) ^= bit << bitIndex;
      }
    }
  }

  DEBUG_LOG("Opprf inputs prepared.");

  RsOpprfSender mOpprfSender;
  mOpprfSender.setTimer(timer);
  macoro::sync_wait(mOpprfSender.send(mSenderSize, Opprf_key, Opprf_val, mPrng,
                                      mNumThreads, chl));
  End = timer.setTimePoint("FuzzyPsiReceiver::run-opprf end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("Opprf done.");

  Begin = timer.setTimePoint("FuzzyPsiReceiver::run-mIMT begin");
  auto output = BitVector();
  macoro::sync_wait(mmIMTReceiver.run(output, chl, mNumThreads));
  End = timer.setTimePoint("FuzzyPsiReceiver::run-mIMT end");

  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("mIMT done.");

  Begin = timer.setTimePoint("FuzzyPsiSender::run PEQT begin");
  auto peqtSend = BitVector(TableSize * peqtLength);
  auto peqtoutput = BitVector(TableSize);
  for (u64 i = 0; i < TableSize; i++) {
    for (u64 j = 0; j < peqtLength; j++) {
      if (j < (peqtLength - mDim)) {
        peqtSend[i * peqtLength + j] = masks[i * (peqtLength - mDim) + j];
      } else {
        peqtSend[i * peqtLength + j] =
            !output[i * mDim + (j - (peqtLength - mDim))];
      }
    }
  }
  macoro::sync_wait(mPeqtReceiver.run(peqtSend, peqtoutput, chl));

  End = timer.setTimePoint("FuzzyPsiSender::run PEQT end");
  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;
  DEBUG_LOG("PEQT done.");
  // output OT begin
  Begin = timer.setTimePoint("FuzzyPsiSender:: Output OT begin");
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

  End = timer.setTimePoint("FuzzyPsiSender:: Output OT end");
  onlineTime +=
      std::chrono::duration_cast<std::chrono::milliseconds>(End - Begin)
          .count();
  onlineComm = chl.bytesSent() + chl.bytesReceived() - offlineComm;

  online_time = onlineTime;
  online_commu = onlineComm;
  offline_commu = offlineComm;
  offline_time = offlineTime;



  co_return;
}
} // namespace CmpFuzzyPSI
