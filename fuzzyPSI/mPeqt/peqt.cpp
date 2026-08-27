#include "peqt.h"
#include "sspeqt_support.h"

namespace CmpFuzzyPSI {
namespace {

u32 peqtModulus(u64 equalityLength) {
  u32 modulus = 1u << oc::log2ceil(equalityLength);
  if (modulus == equalityLength)
    modulus <<= 1;
  return modulus;
}

} // namespace

Proto PeqtSender::setUp(u64 dataSize, u64 equalityLength, PRNG &prng,
                        Socket &chl) {
  mDataSize = dataSize;
  mEqLength = equalityLength;
  const u32 modulus = peqtModulus(mEqLength);

  B2AOfflineSender(mDataSize * mEqLength, modulus, chl, prng, rShare,
                    tShare);
  VoseOfflineSender(modulus, mDataSize, prng, chl, U, V);

  vose_val.resize(mDataSize);
  for (u32 i = 0; i < mDataSize; ++i)
    vose_val[i] = prng.get<u32>() % modulus;
  VoseOnlineSender(modulus, mDataSize, U, V, vose_val, chl, vose_table);
  co_return;
}

Proto PeqtReceiver::setUp(u64 dataSize, u64 equalityLength, PRNG &prng,
                          Socket &chl) {
  mDataSize = dataSize;
  mEqLength = equalityLength;
  const u32 modulus = peqtModulus(mEqLength);

  B2AOfflineReceiver(mDataSize * mEqLength, modulus, chl, prng, rShare,
                      tShare);
  VoseOfflineReceiver(modulus, mDataSize, prng, chl, eps, W);
  VoseOnlineReceiver(modulus, mDataSize, W, eps, chl, vose_table);
  co_return;
}

Proto PeqtSender::run(BitVector &data, BitVector &output, Socket &chl,
                      u64 numThreads) {
  (void)numThreads;
  const u32 modulus = peqtModulus(mEqLength);
  const u32 mask = modulus - 1;
  std::vector<u64> bitShares;
  B2AOnline(mDataSize * mEqLength, modulus, true, data, rShare, tShare, chl,
            bitShares);

  std::vector<u8> local(mDataSize, 0), remote(mDataSize);
  for (u64 i = 0; i < mDataSize; ++i) {
    for (u64 j = 0; j < mEqLength; ++j)
      local[i] = static_cast<u8>((local[i] + bitShares[i * mEqLength + j]) & mask);
    local[i] = static_cast<u8>((local[i] + vose_val[i]) & mask);
  }
  coproto::sync_wait(chl.send(coproto::span<u8>(local)));
  coproto::sync_wait(chl.recv(coproto::span<u8>(remote)));
  coproto::sync_wait(chl.flush());

  output.resize(mDataSize);
  for (u64 i = 0; i < mDataSize; ++i)
    output[i] = vose_table[i * modulus + ((local[i] + remote[i]) & mask)];
  co_return;
}

Proto PeqtReceiver::run(BitVector &data, BitVector &output, Socket &chl,
                        u64 numThreads) {
  (void)numThreads;
  const u32 modulus = peqtModulus(mEqLength);
  const u32 mask = modulus - 1;
  std::vector<u64> bitShares;
  B2AOnline(mDataSize * mEqLength, modulus, false, data, rShare, tShare, chl,
            bitShares);

  std::vector<u8> local(mDataSize, 0), remote(mDataSize);
  for (u64 i = 0; i < mDataSize; ++i) {
    for (u64 j = 0; j < mEqLength; ++j)
      local[i] = static_cast<u8>((local[i] + bitShares[i * mEqLength + j]) & mask);
    local[i] = static_cast<u8>((local[i] + eps[i]) & mask);
  }
  coproto::sync_wait(chl.recv(coproto::span<u8>(remote)));
  coproto::sync_wait(chl.send(coproto::span<u8>(local)));
  coproto::sync_wait(chl.flush());

  output.resize(mDataSize);
  for (u64 i = 0; i < mDataSize; ++i)
    output[i] = vose_table[i * modulus + ((local[i] + remote[i]) & mask)];
  co_return;
}

} // namespace CmpFuzzyPSI
