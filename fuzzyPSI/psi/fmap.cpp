#include "fmap.h"

#include <cstring>

namespace CmpFuzzyPSI {

Proto FmapSender::setUp(u64 senderSize, u64 receiverSize, u64 dim, u64 delta,
                        PRNG &prng, Socket &chl, u64 mNumThreads) {
  (void)prng;
  (void)chl;
  (void)mNumThreads;
  mSenderSize = senderSize;
  mRecverSize = receiverSize;
  mDim = dim;
  mDelta = delta;
  myExpansionRate = mDim + 1;
  anotherExpansionRate = mDim + 1;
  orgSize = osuCrypto::log2ceil(mDelta * (mDim + 1));
  co_return;
}

Proto FmapReceiver::setUp(u64 senderSize, u64 receiverSize, u64 dim,
                          u64 delta, PRNG &prng, Socket &chl,
                          u64 mNumThreads) {
  (void)prng;
  (void)chl;
  (void)mNumThreads;
  mSenderSize = senderSize;
  mRecverSize = receiverSize;
  mDim = dim;
  mDelta = delta;
  myExpansionRate = mDim + 1;
  anotherExpansionRate = mDim + 1;
  orgSize = osuCrypto::log2ceil(mDelta * (mDim + 1));
  co_return;
}

Proto FmapSender::fuzzyMap(span<block> inputs, span<block> identifiers,
                           span<block> origins, PRNG &prng, Socket &chl,
                           u64 mNumThreads) {
  (void)prng;
  (void)mNumThreads;
  const __uint128_t sideLength = (mDim + 1) * mDelta;
  oc::AES hasher(oc::toBlock(12345));
  for (u64 i = 0; i < mSenderSize; ++i) {
    for (u64 j = 0; j < mDim + 1; ++j) {
      for (u64 k = 0; k < mDim; ++k) {
        const __uint128_t input = *reinterpret_cast<const __uint128_t *>(
            &inputs[i * mDim + k]);
        const __uint128_t cell =
            sideLength * ((input - j * mDim) / sideLength) + j * mDim;
        std::memcpy(&origins[i * (mDim + 1) * mDim + j * mDim + k], &cell,
                    sizeof(block));
        block cellBlock;
        std::memcpy(&cellBlock, &cell, sizeof(block));
        identifiers[i * (mDim + 1) + j] ^= hasher.ecbEncBlock(cellBlock);
      }
    }
  }
  co_await chl.flush();
}

Proto FmapReceiver::fuzzyMap(span<block> inputs, span<block> identifiers,
                             PRNG &prng, Socket &chl, u64 mNumThreads) {
  (void)prng;
  (void)mNumThreads;
  const __uint128_t sideLength = (mDim + 1) * mDelta;
  oc::AES hasher(oc::toBlock(12345));
  for (u64 i = 0; i < mSenderSize; ++i) {
    for (u64 j = 0; j < mDim + 1; ++j) {
      for (u64 k = 0; k < mDim; ++k) {
        const __uint128_t input = *reinterpret_cast<const __uint128_t *>(
            &inputs[i * mDim + k]);
        const __uint128_t cell =
            sideLength * ((input - j * mDim) / sideLength) + j * mDim;
        block cellBlock;
        std::memcpy(&cellBlock, &cell, sizeof(block));
        identifiers[i * (mDim + 1) + j] ^= hasher.ecbEncBlock(cellBlock);
      }
    }
  }
  co_await chl.flush();
}

} // namespace CmpFuzzyPSI
