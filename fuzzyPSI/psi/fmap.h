#pragma once

#include "psi/Defines.h"
#include "cryptoTools/Common/Timer.h"

using namespace oc;

namespace CmpFuzzyPSI {

class FmapSender : public oc::TimerAdapter {
public:
  u64 mSenderSize;
  u64 mRecverSize;
  u64 mDim;
  u64 mDelta;
  u64 orgSize;
  u64 myExpansionRate;
  u64 anotherExpansionRate;

  Proto setUp(u64 senderSize, u64 receiverSize, u64 dim, u64 delta,
              PRNG &prng, Socket &chl, u64 mNumThreads = 1);
  Proto fuzzyMap(span<block> inputs, span<block> identifiers,
                 span<block> origins, PRNG &prng, Socket &chl,
                 u64 mNumThreads = 1);
};

class FmapReceiver : public oc::TimerAdapter {
public:
  u64 mSenderSize;
  u64 mRecverSize;
  u64 mDim;
  u64 mDelta;
  u64 orgSize;
  u64 myExpansionRate;
  u64 anotherExpansionRate;

  Proto setUp(u64 senderSize, u64 receiverSize, u64 dim, u64 delta,
              PRNG &prng, Socket &chl, u64 mNumThreads = 1);
  Proto fuzzyMap(span<block> inputs, span<block> identifiers, PRNG &prng,
                 Socket &chl, u64 mNumThreads = 1);
};

} // namespace CmpFuzzyPSI
