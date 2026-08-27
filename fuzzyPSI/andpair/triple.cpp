#include "triple.h"

#include <algorithm>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtSender.h>
#include <stdexcept>

namespace {

// Keep the Silent OT working buffers bounded for the large mIMT batches. This
// is byte-aligned, so the BitVector copies do not need bit shifting.
constexpr u64 kOtChunkSize = 1ull << 22;

// P0 is the sender of one random 1-out-of-2 OT per output bit. If the OT
// outputs are (M0, M1) and P1's natural random choice is a1, define
//   c0 = lsb(M0), b0 = lsb(M0) ^ lsb(M1), c1 = lsb(M_a1).
// Then c0 ^ c1 = b0 & a1. This is exactly the cross correlation consumed by
// the existing mIMT online path; no Beaver conversion or correction exchange
// is needed.
coproto::task<> generateCrossP0(coproto::Socket &chl, BitVector &a0,
                                BitVector &b0, BitVector &c0, bool silent) {
  if (!silent)
    throw std::invalid_argument("mIMT cross triples require Silent OT");

  const u64 count = a0.size();
  PRNG prng(sysRandomSeed());
  a0.randomize(prng); // Not consumed by mIMT's P0 path.
  b0.resize(count);
  c0.resize(count);

  for (u64 offset = 0; offset < count; offset += kOtChunkSize) {
    const u64 chunk = std::min(kOtChunkSize, count - offset);
    SilentOtExtSender sender;
    sender.configure(chunk, 2, 1, SilentSecType::SemiHonest);
    co_await sender.genSilentBaseOts(prng, chl, true);

    std::vector<std::array<block, 2>> messages(chunk);
    co_await sender.silentSend(messages, prng, chl);
    for (u64 i = 0; i < chunk; ++i) {
      c0[offset + i] = block_to_bool(messages[i][0]);
      b0[offset + i] = c0[offset + i] ^ block_to_bool(messages[i][1]);
    }
  }
}

coproto::task<> generateCrossP1(coproto::Socket &chl, BitVector &a1,
                                BitVector &b1, BitVector &c1, bool silent) {
  if (!silent)
    throw std::invalid_argument("mIMT cross triples require Silent OT");

  const u64 count = a1.size();
  PRNG prng(sysRandomSeed());
  b1.resize(count);
  b1.randomize(prng); // Not consumed by mIMT's P1 path.
  c1.resize(count);

  for (u64 offset = 0; offset < count; offset += kOtChunkSize) {
    const u64 chunk = std::min(kOtChunkSize, count - offset);
    SilentOtExtReceiver receiver;
    receiver.configure(chunk, 2, 1, SilentSecType::SemiHonest);
    co_await receiver.genSilentBaseOts(prng, chl, true);

    BitVector choices(chunk);
    std::vector<block> messages(chunk);
    co_await receiver.silentReceive(choices, messages, prng, chl);
    for (u64 i = 0; i < chunk; ++i) {
      a1[offset + i] = choices[i];
      c1[offset + i] = block_to_bool(messages[i]);
    }
  }
}

} // namespace

coproto::task<> triple0(coproto::Socket &chl, BitVector &a0, BitVector &b0,
                        BitVector &c0, bool silent) {
  co_await generateCrossP0(chl, a0, b0, c0, silent);
}

coproto::task<> triple1(coproto::Socket &chl, BitVector &a1, BitVector &b1,
                        BitVector &c1, bool silent) {
  co_await generateCrossP1(chl, a1, b1, c1, silent);
}

// The OT directly produced the cross correlation, so these compatibility
// hooks intentionally have no communication.
coproto::task<> trans_andpair0(coproto::Socket &chl, Triples &triple) {
  (void)chl;
  (void)triple;
  co_return;
}

coproto::task<> trans_andpair1(coproto::Socket &chl, Triples &triple) {
  (void)chl;
  (void)triple;
  co_return;
}
