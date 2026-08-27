#include "ot_ole.h"

#include <algorithm>
#include <libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h>

namespace CmpFuzzyPSI {
namespace {

u64 bitLength(u64 modulus) {
  if (modulus < 2)
    throw std::invalid_argument("OLE modulus must be at least two");
  return oc::log2ceil(modulus);
}

u64 addMod(u64 lhs, u64 rhs, u64 modulus) {
  return static_cast<u64>((static_cast<__uint128_t>(lhs) + rhs) % modulus);
}

u64 mulMod(u64 lhs, u64 rhs, u64 modulus) {
  return static_cast<u64>((static_cast<__uint128_t>(lhs) * rhs) % modulus);
}

constexpr u64 kMaxOtsPerBatch = 1ull << 20;

Proto otOleSenderRange(u64 begin, u64 end, u64 modulus, PRNG &prng,
                       Socket &chl, std::vector<u64> &b,
                       std::vector<u64> &d) {
  const u64 bits = bitLength(modulus);
  const u64 maxElements = std::max<u64>(1, kMaxOtsPerBatch / bits);
  osuCrypto::IknpOtExtSender sender;

  for (u64 batchBegin = begin; batchBegin < end;
       batchBegin += maxElements) {
    const u64 elements = std::min(maxElements, end - batchBegin);
    const u64 otCount = elements * bits;
    std::vector<std::array<block, 2>> pads(otCount);
    macoro::sync_wait(sender.send(pads, prng, chl));

    // Only the low 64 bits of each random-OT pad are needed to mask a u64
    // ring element. This retains more than the protocol's 40-bit security
    // target while halving the correction payload from 32 to 16 bytes/OT.
    std::vector<u64> corrections(2 * otCount);
    for (u64 i = 0; i < elements; ++i) {
      const u64 index = batchBegin + i;
      b[index] = prng.get<u64>() % modulus;
      u64 maskSum = 0;
      for (u64 bit = 0; bit < bits; ++bit) {
        const u64 otIndex = i * bits + bit;
        const u64 mask = prng.get<u64>() % modulus;
        const u64 product = mulMod(b[index], 1ull << bit, modulus);
        maskSum = addMod(maskSum, mask, modulus);
        corrections[2 * otIndex] = pads[otIndex][0].get<u64>(0) ^ mask;
        corrections[2 * otIndex + 1] =
            pads[otIndex][1].get<u64>(0) ^ addMod(mask, product, modulus);
      }
      d[index] = maskSum == 0 ? 0 : modulus - maskSum;
    }
    macoro::sync_wait(chl.send(std::move(corrections)));
    macoro::sync_wait(chl.flush());
  }
  co_return;
}

Proto otOleReceiverRange(u64 begin, u64 end, u64 modulus, PRNG &prng,
                         Socket &chl, std::vector<u64> &a,
                         std::vector<u64> &c) {
  const u64 bits = bitLength(modulus);
  const u64 maxElements = std::max<u64>(1, kMaxOtsPerBatch / bits);
  osuCrypto::IknpOtExtReceiver receiver;

  for (u64 batchBegin = begin; batchBegin < end;
       batchBegin += maxElements) {
    const u64 elements = std::min(maxElements, end - batchBegin);
    const u64 otCount = elements * bits;
    BitVector choices(otCount);
    for (u64 i = 0; i < elements; ++i) {
      const u64 index = batchBegin + i;
      a[index] = prng.get<u64>() % modulus;
      for (u64 bit = 0; bit < bits; ++bit)
        choices[i * bits + bit] = (a[index] >> bit) & 1;
    }

    std::vector<block> pads(otCount);
    macoro::sync_wait(receiver.receive(choices, pads, prng, chl));
    std::vector<u64> corrections(2 * otCount);
    macoro::sync_wait(chl.recv(corrections));
    macoro::sync_wait(chl.flush());

    for (u64 i = 0; i < elements; ++i) {
      u64 value = 0;
      for (u64 bit = 0; bit < bits; ++bit) {
        const u64 otIndex = i * bits + bit;
        const u64 selected = choices[otIndex] ? 1 : 0;
        const u64 term =
            (pads[otIndex].get<u64>(0) ^ corrections[2 * otIndex + selected]) %
            modulus;
        value = addMod(value, term, modulus);
      }
      c[batchBegin + i] = value;
    }
  }
  co_return;
}

} // namespace

Proto otOleSender(u64 count, u64 modulus, PRNG &prng, Socket &chl,
                  std::vector<u64> &b, std::vector<u64> &d) {
  b.resize(count);
  d.assign(count, 0);
  if (count == 0)
    co_return;
  macoro::sync_wait(otOleSenderRange(0, count, modulus, prng, chl, b, d));
  co_return;
}

Proto otOleReceiver(u64 count, u64 modulus, PRNG &prng, Socket &chl,
                    std::vector<u64> &a, std::vector<u64> &c) {
  a.resize(count);
  c.assign(count, 0);
  if (count == 0)
    co_return;
  macoro::sync_wait(otOleReceiverRange(0, count, modulus, prng, chl, a, c));
  co_return;
}

} // namespace CmpFuzzyPSI
