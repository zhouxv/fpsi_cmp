#include "sspeqt_support.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <deque>

namespace CmpFuzzyPSI {
using namespace oc;
namespace {

constexpr u64 kSoftSpokenFieldBits = 5;

void softSend(u32 count, Socket &chl, PRNG &prng,
              AlignedUnVector<std::array<block, 2>> &messages) {
  SoftSpokenShOtSender<> sender;
  sender.init(kSoftSpokenFieldBits, true);
  PRNG otPrng(prng.get<block>());
  AlignedUnVector<block> baseMessages(sender.baseOtCount());
  BitVector baseChoices(sender.baseOtCount());
  baseChoices.randomize(otPrng);
  DefaultBaseOT base;
  coproto::sync_wait(base.receive(baseChoices, baseMessages, otPrng, chl));
  sender.setBaseOts(baseMessages, baseChoices);
  messages.resize(count);
  coproto::sync_wait(sender.send(messages, otPrng, chl));
}

void softRecv(u32 count, const BitVector &choices, Socket &chl, PRNG &prng,
              AlignedUnVector<block> &messages) {
  SoftSpokenShOtReceiver<> receiver;
  receiver.init(kSoftSpokenFieldBits, true);
  PRNG otPrng(prng.get<block>());
  AlignedUnVector<std::array<block, 2>> baseMessages(receiver.baseOtCount());
  DefaultBaseOT base;
  coproto::sync_wait(base.send(baseMessages, otPrng, chl));
  receiver.setBaseOts(baseMessages);
  messages.resize(count);
  coproto::sync_wait(receiver.receive(choices, messages, otPrng, chl));
}

u32 ceilLog2(u32 value) {
  u32 result = 0;
  while ((1u << result) < value)
    ++result;
  return result;
}

u8 blockBit(const block &value, u32 bit) {
  return static_cast<u8>((value.get<u64>(bit / 64) >> (bit % 64)) & 1);
}

void expandSeed(block seed, block &left, block &right) {
  AES aes(seed);
  aes.ecbEncBlock(AllOneBlock, left);
  aes.ecbEncBlock(AllOneBlock ^ OneBlock, right);
}

void prepareCorrelation(u32 depth, u32 level, block seed, block *zero,
                        block *one) {
  if (level == depth)
    return;
  block left, right;
  expandSeed(seed, left, right);
  zero[level] ^= left;
  one[level] ^= right;
  prepareCorrelation(depth, level + 1, left, zero, one);
  prepareCorrelation(depth, level + 1, right, zero, one);
}

void expandTree(u32 depth, u32 level, block seed, block *leaves, u32 &index) {
  if (level == depth) {
    leaves[index++] = seed;
    return;
  }
  block left, right;
  expandSeed(seed, left, right);
  expandTree(depth, level + 1, left, leaves, index);
  expandTree(depth, level + 1, right, leaves, index);
}

void punctureExpand(u32 depth, u32 puncture, const block *seeds,
                    block *leaves) {
  std::vector<u8> bits(depth);
  for (i32 bit = static_cast<i32>(depth) - 1; bit >= 0; --bit) {
    bits[bit] = (puncture & 1) ^ 1;
    puncture >>= 1;
  }

  std::deque<block> queue;
  if (bits[0]) {
    queue.push_back(ZeroBlock);
    queue.push_back(seeds[0]);
  } else {
    queue.push_back(seeds[0]);
    queue.push_back(ZeroBlock);
  }

  u32 path = bits[0];
  for (u32 level = 1; level < depth; ++level) {
    const auto size = queue.size();
    block aggregate = ZeroBlock;
    for (size_t idx = 0; idx < size; ++idx) {
      block seed = queue.front();
      queue.pop_front();
      if (seed == ZeroBlock) {
        if (bits[level]) {
          queue.push_back(ZeroBlock);
          queue.push_back(seeds[level]);
        } else {
          queue.push_back(seeds[level]);
          queue.push_back(ZeroBlock);
        }
      } else {
        block left, right;
        expandSeed(seed, left, right);
        queue.push_back(left);
        queue.push_back(right);
        aggregate ^= bits[level] ? right : left;
      }
    }
    path = ((path ^ 1) << 1) ^ bits[level];
    queue[path] ^= aggregate;
  }
  u32 index = 0;
  for (const auto &leaf : queue)
    leaves[index++] = leaf;
}

void l1lRotSend(u32 domainSize, u32 count, PRNG &prng, Socket &chl,
                std::vector<std::vector<block>> &columns) {
  const u32 depth = ceilLog2(domainSize);
  const u32 leafCount = 1u << depth;
  const u32 otCount = std::max(depth * count, 128u);
  std::vector<block> roots(count);
  std::vector<std::vector<block>> zero(count,
                                       std::vector<block>(depth, ZeroBlock));
  std::vector<std::vector<block>> one(count,
                                      std::vector<block>(depth, ZeroBlock));
  for (u32 batch = 0; batch < count; ++batch) {
    roots[batch] = prng.get<block>();
    prepareCorrelation(depth, 0, roots[batch], zero[batch].data(),
                       one[batch].data());
  }

  AlignedUnVector<std::array<block, 2>> otMessages;
  softSend(otCount, chl, prng, otMessages);
  std::vector<block> pads(2 * depth * count);
  for (u32 batch = 0; batch < count; ++batch) {
    for (u32 level = 0; level < depth; ++level) {
      const u32 padIndex = 2 * (batch * depth + level);
      const u32 otIndex = batch * depth + level;
      pads[padIndex] = otMessages[otIndex][0] ^ zero[batch][level];
      pads[padIndex + 1] = otMessages[otIndex][1] ^ one[batch][level];
    }
  }
  coproto::sync_wait(chl.send(coproto::span<block>(pads)));
  coproto::sync_wait(chl.flush());

  columns.resize(count);
  for (u32 batch = 0; batch < count; ++batch) {
    columns[batch].resize(leafCount);
    u32 index = 0;
    expandTree(depth, 0, roots[batch], columns[batch].data(), index);
  }
}

std::vector<u32> l1lRotRecv(u32 domainSize, u32 count, PRNG &prng,
                            Socket &chl,
                            std::vector<std::vector<block>> &columns) {
  const u32 depth = ceilLog2(domainSize);
  const u32 leafCount = 1u << depth;
  const u32 otCount = std::max(depth * count, 128u);
  std::vector<u32> punctures(count);
  BitVector choices(otCount);
  for (u32 batch = 0; batch < count; ++batch) {
    punctures[batch] = prng.get<u32>() % domainSize;
    u32 value = punctures[batch];
    for (i32 bit = static_cast<i32>(depth) - 1; bit >= 0; --bit) {
      choices[batch * depth + bit] = (value & 1) ^ 1;
      value >>= 1;
    }
  }

  AlignedUnVector<block> otMessages;
  softRecv(otCount, choices, chl, prng, otMessages);
  std::vector<block> pads(2 * depth * count);
  coproto::sync_wait(chl.recv(coproto::span<block>(pads)));
  coproto::sync_wait(chl.flush());

  columns.resize(count);
  for (u32 batch = 0; batch < count; ++batch) {
    std::vector<block> seeds(depth);
    for (u32 level = 0; level < depth; ++level) {
      const u32 padIndex = 2 * (batch * depth + level);
      seeds[level] = pads[padIndex + static_cast<u32>(choices[batch * depth + level])] ^
                     otMessages[batch * depth + level];
    }
    columns[batch].resize(leafCount);
    punctureExpand(depth, punctures[batch], seeds.data(), columns[batch].data());
  }
  return punctures;
}

} // namespace

void B2AOfflineSender(u32 count, u32 modulus, Socket &chl, PRNG &prng,
                      BitVector &rShare, std::vector<u64> &tShare) {
  SoftSpokenShOtSender<> sender;
  sender.init(kSoftSpokenFieldBits, true);
  PRNG otPrng(prng.get<block>());
  AlignedUnVector<block> baseMessages(sender.baseOtCount());
  BitVector baseChoices(sender.baseOtCount());
  baseChoices.randomize(otPrng);
  DefaultBaseOT base;
  coproto::sync_wait(base.receive(baseChoices, baseMessages, otPrng, chl));
  sender.setBaseOts(baseMessages, baseChoices);

  AlignedUnVector<std::array<block, 2>> otMessages(count);
  coproto::sync_wait(sender.send(otMessages, otPrng, chl));
  rShare.resize(count);
  tShare.resize(count);
  std::vector<block> pads(2 * count);
  for (u32 i = 0; i < count; ++i) {
    rShare[i] = prng.getBit();
    tShare[i] = prng.get<u64>() % modulus;
    const u64 neg = (modulus - tShare[i]) % modulus;
    const u64 oneMinus = (1 + modulus - tShare[i]) % modulus;
    pads[2 * i] = otMessages[i][0] ^ toBlock(rShare[i] ? oneMinus : neg);
    pads[2 * i + 1] = otMessages[i][1] ^ toBlock(rShare[i] ? neg : oneMinus);
  }
  coproto::sync_wait(chl.send(coproto::span<block>(pads)));
  coproto::sync_wait(chl.flush());
}

void B2AOfflineReceiver(u32 count, u32 modulus, Socket &chl, PRNG &prng,
                        BitVector &rShare, std::vector<u64> &tShare) {
  rShare.resize(count);
  for (u32 i = 0; i < count; ++i)
    rShare[i] = prng.getBit();
  SoftSpokenShOtReceiver<> receiver;
  receiver.init(kSoftSpokenFieldBits, true);
  PRNG otPrng(prng.get<block>());
  AlignedUnVector<std::array<block, 2>> baseMessages(receiver.baseOtCount());
  DefaultBaseOT base;
  coproto::sync_wait(base.send(baseMessages, otPrng, chl));
  receiver.setBaseOts(baseMessages);

  AlignedUnVector<block> otMessages(count);
  coproto::sync_wait(receiver.receive(rShare, otMessages, otPrng, chl));
  std::vector<block> pads(2 * count);
  coproto::sync_wait(chl.recv(coproto::span<block>(pads)));
  tShare.resize(count);
  for (u32 i = 0; i < count; ++i)
    tShare[i] = (pads[2 * i + (rShare[i] ? 1 : 0)] ^ otMessages[i]).get<u64>(0) % modulus;
  coproto::sync_wait(chl.flush());
}

void B2AOnline(u32 count, u32 modulus, bool isSender,
               const BitVector &inputShare, const BitVector &rShare,
               const std::vector<u64> &tShare, Socket &chl,
               std::vector<u64> &outputShare) {
  BitVector local(count), remote(count);
  for (u32 i = 0; i < count; ++i)
    local[i] = inputShare[i] ^ rShare[i];
  if (isSender) {
    coproto::sync_wait(chl.send(local));
    coproto::sync_wait(chl.recv(remote));
  } else {
    coproto::sync_wait(chl.recv(remote));
    coproto::sync_wait(chl.send(local));
  }
  coproto::sync_wait(chl.flush());
  outputShare.resize(count);
  for (u32 i = 0; i < count; ++i) {
    const bool value = local[i] ^ remote[i];
    const u64 t = tShare[i] % modulus;
    outputShare[i] = isSender ? (value ? (1 + modulus - t) % modulus : t)
                              : (value ? (modulus - t) % modulus : t);
  }
}

void VoseOfflineSender(u32 domainSize, u32 count, PRNG &prng, Socket &chl,
                       std::vector<std::vector<u8>> &U,
                       std::vector<std::vector<u8>> &V) {
  std::vector<std::vector<block>> rawColumns;
  l1lRotSend(domainSize, count, prng, chl, rawColumns);
  U.assign(count, std::vector<u8>(domainSize));
  V.assign(count, std::vector<u8>(domainSize));
  for (u32 batch = 0; batch < count; ++batch) {
    for (u32 i = 0; i < domainSize; ++i) {
      for (u32 j = 0; j < domainSize; ++j) {
        U[batch][i] ^= blockBit(rawColumns[batch][(i + domainSize - j) % domainSize], j);
        V[batch][i] ^= blockBit(rawColumns[batch][j], i);
      }
    }
  }
  coproto::sync_wait(chl.flush());
}

void VoseOfflineReceiver(u32 domainSize, u32 count, PRNG &prng, Socket &chl,
                         std::vector<u32> &eps,
                         std::vector<std::vector<u8>> &W) {
  std::vector<std::vector<block>> rawColumns;
  eps = l1lRotRecv(domainSize, count, prng, chl, rawColumns);
  W.assign(count, std::vector<u8>(domainSize));
  for (u32 batch = 0; batch < count; ++batch) {
    eps[batch] = (domainSize - eps[batch]) % domainSize;
    std::vector<block> columns(domainSize, ZeroBlock);
    for (u32 j = 0; j < domainSize; ++j)
      columns[j] = rawColumns[batch][j];
    for (u32 i = 0; i < domainSize; ++i) {
      u8 value = 0;
      for (u32 j = 0; j < domainSize; ++j)
        value ^= blockBit(columns[j], i);
      for (u32 j = 0; j < domainSize; ++j) {
        if (j != i)
          value ^= blockBit(columns[(i + domainSize - ((eps[batch] + j) % domainSize)) % domainSize], j);
      }
      W[batch][i] = value;
    }
  }
  coproto::sync_wait(chl.flush());
}

void VoseOnlineSender(u32 domainSize, u32 count,
                      const std::vector<std::vector<u8>> &U,
                      const std::vector<std::vector<u8>> &V,
                      const std::vector<u32> &voseValue, Socket &chl,
                      BitVector &tableShare) {
  BitVector correction(count * domainSize);
  for (u32 batch = 0; batch < count; ++batch)
    for (u32 j = 0; j < domainSize; ++j)
      correction[batch * domainSize + j] = U[batch][j] ^
                                           (j == voseValue[batch] % domainSize);
  coproto::sync_wait(chl.send(std::move(correction)));
  coproto::sync_wait(chl.flush());
  tableShare.resize(count * domainSize);
  for (u32 batch = 0; batch < count; ++batch)
    for (u32 i = 0; i < domainSize; ++i)
      tableShare[batch * domainSize + i] = V[batch][i];
}

void VoseOnlineReceiver(u32 domainSize, u32 count,
                        const std::vector<std::vector<u8>> &W,
                        const std::vector<u32> &eps, Socket &chl,
                        BitVector &tableShare) {
  BitVector correction(count * domainSize);
  coproto::sync_wait(chl.recv(correction));
  coproto::sync_wait(chl.flush());
  tableShare.resize(count * domainSize);
  for (u32 batch = 0; batch < count; ++batch)
    for (u32 i = 0; i < domainSize; ++i)
      tableShare[batch * domainSize + i] =
          W[batch][i] ^ correction[batch * domainSize +
                                    ((i + domainSize - eps[batch] % domainSize) % domainSize)];
}

} // namespace CmpFuzzyPSI
