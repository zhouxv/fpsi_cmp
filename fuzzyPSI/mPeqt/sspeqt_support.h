#pragma once

#include "psi/Defines.h"
#include <libOTe/Base/BaseOT.h>
#include <libOTe/TwoChooseOne/SoftSpokenOT/SoftSpokenShOtExt.h>

namespace CmpFuzzyPSI {

// The helpers below implement the offline preprocessing used by ssPEQT:
// Boolean-to-arithmetic conversion (B2A) and the two-party VOSE table.
void B2AOfflineSender(u32 count, u32 modulus, Socket &chl, PRNG &prng,
                      BitVector &rShare, std::vector<u64> &tShare);
void B2AOfflineReceiver(u32 count, u32 modulus, Socket &chl, PRNG &prng,
                        BitVector &rShare, std::vector<u64> &tShare);
void B2AOnline(u32 count, u32 modulus, bool isSender,
               const BitVector &inputShare, const BitVector &rShare,
               const std::vector<u64> &tShare, Socket &chl,
               std::vector<u64> &outputShare);

void VoseOfflineSender(u32 domainSize, u32 count, PRNG &prng, Socket &chl,
                       std::vector<std::vector<u8>> &U,
                       std::vector<std::vector<u8>> &V);
void VoseOfflineReceiver(u32 domainSize, u32 count, PRNG &prng, Socket &chl,
                         std::vector<u32> &eps,
                         std::vector<std::vector<u8>> &W);
void VoseOnlineSender(u32 domainSize, u32 count,
                      const std::vector<std::vector<u8>> &U,
                      const std::vector<std::vector<u8>> &V,
                      const std::vector<u32> &voseValue, Socket &chl,
                      BitVector &tableShare);
void VoseOnlineReceiver(u32 domainSize, u32 count,
                        const std::vector<std::vector<u8>> &W,
                        const std::vector<u32> &eps, Socket &chl,
                        BitVector &tableShare);

} // namespace CmpFuzzyPSI
