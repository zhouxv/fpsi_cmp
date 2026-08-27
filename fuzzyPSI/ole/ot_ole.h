#pragma once

#include "psi/Defines.h"

namespace CmpFuzzyPSI {

// Batched semi-honest OLE over Z_modulus, built directly from IKNP 1-out-of-2
// OT extension. The sender receives (b, d), the receiver receives (a, c), and
// the outputs satisfy a * b = c + d (mod modulus).
Proto otOleSender(u64 count, u64 modulus, PRNG &prng, Socket &chl,
                  std::vector<u64> &b, std::vector<u64> &d);
Proto otOleReceiver(u64 count, u64 modulus, PRNG &prng, Socket &chl,
                    std::vector<u64> &a, std::vector<u64> &c);

} // namespace CmpFuzzyPSI
