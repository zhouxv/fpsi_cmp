#include <iostream>
#include <vector>
#include <iomanip>
#include "debug.h"
#include "sparsehash/dense_hash_map"
#include "cryptoTools/Common/Timer.h" 
#include "cryptoTools/Common/CuckooIndex.h"
#include "volePSI/SimpleIndex.h"

using namespace oc;
using namespace volePSI;

int main() {
    PRNG prng(_mm_set_epi32(42, 43, 44, 45));

    // Parameters (mimic RsCpsi)
    u64 senderSize = 1<<12;
    u64 recverSize = senderSize;
    u64 ssp = 40;

    // Generate overlapping inputs: X (Receiver), Y (Sender)
    std::vector<block> X(recverSize), Y(senderSize);
    prng.get(X.data(), recverSize);
    // Let first 10 of Y = first 10 of X (common elements)
    for (u64 i = 0; i < 10; ++i) Y[i] = X[i];
    // Rest random
    prng.get(Y.data() + 10, senderSize - 10);

    // ==== STEP 1: Receiver builds Cuckoo Table (final placement) ====
    block cuckooSeed = prng.get();  // Receiver chooses seed

    CuckooIndex<> cuckoo;
    cuckoo.init(recverSize, ssp, 0, 3);
    cuckoo.insert(X, cuckooSeed);

    u64 numBinsRecv = cuckoo.mBins.size();

    // CuckooTable[b] = bin index where X[b] finally resides
    std::vector<u64> CuckooTable(recverSize, u64(-1));
    for (u64 i = 0; i < numBinsRecv; ++i) {
        if (!cuckoo.mBins[i].isEmpty()) {
            u64 b = cuckoo.mBins[i].idx();
            CuckooTable[b] = i;
        }
    }

    // Sanity: all placed
    for (u64 b = 0; b < recverSize; ++b) {
        if (CuckooTable[b] == u64(-1)) {
            std::cerr << "❌ Cuckoo insertion failed for X[" << b << "]\n";
            return 1;
        }
    }

    auto params = CuckooIndex<>::selectParams(recverSize, ssp, 0, 3);
    u64 TableSize = params.numBins();
    std::cout << "NumBins = " << TableSize << "\n";

    SimpleIndex sIdx;
    sIdx.init(TableSize, senderSize, ssp, 3);
    sIdx.insertItems(Y, cuckooSeed);

    std::vector<u64> SimpleTable(Y.size() * 3, u64(-1));

    for (u64 i = 0; i < TableSize; ++i) {
        const auto& bin = sIdx.mBins[i];
        u64 size = sIdx.mBinSizes[i];
        for (u64 p = 0; p < size; ++p) {
            u64 j = bin[p].hashIdx();
            u64 b = bin[p].idx();
            assert(b < Y.size());
            assert(j < 3);
            SimpleTable[b * 3 + j] = i;
        }
    }

    std::cout << "\n=== Main Check: For each common X[b], is CuckooTable[b] in {H0,H1,H2} bins? ===\n";
    bool all_ok = true;

    // Build map: X[b] → b (for fast lookup in Y)
    std::unordered_map<block, u64> Y_to_idx;
    for (u64 b = 0; b < senderSize; ++b) {
        Y_to_idx[Y[b]] = b;
    }

    for (u64 b = 0; b < recverSize; ++b) {
        auto it = Y_to_idx.find(X[b]);
        if (it != Y_to_idx.end()) {
            u64 y_idx = it->second;
            u64 final_bin = CuckooTable[b];
            std::unordered_set<u64> candidates = {
                SimpleTable[y_idx * 3 + 0],
                SimpleTable[y_idx * 3 + 1],
                SimpleTable[y_idx * 3 + 2]
            };

            if (candidates.count(final_bin) == 0) {
                std::cout << "❌ X[" << b << "] (Y[" << y_idx << "]) placed in bin " << final_bin
                          << ", but candidates are {" << *candidates.begin();
                auto iter = candidates.begin(); ++iter;
                for (; iter != candidates.end(); ++iter) std::cout << ", " << *iter;
                std::cout << "}\n";
                all_ok = false;
            }
        }
    }
    if (all_ok) {
        std::cout << "✅ All common elements: final bin ∈ {3 candidate bins}.\n";
    } else {
        std::cout << "❌ Some common elements violate cuckoo candidate rule.\n";
        return 1;
}

    // === Validation & Output ===
    // std::cout << "\nBin occupancy:\n";
    // for (u64 i = 0; i < TableSize && i < 20; ++i) { // print first 20 bins
    //     if (sIdx.mBinSizes[i] > 0) {
    //         std::cout << "Bin[" << std::setw(3) << i << "]: ";
    //         for (u64 p = 0; p < sIdx.mBinSizes[i]; ++p) {
    //             u64 b = sIdx.mBins[i][p].idx();
    //             u64 j = sIdx.mBins[i][p].hashIdx();
    //             std::cout << "Y[" << b << "]@H" << j << " ";
    //         }
    //         std::cout << "\n";
    //     }
    // }

    // std::cout << "\nSimpleTable (b*3 + j -> bin):\n";
    // for (u64 b = 0; b < Identifiers.size(); ++b) {
    //     std::cout << "Y[" << b << "]: ";
    //     for (u64 j = 0; j < 3; ++j) {
    //         u64 binIdx = SimpleTable[b * 3 + j];
    //         if (binIdx == u64(-1)) {
    //             // This can happen if (b,j) pair never landed in any bin? 
    //             // Actually, SimpleIndex inserts ALL 3*mSenderSize entries, so no -1 expected.
    //             std::cout << "??? ";
    //         } else {
    //             std::cout << "H" << j << "->" << binIdx << "  ";
    //         }
    //     }
    //     std::cout << "\n";
    // }

    // Cross-check: recompute bin index via hashers and mod, compare with SimpleTable
    // std::cout << "\nCross-check (recompute bin = H_j(Y[b]) mod TableSize):\n";
    // bool all_match = true;
    // for (u64 b = 0; b < Identifiers.size(); ++b) {
    //     for (u64 j = 0; j < 3; ++j) {
    //         block h = hashers[j].hashBlock(Identifiers[b]);
    //         u64 expected_bin = ((u64*)&h)[0] % TableSize; // low 64-bit mod TableSize
    //         u64 recorded_bin = SimpleTable[b * 3 + j];
    //         if (expected_bin != recorded_bin) {
    //             std::cout << "❌ mismatch: Y[" << b << "], H" << j 
    //                       << ": expected " << expected_bin << ", got " << recorded_bin << "\n";
    //             all_match = false;
    //         }
    //     }
    // }
    // if (all_match) {
    //     std::cout << "✅ All (b,j) -> bin mappings match recomputation!\n";
    // } else {
    //     std::cout << "❌ Some mismatches found.\n";
    //     return 1;
    // }

    return 0;
}