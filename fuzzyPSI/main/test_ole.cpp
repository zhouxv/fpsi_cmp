#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Common/Defines.h>
#include <libOTe/TwoChooseOne/TcoOtDefines.h>
#include <libOTe/Vole/Silent/SilentVoleReceiver.h>
#include <libOTe/Vole/Silent/SilentVoleSender.h>

#include "ole/ole.h"

using namespace oc;

typedef std::chrono::high_resolution_clock::time_point tVar;
#define tNow() std::chrono::high_resolution_clock::now()
#define tStart(t) t = tNow()
#define tEnd(t)                                                                \
  std::chrono::duration_cast<std::chrono::milliseconds>(tNow() - t).count()

int main(int argc, char **argv) {
  CLP cmd(argc, argv);

  auto numVole = cmd.getOr<u64>("n", 1 << 20);
  auto mod_val = cmd.getOr<u32>("m", 65536);

  // get up the networking
  auto chl = cp::LocalAsyncSocket::makePair();
  // get a random number generator seeded from the system
  PRNG prng(sysRandomSeed());

  // set up the VOLE sender and receiver
  CoeffCtxIntegerMod<u32> mod_ops(mod_val);
  SilentVoleSenderMod<u32, u32> sender(mod_ops);
  SilentVoleReceiverMod<u32, u32> receiver(mod_ops);

  tVar t;
  tStart(t);

  // configure them
  sender.configure(numVole, SilentBaseType::BaseExtend, 128, mod_ops);
  receiver.configure(numVole, SilentBaseType::BaseExtend, 128, mod_ops);

  // 𝔽   𝔽   𝔾    𝔽
  // A = B + C * DELTA
  // ℛ  𝒮   ℛ   𝒮
  AlignedUnVector<u32> A(numVole);
  AlignedUnVector<u32> B(numVole);
  AlignedUnVector<u32> C(numVole);
  u32 DELTA = prng.get<u32>() % mod_val;

  // execute the VOLE
  auto t_s = sender.silentSend(DELTA, B, prng, chl[0]);
  auto t_rs = receiver.silentReceive(C, A, prng, chl[1]);

  cp::sync_wait(cp::when_all_ready(std::move(t_s), std::move(t_rs)));

  // check the results
  for (u64 i = 0; i < 10; i++) {
    u32 minus, mul, sum;
    mod_ops.mul(mul, C[i], DELTA);
    mod_ops.plus(sum, B[i], mul);
    mod_ops.minus(minus, A[i], B[i]);

    std::cout
        << std::format(
               "VOLE[{}]: A={:12} B={:12} C={:12} DELTA={:12}, A-B={:12}, "
               "C*DELTA={:12}, B+C*DELTA:{:12}",
               i, A[i], B[i], C[i], DELTA, minus, mul, sum)
        << std::endl;
  }

  auto duration = tEnd(t);
  std::cout << std::format(
                   "✓ Silent VOLE test passed! ({} elements), elapsed {} s",
                   numVole, duration / 1000.0)
            << std::endl;

  return 0;
}
