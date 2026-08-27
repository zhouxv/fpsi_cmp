#include "cryptoTools/Common/CLP.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "ole/ot_ole.h"

#include <iostream>
#include <chrono>
#include <thread>

using namespace CmpFuzzyPSI;

int main(int argc, char **argv) {
  oc::CLP cmd(argc, argv);
  const u64 count = cmd.getOr<u64>("n", 128);
  const u64 modulus = cmd.getOr<u64>("m", 65536);
  auto [senderSocket, receiverSocket] = coproto::LocalAsyncSocket::makePair();
  PRNG senderPrng(oc::toBlock(123));
  PRNG receiverPrng(oc::toBlock(456));
  std::vector<u64> b, d, a, c;

  const auto begin = std::chrono::steady_clock::now();

  std::thread sender([&] {
    coproto::sync_wait(otOleSender(count, modulus, senderPrng, senderSocket,
                                    b, d));
  });
  std::thread receiver([&] {
    coproto::sync_wait(otOleReceiver(count, modulus, receiverPrng,
                                      receiverSocket, a, c));
  });
  sender.join();
  receiver.join();
  const auto end = std::chrono::steady_clock::now();

  bool passed = true;
  for (u64 i = 0; i < count; ++i) {
    const u64 product = static_cast<u64>((static_cast<__uint128_t>(a[i]) * b[i]) % modulus);
    const u64 shares = static_cast<u64>((static_cast<__uint128_t>(c[i]) + d[i]) % modulus);
    if (product != shares) {
      std::cerr << "FAIL at " << i << ": " << a[i] << " * " << b[i]
                << " != " << c[i] << " + " << d[i] << " (mod "
                << modulus << ")\n";
      passed = false;
      break;
    }
  }
  const u64 communication = receiverSocket.bytesSent() + receiverSocket.bytesReceived();
  const double elapsed = std::chrono::duration<double>(end - begin).count();
  std::cout << (passed ? "OT-OLE PASS" : "OT-OLE FAIL")
            << "  n=" << count
            << "  time=" << elapsed << " s"
            << "  communication=" << communication / 1024.0 / 1024.0
            << " MiB\n";
  return passed ? 0 : 1;
}
