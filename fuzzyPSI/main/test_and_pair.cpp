#include "andpair/triple.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "cryptoTools/Common/CLP.h"

#include <iostream>
#include <thread>

int main(int argc, char **argv) {
  CLP cmd(argc, argv);
  const u64 count = cmd.getOr<u64>("n", 128);
  if (count == 0) {
    std::cerr << "-n must be positive\n";
    return 1;
  }

  // mIMT's existing setup calls these functions from the two parties.
  auto [socket0, socket1] = coproto::LocalAsyncSocket::makePair();
  Triples cross0(count, true), cross1(count, true);
  std::thread crossParty0([&] {
    macoro::sync_wait(cross0.gen0(socket0));
    macoro::sync_wait(trans_andpair0(socket0, cross0));
  });
  std::thread crossParty1([&] {
    macoro::sync_wait(cross1.gen1(socket1));
    macoro::sync_wait(trans_andpair1(socket1, cross1));
  });
  crossParty0.join();
  crossParty1.join();
  bool adapterPass = true;
  for (u64 i = 0; i < count; ++i) {
    if ((cross0.c[i] ^ cross1.c[i]) != (cross0.b[i] & cross1.a[i])) {
      adapterPass = false;
      break;
    }
  }

  std::cout << "mIMT one-OT cross triple: "
            << (adapterPass ? "PASS" : "FAIL") << "\n";
  return adapterPass ? 0 : 1;
}
