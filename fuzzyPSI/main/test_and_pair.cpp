#include "andpair/triple.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "cryptoTools/Common/CLP.h"

typedef std::chrono::high_resolution_clock::time_point tVar;
#define tNow() std::chrono::high_resolution_clock::now()
#define tStart(t) t = tNow()
#define tEnd(t)                                                                \
  std::chrono::duration_cast<std::chrono::milliseconds>(tNow() - t).count()

int main(int argc, char **argv) {
  CLP cmd(argc, argv);

  auto silent = cmd.isSet("silent");
  auto num_triples = cmd.getOr("n", 128);

  auto [socket0, socket1] = coproto::LocalAsyncSocket::makePair();

  Triples triples_0(num_triples, silent);
  Triples triples_1(num_triples, silent);

  tVar t;
  tStart(t);
  std::thread recv_thread([&]() {
    macoro::sync_wait(triples_0.gen0(socket0));
    macoro::sync_wait(trans_andpair0(socket0, triples_0));
  });

  std::thread send_thread([&]() {
    macoro::sync_wait(triples_1.gen1(socket1));
    macoro::sync_wait(trans_andpair1(socket1, triples_1));
  });

  // wait for both sides to finish
  recv_thread.join();
  send_thread.join();

  auto duration = tEnd(t);

  for (u64 i = 0; i < 128; i++) {
    auto e_s = triples_0.b[i];
    auto and_s = triples_0.c[i];

    auto e_r = triples_1.a[i];
    auto and_r = triples_1.c[i];

    std::cout << "Triple " << i << ": " << "e_s = " << e_s << ", e_r = " << e_r
              << ", A = " << (e_s & e_r) << ", and_s= " << and_s
              << ", and_r = " << and_r << ", AND = " << (and_s ^ and_r) << " "
              << ((e_s & e_r) == (and_s ^ and_r)) << std::endl;
  }

  std::cout << "Generated " << num_triples << " triples in " << duration
            << " ms." << std::endl;

  return 0;
}