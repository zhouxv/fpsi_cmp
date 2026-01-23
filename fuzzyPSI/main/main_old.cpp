
#include "coproto/Socket/AsioSocket.h"
#include "coproto/Socket/LocalAsyncSock.h" // for localAsyncSocketPair
#include "cryptoTools/Common/CLP.h"
#include "cryptoTools/Common/Timer.h"

#include "debug.h"
#include "psi/psi.h"
#include "volePSI/RsOpprf.h"
#include <iomanip>
#include <iostream>
#include <vector>

using namespace volePSI;
using namespace oc;

void printUsage(const char *prog) {
  std::cout << "Usage: " << prog << "[options]\n"
            << "  Options:\n"
            << "    -n <N>          : input size (default: 1024)\n"
            << "    -nn <N>          : input size 2^n (default: 12)\n"
            << "    -dim <N>        : diemnsion (default: 5)\n"
            << "    -delta <N>      : distamce threshold (default: 30)\n"
            << "    -metric <N>     : which p for L_p distance, 0 is for "
               "infinate case (default: 0)\n"
            //   << "    -t <T>          : Thread count (default: 2)\n"
            << "    -a <addr>       : Address for network mode, e.g., "
               "'localhost:1212'\n"
            << "    -LorH           : Low diemsnional protocol (1) or High "
               "dimensional protocol (0)\n"
            << "\n"
            << "Examples:\n"
            << "  " << prog << " -m local\n"
            << "  " << prog << " -m sender -a localhost:1212 &\n"
            << "  " << prog << " -m receiver -a localhost:1212\n";
}

int main(int argc, char **argv) {
  CLP cmd;
  cmd.parse(argc, argv);

  // print help message
  if (cmd.isSet("h") || cmd.isSet("help") || argc < 2) {
    printUsage(argv[0]);
    return 0;
  }

  // obtain parameters
  u64 n = cmd.getOr<u64>("n", 1 << 12); // default 4096
  if (cmd.isSet("nn")) {
    n = 1ULL << cmd.get<u64>("nn"); // if -nn is set, n = 2^nn
  }

  const u64 dim = cmd.getOr<u64>("dim", 6);
  const u64 delta = cmd.getOr<u64>("delta", 60);
  const u64 metric = cmd.getOr<u64>("metric", 0);
  const u64 numThreads = cmd.getOr<u64>("t", 1);
  const u64 LorH = cmd.getOr<u64>("LorH", 0);
  const std::string ip = cmd.getOr<std::string>("ip", "localhost");
  const u64 port = cmd.getOr<u64>("port", 1212);
  const u64 trait = cmd.getOr<u64>("trait", 5);
  std::string addr = ip + ":" + std::to_string(port);

  const int role = cmd.getOr<int>("r", 0); // default receiver (0)

  if (role == 1) {
    for (u64 i = 0; i < trait; i++) {
      // sender side
      CmpFuzzyPSI::FuzzyPsiSender sender;
      block seed = oc::toBlock(123);
      sender.init(n, n, 40, dim, metric, delta, LorH, seed, numThreads, false);

      // connect to receiver
      coproto::Socket chl = coproto::asioConnect(addr, true);

      // generate input data
      std::vector<block> inputs(n * dim);
      PRNG prng(oc::toBlock(456));
      for (u64 i = 0; i < n; ++i) {
        for (u64 j = 0; j < dim; ++j) {
          inputs[i * dim + j] = prng.get<block>();
          if (i < 5) {
            inputs[i * dim + j] = block(i * dim + j, j + 10);
          }
        }
      }

      // for (u64 i=0; i < 5; i++){
      //     inputs[i] = i;
      // }

      // run the protocol
      macoro::sync_wait(sender.run(inputs, chl));
    }

  } else {
    vector<double> times(trait), commus(trait);
    for (u64 i = 0; i < trait; i++) {
      // receiver side
      CmpFuzzyPSI::FuzzyPsiReceiver receiver;
      block seed = oc::toBlock(123);
      receiver.init(n, n, 40, dim, metric, delta, LorH, seed, numThreads,
                    false);

      // connect to sender
      coproto::Socket chl = coproto::asioConnect(addr, false);

      // generate input data
      std::vector<block> inputs(n * dim);
      PRNG prng(oc::toBlock(789));
      for (u64 i = 0; i < n; ++i) {
        for (u64 j = 0; j < dim; ++j) {
          inputs[i * dim + j] = prng.get<block>();
          if (i < 5) {
            inputs[i * dim + j] = block(i * dim + j, j + 15);
          }
        }
      }

      // run the protocol
      macoro::sync_wait(receiver.run(inputs, chl));

      times[i] = receiver.online_time / 1000.0;
      commus[i] = receiver.online_commu / 1024.0 / 1024.0;
    }

    double avg_online_time =
        accumulate(times.begin(), times.end(), 0.0) / trait;

    double avg_com = accumulate(commus.begin(), commus.end(), 0.0) / trait;

    string mertric_str = (metric == 0) ? "inf" : std::to_string(metric);

    cout << std::format("{:^5}  𝐿{}  {:^5}  {:^5}  "
                        "{:^10.3f}  {:^10.3f}",
                        n, mertric_str, dim, delta, avg_com, avg_online_time)
         << endl;
  }

  return 0;
}