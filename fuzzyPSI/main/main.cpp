
#include "coproto/Socket/AsioSocket.h"
#include "coproto/Socket/LocalAsyncSock.h" // for localAsyncSocketPair
#include "cryptoTools/Common/CLP.h"
#include "cryptoTools/Common/Timer.h"

#include "debug.h"
#include "psi/psi.h"
#include "volePSI/RsOpprf.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace volePSI;
using namespace oc;

void printUsage(const char *prog) {
  std::cout
      << "Usage: " << prog << " [options]\n"
      << "  Options:\n"
      << "    -n <N>          : Set size (direct), default: 4096\n"
      << "    -nn <N>         : Set size (logarithm), input size = 2^nn "
         "(overrides -n), default: 12\n"
      << "    -dim <N>        : Dimension of the points, default: 6\n"
      << "    -delta <N>      : Distance threshold δ for fuzzy matching, "
         "default: 60\n"
      << "    -metric <N>     : Distance metric (0: L∞, 1: L₁, 2: L₂), "
         "default: 0\n"
      << "    -ip <addr>      : Server IP address, default: localhost\n"
      << "    -port <N>       : Server port number, default: 1212\n"
      << "    -trait <N>      : Number of trials for averaging results, "
         "default: 5\n"
      << "    -out <file>     : Append the averaged result to a CSV file\n"
      << "    -offlineCache <dir> : One-time L2 offline-material cache\n"
      << "    -offlineOnly    : Generate L2 cache then exit (requires "
         "-offlineCache)\n"
      << "    -h/--help       : Print this help message\n";
}

int main(int argc, char **argv) {
  CLP cmd;
  cmd.parse(argc, argv);

  // print help message
  if (cmd.isSet("h") || cmd.isSet("help")) {
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
  const std::string ip = cmd.getOr<std::string>("ip", "localhost");
  const u64 port = cmd.getOr<u64>("port", 1212);
  const u64 trait = cmd.getOr<u64>("trait", 5);
  const std::string offlineCache = cmd.getOr<std::string>("offlineCache", "");
  const bool offlineOnly = cmd.isSet("offlineOnly");
  if (offlineOnly && offlineCache.empty()) {
    std::cerr << "-offlineOnly requires -offlineCache <dir>\n";
    return 1;
  }
  if (!offlineCache.empty() && metric != 2) {
    std::cerr << "offline cache currently supports only L2\n";
    return 1;
  }
  if (!offlineCache.empty() && trait != 1) {
    std::cerr << "offline cache is one-time material; use -trait 1\n";
    return 1;
  }

  vector<double> online_times(trait), online_commus(trait),
      offline_times(trait), offline_commus(trait);
  for (u64 i = 0; i < trait; i++) {
    // sender side
    CmpFuzzyPSI::FuzzyPsiSender sender;
    CmpFuzzyPSI::FuzzyPsiReceiver receiver;
    block seed = oc::toBlock(123);
    sender.init(n, n, 40, dim, metric, delta, seed, numThreads, false);
    receiver.init(n, n, 40, dim, metric, delta, seed, numThreads, false);
    if (!offlineCache.empty()) {
      sender.setL2OfflineCache(offlineCache, offlineOnly);
      receiver.setL2OfflineCache(offlineCache, offlineOnly);
    }

    // generateinput data
    std::vector<block> sender_inputs(n * dim);
    PRNG sender_prng(oc::toBlock(456));
    for (u64 i = 0; i < n; ++i) {
      for (u64 j = 0; j < dim; ++j) {
        sender_inputs[i * dim + j] = sender_prng.get<block>();
        if (i < 5) {
          sender_inputs[i * dim + j] = block(i * dim + j, j + 10);
        }
      }
    }

    std::vector<block> recv_inputs(n * dim);
    PRNG recv_prng(oc::toBlock(789));
    for (u64 i = 0; i < n; ++i) {
      for (u64 j = 0; j < dim; ++j) {
        recv_inputs[i * dim + j] = recv_prng.get<block>();
        if (i < 5) {
          recv_inputs[i * dim + j] = block(i * dim + j, j + 15);
        }
      }
    }

    // setup channel

    // connect to receiver
    coproto::Socket send_chl, recv_chl;
    auto init_chl = [&](bool is_server) {
      std::string addr = ip + ":" + std::to_string(port + i);
      if (is_server) {
        send_chl = coproto::asioConnect(addr, true);
      } else {
        recv_chl = coproto::asioConnect(addr, false);
      }
    };

    std::thread sender_sock(init_chl, true);
    std::thread recv_sock(init_chl, false);
    sender_sock.join();
    recv_sock.join();

    // run the protocol
    auto send_run = [&]() {
      macoro::sync_wait(sender.run(sender_inputs, send_chl));
    };
    auto recv_run = [&]() {
      macoro::sync_wait(receiver.run(recv_inputs, recv_chl));
    };

    std::thread sender_run_th(send_run);
    std::thread recv_run_th(recv_run);
    sender_run_th.join();
    recv_run_th.join();

    if (offlineOnly) {
      std::cout << "L2 offline material generated in " << offlineCache
                << " (one-time use)\n";
      return 0;
    }

    online_times[i] = receiver.online_time / 1000.0;
    online_commus[i] = receiver.online_commu / 1024.0 / 1024.0;
    offline_times[i] = receiver.offline_time / 1000.0;
    offline_commus[i] = receiver.offline_commu / 1024.0 / 1024.0;
  }

  double avg_online_time =
      accumulate(online_times.begin(), online_times.end(), 0.0) / trait;
  double avg_online_com =
      accumulate(online_commus.begin(), online_commus.end(), 0.0) / trait;

  double avg_offline_time =
      accumulate(offline_times.begin(), offline_times.end(), 0.0) / trait;
  double avg_offline_com =
      accumulate(offline_commus.begin(), offline_commus.end(), 0.0) / trait;

  double total_time = avg_online_time + avg_offline_time;
  double total_com = avg_online_com + avg_offline_com;

  string metric_str = (metric == 0) ? "inf" : std::to_string(metric);

  cout
      << std::format(
             "{:^5}  𝐿{}  {:^5}  {:^5}  "
             "{:^10.3f}  {:^10.3f}  {:^10.3f}  {:^10.3f}  {:^10.3f}  {:^10.3f}",
             n, metric_str, dim, delta, avg_online_com, avg_online_time,
             avg_offline_com, avg_offline_time, total_com, total_time)
      << endl;

  if (cmd.isSet("out")) {
    const auto outputPath = cmd.get<std::string>("out");
    std::ifstream existing(outputPath, std::ios::binary | std::ios::ate);
    const bool writeHeader = !existing || existing.tellg() == 0;
    std::ofstream output(outputPath, std::ios::app);
    if (!output) {
      throw std::runtime_error("failed to open result file: " + outputPath);
    }
    if (writeHeader) {
      output << "Protocol,Metric,Dim,Delta,Size,Online_Com.(MB),Online(s),"
                "Offline_Com.(MB),Offline(s)\n";
    }

    const string csv_metric = (metric == 0) ? "Linf" : "L" + metric_str;
    output << "fpsi_cmp," << csv_metric << ',' << dim << ',' << delta << ','
           << n << ',' << std::fixed << std::setprecision(3) << avg_online_com
           << ',' << avg_online_time << ',' << avg_offline_com << ','
           << avg_offline_time << ',' << total_com << ',' << total_time << '\n';
  }

  return 0;
}
