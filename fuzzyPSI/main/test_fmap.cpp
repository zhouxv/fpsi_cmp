#include "coproto/Socket/AsioSocket.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "cryptoTools/Common/CLP.h"
#include "macoro/async_scope.h"
#include "psi/Defines.h"
#include "psi/fmap.h"
#include <macoro/when_all.h>
#include <ranges>

using namespace oc;
using namespace CmpFuzzyPSI;
using namespace secJoin;

void printUsage(const char *prog) {
  std::cout << "Usage: " << prog << " [options]\n"
            << "  Options:\n"
            << "    -n <N>          : Set size (direct), default: 4096\n"
            << "    -nn <N>         : Set size (logarithm), input size = 2^nn "
               "(overrides -n), default: 12\n"
            << "    -dim <N>        : Dimension of the points, default: 6\n"
            << "    -delta <N>      : Distance threshold δ for fuzzy matching, "
               "default: 60\n"
            << "    -ip <addr>      : Server IP address, default: localhost\n"
            << "    -port <N>       : Server port number, default: 1212\n"
            << "    -trait <N>      : Number of trials for averaging results, "
               "default: 5\n"
            << "    -h/--help       : Print this help message\n";
}

int main(int argc, char **argv) {
  oc::CLP cmd(argc, argv);

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
  const u64 LorH = cmd.getOr<u64>("LorH", 0);
  const u64 trait = cmd.getOr<u64>("trait", 5);
  const std::string ip = cmd.getOr<std::string>("ip", "localhost");
  const u64 port = cmd.getOr<u64>("port", 1213);
  std::string addr = ip + ":" + std::to_string(port);

  std::vector<double> offline_times(trait), offline_commus(trait),
      online_times(trait), online_commus(trait);

  for (u64 i = 0; i < trait; i++) {
    FmapSender sender;
    FmapReceiver recver;
    // set timer
    oc::Timer timer;

    PRNG send_prng(oc::ZeroBlock);
    PRNG recv_prng(oc::OneBlock);

    // generate data
    std::vector<oc::block> sender_data(n * dim);
    send_prng.get(sender_data.data(), sender_data.size());
    for (u64 i = 0; i < 5 * dim; i++) {
      sender_data[i] = oc::block(i, i / dim);
    }

    std::vector<oc::block> recv_data(n * dim);
    recv_prng.get(recv_data.data(), recv_data.size());
    for (u64 i = 0; i < 5 * dim; i++) {
      recv_data[i] = oc::block(i, i / dim);
    }

    // setup sockets
    coproto::Socket send_chl, recv_chl;
    auto init_chl = [&](bool is_server) {
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

    std::vector<oc::block> send_ID(n), recv_ID(n);
    std::vector<oc::block> send_oringins(n * dim), recv_oringins(n * dim);

    // setup phase
    auto begin0 = timer.setTimePoint("begin");

    auto send_setup = [&]() {
      macoro::sync_wait(
          sender.setUp(n, n, dim, delta, LorH, send_prng, send_chl, 1));
    };
    auto recv_setup = [&]() {
      macoro::sync_wait(
          recver.setUp(n, n, dim, delta, LorH, recv_prng, recv_chl, 1));
    };
    std::thread send_setup_th(send_setup);
    std::thread recv_setup_th(recv_setup);
    send_setup_th.join();
    recv_setup_th.join();

    auto begin = timer.setTimePoint("begin");
    auto setup_comm = recv_chl.bytesSent() + recv_chl.bytesReceived();

    // online phase

    auto send_run = [&]() {
      macoro::sync_wait(sender.fuzzyMap(sender_data, send_ID, send_oringins,
                                        send_prng, send_chl, 1));
    };
    auto recv_run = [&]() {
      macoro::sync_wait(
          recver.fuzzyMap(recv_data, recv_ID, recv_prng, recv_chl, 1));
    };
    std::thread send_run_th(send_run);
    std::thread recv_run_th(recv_run);
    send_run_th.join();
    recv_run_th.join();

    auto end = timer.setTimePoint("end");

    auto offlineTime =
        std::chrono::duration_cast<std::chrono::milliseconds>(begin - begin0)
            .count();
    auto onlineTime =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
            .count();
    auto onlineComm =
        recv_chl.bytesSent() + recv_chl.bytesReceived() - setup_comm;

    offline_times[i] = offlineTime / 1000.0;
    online_times[i] = onlineTime / 1000.0;
    offline_commus[i] = setup_comm / 1024.0 / 1024.0;
    online_commus[i] = onlineComm / 1024.0 / 1024.0;
  }

  double avg_offline_commu =
      accumulate(offline_commus.begin(), offline_commus.end(), 0.0) / trait;
  double avg_offline_time =
      accumulate(offline_times.begin(), offline_times.end(), 0.0) / trait;
  double avg_online_time =
      accumulate(online_times.begin(), online_times.end(), 0.0) / trait;
  double avg_online_commu =
      accumulate(online_commus.begin(), online_commus.end(), 0.0) / trait;

  std::cout << std::format("{:^5}  {:^5}  {:^5}  "
                           "{:^10.3f}  {:^10.3f}  {:^10.3f}  {:^10.3f}",
                           n, dim, delta, avg_online_commu, avg_online_time,
                           avg_offline_commu, avg_offline_time)
            << std::endl;

  return 0;
}