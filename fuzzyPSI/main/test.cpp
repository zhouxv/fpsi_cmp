#include <iostream>
#include <vector>
#include <iomanip>
#include <future>
#include "coproto/Socket/LocalAsyncSock.h"   // for localAsyncSocketPair
#include "coproto/Socket/AsioSocket.h"
#include "cryptoTools/Common/Timer.h"
#include "psi/mImt.h"
#include "psi/Defines.h"
// #include "coproto/Proto/task.h" 

// using namespace volePSI;
using namespace oc;
using namespace CmpFuzzyPSI;

// Helper: 生成随机 blocks
std::vector<block> randomBlocks(u64 n, PRNG& prng) {
    std::vector<block> v(n);
    prng.get(oc::span<block>(v));
    return v;
}

// Helper: 生成随机 u64 向量
std::vector<u64> randomBytes(u64 n, PRNG& prng) {
    std::vector<u64> v(n);
    prng.get(oc::span<u64>(v));
    return v;
}

// ───────────────────────────────────────
// Standalone Receiver (connect to sender)
// ───────────────────────────────────────
void runStandaloneReceiver(
    const char* addr,
    u64 cmpsize, u64 dim, u64 delta, u64 metric, u64 Cmp_len,
    u64 numThreads) 
{
    std::cout << "[Receiver] Connecting to " << addr << "...\n";
    PRNG prng(_mm_set_epi32(0xfedc, 0xba98, 0x7654, 0x3210));
    Timer timer;
    timer.setTimePoint("start");

    coproto::Socket chl = coproto::asioConnect(addr, /*isServer=*/0);
    std::cout << "[Receiver] Connected.\n";

    auto hashes = randomBlocks(cmpsize, prng);
    u64 valSize = cmpsize * dim * Cmp_len;
    // auto Val = randomBytes(valSize, prng);
    std::vector<u64> Val(cmpsize*dim);

    for (u64 i=0; i<5; i++){
        Val[i] = i+10;
    }
    Val[0] = 425029;

    mIMTReceiver receiver;
    receiver.setTimer(timer);
    macoro::sync_wait(receiver.setUp(cmpsize, dim, delta, metric, Cmp_len, prng, chl, numThreads));
    auto recvBegin = timer.setTimePoint("run begin");

    auto setup_comm = chl.bytesSent() + chl.bytesReceived();
    DEBUG_LOG("Setup communication: " << setup_comm << " bytes");

    // BitVector v_s(receiver.mShareSize);
    // for (u64 i=0; i<cmpsize*dim; i++){
    //     for (u64 j=0; j<Cmp_len; j++){
    //         bool bit = (Val[i]>> j) & 1;
    //         bit ^= receiver.e[2*i*Cmp_len + 2*j];
    //         v_s[3*i*Cmp_len + j] = bit;
    //         v_s[3*i*Cmp_len + Cmp_len + j] = bit;
    //         v_s[3*i*Cmp_len + 2*Cmp_len + j] = bit;
    //     }
    // }
    BitVector v_s(receiver.mShareSize);
    for (u64 i=0; i<cmpsize*dim; i++){
        for (u64 j=0; j<Cmp_len; j++){
            bool bit = (Val[i]>> j) & 1;
            bit ^= receiver.e[2*i*Cmp_len + 2*j];
            v_s[2*i*Cmp_len + j] = bit;
            v_s[2*i*Cmp_len + Cmp_len + j] = bit;
        }
    }

    // v_s no problem
    // for (u64 i=0; i < 5; i++){
    //     std::cout << "v_s " << i << ": ";
    //     for (u64 j=0; j<Cmp_len; j++)
    //         std::cout << v_s[2*i*Cmp_len+j] << " ";
    //     std::cout << std::endl;
    // }

    macoro::sync_wait(chl.send(std::move(v_s)));

    BitVector output;
    macoro::sync_wait(receiver.run(output, chl, numThreads));
    macoro::sync_wait(chl.flush());
    // co_await receiver.run(Val, chl, numThreads);
    auto recvEnd = timer.setTimePoint("run done");

    for (u64 i=0; i<5; i++){
        std::cout << "output: " << output[i] << " ";
    }
    std::cout << std::endl;

    DEBUG_LOG("Receiver done.");
    DEBUG_LOG("Computation time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(recvEnd - recvBegin).count()
        << " ms");
    DEBUG_LOG("Communication: " << chl.bytesSent() + chl.bytesReceived() - setup_comm << " bytes");
    // chl.close();
}

// ───────────────────────────────────────
// Standalone Sender (listen for receiver)
// ───────────────────────────────────────
void runStandaloneSender(
    const char* addr,
    u64 cmpsize, u64 dim, u64 delta, u64 metric, u64 Cmp_len,
    u64 numThreads) 
{
    std::cout << "[Sender] Listening on " << addr << "...\n";
    PRNG prng(_mm_set_epi32(0x1234, 0x5678, 0x9abc, 0xdef0));
    Timer timer;
    timer.setTimePoint("start");

    coproto::Socket chl = coproto::asioConnect(addr, /*isServer=*/1);
    std::cout << "[Sender] Client connected.\n";

    auto hashes = randomBlocks(cmpsize, prng);
    u64 valSize = cmpsize * dim * Cmp_len;
    std::vector<u64> Val(cmpsize*dim*2);

    // for (u64 i=0; i<5; i++){
    //     Val[3*i] = i+9;
    //     Val[3*i+1] = i+11;
    //     Val[3*i+2] = i+15;
    // }
    // Val[0] = 11;
    // Val[1] = 13;
    // Val[2] = 3;

    for (u64 i=0; i<5; i++){
        Val[2*i] = i+11;
        Val[2*i+1] = i+16;
    }
    Val[0] = 441841;
    Val[1] = 441961;
    

    mIMTSender sender;
    sender.setTimer(timer);
    // sender.setUp(cmpsize, dim, delta, metric, Cmp_len, hashes, prng, numThreads);
    macoro::sync_wait(sender.setUp(cmpsize, dim, delta, metric, Cmp_len, prng, chl, numThreads));
    auto sendBegin = timer.setTimePoint("run begin");

    auto setup_comm = chl.bytesSent() + chl.bytesReceived();
    DEBUG_LOG("Setup communication: " << setup_comm << " bytes");

    BitVector output;
    BitVector v_s(sender.mShareSize);
    macoro::sync_wait(chl.recv(v_s));

    // std::cout << sender.mShareSize << std::endl;

    // for (u64 i=0; i < 5; i++){
    //     std::cout << "v_s " << i << ": ";
    //     for (u64 j=0; j<Cmp_len; j++)
    //         std::cout << v_s[2*i*Cmp_len+j] << " ";
    //     std::cout << std::endl;
    // }

    macoro::sync_wait(sender.run(Val, v_s, output, chl, numThreads));
    macoro::sync_wait(chl.flush());
    auto sendEnd = timer.setTimePoint("run done");

    for (u64 i=0; i<5; i++){
        std::cout << "output: " << output[i] << " ";
    }
    std::cout << std::endl;

    // for (int i=0; i < sender.e.size(); i++){
    //     if (sender.e[i] != sender.d[i]){
    //         throw std::runtime_error("mIMT Sender::e and d mismatch at index "+std::to_string(i));
    //     }
    // }
    DEBUG_LOG("Sender done.");

    DEBUG_LOG("Sender done.");
    DEBUG_LOG("Computation time: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(sendEnd - sendBegin).count()
        << " ms");
    DEBUG_LOG("Communication: " << chl.bytesSent() + chl.bytesReceived() - setup_comm << " bytes");
    // chl.close();
}

// ───────────────────────────────────────
// Usage & main
// ───────────────────────────────────────
void printUsage(const char* prog) {
    std::cout << "Usage: " << prog << " [mode] [options]\n"
        << "Modes:\n"
        << "  local      (default) run sender+receiver in one process\n"
        << "  sender     wait for remote receiver (use with --addr)\n"
        << "  receiver   connect to remote sender (use with --addr)\n"
        << "Options:\n"
        << "  --addr IP:PORT          (default: localhost:1212)\n"
        << "  --N TABLE_SIZE          (default: 1024)\n"
        << "  --dim DIM               (default: 2)\n"
        << "  --delta DELTA           (default: 1)\n"
        << "  --metric METRIC         (0=L1, 1=L2, 2=Hamming; default: 0)\n"
        << "  --cmp-len CMP_LEN       (default: 8)\n"
        << "  --threads NUM           (default: 2)\n";
}

int main(int argc, char** argv) {
    // Defaults
    std::string mode = "local";
    std::string addr = "localhost:1212";
    u64 cmpsize = 1 << 12;   // N = 1024
    u64 dim = 10;
    u64 delta = 30;
    u64 metric = 0;            // L1
    u64 Cmp_len = 19;
    u64 numThreads = 2;

    // Parse args (simple kv)
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mode" && i + 1 < argc) mode = argv[++i];
        else if (arg == "--addr" && i + 1 < argc) addr = argv[++i];
        else if (arg == "--N" && i + 1 < argc) cmpsize = std::stoull(argv[++i]);
        else if (arg == "--dim" && i + 1 < argc) dim = std::stoull(argv[++i]);
        else if (arg == "--delta" && i + 1 < argc) delta = std::stoull(argv[++i]);
        else if (arg == "--metric" && i + 1 < argc) metric = std::stoull(argv[++i]);
        else if (arg == "--cmp-len" && i + 1 < argc) Cmp_len = std::stoull(argv[++i]);
        else if (arg == "--threads" && i + 1 < argc) numThreads = std::stoull(argv[++i]);
        else if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    std::cout << "mIMT Test:\n"
              << "  mode       = " << mode << "\n"
            //   << "  addr       = " << addr << "\n"
              << "  N          = " << cmpsize << "\n"
              << "  dim        = " << dim << "\n"
              << "  delta      = " << delta << "\n"
              << "  metric     = " << metric << " (0=L_infty,1=L1,2=L2)\n"
              << "  cmp-len    = " << Cmp_len << " bits/comp\n";
            //   << "  threads    = " << numThreads << "\n";

    try {
        if (mode == "local") {
            // runLocal(cmpsize, dim, delta, metric, Cmp_len, numThreads);
        } else if (mode == "sender") {
            runStandaloneSender(addr.c_str(), cmpsize, dim, delta, metric, Cmp_len, numThreads);
        } else if (mode == "receiver") {
            runStandaloneReceiver(addr.c_str(), cmpsize, dim, delta, metric, Cmp_len, numThreads);
        } else {
            std::cerr << "Invalid mode: " << mode << "\n";
            printUsage(argv[0]);
            return 1;
        }
        std::cout << "✅ mIMT test completed.\n";

    } catch (const std::exception& e) {
        std::cerr << "❌ Exception: " << e.what() << "\n";
        return -1;
    }

    return 0;
}