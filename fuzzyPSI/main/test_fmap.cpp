#include <ranges>
#include "macoro/async_scope.h"
#include "cryptoTools/Common/CLP.h"
#include "psi/Defines.h"
#include "psi/fmap.h"
#include "coproto/Socket/LocalAsyncSock.h"
#include "coproto/Socket/AsioSocket.h"


using namespace oc;
using namespace CmpFuzzyPSI;
using namespace secJoin;

int main(int argc, char** argv) {

    oc::CLP cmd(argc, argv);

    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 12));
    u64 trials = cmd.getOr("trials", 1);
    bool nt = cmd.getOr("nt", 1);
    auto useOle = cmd.isSet("ole");
    u64 role = cmd.getOr("r", 0);
    u64 dim = cmd.getOr("dim", 6);
    u64 delta = cmd.getOr("delta", 60);
    u64 LorH = 0;

    oc::Timer timer;

    if (role==0){
        oc::Timer timer;

        FmapSender sender;
        sender.setTimer(timer);

        coproto::Socket chl = coproto::asioConnect("localhost:1212", 0);
        PRNG prng(oc::ZeroBlock);

        std::vector<oc::block> data(n*dim);
        prng.get(data.data(), data.size());
        for (u64 i=0; i<5*dim; i++){
            data[i] = oc::block(i, i/dim);
        } 

        auto begin0 = timer.setTimePoint("begin");
        coproto::sync_wait(sender.setUp(n, n, dim, delta, LorH, prng, chl, 1));
        auto setup_comm = chl.bytesSent() + chl.bytesReceived();
        auto begin = timer.setTimePoint("begin");
        std::vector<oc::block> ID(n);
        std::vector<oc::block> oringins(n*dim);
        coproto::sync_wait(sender.fuzzyMap(data, ID, oringins, prng, chl, 1));
        auto end = timer.setTimePoint("end");

        auto offlineTime = std::chrono::duration_cast<std::chrono::milliseconds>(begin - begin0).count();
        auto onlineTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        auto onlineComm = chl.bytesSent() + chl.bytesReceived() - setup_comm;
        
        std::cout <<"Offline time: " << offlineTime << " ms " << "Online time: " << onlineTime << " ms " << std::endl;
        std::cout << "Offline comm: " << setup_comm << " Bytes " << "Online comm: " << onlineComm << " Bytes " << std::endl;


        DEBUG_LOG("offline: " << std::chrono::duration_cast<std::chrono::milliseconds>(begin - begin0).count() << " ms ");
        DEBUG_LOG("Online: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms");
        DEBUG_LOG("Communication: " << chl.bytesSent() + chl.bytesReceived() - setup_comm << " bytes");

    } else {
        oc::Timer timer;

        FmapReceiver recver;
        recver.setTimer(timer);

        coproto::Socket chl = coproto::asioConnect("localhost:1212", 1);
        PRNG prng(oc::OneBlock);

        std::vector<oc::block> data(n*dim);
        prng.get(data.data(), data.size());
        for (u64 i=0; i<5*dim; i++){
            data[i] = oc::block(i,i/dim);
        } 

        auto begin0 = timer.setTimePoint("begin");
        coproto::sync_wait(recver.setUp(n, n, dim, delta, LorH, prng, chl, 1));
        auto setup_comm = chl.bytesSent() + chl.bytesReceived();
        auto begin = timer.setTimePoint("begin");
        std::vector<oc::block> ID(n);
        std::vector<oc::block> oringins(n*dim);
        coproto::sync_wait(recver.fuzzyMap(data, ID, prng, chl, 1));
        auto end = timer.setTimePoint("end");

        auto offlineTime = std::chrono::duration_cast<std::chrono::milliseconds>(begin - begin0).count();
        auto onlineTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        auto onlineComm = chl.bytesSent() + chl.bytesReceived() - setup_comm;
        
        std::cout <<"Offline time: " << offlineTime << " ms " << "Online time: " << onlineTime << " ms " << std::endl;
        std::cout << "Offline comm: " << setup_comm << " Bytes " << "Online comm: " << onlineComm << " Bytes " << std::endl;

        DEBUG_LOG("offline: " << std::chrono::duration_cast<std::chrono::milliseconds>(begin - begin0).count() << " ms ");
        DEBUG_LOG("Online: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms");
        DEBUG_LOG("Communication: " << chl.bytesSent() + chl.bytesReceived() - setup_comm << " bytes");
    }

    return 0;

}