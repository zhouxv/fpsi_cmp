#include <ranges>
#include "macoro/async_scope.h"
#include "cryptoTools/Common/CLP.h"
#include "psi/Defines.h"
#include "oprf_so/oprfso.h"
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

    oc::Timer timer;

    std::cout << n <<  std::endl;

    if (role==0){
        oc::Timer timer;

        OprfsoSender sender;
        sender.setTimer(timer);

        std::vector<oc::block> share(n);
        
        coproto::Socket chl = coproto::asioConnect("localhost:1212", 0);
        PRNG prng(oc::ZeroBlock);

        macoro::sync_wait(sender.setUp(n, prng, chl, 1));
        auto begin = timer.setTimePoint("begin");
        macoro::sync_wait(sender.oprfSo(share, prng, chl));
        auto end = timer.setTimePoint("end");

        // macoro::sync_wait(chl.flush());
        // chl.close();
        AltModPrf::KeyType kk = sender.mAltModWPrfSender.getKey();

        macoro::sync_wait(chl.send(kk));
        macoro::sync_wait(chl.send(share));

        std::cout << "AltModWPrf n:" << n << ", " <<
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << " ns " << std::endl;

    } else {
        oc::Timer timer;

        OprfsoReceiver recver;
        recver.setTimer(timer);

        std::vector<oc::block> data(n);
        std::vector<oc::block> share(n);

        coproto::Socket chl = coproto::asioConnect("localhost:1212", 1);
        PRNG prng(oc::OneBlock);

        prng.get(data.data(), n);

        macoro::sync_wait(recver.setUp(n, prng, chl));
        auto begin = timer.setTimePoint("begin");
        macoro::sync_wait(recver.oprfSo(data, share, prng, chl));
        auto end = timer.setTimePoint("end");

        // macoro::sync_wait(chl.flush());
        // chl.close();

        AltModPrf::KeyType kk;
        std::vector<oc::block> sShare(n);
        macoro::sync_wait(chl.recv(kk));
        macoro::sync_wait(chl.recv(sShare));
        AltModPrf prf(kk);
        std::vector<oc::block> y(n);
        prf.eval(data, y);

        for (u64 i=0; i<10; i++){
            if ((sShare[i] ^ share[i]) != prf.eval(data[i]))
                std::cout << i << " wrong! " <<  std::endl;
            else    
                std::cout << prf.eval(data[i]) << " not wrong! " <<  std::endl;
        }

        std::cout << "AltModWPrf n:" << n << ", " <<
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "ns" << std::endl;
    }

    return 0;
}