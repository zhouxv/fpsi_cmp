
#include "oprfso.h"

namespace CmpFuzzyPSI{

    
    Proto OprfsoSender::setUp(u64 receiverSize, PRNG& prng, Socket& chl, u64 mNumThreads){
        mAltModWPrfSender.mUseMod2F4Ot = 1;
        ole.init(chl.fork(), prng, 0, 1, 1 << 18, 1);

        AltModPrf dm(prng.get());
        std::vector<oc::block> rk(AltModPrf::KeySize);
		for (u64 i = 0; i < AltModPrf::KeySize; ++i)
		{
			rk[i] = oc::block(i, *oc::BitIterator((u8*)&dm.mExpandedKey, i));
		}

        mAltModWPrfSender.init(receiverSize, ole, AltModPrfKeyMode::SenderOnly, AltModPrfInputMode::ReceiverOnly, dm.getKey(), rk);

        co_return;
    }

    Proto OprfsoReceiver::setUp(u64 receiverSize, PRNG& prng, Socket& chl, u64 mNumThreads){
        mAltModWPrfReceiver.mUseMod2F4Ot = 1;
        ole.init(chl.fork(), prng, 1, 1, 1 << 18, 1);

		std::vector<std::array<oc::block, 2>> sk(AltModPrf::KeySize);
		for (u64 i = 0; i < AltModPrf::KeySize; ++i)
		{
			sk[i][0] = oc::block(i, 0);
			sk[i][1] = oc::block(i, 1);
		}     

        mAltModWPrfReceiver.init(receiverSize, ole, AltModPrfKeyMode::SenderOnly, AltModPrfInputMode::ReceiverOnly, {}, sk);

        co_return;
    }

    Proto OprfsoSender::oprfSo(span<oc::block> oprfShare, PRNG& prng, Socket& chl, u64 mNumThreads){
        auto r = coproto::sync_wait(coproto::when_all_ready(
            mAltModWPrfSender.evaluate({}, oprfShare, chl, prng),
            ole.start()
        ));

        std::get<0>(r).result();
        std::get<1>(r).result();

        co_return;

    }

    Proto OprfsoReceiver::oprfSo(span<oc::block> data, span<oc::block> oprfShare, PRNG& prng, Socket& chl, u64 mNumThreads){
        auto r = coproto::sync_wait(coproto::when_all_ready(
            mAltModWPrfReceiver.evaluate(data, oprfShare, chl, prng),
            ole.start()
        ));

        std::get<0>(r).result();
        std::get<1>(r).result();
        
        co_return;
    }

}