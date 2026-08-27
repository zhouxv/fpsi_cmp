#include "coproto/Socket/LocalAsyncSock.h"
#include "mPeqt/peqt.h"

#include <iostream>
#include <thread>

using namespace CmpFuzzyPSI;

int main() {
  bool allPassed = true;
  for (u64 equalityLength : {4ull, 8ull, 16ull, 32ull, 64ull}) {
    for (u64 dataSize : {1ull, 16ull, 64ull}) {
      auto [senderSocket, receiverSocket] = coproto::LocalAsyncSocket::makePair();
      PRNG senderPrng(oc::toBlock(100 + equalityLength));
      PRNG receiverPrng(oc::toBlock(200 + equalityLength));
      PeqtSender sender;
      PeqtReceiver receiver;

      std::thread offlineSender([&] {
        coproto::sync_wait(sender.setUp(dataSize, equalityLength, senderPrng,
                                        senderSocket));
      });
      std::thread offlineReceiver([&] {
        coproto::sync_wait(receiver.setUp(dataSize, equalityLength, receiverPrng,
                                          receiverSocket));
      });
      offlineSender.join();
      offlineReceiver.join();

      BitVector senderData(dataSize * equalityLength);
      BitVector receiverData(dataSize * equalityLength);
      PRNG dataPrng(oc::toBlock(300 + equalityLength));
      dataPrng.get(senderData.data(), senderData.sizeBytes());
      dataPrng.get(receiverData.data(), receiverData.sizeBytes());
      for (u64 i = 0; i < dataSize / 2; ++i)
        for (u64 bit = 0; bit < equalityLength; ++bit)
          receiverData[i * equalityLength + bit] =
              senderData[i * equalityLength + bit];

      BitVector senderOutput, receiverOutput;
      std::thread onlineSender([&] {
        coproto::sync_wait(sender.run(senderData, senderOutput, senderSocket));
      });
      std::thread onlineReceiver([&] {
        coproto::sync_wait(
            receiver.run(receiverData, receiverOutput, receiverSocket));
      });
      onlineSender.join();
      onlineReceiver.join();

      for (u64 i = 0; i < dataSize; ++i) {
        bool equal = true;
        for (u64 bit = 0; bit < equalityLength; ++bit)
          equal &= senderData[i * equalityLength + bit] ==
                   receiverData[i * equalityLength + bit];
        if ((senderOutput[i] ^ receiverOutput[i]) != equal) {
          std::cerr << "FAIL: length=" << equalityLength << ", size=" << dataSize
                    << ", item=" << i << '\n';
          allPassed = false;
        }
      }
    }
  }
  std::cout << (allPassed ? "ssPEQT PASS" : "ssPEQT FAIL") << '\n';
  return allPassed ? 0 : 1;
}
