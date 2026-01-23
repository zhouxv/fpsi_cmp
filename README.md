debug模式：-DENABLE_DEBUG_LOG=1 -DCMAKE_BUILD_TYPE=Debug（会慢很多）

mImt:

cmake -S . -B ./build && cmake --build ./build

./build/mtest --mode sender --addr localhost:1213 --N 1000 --metric 0

./build/mtest --mode receiver --addr localhost:1213 --N 1000 --metric 0

main:

cmake -S . -B ./build && cmake --build ./build

./build/main -r 0 -metric 1

./build/main -r 1 -metric 1

m_peqt_test (cmake --build ./build -t m_peqt_test)

./build/m_peqt_test -sender -n 4096 -l 40
./build/m_peqt_test -receiver -n 4096 -l 40

m_oprf_test

./build/m_oprf_test -r 0 -nn 12
./build/m_oprf_test -r 1 -nn 12

m_fmap_test

./build/m_fmap_test -r 0 -nn 12
./build/m_fmap_test -r 1 -nn 12

需要补充offline部分:
psiL2(mod L OLE), mImt(AND pairs)