nvcc -c GPU-wrapper.cu kernel.cu includes/uint256_t-master/uint128_t.cpp includes/uint256_t-master/uint256_t.cpp csprng.cpp
g++ -o tests -pthread -I/usr/local/cuda/include -L/usr/local/cuda/lib64 tests.cpp matrix.cpp GSW.cpp somewhat.cpp csprng.o GPU-wrapper.o kernel.o  uint256_t.o uint128_t.o -lntl -lgmp -lm -lcuda -lcudart 
rm kernel.o uint256_t.o uint128_t.o GPU-wrapper.o csprng.o