#ifndef KER_H
#define KER_H
#include "GPU_wrapper.h"

typedef uint4 u128;
typedef u128 big;

__device__ void add_uint128(u128 *a, u128 augend);
__device__ u128 mul_uint128(u128 a, u128 b);
__device__ void sub_uint128(u128 *a, u128 augend);
__global__ void gpu_sparse_mult(big *a, big *b, big *c, int m, int k, uint bits, big q);
__global__ void gpu_add(uint size, big q, big *x, big *y);

#endif