#ifndef U128_H
#define U128_H

#include <iostream>

#include <bits/stdc++.h>
#include <assert.h>
#include <malloc.h>
#include "csprng.hpp"


using namespace std;

#define BLOCK_SIZE 16

#include "uint256_t-master/uint256_t.h"


typedef uint128_t bigH;
typedef uint256_t biggerH;

void MAIN_TEST_GPU(bigH* A_h,bigH* R_h,bigH* result_h,bigH g,uint bits,int n,int m);
bigH* encrypt(bigH* pk_h,bigH* R_h,bigH* G_h,bigH q_h,uint n,uint m,uint bits,unsigned char bit);
bigH* encryptFast(bigH* pk_h, bigH* G_h,bigH q_h,uint n,uint m,uint bits,unsigned char bit);

#endif