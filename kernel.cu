#include "includes/kernel.cuh"

__device__ void add_uint128 (u128* a,u128 augend)
{
    u128 addend = a[0];
    u128 res;
    asm ("add.cc.u32      %0, %4, %8;\n\t"
         "addc.cc.u32     %1, %5, %9;\n\t"
         "addc.cc.u32     %2, %6, %10;\n\t"
         "addc.u32        %3, %7, %11;\n\t"
         : "=r"(res.x), "=r"(res.y), "=r"(res.z), "=r"(res.w)
         : "r"(addend.x), "r"(addend.y), "r"(addend.z), "r"(addend.w),
           "r"(augend.x), "r"(augend.y), "r"(augend.z), "r"(augend.w));
    a[0] = res;
}


__device__ u128 mul_uint128 (u128 a, u128 b)
{
    u128 res;
    asm ("{\n\t"
         "mul.lo.u32      %0, %4, %8;    \n\t"
         "mul.hi.u32      %1, %4, %8;    \n\t"
         "mad.lo.cc.u32   %1, %4, %9, %1;\n\t"
         "madc.hi.u32     %2, %4, %9,  0;\n\t"
         "mad.lo.cc.u32   %1, %5, %8, %1;\n\t"
         "madc.hi.cc.u32  %2, %5, %8, %2;\n\t"
         "madc.hi.u32     %3, %4,%10,  0;\n\t"
         "mad.lo.cc.u32   %2, %4,%10, %2;\n\t"
         "madc.hi.u32     %3, %5, %9, %3;\n\t"
         "mad.lo.cc.u32   %2, %5, %9, %2;\n\t"
         "madc.hi.u32     %3, %6, %8, %3;\n\t"
         "mad.lo.cc.u32   %2, %6, %8, %2;\n\t"
         "madc.lo.u32     %3, %4,%11, %3;\n\t"
         "mad.lo.u32      %3, %5,%10, %3;\n\t"
         "mad.lo.u32      %3, %6, %9, %3;\n\t"
         "mad.lo.u32      %3, %7, %8, %3;\n\t"
         "}"
         : "=r"(res.x), "=r"(res.y), "=r"(res.z), "=r"(res.w)
         : "r"(a.x), "r"(a.y), "r"(a.z), "r"(a.w),
           "r"(b.x), "r"(b.y), "r"(b.z), "r"(b.w));
    return res;
}

__device__  void sub_uint128 (u128* a,u128 augend)
{
    u128 subend = *a;
    u128 res;
    uint carry = 0;
    asm ("sub.cc.u32      %0, %5, %9;\n\t"
         "subc.cc.u32     %1, %6, %10;\n\t"
         "subc.cc.u32     %2, %7, %11;\n\t"
         "subc.cc.u32     %3, %8, %12;\n\t"
         "addc.u32        %4, 0, 0;\n\t"
         : "=r"(res.x), "=r"(res.y), "=r"(res.z), "=r"(res.w),"=r"(carry)
         : "r"(subend.x), "r"(subend.y), "r"(subend.z), "r"(subend.w), "r"(augend.x), "r"(augend.y), "r"(augend.z), "r"(augend.w));
    *a = (carry) ? res : subend;
}

__global__ void gpu_sparse_mult(big *a,big *b, big *c, int m, int k,uint bits,big q)
{ 
    assert(k/m == bits);
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    big sum;
    sum.x = 0;sum.y = 0;sum.z= 0;sum.w = 0;
    if( col < k && row < m) 
    {
        for(int i = 0; i < m; i++) 
        {
            big B = b[i*k + col];
            uint loc = bits,word = B.x,counter = 0;
            while(loc != 0){
                if(word&1){
                    add_uint128(&sum,a[row*k + i*bits + counter]);
                    sub_uint128(&sum,q);
                }
                word >>= 1;
                loc--;
                counter++;
                if(counter == 32){ word = B.y;}
                else if(counter == 64){ word = B.z;}
                else if(counter == 96){ word = B.w;}
            }
        }
        c[row * k + col] = sum;
    }
} 

__global__
void gpu_add(uint size,big q, big *x, big *y) // adds x to y mod q
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < size; i += stride){
    // y[i] = x[i] + y[i];
    add_uint128(&y[i],x[i]);
    sub_uint128(&y[i],q);
  }
}