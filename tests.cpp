#include "includes/GPU_wrapper.h"
#include <cuda.h>
#include "includes/matrix.h"
#include "includes/GSW.h"

// nvcc -c GPU_wrapper.cu kernel.cu includes/uint256_t-master/uint128_t.cpp includes/uint256_t-master/uint256_t.cpp && g++ -o tests -I/usr/local/cuda/include -L/usr/local/cuda/lib64 tests.cpp GPU_wrapper.o kernel.o  uint256_t.o uint128_t.o -lcuda -lcudart && rm kernel.o uint256_t.o uint128_t.o GPU_wrapper.o

unsigned char* inverseG(bigH* matrix,int n,int m,int bits){
    unsigned char* ret = (unsigned char* )malloc(sizeof(unsigned char)*n*bits*m);
    bigH temp;
    for(int i = 0;i<m;i++){
        for(int j = 0;j<n;j++){
            temp = matrix[j*m+i];
            // cout << temp << endl;
            for(int k = 0;k<bits;k++){
                ret[(j*bits + k)*m + i] = temp&1;
                temp >>= 1;
            }
        }
    }
    return ret;
}

void test_inverseG(int m,int n,int bits,bigH q){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;
    bigH* test = (bigH*)malloc(sizeof(bigH)*m*n);
    for(int i = 0;i<m*n;i++){
        test[i] = dis(gen);
        test[i] <<= 64;
        test[i] |= dis(gen);
        test[i] %= q;
    }
    print_martix(test,m,n);
    cout << endl;
    unsigned char* result = inverseG(test,m,n,bits);
    print_martix(result,n,n);
    free(result);
    free(test);
}

void test_mul(int m,int n,int bits,bigH q){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;
    bigH* test = (bigH*)malloc(sizeof(bigH)*m*n);
    for(int i = 0;i<m*n;i++){
        test[i] = dis(gen);
        test[i] <<= 64;
        test[i] |= dis(gen);
        test[i] %= q;
    }
    print_martix(test,m,n);
    unsigned char* sparse = inverseG(test,m,n,bits);
    print_martix(sparse,n,n);
    
    bigH* result = sparse_mul(test,sparse,m,n,q);
    print_martix(test,m,n);

    free(sparse);
    free(test);
    free(result);
    
}

void test_main(int m,int n,int bits,bigH q)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;
    bigH* test = (bigH*)malloc(sizeof(bigH)*m*n);
    for(int i = 0;i<m*n;i++){
        test[i] = dis(gen);
        test[i] <<= 64;
        test[i] |= dis(gen);
        test[i] %= q;
    }
    unsigned char* sparse = inverseG(test,m,n,bits);
    bigH* result = sparse_mul(test,sparse,m,n,q);

    cout << "moving to GPU...\n";
    MAIN_TEST_GPU(test,test,result,q,bits,m,n);
    free(result);

}

void test_matmul(int m,int n,int k,bigH q){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;
    

    matrix test1(m,n,q);
    matrix test2(n,k,q);
    matrix result(m,k,q);

    for(int i = 0;i<m*n;i++){
        test1.vec[i] = dis(gen);
        test1.vec[i] <<= 64;
        test1.vec[i] |= dis(gen);
        if(!(test1.vec[i] < q))
        test1.vec[i] %= q;
    }

    cout << "Generated 1\n";

    for(int i = 0;i<k*n;i++){
        test2.vec[i] = dis(gen);
        test2.vec[i] <<= 64;
        test2.vec[i] += dis(gen);
        if(!(test2.vec[i] < q))
        test2.vec[i] %= q;
    }
    
    // print_martix(test1.vec,m,n);
    // print_martix(test2.vec,n,k);

    cout << "Generated 2\n";

    matmul(result,test1,test2);

    cout << "done\n";
    // print_martix(result.vec,m,k);

}

void test_invG(int n,int bits,bigH q){
    matrix G(n,n*bits,q);
    genGadget(bits,G);
    print_martix(G.vec,G.rows,G.cols);
}

void test_gaussian(int n,float b,bigH q){
    matrix test(1,n,q);
    for(int i = 0;i<n;i++){
        test.vec[i] = 0;
    }
    gaussian(test,b);
    print_martix(test.vec,1,n);
}

int main(){
    cout << "Enter n and bits:";
    uint n;
    float b;
    cin >> n>>b;
    bigH q = 100;
    test_gaussian(n,b,q);
}