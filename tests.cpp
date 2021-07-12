#include "includes/GPU_wrapper.h"
#include <cuda.h>
#include "includes/matrix.h"
#include "includes/GSW.h"
#include <chrono>

using namespace std::chrono;
// nvcc -c GPU_wrapper.cu kernel.cu includes/uint256_t-master/uint128_t.cpp includes/uint256_t-master/uint256_t.cpp && g++ -o tests -I/usr/local/cuda/include -L/usr/local/cuda/lib64 tests.cpp GPU_wrapper.o kernel.o  uint256_t.o uint128_t.o -lcuda -lcudart && rm kernel.o uint256_t.o uint128_t.o GPU_wrapper.o

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

void testKeygen(long n,long bits,float b,bigH q){
    parameters *p = new parameters(n,(n+1)*bits,q,b,bits);
    GSW test(p);
    test.keygen();

    matrix B(1,p->m,q);
    matmul(B,test.sk,test.pk);

    print_martix(B.vec,B.rows,B.cols);
}

void testEnc(long n,long bits,float b,bigH q,unsigned char bit){
    parameters *p = new parameters(n,(n+1)*bits,q,b,bits);
    GSW test(p);    
    test.keygen();

    cout << "Keygen Done\n";

    matrix R(test.params.n + 1,test.params.m,test.params.q);
    fillRand(R);
    unsigned char* spar = inverseG(R.vec,R.rows,R.cols,bits);
    bigH* check = sparse_mul(test.pk.vec,spar,test.pk.rows,test.pk.cols,q);
    if(bit == 1){
        cout << "adding G\n";
        add_cpu(check,test.G.vec,test.G.rows*test.G.cols,q);
    }
    cout << "moving to GPU\n";
    bigH* ret = encrypt(test.pk.vec,R.vec,test.G.vec,test.params.q,test.pk.rows,test.pk.cols,bits,bit);
    for(int i = 0;i<R.rows*R.cols;i++){
        assert(check[i] == ret[i]);
    }
    cout << "test passed\n";
}

void test_speed(long n,long bits,float b,bigH q,unsigned char bit){
    parameters *p = new parameters(n,(n+1)*bits,q,b,bits);
    GSW test(p);

    auto start = high_resolution_clock::now();

    test.keygen();

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<seconds>(stop - start);

    cout << "keygen done, time taken: " << duration.count() << "\n";
    
    start = high_resolution_clock::now();
    matrix m;
    test.encryptBit(bit,m);
    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);

    cout << "Enc. done, time taken: " << duration.count() << "\n";
}

void test_w(long n,long bits,float b,bigH q,unsigned char bit){
    parameters *p = new parameters(n,(n+1)*bits,q,b,bits);
    GSW test(p);

    print_martix(test.W,test.params.bits,1);
}

void test_d(){
    int k, l;
    cout << "Enter k and L: ";
    cin >> k >> l;
    parameters *p = setup(k,l);
    // parameters *p = new parameters(n,(n+1)*bits,q,b,bits);
    GSW test(p);
    test.keygen();
    cout << "Keygen Done\n";

    matrix m;
    int bit1,bit2;

    for (int i = 0;i<10;i++){
        bit1 = rand()%2;
        cout << "Encrypting.. ";
        test.encryptBit(bit1,m);
        cout << "Done \nDecrypting..\n";
        bit2 = test.decryptBit(m);
        assert(bit1 == bit2);
        if(bit1 != bit2){ cout << "Test"  << i << " failed!!\n"; return;} 
        else cout << "Test " << i + 1 << " Passed\n\n";
    }
    cout << "All tests passed\n";
}

void test_saving(){
    int k, l;
    cout << "Enter k and L: ";
    cin >> k >> l;
    parameters *p = setup(k,l);
    GSW test(p);
    test.keygen();
    test.saveState();
}

int main(){
    test_d();
}