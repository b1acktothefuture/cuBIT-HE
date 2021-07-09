#include "includes/matrix.h"

matrix::matrix(long n,long m,bigH mod){
    vec = (bigH*)malloc(sizeof(bigH)*m*n);
    rows = n;
    cols = m;
    q = mod;
}

bigH matrix::get(long i,long j){
    assert(i>=0 && i<rows && j>=0 && j<cols);
    return vec[cols*i + j];
}

void matrix::set(long i,long j,bigH m){
    assert(i>=0 && i<rows && j>=0 && j<cols);
    if(m>q) m = m%q;
    vec[cols*i + j] = m;
}

void matrix::reshape(long m,long n){
    free(vec);
    vec = (bigH*)malloc(sizeof(bigH)*m*n);
    rows = m;
    cols = n;
}

void matmul(matrix& x,const matrix& m,const matrix& n){
    assert(m.cols == n.rows);
    biggerH holder,h1,h2;
    bigH sum,q;
    q = m.q;
    bigH* one = m.vec,*two = n.vec,*res = x.vec;
    long M = m.rows,N = n.cols,K = n.cols;
    for(long i = 0;i<m.rows;i++)
    for(long j = 0;j<n.cols;j++){
        sum = 0;
        for(int k = 0;k<m.cols;k++){
            h1 = biggerH(one[N*i + k]);
            h2 = biggerH(two[k*K + j]);
            holder = (h1*h2);
            if(!(holder < q)) holder %= q;
            sum += holder.lower();
            if(!(sum < q)) sum -= q;
        }
        res[K*i + j] = sum;
    }
}

bigH* sparse_mul(bigH* A,unsigned char* B,int m,int n,bigH q){  // m*n & n*n
    bigH* result = (bigH* )malloc(sizeof(bigH)*m*n);
    bigH sum;
    for(int i = 0;i<m;i++)
    for(int j = 0;j<n;j++){
        sum = 0;
        for(int k = 0;k<n;k++){
            if(B[k*n + j]){
                sum += A[i*n + k];
                if(sum >= q) sum -= q;
            }
        }
        result[n*i + j] = sum;
    }
    return result;
}

void print_martix(bigH* b,int rows,int cols){
    unsigned int u = UINT32_MAX;
    for(int i = 0;i<rows;i++){
    for(int j = 0;j<cols;j++)
        std::cout << (b[i*cols + j]) << " ";
    std::cout << std::endl;
    }
    cout << endl;
}


void print_martix(unsigned char* b,int rows,int cols){
    for(int i = 0;i<rows;i++){
    for(int j = 0;j<cols;j++)
        std::cout << int(b[i*cols + j]) << " ";
    std::cout << std::endl;
    }
    cout << endl;
}