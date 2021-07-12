#include "includes/GPU_wrapper.h"
#include "includes/kernel.cuh"

/******************************************************************************/
// helper functions

void print_martix(big* b,int rows,int cols){
    for(int i = 0;i<rows;i++){
    for(int j = 0;j<cols;j++)
        std::cout << b[i*cols + j].x << " ";
    std::cout << std::endl;
    }
    cout << endl;
}

big bighToBig(bigH g){
    uint AND = UINT32_MAX;
    big q;
    uint64_t word = g.lower();
    q.x = word&AND;
    word >>= 32;
    q.y = word&AND;
    word = g.upper();
    q.z = word&AND;
    word >>= 32;
    q.w = word&AND;
    return q;
}

bigH bigTobigH(big q){
    bigH g(0);
    g += q.w;
    g <<= 32;
    g += q.z;
    g <<= 32;
    g += q.y;
    g <<= 32;
    g += q.x;
    return g;
}

big* convert(bigH* matrix,long size){
    big* ret = (big *)malloc(sizeof(big)*size);
    for(long i =0;i<size;i++){
        ret[i] = bighToBig(matrix[i]);
    }
    return ret;
}

bigH* convertBack(big* matrix,long size){
    bigH* ret = (bigH *)malloc(sizeof(bigH)*size);
    for(long i =0;i<size;i++){
        ret[i] = bigTobigH(matrix[i]);
    }
    return ret;
}


void encryptHelper(big* A,big* R,big* G,big* result,big q,uint bits,unsigned char bit,int n,int m){

    unsigned int grid_rows = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (m + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

    gpu_sparse_mult<<<dimGrid,dimBlock>>>(A,R,result,n,m,bits,q);
    cudaThreadSynchronize();
    if(bit == 1){
        int blockSize = 256;
        int numBlocks = (m*n + blockSize - 1) / blockSize;
        gpu_add<<<numBlocks, blockSize>>>(n*m,q, G, result);
    }
}

void fillRandom(big* R,long size,long len){ // will work only for modulus size strictly less than 128
    long words = size/32;
    long rem = size%32;
    uint t = 1,arr[4];
    t <<= rem;
    duthomhas::csprng rng;
    arr[0] = 0;arr[1] = 0;arr[2] = 0;arr[3] = 0;
    uint AND = 0xFFFFFFFF;
    for(long i = 0;i<len;i++){
        for(long j = 0;j<words;j++){
            arr[j] = rng(uint());
        }
        arr[words] = rng(uint())%t; // pathological case, when rem = 31
        R[i].x = arr[0];R[i].y = arr[1];R[i].z = arr[2];R[i].w = arr[3];
    }
    ~rng();
}
/******************************************************************************/
// Tests

void test(big* A,big* R,big* result,big q,uint bits,uint n,uint m){
    big* d_matrix1,*d_matrix2,*d_result;
    cudaMalloc((void **)&d_matrix1,sizeof(u128)*n*m);
    cudaMalloc((void **)&d_matrix2,sizeof(u128)*n*m);
    cudaMalloc((void **)&d_result,sizeof(u128)*n*m);

    cudaMemcpy(d_matrix1,A,sizeof(u128)*m*n,cudaMemcpyHostToDevice);
    cudaMemcpy(d_matrix2,R,sizeof(u128)*m*n,cudaMemcpyHostToDevice);

    unsigned int grid_rows = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (m + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

    gpu_sparse_mult<<<dimGrid,dimBlock>>>(d_matrix1,d_matrix2,d_result,n,m,bits,q);

    big* chck = (big* )malloc(sizeof(big)*m*n);
    cudaMemcpy(chck,d_result,sizeof(u128)*m*n,cudaMemcpyDeviceToHost);

    bool pass = 1;
    for(int i = 0;i<m*n;i++){
        if(chck[i].x != result[i].x || chck[i].y != result[i].y || chck[i].z != result[i].z || chck[i].w != result[i].w){
            cout << "test failed\n\n";
            pass = 0;
            break;
        }
    }

    if(pass) cout << "test passed\n";

    cudaFree(d_result);
    cudaFree(d_matrix1);
    cudaFree(d_matrix2);

    free(chck);
}

void MAIN_TEST_GPU(bigH* A_h,bigH* R_h,bigH* result_h,bigH g,uint bits,int n,int m){
    big* A = convert(A_h,m*n);
    big* R = convert(R_h,m*n);
    big* result = convert(result_h,m*n);

    big q = bighToBig(g);

    test(A,R,result,q,bits,n,m);
    
    free(A);
    free(R);
    free(result);
    
}


/******************************************************************************/


bigH* encrypt(bigH* pk_h,bigH* R_h,bigH* G_h,bigH q_h,uint n,uint m,uint bits,unsigned char bit){
    long size = m*n;

// converting between u128 representation on cpu and GPU

    big* PK = convert(pk_h,size);
    big* R = convert(R_h,size);
    big* G = convert(G_h,size);
    big g = bighToBig(q_h);
    
    big* pk_d,*R_d,*G_d,*result_d;
/*
copying keys,random matrix R, and Gadget matrix G from host to device
*/
    cudaMalloc((void **)&pk_d,sizeof(big)*size); 
    cudaMalloc((void **)&R_d,sizeof(big)*size);
    cudaMalloc((void **)&G_d,sizeof(big)*size);
    cudaMalloc((void **)&result_d,sizeof(big)*size);

    cudaMemcpy(pk_d,PK,sizeof(u128)*size,cudaMemcpyHostToDevice);
    cudaMemcpy(R_d,R,sizeof(u128)*size,cudaMemcpyHostToDevice);
    cudaMemcpy(G_d,G,sizeof(u128)*size,cudaMemcpyHostToDevice);

// freeing unncessary representation of pk and G on CPU

    free(PK);
    free(G);
// actual encryption function
    encryptHelper(pk_d,R_d,G_d,result_d,g,bits,bit,n,m);
// copying the result back to host
    cudaMemcpy(R,result_d,sizeof(u128)*size,cudaMemcpyDeviceToHost);;
    bigH* cipherText = convertBack(R,size);

// free memory from GPU
    cudaFree(pk_d);
    cudaFree(R_d);
    cudaFree(G_d);
    cudaFree(result_d);

    return cipherText;

}

bigH* encryptFast(bigH* pk_h, bigH* G_h,bigH q_h,uint n,uint m,uint bits,unsigned char bit){
    long size = m*n;
    big* PK = convert(pk_h,size);
    big* G = convert(G_h,size);
    big g = bighToBig(q_h);

    big* R = (big*)malloc(sizeof(big)*size);
    fillRandom(R,bits,size);

    big* pk_d,*R_d,*G_d,*result_d;
    cudaMalloc((void **)&pk_d,sizeof(big)*size);
    cudaMalloc((void **)&R_d,sizeof(big)*size);
    cudaMalloc((void **)&G_d,sizeof(big)*size);
    cudaMalloc((void **)&result_d,sizeof(big)*size);

    cudaMemcpy(pk_d,PK,sizeof(u128)*size,cudaMemcpyHostToDevice);
    cudaMemcpy(R_d,R,sizeof(u128)*size,cudaMemcpyHostToDevice);
    cudaMemcpy(G_d,G,sizeof(u128)*size,cudaMemcpyHostToDevice);

    free(PK);
    free(G);

    encryptHelper(pk_d,R_d,G_d,result_d,g,bits,bit,n,m);

    cudaMemcpy(R,result_d,sizeof(u128)*size,cudaMemcpyDeviceToHost);
    bigH* cipherText = convertBack(R,size);


    cudaFree(pk_d);
    cudaFree(R_d);
    cudaFree(G_d);
    cudaFree(result_d);
    free(R);
    return cipherText;
}

void encryptBatch(bigH* pk_h, bigH* G_h,bigH q_h,uint n,uint m,uint bits,unsigned char* bit,bigH** cipherTexts,int len){
    long size = m*n;
    big* PK = convert(pk_h,size);
    big* G = convert(G_h,size);
    big g = bighToBig(q_h);

    big* R = (big*)malloc(sizeof(big)*size);

    big* pk_d,*R_d,*G_d,*result_d;
    cudaMalloc((void **)&pk_d,sizeof(big)*size);
    cudaMalloc((void **)&R_d,sizeof(big)*size);
    cudaMalloc((void **)&G_d,sizeof(big)*size);
    cudaMalloc((void **)&result_d,sizeof(big)*size);

    cudaMemcpy(pk_d,PK,sizeof(u128)*size,cudaMemcpyHostToDevice);
    cudaMemcpy(G_d,G,sizeof(u128)*size,cudaMemcpyHostToDevice);

    free(PK);
    free(G);

    for(int i = 0;i<len;i++){
        fillRandom(R,bits,size);
        cudaMemcpy(R_d,R,sizeof(u128)*size,cudaMemcpyHostToDevice);
        encryptHelper(pk_d,R_d,G_d,result_d,g,bits,bit[i],n,m);
        cudaMemcpy(R,result_d,sizeof(u128)*size,cudaMemcpyDeviceToHost);
        cipherTexts[i] = convertBack(R,size);
    }

    cudaFree(pk_d);
    cudaFree(R_d);
    cudaFree(G_d);
    cudaFree(result_d);
    free(R);
}