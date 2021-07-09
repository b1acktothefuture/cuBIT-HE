#include "includes/GPU_wrapper.h"
#include "includes/kernel.cuh"

void print_martix(big* b,int rows,int cols){
    for(int i = 0;i<rows;i++){
    for(int j = 0;j<cols;j++)
        std::cout << b[i*cols + j].w << " ";
    std::cout << std::endl;
    }
    cout << endl;
}

big bighToBig(bigH g){
    uint AND = UINT32_MAX;
    big q;
    q.x = g&AND;
    g >>= 32;
    q.y = g&AND;
    g >>= 32;
    q.z = g&AND;
    g >>= 32;
    q.w = g&AND;
    return q;
}

bigH bigTobigH(big q){
    bigH g(0);
    g += q.x;
    g <<= 32;
    g += q.y;
    g <<= 32;
    g += q.z;
    g <<= 32;
    g += q.w;
    return g;
}

big* convert(bigH* matrix,int size){
    big* ret = (big *)malloc(sizeof(big)*size);
    for(int i =0;i<size;i++){
        ret[i] = bighToBig(matrix[i]);
    }
    return ret;
}

bigH* convertBack(big* matrix,int size){
    bigH* ret = (bigH *)malloc(sizeof(bigH)*size);
    for(int i =0;i<size;i++){
        ret[i] = bigTobigH(matrix[i]);
    }
    return ret;
}


void encryptH(big* A,big* R,big* G,big* result,big q,uint bits,unsigned char bit,int n,int m){

    unsigned int grid_rows = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int grid_cols = (m + BLOCK_SIZE - 1) / BLOCK_SIZE;
    dim3 dimGrid(grid_cols, grid_rows);
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

    gpu_sparse_mult<<<dimGrid,dimBlock>>>(A,R,result,n,m,bits,q);

    if(bit == 1){
        int blockSize = 256;
        int numBlocks = (m*n + blockSize - 1) / blockSize;
        gpu_add<<<numBlocks, blockSize>>>(n*m,q, G, result);
    }
}




/**************************************************************************/


void test(big* A,big* R,big* result,big q,uint bits,int n,int m){
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


/******************************************************************************/


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