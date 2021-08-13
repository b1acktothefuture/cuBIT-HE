#ifndef MATRIX_H
#define MATRIX_H

#include "GPU_wrapper.h"
#include <cassert>
#include <omp.h>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

class matrix;

class matrix
{
public:
    bigH *vec;
    long rows, cols;
    bigH q;

public:
    matrix(long n = 1, long m = 1, bigH mod = bigH(1)); // bigH mod = 1
    matrix(long n, long m, bigH *v, bigH mod)
    {
        vec = v;
        rows = n;
        cols = m;
        q = mod;
    }
    bigH get(long i, long j);
    void set(long i, long j, bigH m);

    void reshape(long m, long n); // erases data and fills with zeros

    ~matrix() { free(vec); };
};

void matmul(matrix &x, const matrix &m, const matrix &n);
bigH *sparse_mul(bigH *A, unsigned char *B, int m, int n, bigH q);
void print_martix(bigH *b, int rows, int cols);
void print_martix(matrix &m);
void print_martix(unsigned char *b, int rows, int cols);
unsigned char *inverseG(bigH *matrix, int n, int m, int bits);
void add_cpu(bigH *a, bigH *B, long size, bigH q);
void writeMatrix(matrix &c, string s);
void readMatrix(matrix &m, string s);

// void add(matrix& x,const matrix& a,const matrix& c);

#endif