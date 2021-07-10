#ifndef UTILS_H
#define UTILS_H

#include "matrix.h"
#include "csprng.hpp"


void genGadget(long n,matrix &G);
void fillRand(matrix &mat);
void gaussian(matrix &m,double b);

struct parameters{
    bigH q;
    double b;
    unsigned long n,m,bits;
    parameters(unsigned long dimension = 0,unsigned long samples = 0,bigH quotient = 1,double bound = 0,long len = 0){
        n = dimension;
        m = samples;
        q = quotient;
        b = bound;
        bits = len;
    }
};

class GSW{
    int z = 0;
    public:
    parameters params;
    matrix sk,pk,G;

    GSW(parameters* p);

    void keygen();

    void encryptBit(int t,matrix& m);
    int decryptBit(matrix& C); // check if decryption is closer to q/2 or {0,q};

};

#endif