#ifndef UTILS_H
#define UTILS_H

#include "matrix.h"
#include "csprng.hpp"
#include <fstream>

struct parameters
{
    bigH q;
    double b;
    unsigned long n, m, bits;
    parameters(unsigned long dimension = 0, unsigned long samples = 0, bigH quotient = 1, double bound = 0, long len = 0)
    {
        n = dimension;
        m = samples;
        q = quotient;
        b = bound;
        bits = len;
    }
};

parameters *setup(int kappa, int L);

class GSW
{
protected:
    int z = 0;

public:
    parameters params;
    matrix sk, pk, G;
    unsigned char *W;

    GSW(parameters *p = setup(10, 1));
    GSW(string s, string p, string par);

    void keygen();
    void saveState();

    void encryptBit(int t, matrix &m, int i = 0);
    int decryptBit(matrix &C); // check if decryption is closer to q/2 or {0,q};

    matrix *encryptBits(unsigned char *C, int len);
    unsigned char *decryptBits(matrix *C, int len);
};

void genGadget(long n, matrix &G);
void fillRand(matrix &mat);
void gaussian(matrix &m, double b);

#endif