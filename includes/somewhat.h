#ifndef SOMEWHAT_H
#define SOMEWHAT_H

#include "GSW.h"

class somewhatGSW : public GSW
{
public:
    somewhatGSW(parameters *p) : GSW(p){};

    void encryptSW(bigH u, matrix &m, int T = 0);
    bigH decryptSW(matrix &m);
    bigH decryptMP(matrix &m);

    void add(matrix &c1, matrix &c2, matrix &res);
    void mul(matrix &c1, bigH c);
};

parameters *setupSW(int kappa);

#endif