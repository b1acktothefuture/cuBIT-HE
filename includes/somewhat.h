#include "GSW.h"

class somewhatGSW : public GSW
{
public:
    somewhatGSW(parameters *p) : GSW(p){};

    void encryptSW(bigH u, matrix &m);
    bigH decryptSW(matrix &m);

    void add(matrix &c1, matrix &c2, matrix &res);
    void mul(matrix &c1, bigH c);
};

parameters *setupSW(int kappa);