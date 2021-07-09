#include "includes/GSW.h"

void genGadget(long bits,matrix &G){
    bigH* b = G.vec;
    int m = G.rows,n = G.cols;
    b[0] = 1;
    for(long i = 1;i<bits;i++)
        b[i] = b[i-1] << 1;
    for(long i = 1;i<G.rows;i++)
        for(long j = 0;j<bits;j++)
            b[n*i + i*bits + j] = b[j];
}

void fillRand(matrix &mat){
    duthomhas::csprng rng;
    bigH mod = mat.q;
    bigH* v = mat.vec;
    for(long i = 0;i<mat.rows*mat.cols;i++){
        v[i] = rng(bigH());
        if(!(v[i] < mod)) v[i] %= mod;
    }
}

void gaussian(matrix &m,double b){ // the distribution is uniform, change it to gaussian
    duthomhas::csprng rng;
    normal_distribution<double> distribution(0.0,b);
    bigH t;
    int err;
    bigH* v = m.vec;
    for(long i = 0;i<m.cols*m.rows;i++){
        //t = rng(int())%m.q; //edit this to gaussian distribution
        err = int(round(distribution(rng)));
        if(err < 0) t = (m.q + err);
        else t = err;
        v[i] += t;
        if(!(v[i] < m.q)) v[i] -= m.q;
    }
    ~rng();
}

GSW::GSW(parameters *p){
    params = *p;
    G.reshape(params.n + 1,params.bits*(params.n + 1));
    G.q = p->q;
    genGadget(params.bits,G);
}

void GSW::keygen(){

    sk.reshape(1,params.n + 1);
    pk.reshape(params.n + 1,params.m);
    sk.q = params.q;
    pk.q = params.q;
    
    fillRand(sk);
    fillRand(pk);
    
    sk.set(1,sk.cols-1,0);
    matrix B(1,params.m,params.q);

    matmul(B,sk,pk);
    gaussian(B,params.b);

}