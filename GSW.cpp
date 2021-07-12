#include "includes/GSW.h"


parameters* setup(int kappa,int L){
    float sigma = 3.6;float sigma6 = 6*sigma;
    ZZ lower_bound;
    unsigned long n = (kappa+110)/7.2;
    ZZ quotient(4);
    unsigned long l = floor(log(quotient)/log(2)) + 1;
    unsigned long N = (n + 1) * l;
    while (true) {
        power(lower_bound, N+1,L);
        lower_bound *= 8 * sigma6;
        if (quotient <= lower_bound) {
            NextPrime(quotient, lower_bound);
        } else {
            break;
        }
        n = log(quotient/ceil(sigma))*(kappa + 110)/(7.2*log(2));
        l = floor(log(quotient)/log(2)) + 1;
        N = (n + 1) * l;
    }
    long t = NumBits(quotient);
    if(t<128) {
        bigH q = 0;
        bigH temp;
        for(long i = 0;i<l;i++){
            temp = bit(quotient,i);
            q |= temp<<i;
        }
        parameters *p = new parameters(n,N,q,sigma,t);
        cout <<"Size of Modulus: " << l << endl;
        cout << "Dimension: "<< n << endl;
        return p;
    }
    else return setup(64,3);
}

void writeMatrix(matrix& c,string s){
    ofstream CT(s+ ".txt");
    CT << c.rows << " " << c.cols << endl;
    for(long i = 0;i<c.rows*c.cols;i++){
        CT << c.vec[i] << endl;
    }
    CT.close();
}

void readMatrix(matrix& c,string s){
    fstream CT(s, std::ios_base::in);
    long b;
    ulong n,m;
    CT >> n;
    CT >> m;
    c.reshape(n,m);
    long i = 0;
    while (CT >> b)
    {
        c.vec[i] = b;
        i++;
    }
    CT.close();
}

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
    ~rng();
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
        if(err < 0) t = (m.q - (-1*err));
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

    matrix w(1,1,params.q);
    bigH t = params.q;
    if(t&1 == 1) t = (t>>1) + 1;
    else t >>= 1;
    w.vec[0] = t;

    W = inverseG(w.vec,w.rows,w.cols,params.bits);

    z = 1;
}

void GSW::keygen(){
    assert(z==1);
    sk.reshape(1,params.n + 1);
    pk.reshape(params.n + 1,params.m);
    sk.q = params.q;
    pk.q = params.q;
    
    fillRand(sk);
    fillRand(pk);
    
    sk.set(0,sk.cols-1,0);
    matrix B(1,params.m,params.q);

    matmul(B,sk,pk);
    gaussian(B,params.b);

    long last = pk.cols*(pk.rows-1);
    for(long i = 0;i<pk.cols;i++)
        pk.vec[last + i] = B.vec[i];
    
    sk.set(0,sk.cols-1,sk.q - 1);
    z = 2;
}

void GSW::encryptBit(int t,matrix& m){
    assert(t==0 || t==1);
    assert(z==2);
    matrix R(params.n+1,params.m,params.q);
    fillRand(R);
    free(m.vec);
    m.vec = encrypt(pk.vec,R.vec,G.vec,params.q,pk.rows,pk.cols,params.bits,t);
    m.q =params.q;
    m.rows = params.n + 1;
    m.cols = params.m;
}

int GSW::decryptBit(matrix& C){
    assert(z==2);
    bigH t = params.q;
    if(t&1 == 1) t = (t>>1) + 1;
    else t >>= 1;
    bigH sum,q;
    q = params.q;
    // bigH* temp = sparse_mul(C.vec,W,C.rows,C.cols,C.q);
    matrix ret(params.n+1,1,q);
    bigH* ct = C.vec;
    long m = params.m;

    for(long i = 0;i<C.rows;i++){
        sum = 0;
    for(long k = 0;k<params.bits;k++){
        if(W[k] == 1){
            sum += ct[(i+1)*m - params.bits + k];
            if(!(sum < q)) sum -= q;
        }
    }
        ret.vec[i] = sum;    
    }

    matrix res;
    res.q = params.q;

    matmul(res,sk,ret);
    bigH bit = res.vec[0];

    bigH zero = min(bit, params.q-bit);
    bigH one;
    if(bit > t) one = bit - t;
    else one = t - bit;
    return (one < zero);

}

void GSW::saveState(){
    cout << "warning!!, this will overwrite existing files (if any), want to proceed?? (0/1): ";
    int i;
    cin >> i;
    if(i==1){
    cout << "writing secret key...\n";
    writeMatrix(sk,"secretKey");
    cout << "writing public key...\n";
    writeMatrix(pk,"publicKey");
    cout << "writing parameters...\n";
    ofstream CT("params.txt");
    CT << params.n << "\n" << params.m << "\n" << params.q << "\n" << params.b << "\n"<<params.bits;
    }
}