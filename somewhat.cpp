#include "includes/somewhat.h"

void getGadget(bigH u, long bits, matrix &G)
{
    bigH *b = G.vec;
    int m = G.rows, n = G.cols;
    for (long i = 0; i < m * n; i++)
        b[i] = 0;
    b[0] = u;
    for (long i = 1; i < bits; i++)
    {
        b[i] = (b[i - 1] << 1);
        if (b[i] > G.q)
            b[i] -= G.q;
    }
    for (long i = 1; i < G.rows; i++)
        for (long j = 0; j < bits; j++)
            b[n * i + i * bits + j] = b[j];
}

parameters *setupSW(int kappa)
{
    ZZ quotient = GenGermainPrime_ZZ(kappa, 80);
    bigH q = 0, temp;
    long l = floor(log(quotient) / log(2)) + 1;
    long t = NumBits(quotient);
    for (long i = 0; i < l; i++)
    {
        temp = bit(quotient, i);
        q |= temp << i;
    }
    long n = kappa;
    long N = (n + 1) * l;
    parameters *p = new parameters(n, N, q, 1.0, t);
    return p;
}

void somewhatGSW::encryptSW(bigH u, matrix &m)
{
    assert(z == 2);
    unsigned char *R = new unsigned char[params.m * params.m];
    for (int i = 0; i < params.m * params.m; i++)
        R[i] = rand() % 2;
    bigH *vec = sparse_mul(pk.vec, R, params.n + 1, params.m, params.q);
    matrix A(params.n + 1, params.m, vec, params.q);
    // matrix b(sk.rows, pk.cols, params.q);
    m.q = params.q;
    m.reshape(G.rows, G.cols);
    getGadget(u, params.bits, m);

    for (long i = 0; i < G.rows * G.cols; i++)
    {
        m.vec[i] = m.vec[i] + A.vec[i];
        if (m.vec[i] > m.q)
            m.vec[i] = m.vec[i] - m.q;
    }
}

bigH somewhatGSW::decryptSW(matrix &C)
{
    assert(z == 2);
    matrix msg(sk.rows, pk.cols, params.q);
    matmul(msg, sk, C);
    matrix sg(sk.rows, pk.cols, params.q);
    matmul(sg, sk, G);
    uint64_t AND = 0xFFFFFFFFFFFFFFFF;
    unordered_map<uint64_t, long> modes; // make this work with uint128_t
    for (long i = 0; i < sk.rows * pk.cols; i++)
    {
        bigH t = msg.vec[i] / sg.vec[i], r = sg.vec[i] / 2, k = msg.vec[i] % sg.vec[i];
        if (k > r)
            t = t + 1;
        uint64_t key = t & AND;
        if (modes.find(key) == modes.end())
            modes[key] = 1;
        else
            modes[key]++;
    }
    uint64_t best_num = 0;
    // infinity
    bigH best_dist = AND;
    best_dist <<= 64;
    best_dist |= AND;

    for (auto i : modes)
    {
        bigH dist = 0;
        bigH container;
        for (long j = 0; j < sg.rows * sg.cols; j++)
        {
            container = sg.vec[j];
            container *= i.first;
            container %= params.q;
            container = (msg.vec[j] + (params.q - container)) % params.q;
            if (2 * container > params.q)
                container = params.q - container;
            container = container * container;
            dist += container;
        }
        if (dist < best_dist)
        {
            best_num = i.first;
            best_dist = dist;
        }
    }
    //cout << best_num << endl;
    return best_num;
}

void somewhatGSW::add(matrix &a, matrix &b, matrix &r)
{
    assert(a.rows == b.rows && a.cols == b.cols);
    r.reshape(a.rows, a.cols);
    r.q = a.q;
    for (long j = 0; j < a.rows * a.cols; j++)
    {
        r.vec[j] = a.vec[j] + b.vec[j];
        if (!(r.vec[j] < params.q))
            r.vec[j] -= params.q;
    }
}

void somewhatGSW::mul(matrix &C1, bigH c)
{
    biggerH container;
    for (long i = 0; i < C1.rows * C1.cols; i++)
    {
        container = biggerH(C1.vec[i]) * biggerH(c);
        if (!(container < C1.q))
            container %= C1.q;
        C1.vec[i] = container.lower();
    }
}