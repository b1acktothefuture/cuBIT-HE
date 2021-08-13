#include "../includes/matrix.h"

matrix::matrix(long n, long m, bigH mod)
{
    vec = (bigH *)malloc(sizeof(bigH) * m * n);
    rows = n;
    cols = m;
    q = mod;
}

bigH matrix::get(long i, long j)
{
    assert(i >= 0 && i < rows && j >= 0 && j < cols);
    return vec[cols * i + j];
}

void matrix::set(long i, long j, bigH m)
{
    assert(i >= 0 && i < rows && j >= 0 && j < cols);
    if (m > q)
        m = m % q;
    vec[cols * i + j] = m;
}

void matrix::reshape(long m, long n)
{
    free(vec);
    vec = (bigH *)malloc(sizeof(bigH) * m * n);
    for (long i = 0; i < m * n; i++)
        vec[i] = 0;
    rows = m;
    cols = n;
}

void matmul(matrix &x, const matrix &m, const matrix &n)
{
    assert(m.cols == n.rows);
    bigH *one = m.vec, *two = n.vec, *res = x.vec;
    long i, j, k;
    long M = m.rows, N = n.cols, K = n.cols;
    // #pragma omp parallel shared(one,two,res) private(i,j,k)

    biggerH holder, h1, h2;
    bigH sum, q;
    q = m.q;

    // #pragma omp for  schedule(static)
    for (i = 0; i < m.rows; i++)
    {
        for (j = 0; j < n.cols; j++)
        {
            // cout << (float((j+1)*100)/float(n.cols)) <<" %" << endl;
            sum = 0;
            for (k = 0; k < m.cols; k++)
            {
                // 256 bit arithmetic
                h1 = biggerH(one[N * i + k]);
                h2 = biggerH(two[k * K + j]);
                holder = (h1 * h2);
                if (!(holder < q))
                    holder %= q;

                sum += holder.lower();
                if (!(sum < q))
                    sum -= q;
            }
            res[K * i + j] = sum;
        }
    }
}

bigH *sparse_mul(bigH *A, unsigned char *B, int m, int n, bigH q)
{ // m*n & n*n
    bigH *result = (bigH *)malloc(sizeof(bigH) * m * n);
    bigH sum;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
        {
            sum = 0;
            for (int k = 0; k < n; k++)
            {
                if (B[k * n + j])
                {
                    sum += A[i * n + k];
                    if (sum >= q)
                        sum -= q;
                }
            }
            result[n * i + j] = sum;
        }
    return result;
}

unsigned char *inverseG(bigH *matrix, int n, int m, int bits)
{
    unsigned char *ret = (unsigned char *)malloc(sizeof(unsigned char) * n * bits * m);
    bigH temp;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            temp = matrix[j * m + i];
            // cout << temp << endl;
            for (int k = 0; k < bits; k++)
            {
                ret[(j * bits + k) * m + i] = temp & 1;
                temp >>= 1;
            }
        }
    }
    return ret;
}

void add_cpu(bigH *a, bigH *b, long size, bigH q)
{
    for (int i = 0; i < size; i++)
    {
        a[i] += b[i];
        if (!(a[i] < q))
            a[i] -= q;
    }
}

void print_martix(bigH *b, int rows, int cols)
{
    unsigned int u = UINT32_MAX;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
            std::cout << (b[i * cols + j]) << " ";
        std::cout << std::endl;
    }
    cout << endl;
}

void print_martix(unsigned char *b, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
            std::cout << int(b[i * cols + j]) << " ";
        std::cout << std::endl;
    }
    cout << endl;
}

void print_martix(matrix &m)
{
    cout << m.q << endl;
    print_martix(m.vec, m.rows, m.cols);
}

void writeMatrix(matrix &c, string s)
{
    ofstream CT(s + ".txt");
    CT << c.rows << " " << c.cols << endl;
    for (long i = 0; i < c.rows * c.cols; i++)
    {
        CT << c.vec[i].upper() << " " << c.vec[i].lower() << endl;
    }
    CT.close();
}

void readMatrix(matrix &c, string s)
{
    fstream CT(s, std::ios_base::in);
    bigH b;
    ulong n, m;
    CT >> n;
    CT >> m;
    c.reshape(n, m);
    uint64_t t;
    long i = 0;
    while (CT >> t)
    {
        b = t;
        b <<= 64;
        CT >> t;
        b += t;
        c.vec[i] = b;
        i++;
    }
    CT.close();
}