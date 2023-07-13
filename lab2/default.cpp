#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <time.h>

using namespace std;

typedef long long ll;

vector<double> Mul(vector<double> x, double *A, ll n) {
    vector<double> b(n);
    for (ll i = 0; i < n; i++) {
        for (ll j = 0; j < n; j++)
            b[i] += A[i * n + j] * x[j];
    }
    return b;
}

vector<double> sub(vector<double> a, vector<double> b, ll n) {
    vector<double> c(n);
    for (ll i = 0; i < n; i++)
        c[i] = a[i] - b[i];
    return c;
}

vector<double> Const_Mul(vector<double> a, double c, ll n) {
    vector<double> x(n);
    for (ll i = 0; i < n; i++)
        x[i] = a[i] * c;
    return x;
}

double Scalar_Mul(vector<double> a, vector<double> b, ll n) {
    double answ = 0;
    for (ll i = 0; i < n; i++)
        answ += a[i] * b[i];
    return answ;
}

double norm(vector<double> x, ll n) {
    double res = 0;
    for (ll i = 0; i < n; i++)
        res += x[i] * x[i];
    return sqrt(res);
}

void fill_matrix(double *A, ll n) {
    for (ll i = 0; i < n; i++)
        for (ll j = 0; j < n; j++) {
            if (i == j)A[i * n + j] = rand() + n;
            else {
                A[i * n + j] = rand();
                A[i * n + j] = A[j * n + i];
            }
        }
}

int main(int argc, char **argv) {
    bool not_end = true;
    srand(time(0));
    clock_t start = clock();

    ll n = 3000, Max_Iter = 3000;
    double eps = 0.0000000001;

    vector<double> x(n), b(n), r(n), z(n);

    double *A = new double[n * n];
    fill_matrix(A, n);

    for (ll i = 0; i < n; i++)b[i] = rand();

    double B_norm = norm(b, n);

    r = sub(b, Mul(x, A, n), n);
    z = r;
    for (ll i = 0; i < Max_Iter && not_end; i++) {
        vector<double> buf(n);
        buf = Mul(z, A, n);

        double alph = Scalar_Mul(r, r, n) / Scalar_Mul(buf, z, n);

        x = sub(x, Const_Mul(z, -1 * alph, n), n);
        vector<double> r_next(n);

        r_next = sub(r, Const_Mul(buf, alph, n), n);
        double bet = Scalar_Mul(r_next, r_next, n) / Scalar_Mul(r, r, n);

        r = r_next;
        z = sub(r, Const_Mul(z, -1 * bet, n), n);
        not_end *= (norm(r, n) / B_norm >= eps);
    }

    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    ofstream fout;

    fout.open("out_default.txt");
    fout<<seconds<<"\n";
    fout << "Norm of B " << norm(b, n) << "\n";
    fout << "Eps diff " << norm(r, n) / B_norm<<"\n";
    fout.close();

    return 0;
}