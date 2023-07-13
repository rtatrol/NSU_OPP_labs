#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <time.h>
#include <omp.h>

using namespace std;

typedef long long ll;

double Mul(vector<double> x, double *A, ll n) {
    double b=0;
    for (ll i = 0; i < n; i++) {
            b+= A[i] * x[i];
    }
    return b;
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

double res1=0,res2=0,res3=0;
double alph=0,bet=0;

int main(int argc, char **argv) 
{
    bool not_end = true;
    srand(time(0));
    double start,end;

    ll n = 3000, Max_Iter = 3000;
    double eps = 0.0000000001;

    vector<double> x(n), b(n), r(n), z(n);

    double *A = new double[n * n];
    fill_matrix(A, n);

    for (ll i = 0; i < n; i++)b[i] = rand();

    double B_norm = norm(b, n);

    r = b;
    z = r;

    vector<double> buf(n);

    #pragma omp parallel num_threads(4)
    {
        #pragma omp single
            start = omp_get_wtime();

        for (ll i = 0; i < Max_Iter && not_end; i++)
        {             
            #pragma omp for schedule(static)
            for(ll j=0;j<n;j++)
                buf[j] = Mul(z, &A[j*n], n);

            #pragma omp single
                res1=0;
            #pragma omp for reduction(+:res1) schedule(static)
            for(ll j=0;j<n;j++)
                res1 += r[j] * r[j];
                
            #pragma omp single
                res2=0;
            #pragma omp for reduction(+:res2) schedule(static)
            for(ll j=0;j<n;j++)
                res2+=buf[j]*z[j];

            #pragma omp single
            alph= res1 / res2;

            #pragma omp for schedule(static)
            for(ll j=0;j<n;j++)
                x[j]+=alph*z[j];

            #pragma omp for schedule(static)
            for(ll j=0;j<n;j++)
                r[j]-=alph*buf[j];

            #pragma omp single 
                res3=0;
            #pragma omp for reduction(+:res3) schedule(static)
            for(ll j=0;j<n;j++)
                res3+=r[j]*r[j];

            #pragma omp single
            bet = res3 / res1;

            #pragma omp for schedule(static)
            for(ll j=0;j<n;j++)
                z[j]=r[j]+bet*z[j];

            #pragma omp single
                not_end *= (sqrt(res3) / B_norm >= eps);
        }

        #pragma omp single
        {
            end = omp_get_wtime();
            ofstream fout;
            fout.open("out_default.txt");
            fout<<x[0]<<' '<<x[100]<<' '<<x[2999]<<endl;
            fout <<"Time on "<<omp_get_num_threads()<<" Threads is "<< end - start << "\n";
            fout << "Norm of B " << norm(b, n) << "\n";
            fout << "Eps diff " << norm(r, n) / B_norm<<"\n";
            fout.close();
        }

    }
    return 0;
}