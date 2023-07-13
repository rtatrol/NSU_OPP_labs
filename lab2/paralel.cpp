#include "mpi.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>

using namespace std;

typedef long long ll;

vector<double> Mul(vector<double> x, double *A, ll N, ll n)
{
    vector<double> b(N);
    for (ll i = 0; i < N; i++)
    {
        for (ll j = 0; j < n; j++)
            b[i] += A[i * n + j] * x[j];
    }
    return b;
}

vector<double> sub(vector<double> a, vector<double> b, ll n)
{
    vector<double> c(n);
    for (ll i = 0; i < n; i++)
        c[i] = a[i] - b[i];
    return c;
}

vector<double> Const_Mul(vector<double> a, double c, ll n)
{
    vector<double> x(n);
    for (ll i = 0; i < n; i++)
        x[i] = a[i] * c;
    return x;
}

double Scalar_Mul(vector<double> a, vector<double> b, ll n)
{
    double answ = 0;
    for (ll i = 0; i < n; i++)
        answ += a[i] * b[i];
    return answ;
}

double norm(vector<double> x, ll n)
{
    double res = 0;
    for (ll i = 0; i < n; i++)
        res += x[i] * x[i];
    return sqrt(res);
}

void fill_matrix(double *A, ll n)
{
    for (ll i = 0; i < n; i++)
        for (ll j = 0; j < n; j++)
        {
            if (i == j)
                A[i * n + j] = rand() + n;
            else
            {
                A[i * n + j] = rand();
                A[i * n + j] = A[j * n + i];
            }
        }
}

int main(int argc, char **argv)
{
    srand(time(0));
    bool not_end = true;

    ll n = 3000, Max_Iter = 3000;
    double eps = 0.0000000001, result_rr = 0;

    int rank, size;
    double start_time, end_time, elapsed_time;
    MPI_Init(&argc, &argv);
    start_time = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t N = n / size, start = N * rank;
    if (rank == size - 1)
        N += n % size;
    double *A1 = new double[n * N];
    vector<double> x(N), x1(n), b(n), r(N), r1(n), z(N), z1(n);
    double *A = NULL;

    if (rank == 0)
    {
        A = new double[n * n];
        fill_matrix(A, n);
        for (ll i = 0; i < n; i++)
            b[i] = rand();
    }

    int send_counts[size], displs[size];
    for (ll i = 0; i < size - 1; i++)
    {
        send_counts[i] = (n / size) * n;
        displs[i] = i * (n / size) * n;
    }
    send_counts[size - 1] = (n / size + n % size) * n;
    displs[size - 1] = (size - 1) * (n / size) * n;

    MPI_Scatterv(A, send_counts, displs, MPI_DOUBLE, A1, send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&b[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (ll i = start; i < start + N; i++)
        z1[i - start] = r[i - start] = b[i]; 
    z = r;
    double B_norm = norm(b, n);

    for (ll j = 0; j < size - 1; j++)
    {
        send_counts[j] = n / size;
        displs[j] = j * (n / size);
    }
    send_counts[size - 1] = n / size + (n % size);
    displs[size - 1] = (size - 1) * (n / size);

    for (ll i = 0; i < Max_Iter and not_end; i++)
    {

        MPI_Allgatherv(&z[0], send_counts[rank], MPI_DOUBLE, &z1[0], send_counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        vector<double> buf(N);
        buf = Mul(z1, A1, N, n);

        double Azz = 0, rr = 0;
        double rr_buf = Scalar_Mul(r, r, N);
        double Azz_buf = Scalar_Mul(buf, z, N);
        MPI_Allreduce(&rr_buf, &rr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Azz_buf, &Azz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double alph = rr / Azz;

        x = sub(x, Const_Mul(z, -1 * alph, N), N);

        vector<double> r_next(N);
        r_next = sub(r, Const_Mul(buf, alph, N), N);

        double rr_next = 0;
        double rr_next_buf = Scalar_Mul(r_next, r_next, N);
        MPI_Allreduce(&rr_next_buf, &rr_next, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double bet = rr_next / rr;

        r = r_next;
        z = sub(r, Const_Mul(z, -1 * bet, N), N);

        result_rr = rr;
        not_end *= (sqrt(rr) / B_norm >= eps);
    }
    delete A1;

    MPI_Allgatherv(&x[0], send_counts[rank], MPI_DOUBLE, &x1[0], send_counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    elapsed_time = end_time - start_time;
    MPI_Finalize();

    if (rank == 0)
    {

        ofstream fout;
        fout.open("out_paralel_new.txt");
        fout << "Time paralel on " << size << " Process is  " << elapsed_time << "\n";
        fout << "Norm of B " << B_norm << "\n";
        fout << "Eps diff " << sqrt(result_rr) / B_norm << endl;
        if (sqrt(result_rr) / B_norm < eps)
            fout << "YES it work\n";
        else
            fout << "NO cant find answer\n";

        // for(int i=0;i<100;i++)
        // fout<<b[i]-x1[i]<<' ';

        delete A;
        fout.close();
    }
    return 0;
}