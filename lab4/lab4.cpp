#include "mpi.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>

using namespace std;
typedef long long ll;

ll to_value(char *buf)
{
    ll res = atoi(buf);
    return res;
}

void fill_matrix(int *M, int y, int x)
{
    for (ll i = 0; i < y; i++)
        for (ll l = 0; l < x; l++)
            M[i * x + l] = rand()%10000;
}

void mul(int *a, int *b, int *c, int k1, int n2, int k2)
{
    for (ll i = 0; i < k1; i++)
    {
        for (ll j = 0; j < n2; j++)
        {
            for (ll k = 0; k < k2; k++)
                c[i * k2 + k] += a[i * n2 + j] * b[j * k2 + k];
        }
    }
}
void pull(int *C, int *c, int n3, int k1, int k2, int P1_cord, int P2_cord)
{
    for (ll i = P1_cord * k1; i < (P1_cord + 1) * k1; i++)
        for (ll j = P2_cord * k2; j < (P2_cord + 1) * k2; j++)
            C[i * n3 + j] = c[(i - P1_cord * k1) * k2 + j - P2_cord * k2];
}

int main(int argc, char **argv)
{
    int rank, size;
    int n1, n2, n3, p1, p2, k1, k2;
    n2 = 2500;
    k1 = k2 = 500;

    p1 = to_value(argv[1]);
    n1 = k1 * p1;

    p2 = to_value(argv[2]);
    n3 = k2 * p2;

    double start_time, elapsed_time;

    MPI_Init(&argc, &argv);
    start_time = MPI_Wtime();
    int dims[2] = {p1, p2};
    int periods[2] = {0, 0};

    MPI_Comm cart_comm, row_comm, col_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    MPI_Comm_rank(cart_comm, &rank);
    int mpi_cords[2];
    MPI_Cart_coords(cart_comm, rank, 2, mpi_cords);

    int *a = new int[k1 * n2];
    int *b = new int[n2 * k2];
    int *c = new int[k1 * k2];
    for (ll i = 0; i < k1 * k2; i++)
        c[i] = 0;

    MPI_Datatype col, coltype;
    MPI_Type_vector(n2, k2, n3, MPI_INT, &col);
    MPI_Type_commit(&col);
    MPI_Type_create_resized(col, 0, k2 * sizeof(int), &coltype);
    MPI_Type_commit(&coltype);

    MPI_Comm_split(cart_comm, mpi_cords[0], rank, &row_comm);
    MPI_Comm_split(cart_comm, mpi_cords[1], rank, &col_comm);

    int *A = NULL;
    int *B = NULL;
    int *C = NULL;

    if (rank == 0)
    {
        A = new int[n1 * n2];
        fill_matrix(A, n1, n2);
        B = new int[n2 * n3];
        fill_matrix(B, n2, n3);
        C = new int[n1 * n3];
    }

    if (mpi_cords[1] == 0)
        MPI_Scatter(A, k1 * n2, MPI_INT, a, k1 * n2, MPI_INT, 0, col_comm);
    if (mpi_cords[0] == 0)
        MPI_Scatter(B, 1, coltype, b, n2 * k2, MPI_INT, 0, row_comm);

    MPI_Bcast(a, k1 * n2, MPI_INT, 0, row_comm);

    MPI_Bcast(b, n2 * k2, MPI_INT, 0, col_comm);

    mul(a, b, c, k1, n2, k2);

    if (rank != 0)
    {
        MPI_Send(c, k1 * k2, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    else
    {
        pull(C, c, n3, k1, k2, 0, 0);
        for (int i = 1; i < p1 * p2; i++)
        {
            MPI_Recv(c, k1 * k2, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            pull(C, c, n3, k1, k2, i / p2, i % p2);
        }
    }

    elapsed_time = MPI_Wtime() - start_time;
    MPI_Comm_free(&cart_comm);
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    MPI_Type_free(&col);
    MPI_Type_free(&coltype);

    MPI_Finalize();

    delete a;
    delete b;
    delete c;

    if (rank == 0)
    {
        ofstream fout;
        fout.open("out_program.txt");
        fout << "Time paralel on  " << p1 << '*' << p2 << "  Process is  " << elapsed_time << "\n";
        fout<< "Matrix size is: A "<<n1<<'*'<<n2<<"    "<<"B: "<<n2<<'*'<<n3<<"\n";

        fout<<C[0]<<' '<<C[10]<<endl;
        // for (ll i = 0; i < n1; i++)
        // {
        //     for (ll j = 0; j < n3; j++)
        //         fout << C[i * n3 + j] << ' ';
        //     fout << "\n";
        // }
        
        delete A;
        delete B;
        delete C;
        fout.close();
    }
    return 0;
}