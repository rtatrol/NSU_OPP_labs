#include "mpi.h"
#include <iostream>
#include <vector>

using namespace std;
typedef long long ll;

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ll n = atoi(argv[1]);
    n += size - n % size;
    ll x = n / size;

    int *a = new int[n];
    int *b = new int[n];
    int *a1 = new int[x];

    if (rank == 0) {
        cout << n << ' ' << x << endl;
        for (ll i = 0; i < atoi(argv[1]); ++i)
            a[i] = b[i] = 1;
        for (ll i = atoi(argv[1]); i < n; i++)
            a[i] = b[i] = 0;
    }

    MPI_Scatter(a, x, MPI_INT, a1, x, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, n, MPI_INT, 0, MPI_COMM_WORLD);

    ll buf = 0, sum = 0;
    for (ll i = 0; i < x; i++)
        for (ll j = 0; j < n; j++)
            buf += a1[i] * b[j];
    MPI_Reduce(&buf, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
        cout << "Answer is " << sum << endl;

    MPI_Finalize();

}