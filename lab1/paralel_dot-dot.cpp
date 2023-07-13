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
    n+=size-n%size;
    ll x=n/size;

    if (rank == 0) {
        cout <<n<<' '<< x << endl;
        int *a = new int[n];
        int *b = new int[n];

        for (ll i = 0; i < atoi(argv[1]); ++i)
            a[i] =b[i] = 1;
        for(ll i=atoi(argv[1]);i<n;i++)
            a[i]=b[i]=0;

        for (ll i = 1; i < size; i++) {
            MPI_Send(a + i* x, x, MPI_INT, i, 123, MPI_COMM_WORLD);
            MPI_Send(b, n, MPI_INT, i, 999, MPI_COMM_WORLD);
        }

        int result = 0;
        for (ll i = 0; i < x; i++)
            for (ll l = 0; l < n; l++)
                result += a[i] * b[l];

        for (int buf = 0, i = 1; i < size; i++, result += buf) {
            MPI_Recv(&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        cout << result << endl;
        delete a;
        delete b;

    } else {
        int *a = new int[x];
        int *b = new int[n];
        MPI_Recv(a, x, MPI_INT, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, n, MPI_INT, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int res = 0;
        for (ll i = 0; i < x; i++)
            for (ll l = 0; l < n; l++)
                res += a[i] * b[l];
        MPI_Send(&res, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
        delete a;
        delete b;
    }
    MPI_Finalize();
}