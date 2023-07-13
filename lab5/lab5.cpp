#include "mpi.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>

using namespace std;
typedef long long ll;

int X, Y;
double global_start = -1, global_end = -1;
int diff = 0;

int count(int **field, int i, int j)
{
    ll answ = 0;
    answ += field[i - 1][j] + field[i + 1][j];
    if (j == 0)
    {
        answ += field[i - 1][Y - 1] + field[i][Y - 1] + field[i + 1][Y - 1];
        answ += field[i - 1][j + 1] + field[i][j + 1] + field[i + 1][j + 1];
    }
    else
    {
        if (j == Y - 1)
        {
            answ += field[i - 1][0] + field[i][0] + field[i + 1][0];
            answ += field[i - 1][j - 1] + field[i][j - 1] + field[i + 1][j - 1];
        }
        else
        {
            answ += field[i - 1][j + 1] + field[i][j + 1] + field[i + 1][j + 1];
            answ += field[i - 1][j - 1] + field[i][j - 1] + field[i + 1][j - 1];
        }
    }
    return answ;
}

void gameInside(int **before, int **next)
{
    for (ll i = 1; i < diff - 1; i++)
    {
        for (ll l = 0; l < Y; l++)
        {
            int near = count(before, i, l);
            if (before[i][l] && (near < 2 || near > 3))
                next[i][l] = 0;
            else
            {
                if (!before[i][l] && near == 3)
                    next[i][l] = 1;
                else
                    next[i][l] = before[i][l];
            }
        }
    }
}

void gameLimitsFirst(int **before, int **next, int *first_line)
{
    int near = 0;

    near += before[0][Y - 1] + before[1][Y - 1] + first_line[Y - 1] + before[0][1] + before[1][1] + first_line[1] + first_line[0] + before[1][0];
    if (before[0][0] && (near < 2 || near > 3))
        next[0][0] = 0;
    else
    {
        if (!before[0][0] && near == 3)
            next[0][0] = 1;
        else
            next[0][0] = before[0][0];
    }

    near = 0;
    near = first_line[Y - 1] + before[1][Y - 1] + before[0][Y - 2] + before[1][Y - 2] + first_line[Y - 2] + before[0][0] + before[1][0] + first_line[0];
    if (before[0][Y - 1] && (near < 2 || near > 3))
        next[0][Y - 1] = 0;
    else
    {
        if (!before[0][Y - 1] && near == 3)
            next[0][Y - 1] = 1;
        else
            next[0][Y - 1] = before[0][Y - 1];
    }

    for (ll j = 1; j < Y - 1; j++)
    {
        near = 0;
        near += first_line[j - 1] + first_line[j] + first_line[j + 1] + before[1][j - 1] + before[1][j] + before[1][j + 1] + before[0][j - 1] + before[0][j + 1];
        if (before[0][j] && (near < 2 || near > 3))
            next[0][j] = 0;
        else
        {
            if (!before[0][j] && near == 3)
                next[0][j] = 1;
            else
                next[0][j] = before[0][j];
        }
    }
}

void gameLimitsEnd(int **before, int **next, int *end_line)
{
    int near = 0;
    near += before[diff - 2][0] + end_line[0] + before[diff - 2][1] + before[diff - 1][1] + end_line[1] + before[diff - 1][Y - 1] + before[diff - 2][Y - 1] + end_line[Y - 1];
    if (before[diff - 1][0] && (near < 2 || near > 3))
        next[diff - 1][0] = 0;
    else
    {
        if (!before[diff - 1][0] && near == 3)
            next[diff - 1][0] = 1;
        else
            next[diff - 1][0] = before[diff - 1][0];
    }

    near = 0;
    near += before[diff - 2][Y - 1] + end_line[Y - 1] + before[diff - 2][Y - 2] + before[diff - 1][Y - 2] + end_line[Y - 2] + before[diff - 1][0] + before[diff - 2][0] + end_line[0];

    if (before[diff - 1][Y - 1] && (near < 2 || near > 3))
        next[diff - 1][Y - 1] = 0;
    else
    {
        if (!before[diff - 1][Y - 1] && near == 3)
            next[diff - 1][Y - 1] = 1;
        else
            next[diff - 1][Y - 1] = before[diff - 1][Y - 1];
    }

    for (ll j = 1; j < Y - 1; j++)
    {
        near = 0;
        near += end_line[j - 1] + end_line[j] + end_line[j + 1] + before[diff - 2][j - 1] + before[diff - 2][j] + before[diff - 2][j + 1] + before[diff - 1][j - 1] + before[diff - 1][j + 1];
        if (before[diff - 1][j] && (near < 2 || near > 3))
            next[diff - 1][j] = 0;
        else
        {
            if (!before[diff - 1][j] && near == 3)
                next[diff - 1][j] = 1;
            else
                next[diff - 1][j] = before[diff - 1][j];
        }
    }
}

void savePrevious(int ***previous, int iter, int **field)
{
    previous[iter] = new int *[diff];

    for (ll i = 0; i < diff; i++)
    {
        previous[iter][i] = new int[Y];

        for (ll j = 0; j < Y; j++)
            previous[iter][i][j] = field[i][j];
    }
}

void findFlags(int *flags, int ***previous, int **field, int iter)
{
    for (ll k = 0; k < iter; k++)
    {
        bool answ = 1;
        for (ll i = 0; i < diff && answ; i++)
        {
            for (ll l = 0; l < Y && answ; l++)
                answ *= previous[k][i][l] == field[i][l];
        }
        flags[k] = answ;
    }
}

int checkFlags(int **flagsMatrix, int iter, int size)
{

    for (ll i = 0; i < iter; i++)
    {
        bool isCorrect = 1;
        for (ll j = 0; j < size && isCorrect; j++)
            isCorrect *= flagsMatrix[j][i];

        if (isCorrect)
            return i;
    }
    return -1;
}

void rememberField(int **field, int **next)
{
    for (ll i = 0; i < diff; i++)
    {
        for (ll l = 0; l < Y; l++)
            field[i][l] = next[i][l];
    }
}

int main(int argc, char **argv)
{

    int flagPrev = -1, flagNow = -1;

    int rank, size, limit = 405;

    X = Y = 100;

    double start_time, elapsed_time;
    int *first_line = new int[Y];
    int *end_line = new int[Y];

    MPI_Init(&argc, &argv);
    start_time = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Request sendFirst, sendEnd, recvFirst, recvEnd;

    global_end = ceil(X * (rank + 1) / size);
    global_start = ceil(X * rank / size) + 1;
    diff = global_end - global_start + 1;

    int **field = new int *[diff];
    int **next = new int *[diff];

    for (ll i = 0; i < diff; i++)
    {
        field[i] = new int[Y];
        next[i] = new int[Y];
        for (ll j = 0; j < Y; j++)
            field[i][j] = next[i][j] = 0;
    }

    // MAKE START PARAMETRS
    if (rank == 0)
    {
        field[0][1] = 1;
        field[1][2] = 1;
        field[2][0] = 1;
        field[2][1] = 1;
        field[2][2] = 1;
    }

    int ***previous = new int **[limit];

    for (ll i = 0; i < limit; i++)
    {
        savePrevious(previous, i, field);

        int dest = (rank - 1 + size) % size;
        MPI_Isend(&field[0][0], Y, MPI_INT, dest, 0, MPI_COMM_WORLD, &sendFirst);

        dest = (rank + 1) % size;
        MPI_Isend(&field[diff - 1][0], Y, MPI_INT, dest, 0, MPI_COMM_WORLD, &sendEnd);

        dest = (rank - 1 + size) % size;
        MPI_Irecv(first_line, Y, MPI_INT, dest, 0, MPI_COMM_WORLD, &recvFirst);

        dest = (rank + 1) % size;
        MPI_Irecv(end_line, Y, MPI_INT, dest, 0, MPI_COMM_WORLD, &recvEnd);

        int *flags = new int[i];
        int **flagsMatrix = new int *[size];

        if (i > 0)
        {
            findFlags(flags, previous, field, i);

            for (ll j = 0; j < size; j++)
            {
                if (j != rank)
                    flagsMatrix[j] = new int[i];
                else
                    flagsMatrix[j] = flags;
            }

            for (ll j = 0; j < size; j++)
            {
                if (j != rank)
                {
                    MPI_Send(flags, i, MPI_INT, j, 0, MPI_COMM_WORLD);
                    MPI_Recv(&flagsMatrix[j][0], i, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        gameInside(field, next);

        MPI_Wait(&sendFirst, MPI_STATUS_IGNORE);
        MPI_Wait(&recvEnd, MPI_STATUS_IGNORE);
        gameLimitsEnd(field, next, end_line);

        MPI_Wait(&sendEnd, MPI_STATUS_IGNORE);
        MPI_Wait(&recvFirst, MPI_STATUS_IGNORE);
        gameLimitsFirst(field, next, first_line);
        

        if (i > 0)
        {
            flagPrev = checkFlags(flagsMatrix, i, size);

            for (ll j = 0; j < size; j++)
                delete flagsMatrix[j];
            delete flagsMatrix;

            if (flagPrev >= 0)
            {
                flagNow = i;
                break;
            }
        }

        rememberField(field, next);
    }

    elapsed_time = MPI_Wtime() - start_time;

    delete first_line;
    delete end_line;

    for (ll i = 0; i < diff; i++)
    {
        delete field[i];
        delete next[i];
    }
    delete field;
    delete next;

    if (flagNow != -1)
        limit = flagNow;

    for (ll i = 0; i < limit; i++)
    {
        for (ll l = 0; l < diff; l++)
            delete previous[i][l];
        delete previous[i];
    }
    delete previous;

    MPI_Finalize();
    if (rank == 0)
    {
        ofstream fout;
        fout.open("out_program.txt");
        fout << "Time program on " << size << " processes is " << elapsed_time << "\n";
        if (flagNow != -1)
            fout << "Game has cycle on iterations " << flagPrev << ' ' << flagNow << "\n";
        else
            fout << "Game reach limit\n";
    }
}