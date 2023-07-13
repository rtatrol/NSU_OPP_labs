#include <iostream>
#include <mpi.h>
#include <cmath>
#include <ctime>
#include <cstdarg>
#include <unistd.h>
#include <vector>
#include <pthread.h>
#include <stdio.h>
#define _OPERATION_NUMB 4
struct Task
{
    int repeatCount;
};

int rank;
int size;
pthread_mutex_t mutex;
pthread_t threads[2];
std::vector<Task> tasks;
int currentIter = 0;
int currentTask = 0;
MPI_Datatype Task_Info;
int L = 1000;

int GetNewTask(int askRank)
{
    int check = 1;
    MPI_Send(&check, 1, MPI_INT, askRank, 1, MPI_COMM_WORLD); // tags   1 = ask   2=answer 3=data
    MPI_Recv(&check, 1, MPI_INT, askRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (!check)
    {
        return 0;
    }

    Task task;
    MPI_Recv(&task, 1, Task_Info, askRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    pthread_mutex_lock(&mutex);
    tasks.push_back(task);
    pthread_mutex_unlock(&mutex);
    return 1;
}

double doTask(Task &task)
{
    double answ = 0;
    for (int i = 0; i < task.repeatCount; i++)
        answ += sin(i);
    return answ;
}

void generateTask(int countTask, int iteration)
{
    int tmp = std::abs(rank - (iteration % size));
    pthread_mutex_lock(&mutex);
    tasks.clear();
    Task task;
    for (int i = 0; i < countTask; i++)
    {
        task.repeatCount = std::abs(50 - i % 100) * tmp * L;
        tasks.push_back(task);
    }
    pthread_mutex_unlock(&mutex);
}

void refreshOthers(int *arr)
{
    for (int i = 0; i < size; i++)
        arr[i] = 1;
    arr[rank] = 0;
}

void *calcThread(void* _args)
{
    int *othersStatus = new int[size];
    for (currentIter = 0; currentIter < _OPERATION_NUMB; currentIter++)
    {

        double startTime = MPI_Wtime();
        int taskToGenerate = 100 * size;
        refreshOthers(othersStatus);
        generateTask(taskToGenerate, currentIter);

        currentTask = 0;

        int askRank;
        double result = 0;
        int taskCount = 0;
        int i = 0;

        while (true)
        {
            pthread_mutex_lock(&mutex);
            while (currentTask < tasks.size())
            {
                Task taskToDo = tasks[currentTask];
                pthread_mutex_unlock(&mutex);
                result += doTask(taskToDo);
                taskCount++;
                currentTask++;
                pthread_mutex_lock(&mutex);
            }
            pthread_mutex_unlock(&mutex);

            for (; i < size;)
            {
                askRank = (rank + i) % size;
                if (!othersStatus[askRank])
                {
                    i++;
                    continue;
                }
                else
                {
                    taskCount++;
                    othersStatus[askRank] = GetNewTask(askRank);
                    break;
                }
            }
            if (i == size)
                break;
        }
        double elapsedTime = MPI_Wtime() - startTime;
        printf("Process %d |||currentIteration is %d |||TasksDone is %d |||Time is %lf |||result is %lf\n",
               rank, currentIter, taskCount, elapsedTime, result);
        MPI_Barrier(MPI_COMM_WORLD);

        double maxTime, minTime;
        MPI_Reduce(&elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&elapsedTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            printf("Iteration %d |||Time of imbalance: %lf\nFraction of imbalance: %lf\n", currentIter, maxTime - minTime,
                   (maxTime - minTime) / minTime * 100.0);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    int progEXIT = 0;
    MPI_Send(&progEXIT, 1, MPI_INT, rank, 1, MPI_COMM_WORLD);
    delete othersStatus;
    return NULL;
}

void *dataThread(void* _args)
{
    while (currentIter < 4)
    {
        MPI_Status status;
        int res;
        MPI_Recv(&res, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);

        if (res == 0)
            break;

        pthread_mutex_lock(&mutex);
        if (currentTask >= tasks.size())
        {
            pthread_mutex_unlock(&mutex);
            int answer = 0;
            MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
        }
        else
        {
            Task taskToSend = tasks.back();
            tasks.pop_back();
            pthread_mutex_unlock(&mutex);
            int answer = 1;
            MPI_Send(&answer, 1, MPI_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
            MPI_Send(&taskToSend, 1, Task_Info, status.MPI_SOURCE, 3, MPI_COMM_WORLD);
        }
    }
    return NULL;
}

int Calculations()
{
    pthread_attr_t attrs;
    if (0 != pthread_attr_init(&attrs))
    {
        perror("Cannot initialize attributes");
        return 0;
    }
    if (0 != pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE))
    {
        perror("Error in setting attributes");
        return 0;
    }

    if (pthread_create(&threads[0], &attrs, calcThread, NULL) != 0 or
        pthread_create(&threads[1], &attrs, dataThread, NULL) != 0)
    {
        perror("Cannot create a thread");
        return 0;
    }

    pthread_attr_destroy(&attrs);

    for (int i = 0; i < 2; i++)
    {
        if (pthread_join(threads[i], NULL)!=0)
        {
            perror("Cannot join a thread");
            return 0;
        }
    }
    return 1;
}

int main(int argc, char **argv)
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE)
    {
        perror("cant make MPI_THREAD_MULTIPLE\n");
        return -1;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int count = 1;
    int array_of_blocklengths[] = {1};
    MPI_Aint array_of_displacements = 0;
    MPI_Datatype array_of_types[] = {MPI_INT};
    MPI_Type_create_struct(count, array_of_blocklengths, &array_of_displacements, array_of_types, &Task_Info);
    MPI_Type_commit(&Task_Info);

    pthread_mutex_init(&mutex, NULL);

    double startTime = MPI_Wtime();

    if (!Calculations())
    {
        MPI_Type_free(&Task_Info);
        pthread_mutex_destroy(&mutex);
        MPI_Finalize();
        perror("Error in calculations");
        return -1;
    }

    double elapsedTime = MPI_Wtime() - startTime;

    printf("Process %d time is %lf\n", rank, elapsedTime);

    MPI_Type_free(&Task_Info);
    pthread_mutex_destroy(&mutex);
    MPI_Finalize();
    return 0;
}