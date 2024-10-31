#include <mpi.h>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <thread>
#include <chrono>
#include <omp.h>

constexpr int N = 100000;

using namespace std;

static vector<int> GenerationVector()
{
    random_device rd;
    mt19937 mersenne(rd());
    vector<int> x = vector<int>(N);
    for (int i = 0; i < N; i++)
        x[i] = mersenne() % 73;
    return x;
}

static double Norma(vector<int> &x)
{
    int norm = 0;
#pragma omp parallel for reduction(+ : norm)
    for (int i = 0; i < N; i++)
        norm += x[i] * x[i];
    return sqrt(norm);
}

/**
 * @brief The main function of this program.
 *
 * This program is an MPI example that computes the norm of a vector
 * on each process and sends it to process 0 which prints the results.
 *
 * MPI_Init() and MPI_Finalize() are called to initialize and finalize the MPI
 * environment. MPI_Comm_rank() and MPI_Comm_size() are used to query the rank of
 * the current process and the total number of processes in the MPI_COMM_WORLD
 * communicator. MPI_Send() and MPI_Irecv() are used to send and receive the
 * norms of the vectors. MPI_Test() is used to check if a non-blocking receive
 * has completed.
 *
 * If the rank is not 0, the program generates a vector of random integers and
 * computes its norm, then sends the norm to process 0. The program sleeps for 10
 * seconds if the rank is 1. If the rank is 0, the program receives the norms
 * from all other processes and prints them to stdout.
 */
int main(int argc, char **argv)
{
    double norm_new;
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    vector<MPI_Request> requests(size - 1);

    if (rank != 0)
    {
        vector<int> x = GenerationVector();
        norm_new = Norma(x);

        if (rank == 1)
        {
            this_thread::sleep_for(chrono::seconds(10));
        }

        MPI_Send(&norm_new, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
        cout << "Sent norm: " << norm_new << " from rank " << rank << endl;
    }
    else
    {
        vector<double> norms(size - 1);
        vector<int> received(size - 1, 0);
        int completed_count = 0;

        for (int i = 1; i < size; i++)
        {
            MPI_Irecv(&norms[i - 1], 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &requests[i - 1]);
        }

        while (completed_count < size - 1)
        {
            for (int i = 0; i < size - 1; i++)
            {
                if (!received[i])
                {
                    int flag = 0;
                    MPI_Test(&requests[i], &flag, &status);

                    if (flag)
                    {
                        received[i] = 1;
                        completed_count++;
                        cout << "Process " << status.MPI_TAG << "; received norm of vector " << norms[i] << endl;
                    }
                }
            }
        }
    }

    MPI_Finalize();
    return 0;
}
