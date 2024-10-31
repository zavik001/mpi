#include <mpi.h>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <omp.h>

constexpr int N = 100000;

using namespace std;

static vector<int> GenerationVector()
{
    random_device rd; mt19937 mersenne(rd());
    vector<int> x = vector<int>(N);
    for (int i = 0; i < N; i++)
        x[i] = mersenne() % 73;
    return x;
}

static double Norma(vector<int>& x)
{
    int norm = 0;
    #pragma omp parallel for reduction (+:norm)
    for (int i = 0; i < N; i++)
        norm += x[i] * x[i];
    return sqrt(norm);
}

/**
 * @brief Main function of this program.
 *
 * This program is an MPI example that computes the norm of a vector on each
 * process and sends it to process 0 which prints the results.
 *
 * MPI_Init() and MPI_Finalize() are called to initialize and finalize the MPI
 * environment. MPI_Comm_rank() and MPI_Comm_size() are used to query the rank of
 * the current process and the total number of processes in the MPI_COMM_WORLD
 * communicator. MPI_Send() and MPI_Recv() are used to send and receive the norms
 * of the vectors. If the rank is not 0, the program generates a vector of random
 * integers and computes its norm, then sends the norm to process 0. If the rank
 * is 0, the program receives the norms from all other processes and prints them
 * to stdout.
 */
int main(int one, char** two)
{
    double norm_new;
    int rank, size;

    MPI_Init(&one, &two);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status st;

    if (rank != 0)
    {
        vector<int> x = GenerationVector();
        norm_new = Norma(x);
        MPI_Send(&norm_new, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
        cout << "Send norm: " << norm_new << endl;
    }
    else
    {
        for (int i = 1; i < size; i++)
        {
            MPI_Recv(&norm_new, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
            cout << "Process " << st.MPI_TAG << "; sent norm of vector " << norm_new << endl;
        }
    }

    MPI_Finalize();
    return 0;
}
