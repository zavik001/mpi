#include <mpi.h>
#include <iostream>
#include <random>
#include <cmath>
#include <vector>

constexpr int N = 100000;

using namespace std;

/**
 * @brief Generate a vector of size N filled with random integers
 *        between 0 and 72 (inclusive) using the Mersenne Twister
 *        algorithm.
 *
 * @return A vector of size N filled with random integers
 */
static vector<int> GenerationVector()
{
    random_device rd;
    mt19937 mersenne(rd());
    vector<int> x(N);
    for (int i = 0; i < N; i++)
        x[i] = mersenne() % 73;
    return x;
}


/**
 * @brief Compute the Euclidean norm of a given vector x.
 *
 * The norm is computed as the square root of the sum of the squares of all
 * elements in the vector.
 *
 * @param[in] x The vector to compute the norm of
 * @return The Euclidean norm of the vector x
 */
static double Norma(vector<int>& x)
{
    int norm = 0;
    for (int i = 0; i < N; i++)
        norm += x[i] * x[i];
    return sqrt(norm);
}

/**
 * @brief Simple MPI program to compute the norm of a vector on each
 *        process and send it to process 0 which prints the results.
 *
 * MPI_Init() and MPI_Finalize() are called to initialize and finalize
 * the MPI environment. MPI_Comm_rank() and MPI_Comm_size() are used to
 * query the rank of the current process and the total number of processes
 * in the MPI_COMM_WORLD communicator. If the rank is not 0, the program
 * generates a vector of random integers and computes its norm, then sends
 * the norm to process 0. If the rank is 0, the program receives the norms
 * from all other processes and prints them to stdout.
 */
int main(int argc, char** argv)
{
    double norm_new;
    int rank, size;

    MPI_Init(&argc, &argv);

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
