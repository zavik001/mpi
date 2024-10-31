#include <iostream>

#include <mpi.h>

using namespace std;

/**
 * @brief Simple MPI program to print rank and size of MPI_COMM_WORLD
 *
 * MPI_Init() and MPI_Finalize() are called to initialize and finalize the MPI
 * environment. MPI_Comm_rank() and MPI_Comm_size() are used to query the rank of
 * the current process and the total number of processes in the MPI_COMM_WORLD
 * communicator. The results are printed to stdout.
 */
int main(int argc, char **argv)
{
    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cout << "rank:" << rank << " size:" << size << endl;

    MPI_Finalize();

    return 0;
}