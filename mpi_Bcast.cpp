#include <mpi.h>
#include <iostream>

using namespace std;

/**
 * @brief A simple MPI program that takes an integer from the user and broadcasts
 *        it to all processes in the MPI_COMM_WORLD communicator.
 *
 * MPI_Init() and MPI_Finalize() are called to initialize and finalize the MPI
 * environment. MPI_Comm_rank() and MPI_Comm_size() are used to query the rank of
 * the current process and the total number of processes in the MPI_COMM_WORLD
 * communicator. The results are printed to stdout. MPI_Bcast() is used to
 * broadcast the integer to all processes in the MPI_COMM_WORLD communicator.
 */
int main(int argc, char** argv) {
    int rank, size;
    int number;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        cout << "Enter an integer: ";
        cin >> number;
    }

    MPI_Bcast(&number, 1, MPI_INT, 0, MPI_COMM_WORLD);

    cout << "Process " << rank << " got the number: " << number << std::endl;

    MPI_Finalize();
    return 0;
}
