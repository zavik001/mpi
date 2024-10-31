# MPI

This project is a collection of MPI (Message Passing Interface) programs that demonstrate various concepts and techniques for parallel computing.

## Introduction

MPI is a standardized and portable message-passing system designed to function on a wide variety of parallel computers. This project provides a set of example programs that illustrate how to use MPI for parallel computing.

## Programs

The following programs are included in this project:

* `mpi_size_rank.cpp`: A simple program that demonstrates how to use MPI to determine the rank and size of a process group.
* `mpi_Bcast.cpp`: A program that demonstrates how to use MPI's broadcast function to send data from one process to all other processes in a group.
* `mpi_recv_send.cpp`: A program that demonstrates how to use MPI's send and receive functions to exchange data between processes.
* `mpi_Irecv_Test.cpp`: A program that demonstrates how to use MPI's non-blocking receive function to receive data from other processes.
* `mpi_openmp.cpp`: A program that demonstrates how to use MPI and OpenMP together to perform parallel computations.
