# Kmeans_MPI
K-means algorithm implemented in C, paralelized with MPI

How to generate random data to be clustered:

`python getinput.py k n > dados.txt`

where k is the number of clusters and n is the number of points. Example for k=100 and n=40000:

`python getinput.py 100 40000 > dados.txt`

How to compile the kmeans parallel C program:

`mpicc -o kmeans_MPI kmeans_MPI.c -lm -Wall -O3`

How to execute:

`mpirun -np NPROCS ./kmeans_MPI < dados.txt`

where NPROCS is the number of processes.
