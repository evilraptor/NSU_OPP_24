////////group.c
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#define N 150000

long long Calculate(int* a, int size_a, int* b, int size_b){
	long long result = 0;
	for (int i = 0; i < size_a; i++) {
		for (int j = 0; j < size_b; j++) {
			result += a[i] * b[j];
		}
	}
	return result;	
}

int main(int argc, char** argv) {
	int rank, size;
	long long int s = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	
	int arr_size;
	if (N % size != 0) arr_size = N + size - (N % size);
	else arr_size = N;
	int* a = malloc(arr_size * sizeof(int));
	int* b = malloc(arr_size * sizeof(int));
	int chunk_size = arr_size / size;
	int* recv_buffer = malloc(chunk_size * sizeof(int));

	double start, finish;
	if (rank == 0) {

		if (a && b) {
			for (int i = 0; i < arr_size; i++) {
				if (i < N) {
					a[i] = i;
					b[i] = i;
				}
				else {
					a[i] = 0;
					b[i] = 0;
				}
			}
		}
		start = MPI_Wtime();
	}

	long long int stmp = 0;

	MPI_Scatter(a, chunk_size, MPI_INT, recv_buffer, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(b, N, MPI_INT, 0, MPI_COMM_WORLD);
	stmp = Calculate(recv_buffer, chunk_size, b, N);
	MPI_Reduce(&stmp, &s, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	//printf("rank %d, send %lld, 1st %d\n", rank,stmp, recv_buffer[0]);
	if (rank == 0) {
		finish = MPI_Wtime();

		//printf("s: %lld\n", s);
		printf("Time: %f seconds\n", finish - start);
		printf("Sp = %f\n", 22 / (finish - start));
		printf("Ep = %f\n", (22 / (finish - start) / size) * 100);

	}

	free(a);
	free(b);
	free(recv_buffer);
	MPI_Finalize();
	return 0;

}
