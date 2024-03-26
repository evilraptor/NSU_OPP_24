//////////with_comm_ptp.c
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
	
	if (rank == 0) {
		for(int i = 0; i < arr_size; i++){
			if(i < N){
				a[i] = i;
				b[i] = i;
			}
			else{
				a[i] = 0;
				b[i] = 0;
			}
		}
		double start, finish;
		start = MPI_Wtime();
		int shift = chunk_size;
		for(int i = 1; i < size; i++){
			MPI_Send(a + shift, chunk_size, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(b, N, MPI_INT, i, 1, MPI_COMM_WORLD);
			shift += chunk_size;
		}
		
		s = Calculate(a, chunk_size, b, N);
		long long stmp = 0;
		for (int i = 1; i < size; i++) {
			MPI_Recv(&stmp, 1, MPI_LONG_LONG, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
			//printf("%d got something\n", rank);
			s += stmp;
		}
		finish = MPI_Wtime();
		
		//printf("S: %lld\n", s);
		printf("Time: %f seconds\n", finish - start);
		printf("Sp: %f\n", 25 / (finish - start));
		printf("Ep: %f\n", (25 / (finish - start) / size) * 100);
	} else {
		//printf("%d sending s\n", rank);
		MPI_Recv(a, chunk_size, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(b, N, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		s = Calculate(a, chunk_size, b, N);
		MPI_Send(&s, 1, MPI_LONG_LONG, 0, 1, MPI_COMM_WORLD);
	}
	
	
	free(a);
	free(b);
	MPI_Finalize();
	return 0;
}
