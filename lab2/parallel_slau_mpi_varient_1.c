#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

// XXX
#define N 800//4

#define MAXITER 15000
#define SEQTIME 106

void print_vector(double *vector, int count, char* name) {
	printf("printing vector %s:\n",name);
	for (int i = 0; i < count; i++) printf("%f ", vector[i]);
	printf("\n");
}

void print_matrix(double *matrix, int lines_count, char* name) {
	printf("printing matrix %s:\n",name);
	for (int i = 0; i < lines_count; i++) {
		for (int j = 0; j < N; j++)
			printf("%f ", matrix[i*N+j]);
		printf("\n");
	}
}

void clear_vector(double *vector, int count) {
	for (int i = 0; i < count; i++) vector[i] = 0.0;
}

void multiply_scalar_and_vector(double *vector, int size, double scalar, double *res_vector,int out_put_flag) {
	for (int i = 0; i < size; i++) res_vector[i] = vector[i] * scalar;
	if(out_put_flag==1) print_vector(res_vector, size, "res_vector");
}

void plus_vectors(double *vector1, double *vector2, int size, double *res_vector,int out_put_flag) {
	for (int i = 0; i < size; i++) res_vector[i] = vector1[i] + vector2[i];
	if(out_put_flag==1) print_vector(res_vector, size, "res_vector");
}

void minus_vectors(double *vector1, double *vector2, int size, double *res_vector, int out_put_flag) {
	for (int i = 0; i < size; i++) res_vector[i] = vector1[i] - vector2[i];
	if(out_put_flag==1) print_vector(res_vector, size, "res_vector");
}

void multiply_matrix_and_vector(double *matrix, int lines_count, double *vector, double *res_vector, int out_put_flag) {
	for (int i = 0; i < lines_count; i++) {
		res_vector[i] = 0;
		for (int j = 0; j < N; j++) 
			res_vector[i] += matrix[i * N + j] * vector[j];
	}
	if(out_put_flag==1) print_vector(res_vector,lines_count,"res_vector");
}

double get_norm(double *vector) {
	double sum = 0;
	for (int i = 0; i < N; i++) sum += vector[i] * vector[i];
	sum = sqrt(sum);
	return sum;
}

double multiply_vectors(double *vector1, double *vector2, int size, int out_put_flag) {
	double res = 0;
	for (int i = 0; i < size; i++)
		res += vector1[i] * vector2[i];
	if(out_put_flag==1) printf("res: %f\n",res);	
	return res;
}

int min(int a, int b) {
	return a < b ? a : b;
}

int check_criteria(double *vector1, double *vector2, double eps) {
	if ((get_norm(vector1) / get_norm(vector2)) < eps) return 1;
	return 0;
}

void set_values(double *A, double *u, double *x, double *b, int out_put_flag) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) A[i * N + j] = 2.0;
			else A[i * N + j] = 1.0;
		}	
	}
	
	for (int i = 0; i < N; i++) {
		u[i] = sin(2 * M_PI * i / N);
		x[i] = 0.0;
		b[i] = 0.0;
	}
	multiply_matrix_and_vector(A, N, u, b, 0);
	
	if(out_put_flag==1) for (int i = 0; i < N; i++) printf("b: %f\n",b[i]);
}

void cut_x(double *x, double *cutted_x, int size, int displs) {
	for (int i = 0; i < size; i++) cutted_x[i] = x[i + displs];
}

/////////////////////////

int main(int argc, char **argv) {
	int rank, size;
	MPI_Status status;
	double start_time, end_time;
	double eps = 1e-5, tao = 1e-3;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int *send_counts = (int *)malloc(size * sizeof(int));
	int *divided_send_counts = (int *)malloc(size * sizeof(int));
	int *displs = (int *)malloc(size * sizeof(int));
	int *divided_displs = (int *)malloc(size * sizeof(int));
	
	double *A = NULL;
	double *u = NULL;
	double *x = (double *)malloc(N * sizeof(double));
	double *b = (double *)malloc(N * sizeof(double));
	double *temp_vector = (double *)malloc(N * sizeof(double));
	if(rank==0){
		A = (double *)malloc(N * N * sizeof(double));
		u = (double *)malloc(N * sizeof(double));
		set_values(A, u, x, b, 0);
		
		for (int i = 0; i < size; i++) {
			send_counts[i] = N / size;
			if (i < N % size) ++send_counts[i];
			send_counts[i] *= N;
		}
		displs[0] = 0;
		for (int i = 1; i < size; i++) {
			displs[i] = displs[i-1]+send_counts[i-1];
		}	
		for (int i = 0; i < size; i++) {
		divided_send_counts[i] = send_counts[i] / N;
		divided_displs[i]= displs[i] / N;
		}
		start_time = MPI_Wtime();
	}
	
	MPI_Bcast(send_counts, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(divided_send_counts, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(divided_displs, size, MPI_INT, 0, MPI_COMM_WORLD);
	
	double *part_A = (double *)malloc(send_counts[rank] * sizeof(double));
	MPI_Scatterv(A, send_counts, displs, MPI_DOUBLE, part_A, send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	double *temp_x = (double *)malloc(N * sizeof(double));
	double *Ax = (double *)malloc(divided_send_counts[rank] * sizeof(double));
	double *cutted_b = (double *)malloc(divided_send_counts[rank] * sizeof(double));
	double *cutted_x = (double *)malloc(divided_send_counts[rank] * sizeof(double));
	for (int i = 0; i < N; i++) temp_x[i] = 0.0;
	for (int i = 0; i < divided_send_counts[rank]; i++) {
		Ax[i] = 0.0;
		cutted_b[i] = b[i + divided_displs[rank]];
		cutted_x[i] = x[i + divided_displs[rank]];
	}
	
	
	int state = 1, counter = 0;
	while (state) {
	
		//AX
		multiply_matrix_and_vector(part_A, divided_send_counts[rank], x, Ax, 0);

		//Ax-b->temp_vector
		minus_vectors(Ax, cutted_b, divided_send_counts[rank], temp_vector, 0);	
		
		//temp_vector*tao->temp_vector
		multiply_scalar_and_vector(temp_vector, divided_send_counts[rank], tao, temp_vector, 0);

		//x-temp_vector*tao->temp_vector
		minus_vectors(cutted_x, temp_vector, divided_send_counts[rank],temp_vector, 0);
		
		MPI_Gatherv(temp_vector, divided_send_counts[rank], MPI_DOUBLE, x, divided_send_counts, divided_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (rank == 0) {
//			printf("%d\n",counter);
//			print_vector(x,N,"x");
			minus_vectors(x, b, N, temp_vector, 0);
			if (check_criteria(temp_vector, x, eps)){
				state = 0;
				printf("\ndone in %d iterations\n",counter);	
			}
			if (counter == MAXITER - 1) {
				printf("\niterations ended\n");
			}
		}
		MPI_Bcast(&state, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if (counter < MAXITER - 1) counter++;
		else  break;
		MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		cut_x(x, cutted_x, divided_send_counts[rank], divided_displs[rank]);
	
	}
	
	if(rank==0){
		end_time = MPI_Wtime();
		printf("time: %f\n",end_time-start_time);
		if(argc>1){
			char *tmp=argv[1];
			printf("%s\n",tmp);
			if (strcmp(tmp,"1")==0){
				print_vector(u,min(6,N),"u");
				print_vector(x,min(6,N),"x");
			}
		}
		free(A);
		free(u);
	}
	free(send_counts);
	free(displs);
	free(divided_send_counts);
	free(divided_displs);
	free(x);
	free(b);
	free(temp_vector);
	free(temp_x);
	free(Ax);
	free(cutted_b);
	free(cutted_x);
	
	MPI_Finalize();
	return 0;
}