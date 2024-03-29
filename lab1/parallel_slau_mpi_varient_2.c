#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

// XXX
#define N 100//35000

//или сделать переменными окружения...
#define MAXITER 30000
#define SETVALVARIENT 3
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

void multiply_vector_and_scalar(double *vector, int size, double scalar, double *res_vector,int out_put_flag) {
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

//скалярное произведение
double dot_product(double *vector1, double *vector2, int size) {
	double result = 0;
	for (int i = 0; i < size; i++) 
		result += vector1[i] * vector2[i];
	return result;
}

int min(int a, int b) {
	return a < b ? a : b;
}

int check_criterial(double *vector1, double *vector2, double eps) {
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
	
	if (SETVALVARIENT == 1) {
		for (int i = 0; i < N; i++) {
			u[i] = 1.0;
			x[i] = 0.0;
			b[i] = N+1;
		}
	} else if (SETVALVARIENT == 2) {
		for (int i = 0; i < N; i++) {
			u[i] = sin(2 * M_PI * i / N);
			x[i] = 0.0;
			b[i] = 0.0;
		}
		multiply_matrix_and_vector(A, N, u, b, 0);
	} else if (SETVALVARIENT == 3) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				A[j*N +i] = (i == j ? i : N);
			}
		}
		for (int i = 0; i < N; i++) {
			u[i] = 0.0;
			x[i] = 0.0;
			b[i] = i * i;
		}
	}

	if(out_put_flag==1) 
		for (int i = 0; i < N; i++) printf("b: %f\n",b[i]);
}

void cut_vector(double *vector, double *cutted_vector, int size, int displs) {
	for (int i = 0; i < size; i++) cutted_vector[i] = vector[displs + i];
}

void set_sending_counts(int *send_counts, int *divided_send_counts, int size) {
	for (int i = 0; i < size; i++) {
		send_counts[i] = N / size;
		if (i < N % size) 
			++send_counts[i];
		divided_send_counts[i] = send_counts[i];
		send_counts[i] *= N;
	}
}

void set_displs(int *send_counts, int *displs, int *divided_displs, int size) {
	displs[0] = 0;
	for (int i = 1; i < size; i++) 
		displs[i] = displs[i-1]+send_counts[i-1];		
	for (int i = 0; i < size; i++) 
		divided_displs[i]= displs[i] / N;
}

void set_zero(double *array, int size) {
	for (int i = 0; i < size; i++) {
		array[i] = 0.0;
	}
}
/////////////////////////

int main(int argc, char **argv) {
	int rank, size;
	double start_time = 0.0, end_time = 0.0;
	double eps = 1e-5, tau = 1e-3;
	double global_dot_product1 = 0.0, global_dot_product2 = 0.0;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int *send_counts = (int *)malloc(size * sizeof(int));
	int *divided_send_counts = (int *)malloc(size * sizeof(int));
	int *displs = (int *)malloc(size * sizeof(int));
	int *divided_displs = (int *)malloc(size * sizeof(int));
	
	double *A = NULL;
	double *u = NULL;
	double *Ax = NULL;
	double *x = (double *)malloc(N * sizeof(double));
	double *b = (double *)malloc(N * sizeof(double));
	double *y = (double *)malloc(N * sizeof(double));
	
	if(rank == 0){
		start_time = MPI_Wtime();
		A = (double *)malloc(N * N * sizeof(double));
		Ax = (double *)malloc(N * sizeof(double));
		u = (double *)malloc(N * sizeof(double));
		set_values(A, u, x, b, 0);
		set_zero(Ax, N);
	}
	set_sending_counts(send_counts, divided_send_counts, size);
	set_displs(send_counts, displs, divided_displs, size);
	
	double *local_Ax = (double *)malloc(divided_send_counts[rank] * sizeof(double));
	double *cutted_b = (double *)malloc(divided_send_counts[rank] * sizeof(double));
	double *cutted_x = (double *)malloc(divided_send_counts[rank] * sizeof(double));
	double *cutted_y = (double *)malloc(divided_send_counts[rank] * sizeof(double));
	double *temp_vector1 = (double *)malloc(N * sizeof(double));
	double *temp_vector2 = (double *)malloc(N * sizeof(double));
	
	double *part_A = (double *)malloc(send_counts[rank] * sizeof(double));
	MPI_Scatterv(A, send_counts, displs, MPI_DOUBLE, part_A, send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	cut_vector(b, cutted_b, divided_send_counts[rank], divided_displs[rank]);	
	set_zero(local_Ax, divided_send_counts[rank]);
	set_zero(cutted_y, divided_send_counts[rank]);
	set_zero(temp_vector1, N);
	set_zero(temp_vector2, N);
	
/////////////////////////////
	int state = 1, counter = 0;
	while (state) {
		cut_vector(x, cutted_x, divided_send_counts[rank], divided_displs[rank]);
		
		//Ax		
		multiply_matrix_and_vector(part_A, divided_send_counts[rank], x, local_Ax, 0);
		//Ax-b->cutted_y
		minus_vectors(local_Ax, cutted_b, divided_send_counts[rank], cutted_y, 0);
		MPI_Allgatherv(cutted_y, divided_send_counts[rank], MPI_DOUBLE, y, divided_send_counts, divided_displs, MPI_DOUBLE, MPI_COMM_WORLD);		
		
		//A*y->Ay (temp_vector1)
		multiply_matrix_and_vector(part_A, divided_send_counts[rank], y, temp_vector1, 0);
		double local_dot_product1 = dot_product(cutted_y, temp_vector1, divided_send_counts[rank]);
		double local_dot_product2 = dot_product(temp_vector1, temp_vector1, divided_send_counts[rank]);
		MPI_Reduce(&local_dot_product1, &global_dot_product1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&local_dot_product2, &global_dot_product2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0)
			tau = global_dot_product1 / global_dot_product2;
		MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		//tau*y -> cutted_y
		multiply_vector_and_scalar(cutted_y, divided_send_counts[rank], tau, cutted_y, 0);
		//x-tau*y -> cutted_x
		minus_vectors(cutted_x, cutted_y, divided_send_counts[rank],cutted_x, 0);
		MPI_Allgatherv(cutted_x, divided_send_counts[rank], MPI_DOUBLE, x, divided_send_counts, divided_displs, MPI_DOUBLE, MPI_COMM_WORLD);

/////////////////////////////

		multiply_matrix_and_vector(part_A, divided_send_counts[rank], x, local_Ax, 0);
		minus_vectors(local_Ax, cutted_b, divided_send_counts[rank], temp_vector1, 0);
		MPI_Gather(temp_vector1, divided_send_counts[rank], MPI_DOUBLE, temp_vector2, divided_send_counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (rank == 0)
			if (check_criterial(temp_vector2, b, eps)) state = 0;
		MPI_Bcast(&state, 1, MPI_INT, 0, MPI_COMM_WORLD);
		counter++;
	}
/////////////////////////////
	
	if(rank==0){
		end_time = MPI_Wtime();
		printf("size: %d\n", size);
		printf("time: %f\n",end_time-start_time);
		printf("Sp: %f\n", SEQTIME / (end_time - start_time));
		printf("Ep: %f\n", (SEQTIME / (end_time - start_time) / size) * 100);
		printf("iterations: %d\n",counter);
		if(argc > 1){
			char *tmp=argv[1];
			if (strcmp(tmp,"1")==0){
				print_vector(u,min(6,N),"u");
				print_vector(x,min(6,N),"x");
			}
		}
		free(A);
		free(u);
		free(Ax);
	}
	free(send_counts);
	free(displs);
	free(divided_send_counts);
	free(divided_displs);
	free(x);
	free(b);
	free(local_Ax);
	free(cutted_b);
	free(cutted_x);
	free(cutted_y);
	free(y);
	free(temp_vector1);
	free(temp_vector2);
	
	MPI_Finalize();
	return 0;
}
