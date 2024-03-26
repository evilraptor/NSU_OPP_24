#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// XXX
#define N 100//35000

#define MAXITER 15000
#define SETVALVARIENT 1
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

/////////////////////////

int main(int argc, char **argv) {
	clock_t start_time, end_time;
	double eps = 1e-5, tau = 1e-3;
	
	double *A = NULL;
	double *u = NULL;
	double *x = (double *)malloc(N * sizeof(double));
	double *b = (double *)malloc(N * sizeof(double));
	double *y = (double *)malloc(N * sizeof(double));
	
	start_time = clock();
	A = (double *)malloc(N * N * sizeof(double));
	u = (double *)malloc(N * sizeof(double));
	set_values(A, u, x, b, 0);
	
	double *Ax = (double *)malloc(N * sizeof(double));
	for (int i = 0; i < N; i++) 
		Ax[i] = 0.0;
	double *temp_vector2 = (double *)malloc(N * sizeof(double));
	
	
	for (int i = 0; i < N; i++) {
		Ax[i] = 0.0;
	}
	
	
/////////////////////////////
	int state = 1, counter = 0;
	while (state) {
		
		//Ax 		
		multiply_matrix_and_vector(A, N, x, Ax, 0);
		//Ax-b-> y
		minus_vectors(Ax, b, N, y, 0);
		
		//A*y->Ay (temp_vector2)			
		multiply_matrix_and_vector(A, N, y, temp_vector2, 0);
		//tau=(y,Ay)/(Ay,Ay)
		tau = dot_product(y, temp_vector2, N) / dot_product(temp_vector2, temp_vector2, N);
		
		//tau*y -> y
		multiply_vector_and_scalar(y, N, tau, y, 0);
		//x-tau*y -> x
		minus_vectors(x, y, N,x, 0);
		
		//check
		multiply_matrix_and_vector(A, N, x, Ax, 0);
		minus_vectors(Ax, b, N, temp_vector2, 0);
		if (check_criteria(temp_vector2, b, eps))
			state = 0;
		counter++;
	}
	
	end_time = clock();
	printf("time: %f\n",end_time-start_time);
	printf("iterations: %d\n",counter);	
	if(argc>1){
		char *tmp=argv[1];
		if (strcmp(tmp,"1")==0){
		if (SETVALVARIENT == 2) print_vector(u,min(6,N),"u");
		print_vector(x,min(6,N),"x");
		}
	}
	free(A);
	free(u);
	free(x);
	free(b);
	free(temp_vector2);
	free(Ax);
	
	return 0;
}
