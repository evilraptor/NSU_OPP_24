#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 150000

int main() {
	int a[N], b[N];
	long long s = 0;
	for (int i = 0; i < N; i++) {
		a[i] = i;
		b[i] = i;
	}
	
	clock_t start_time = clock();	
	for (int i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			s += a[i] * b[j];
		}
	}
	clock_t end_time = clock();
	double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

	printf("result: %lld\n", s);
	printf("time: %f\n", total_time);
	return 0;
}
