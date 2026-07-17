#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <malloc.h>

void fileCreation1(int);//ћодельна€ задача с заданным решением
void fileCreation2(int);//ћодельна€ задача с произвольным решением

int main() {
	int N = 10000;
	fileCreation2(N);
	return 0;
}

void fileCreation1(int N) {
	FILE* out = fopen("in.txt", "w");

	//ћодельна€ задача с заданным решением
	fprintf(out, "%d\n", N);
	for (int i = 0; i < N; i++)
		fprintf(out, "%lf\n", N+1);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			if (i == j) fprintf(out, "%lf ", 2);
			else fprintf(out, "%lf", 1);
		fputc('\n', out);
	}

	fclose(out);
}

void fileCreation2(int N) {
	FILE* out = fopen("in.txt", "w");

	double* u = (double*)malloc(sizeof(double) * N);
	for (int i = 0; i < N; i++) 
		u[i] = sin((2 * M_PI * i) / N);

	double* b = (double*)malloc(sizeof(double) * N);
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum = 0;
		//printf_s("%d row\n", i);
		for (int j = 0; j < N; j++) {
			sum += u[j];
			if (i == j) sum += u[j];
		}
		b[i] = sum;
	}

	//ћодельна€ задача с заданным решением
	fprintf(out, "%d\n", N);
	for (int i = 0; i < N; i++)
		fprintf(out, "%lf\n", b[i]);
	for (int i = 0; i < N; i++) {
		//printf("output row %d\n", i);
		for (int j = 0; j < N; j++)
			if (i == j) fprintf(out, "%d ", 2);
			else fprintf(out, "%d ", 1);
		fputc('\n', out);
	}


	free(b);
	free(u);
	fclose(out);
}