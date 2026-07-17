#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>

#define EPSILON 0.00001
#define N 65000

double norm(double*);
double scalarMult(double*, double*);
double* matrixOnVector(double*, double*);
double* scalarOnVector(double*, double);
double* vectorSub(double*, double*);

int main(int argc, const char** argv) {
	clock_t start, end;
	start = clock();

	double* tmp;
	double* b = (double*)malloc(sizeof(double) * N);
	double* A = (double*)malloc(sizeof(double) * N * N);
	double* x = (double*)malloc(sizeof(double) * N);

	//Один из вариантов заполнения
	for (long long int i = 0; i < N; i++) {
		b[i] = N + 1;
		for (long long int j = 0; j < N; j++) {
			if (i == j) A[i * N + j] = 2;
			else A[i * N + j] = 1;
		}
	}

	double* y = vectorSub(matrixOnVector(A, x), b);
	tmp = matrixOnVector(A, y);
	double t = scalarMult(y, tmp) / scalarMult(tmp, tmp);

	int counter = 0;

	while (norm(y) / norm(b) >= EPSILON) {
		x = vectorSub(x, scalarOnVector(y, t));
		
		free(tmp);
		tmp = scalarOnVector(y, t);

		free(tmp);
		free(y);
		tmp = matrixOnVector(A, x);
		y = vectorSub(tmp, b);

		free(tmp);
		tmp = matrixOnVector(A, y);
		t = scalarMult(y, tmp) / scalarMult(tmp, tmp);
	}

	end = clock();

	printf("time taken: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);

	free(tmp);
	free(A);
	free(b);
	free(x);
	free(y);

	return 0;
}

double norm(double* vec) {
	double output = 0;
	for (long long int i = 0; i < N; i++) {
		output += vec[i] * vec[i];
	}
	output = sqrt(output);
	return output;
}

double scalarMult(double* vec1, double* vec2) {
	double output = 0;
	for (long long int i = 0; i < N; i++)
		output += vec1[i] * vec2[i];
	return output;
}

double* matrixOnVector(double* matrix, double* vec){
	double* output = (double*)malloc(sizeof(double) * N);
	for (long long int row = 0; row < N; row++) {
		output[row] = 0;
		for (long long int column = 0; column < N; column++) {
			output[row] += matrix[row * N  + column] * vec[column];
		}
	}
	return output;
}

double* scalarOnVector(double* vec, double a) {
	double* output = (double*)malloc(sizeof(double) * N);
	for (long long int i = 0; i < N; i++)
		output[i] = a * vec[i];
	return output;
}

double* vectorSub(double* vec1, double* vec2) {
	double* output = (double*)malloc(sizeof(double) * N);
	for (long long int i = 0; i < N; i++)
		output[i] = vec1[i] - vec2[i];
	return output;
}