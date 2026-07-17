#define _USE_MATH_DEFINES

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define EPSILON 0.00001
#define N 50000

void fill1(double* AChunk, double* b, int blockSize, int blockLow) {
	for (int i = 0; i < N; i++) {
		b[i] = N + 1;
		if (i < blockSize) {
			for (int j = 0; j < N; j++) {
				if (i + blockLow == j) AChunk[(long long int)i * N + j] = 2;
				else AChunk[(long long int)i * N + j] = 1;
			}
		}
	}
}

void fill2(double* AChunk, double* b, int blockSize, int blockLow, int rank) {

	double* u = (double*)malloc(sizeof(double) * N);

	for (int i = 0; i < N; i++) {
		u[i] = sin(2 * M_PI * i / N);
		if (i < blockSize) {
			for (int j = 0; j < N; j++) {
				if (i + blockLow == j) AChunk[(long long int)i * N + j] = 2;
				else AChunk[(long long int)i * N + j] = 1;
			}
		}
	}

	for (int i = 0; i < N; i++) {
		b[i] = 0;
		for (int j = 0; j < N; j++) {
			if (i == j) b[i] += 2 * u[j];
			else b[i] += u[j];
		}
	}

	free(u);
}

int main(int argc, char** argv) {
	int size, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	clock_t start, end;

	int blockSize = N / size;
	int blockLow = rank * N / size;

	double* AChunk = (double*)malloc(sizeof(double) * blockSize * N);
	double* localTmp = (double*)malloc(sizeof(double) * blockSize);
	double* b = (double*)malloc(sizeof(double) * N);//Дублируются в каждом процессе
	double* y = (double*)malloc(sizeof(double) * N);
	double* Ay = (double*)malloc(sizeof(double) * N);
	double* x = (double*)calloc(N, sizeof(double));//Дублируются в каждом процессе

	double bNorm = 0;
	double criteria = 0;
	double tau = 0;
	double divisible = 0; // (y^n, Ay^n)
	double divisor = 0; // (Ay^n, Ay^n)
	double localBNorm = 0;
	double localCriteria = 0;
	double localDivisible = 0;
	double localDivisor = 0;
	double tmp;

	//Один из вариантов заполнения
	fill2(AChunk, b, blockSize, blockLow, rank);

	start = clock();

	//Расчёт нормы вектора b
	for (int i = 0; i < blockSize; i++)
		localBNorm += b[blockLow + i] * b[blockLow + i];
	MPI_Allreduce(&localBNorm, &bNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	bNorm = sqrt(bNorm);

	do {
		//y^n = Ax^n - b
		for (long long int i = 0; i < blockSize; i++) {
			localTmp[i] = 0;
			for (long long int j = 0; j < N; j++)
				localTmp[i] += AChunk[i * N + j] * x[j];
			localTmp[i] = localTmp[i] - b[i + blockLow];

		}
		MPI_Allgather(localTmp, blockSize, MPI_DOUBLE, y, blockSize, MPI_DOUBLE, MPI_COMM_WORLD);

		//t^n = (y^n, Ay^n)/(Ay^n, Ay^n)
		localDivisible = 0;
		localDivisor = 0;
		for (int i = blockLow; i < blockLow + blockSize; i++) {
			Ay[i] = 0;
			for (int j = 0; j < N; j++) {
				Ay[i] += AChunk[(long long int)(i - blockLow) * N + j] * y[j];
			}
			localDivisible += y[i] * Ay[i];
			localDivisor += Ay[i] * Ay[i];
		}

		MPI_Allreduce(&localDivisible, &divisible, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&localDivisor, &divisor, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tau = divisible / divisor;

		//X^(n+1) = x^n - t^n y^n
		for (int i = 0; i < blockSize; i++)
			localTmp[i] = x[i + blockLow] - tau * y[i + blockLow];
		MPI_Allgather(localTmp, blockSize, MPI_DOUBLE, x, blockSize, MPI_DOUBLE, MPI_COMM_WORLD);

		//Пересчёт критерия
		localCriteria = 0;
		for (int i = 0; i < blockSize; i++) {
			tmp = 0;
			for (int j = 0; j < N; j++)
				tmp += AChunk[(long long int)i * N + j] * x[j];
			tmp -= b[i + blockLow];
			localCriteria += tmp * tmp;
		}
		MPI_Allreduce(&localCriteria, &criteria, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		criteria = sqrt(criteria);
	} while (criteria / bNorm >= EPSILON);

	end = clock();
	if (rank == 0) printf("Time taken: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);

	free(x);
	free(b);
	free(AChunk);
	free(Ay);
	free(y);
	free(localTmp);

	MPI_Finalize();

	return 0;
}