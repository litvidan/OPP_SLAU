#define _USE_MATH_DEFINES

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define EPSILON 0.00001
#define N 70000

void fill1(double* AChunk, double* bChunk, int blockSize, int blockLow) {

	for (long long int i = 0; i < blockSize; i++) {
		bChunk[i] = N + 1;
		for (long long int j = 0; j < N; j++) {
			if (i + blockLow == j) AChunk[i * N + j] = 2;
			else AChunk[i * N + j] = 1;
		}
	}
}

void fill2(double* AChunk, double* b, int blockSize, int blockLow, int rank) {
	double* u = (double*)malloc(sizeof(double) * N);

	for (long long int i = 0; i < N; i++) {
		u[i] = sin(2 * M_PI * i / N);
	}

	for (long long int i = 0; i < N; i++) {
		if (i < blockSize) {
			b[i] = 0;
			for (int j = 0; j < N; j++) {
				if (i + blockLow == j) AChunk[i * N + j] = 2;
				else AChunk[i * N + j] = 1;
				b[i] += AChunk[i * N + j] * u[j];
			}
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

	double* AChunk   = (double*)malloc(sizeof(double) * blockSize * N);//В данной версии программы у каждого процесса столбцы матрицы, а не строки
	double* bChunk   = (double*)malloc(sizeof(double) * blockSize);
	double* xChunk   = (double*)calloc(blockSize, sizeof(double));
	double* yChunk   = (double*)malloc(sizeof(double) * blockSize);
	double* AyChunk  = (double*)malloc(sizeof(double) * N);
	double* AxChunk  = (double*)malloc(sizeof(double) * N);
	double* AxbChunk = (double*)malloc(sizeof(double) * blockSize);

	double bNorm          = 0;
	double criteria       = 0;
	double tau            = 0;
	double divisible      = 0; // (y^n, Ay^n)
	double divisor        = 0; // (Ay^n, Ay^n)
	double localBNorm     = 0;
	double localCriteria  = 0;
	double localDivisible = 0;
	double localDivisor   = 0;

	fill2(AChunk, bChunk, blockSize, blockLow, rank);

	start = clock();

	//Расчёт нормы вектора b
	for (int i = 0; i < blockSize; i++)
		localBNorm += bChunk[i] * bChunk[i];
	MPI_Allreduce(&localBNorm, &bNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	bNorm = sqrt(bNorm);

	do {
		//y^n = Ax^n - b
		for (long long int i = 0; i < blockSize; i++) {
			for (long long int j = 0; j < N; j++) {
				if (i == 0) AxChunk[j] = 0;
				AxChunk[j] += AChunk[i * N + j] * xChunk[i];
			}
			AxChunk[blockLow + i] -= bChunk[i];
		}

		for (long long int i = 0; i < size; i++) {
			MPI_Reduce((AxChunk + i * N / size), yChunk, N / size, MPI_DOUBLE, MPI_SUM, i, MPI_COMM_WORLD);
		}

		//t^n = (y^n, Ay^n)/(Ay^n, Ay^n)
		// Расчёт Ay^n
		for (long long int i = 0; i < blockSize; i++) {
			for (long long int j = 0; j < N; j++) {
				if (i == 0) AyChunk[j] = 0;
				AyChunk[j] += AChunk[i * N + j] * yChunk[i];
			}
		}
		for (long long int i = 0; i < size; i++) {
			MPI_Reduce((AyChunk + i * N / size), AxbChunk, N / size, MPI_DOUBLE, MPI_SUM, i, MPI_COMM_WORLD);//Axb используется, чтобы сэкономить память
		}
		// Расчёт делимого и делителя
		localDivisible = 0;
		localDivisor = 0;
		for (long long int i = 0; i < blockSize; i++) {
			localDivisible += yChunk[i] * AxbChunk[i];
			localDivisor   += AxbChunk[i] * AxbChunk[i];
		}
		MPI_Allreduce(&localDivisible, &divisible, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&localDivisor, &divisor, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tau = divisible / divisor;

		//X^(n+1) = x^n - t^n y^n
		for (long long int i = 0; i < blockSize; i++) {
			xChunk[i] = xChunk[i] - tau * yChunk[i];
		}

		//Расчёт критерия: ||Ax^n-b||
		localCriteria = 0;
		for (long long int i = 0; i < blockSize; i++) {
			for (long long j = 0; j < N; j++) {
				if (i == 0) AxChunk[j] = 0;
				AxChunk[j] += AChunk[i * N + j] * xChunk[i];
			}
			AxChunk[blockLow + i] -= bChunk[i];
		}

		for (int i = 0; i < size; i++) {
			MPI_Reduce(AxChunk + i * N / size, AxbChunk, N / size, MPI_DOUBLE, MPI_SUM, i, MPI_COMM_WORLD);
		}

		for (long long int i = 0; i < blockSize; i++) {
			AxbChunk[i] = AxbChunk[i] * AxbChunk[i];
			localCriteria += AxbChunk[i];
		}
		MPI_Allreduce(&localCriteria, &criteria, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		criteria = sqrt(criteria);

	} while (criteria / bNorm >= EPSILON);

	end = clock();
	if (rank == 0) {
		printf("Time taken: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);

		double* x = (double*)malloc(sizeof(double) * N);
		MPI_Gather(xChunk, blockSize, MPI_DOUBLE, x, blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	}
	else {
		MPI_Gather(xChunk, blockSize, MPI_DOUBLE, NULL, blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	free(AChunk);
	free(bChunk);
	free(xChunk);
	free(yChunk);
	free(AyChunk);
	free(AxChunk);
	free(AxbChunk);

	MPI_Finalize();

	return 0;
}