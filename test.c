#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MATRIX_SIZE 4

void matrix_vector_multiply(int* local_matrix, int* local_vector, int* local_result, int local_size) {
    for (int i = 0; i < local_size; i++) {
        local_result[i] = 0;
        for (int j = 0; j < MATRIX_SIZE; j++) {
            local_result[i] += local_matrix[i * MATRIX_SIZE + j] * local_vector[j];
        }
    }
}

int main(int argc, char* argv[]) {
    int rank, size;
    int matrix[MATRIX_SIZE][MATRIX_SIZE] = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12},
        {13, 14, 15, 16}
    };
    int vector[MATRIX_SIZE] = { 1, 2, 3, 4 };
    int result[MATRIX_SIZE];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_size = MATRIX_SIZE / size;
    int local_matrix[local_size * MATRIX_SIZE];
    int local_vector[MATRIX_SIZE];
    int local_result[local_size];

    MPI_Scatter(matrix, local_size * MATRIX_SIZE, MPI_INT, local_matrix, local_size * MATRIX_SIZE, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(vector, MATRIX_SIZE, MPI_INT, local_vector, MATRIX_SIZE, MPI_INT, 0, MPI_COMM_WORLD);

    matrix_vector_multiply(local_matrix, local_vector, local_result, local_size);

    MPI_Gather(local_result, local_size, MPI_INT, result, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Matrix:\n");
        for (int i = 0; i < MATRIX_SIZE; i++) {
            for (int j = 0; j < MATRIX_SIZE; j++) {
                printf("%d ", matrix[i][j]);
            }
            printf("\n");
        }

        printf("Vector:\n");
        for (int i = 0; i < MATRIX_SIZE; i++) {
            printf("%d ", vector[i]);
        }
        printf("\n");

        printf("Result:\n");
        for (int i = 0; i < MATRIX_SIZE; i++) {
            printf("%d ", result[i]);
        }
        printf("\n");
    }

    MPI_Finalize();

    return 0;
}
