#include <iostream>
#include <mpi.h>
#include <iostream>
#include <cmath>

#define N 10

typedef struct {
    double real;
    double imag;
} complex;

typedef struct {
    int rows;
    int cols;
    complex *data;
} matrix;

void init_matrix(matrix *mat, int rows, int cols) {
    mat->rows = rows;
    mat->cols = cols;
    mat->data = (complex*) malloc(rows * cols * sizeof(complex));
}

void free_matrix(matrix *mat) {
    free(mat->data);
}

void print_matrix(matrix *mat) {
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            std::cout << "(" << mat->data[i * mat->cols + j].real << "," << mat->data[i * mat->cols + j].imag << ") ";
        }
        std::cout << std::endl;
    }
}

bool compare_matrices(matrix *mat1, matrix *mat2, double tolerance) {
    if (mat1->rows != mat2->rows || mat1->cols != mat2->cols) {

        std::cout << "Sizes do not match!" << std::endl;

        return false;
    }
    for (int i = 0; i < mat1->rows; i++) {
        for (int j = 0; j < mat1->cols; j++) {
            double diff_real = std::abs(mat1->data[i * mat1->cols + j].real - mat2->data[i * mat2->cols + j].real);
            double diff_imag = std::abs(mat1->data[i * mat1->cols + j].imag - mat2->data[i * mat2->cols + j].imag);
            if (diff_real > tolerance || diff_imag > tolerance) {
                return false;
            }
        }
    }
    return true;
}

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialize matrices
    matrix A, B, C, C_seq;
    init_matrix(&A, N, N);
    init_matrix(&B, N, N);
    init_matrix(&C, N, N);
    init_matrix(&C_seq, N, N);

    // Initialize random values for matrices A and B
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A.data[i * N + j].real = rand() % 100;
            A.data[i * N + j].imag = rand() % 100;
            B.data[i * N + j].real = rand() % 100;
            B.data[i * N + j].imag = rand() % 100;
        }
    }

    // Perform matrix multiplication sequentially
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C_seq.data[i * N + j].real = 0;
            C_seq.data[i * N + j].imag = 0;
            for (int k = 0; k < N; k++) {
                C_seq.data[i * N + j].real += A.data[i * N + k].real * B.data[k * N + j].real - A.data[i * N + k].imag * B.data[k * N + j].imag;
                C_seq.data[i * N + j].imag += A.data[i * N + k].real * B.data[k * N + j].imag + A.data[i * N + k].imag * B.data[k * N + j].real;
            }
        }
    }

    // Divide rows of matrix A among processes
    int rows_per_process = N / size;
    int rows_remaining = N % size;
    int start_row = rank * rows_per_process;
    int end_row = start_row + rows_per_process - 1;
    if (rank == size - 1) {
        end_row += rows_remaining;
    }

    // Allocate memory for local matrices
    matrix local_A, local_C;
    init_matrix(&local_A, rows_per_process, N);
    init_matrix(&local_C, rows_per_process, N);

    // Scatter rows of matrix A to processes
    MPI_Scatter(A.data, rows_per_process * N * sizeof(complex), MPI_BYTE, local_A.data, rows_per_process * N * sizeof(complex), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Broadcast matrix B to all processes
    MPI_Bcast(B.data, N * N * sizeof(complex), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Perform matrix multiplication
    for (int i = 0; i < rows_per_process; i++) {
        for (int j = 0; j < N; j++) {
            local_C.data[i * N + j].real = 0;
            local_C.data[i * N + j].imag = 0;
            for (int k = 0; k < N; k++) {
                local_C.data[i * N + j].real += local_A.data[i * N + k].real * B.data[k * N + j].real - local_A.data[i * N + k].imag * B.data[k * N + j].imag;
                local_C.data[i * N + j].imag += local_A.data[i * N + k].real * B.data[k * N + j].imag + local_A.data[i * N + k].imag * B.data[k * N + j].real;
            }
        }
    }

    // Gather rows of matrix C from processes
    MPI_Gather(local_C.data, rows_per_process * N * sizeof(complex), MPI_BYTE, C.data, rows_per_process * N * sizeof(complex), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Compare parallel and sequential results
    if (rank == 0) {
        double tolerance = 1e-4;
        if (compare_matrices(&C, &C_seq, tolerance)) {
            std::cout << "Results match!" << std::endl;
        } else {
            std::cout << "Results do not match!" << std::endl;
        }
        print_matrix(&C);
        std::cout << "----------------" << std::endl;

        print_matrix(&C_seq);
    }

    // Free memory
    free_matrix(&A);
    free_matrix(&B);
    free_matrix(&C);
    free_matrix(&C_seq);
    free_matrix(&local_A);
    free_matrix(&local_C);

    MPI_Finalize();
    return 0;
}
