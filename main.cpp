#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Создание первой группы процессов
    if (rank > 0 && rank < size/2) {
        // Обмен информацией по кольцу при помощи сдвига в одномерной декартовой топологии
        int left = (rank - 1 + size/2) % (size/2);
        int right = (rank + 1) % (size/2);
        int data = rank;
        for (int i = 0; i < size/2; i++) {
            MPI_Sendrecv(&data, 1, MPI_INT, right, 0, &data, 1, MPI_INT, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "Rank " << rank << " received data " << data << " from rank " << left << std::endl;
            MPI_Send(&data, 1, MPI_INT, right, 0, MPI_COMM_WORLD);
        }
    }
        // Создание второй группы процессов
    else {

        // Общение осуществляется по схеме master-slave (звезда) на базе топологии графа
        int master = 0;
        int data;
        if (rank == master) {

            data = 42;
        }
        MPI_Bcast(&data, 1, MPI_INT, master, MPI_COMM_WORLD);
        std::cout << "Rank " << rank << " received data " << data << " from master" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
