#include <mpi.h>
#include <iostream>
#include <sstream>
#include <vector>

void star(int N, int M) {

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Создаем массив данных, который будет разделен между процессами
    std::vector<int> data(N * M);
    if (rank == 0) {
        for (int i = 0; i < N * M; i++) {
            data[i] = i;
        }
    }

    // Определяем размеры блоков для каждого процесса
    std::vector<int> counts(size);
    std::vector<int> displs(size);
    int block_size = N * M / size;
    for (int i = 0; i < size; i++) {
        counts[i] = block_size;
        displs[i] = i * block_size;
    }
    counts[size - 1] = N * M - (size - 1) * block_size;

    // Разделяем массив данных между процессами
    std::vector<int> local_data(block_size);
    MPI_Scatterv(data.data(), counts.data(), displs.data(), MPI_INT,
                 local_data.data(), block_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Обрабатываем блок данных
    for (int i = 0; i < block_size; i++) {
        local_data[i] = rank;
    }

    // Отправляем результаты обратно процессу с рангом 0
    std::vector<int> result(N * M);
    MPI_Gatherv(local_data.data(), block_size, MPI_INT,
                result.data(), counts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Выводим результаты
        for (int i = 0; i < N * M; i++) {
            std::cout << "Process 0 received message " << result[i] << std::endl;
        }
    }

    MPI_Finalize();
}

int main(int argc, char** argv) {
    std::stringstream s1(argv[1]);
    int N;
    int M;
    s1 >> N;
    std::stringstream s2(argv[2]);
    s2 >> M;

    MPI_Init(&argc, &argv);
    star(N, M);
    return 0;
}
