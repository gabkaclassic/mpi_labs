#include <mpi.h>
#include <iostream>
#include <sstream>
#include <vector>

void star(int N, int M) {

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != N + 1) {
        std::cerr << "Error: expected " << N + 1 << " processes, but got " << size << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::vector<int> messages(M);

    if (rank == 0) {
        for (int i = 0; i < M; i++) {
            messages[i] = i;
        }
        MPI_Bcast(messages.data(), M, MPI_INT, 0, MPI_COMM_WORLD);

    }


    if (rank != 0) {
        for (int i = 0; i < M; i++) {
            std::cout << "Process " << rank << " received message " << messages[i] << std::endl;
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
    star(N-1, M);
    return 0;
}
