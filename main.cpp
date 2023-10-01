#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <sstream>

using namespace std;

void merge(vector<int>& arr, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;

    vector<int> L(n1), R(n2);

    for (int i = 0; i < n1; i++) {
        L[i] = arr[l + i];
    }
    for (int j = 0; j < n2; j++) {
        R[j] = arr[m + 1 + j];
    }

    int i = 0, j = 0, k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void mergeSort(vector<int>& arr, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;

        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        merge(arr, l, m, r);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::stringstream s1(argv[1]);
    int N;
    s1 >> N;
    const int chunkSize = N / size;

    vector<int> arr(N);

    if (rank == 0) {
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            arr[i] = rand() % 100000;
        }
    }

    vector<int> chunk(chunkSize);
    MPI_Scatter(&arr[0], chunkSize, MPI_INT, &chunk[0], chunkSize, MPI_INT, 0, MPI_COMM_WORLD);

    mergeSort(chunk, 0, chunkSize - 1);

    vector<int> sorted(N);
    MPI_Gather(&chunk[0], chunkSize, MPI_INT, &sorted[0], chunkSize, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        mergeSort(sorted, 0, N - 1);

        cout << "Sorted array: ";
        for (int i = 0; i < N; i++) {
            cout << sorted[i] << " ";
        }
        cout << endl;
    }

    MPI_Finalize();
    return 0;
}
