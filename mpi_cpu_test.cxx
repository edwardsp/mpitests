#include <mpi.h>
#include <iostream>
#include <vector>

const int ITERATIONS = 1000000;
const int DATASIZE = 28000;

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::vector<double> data_a(DATASIZE, 1.0);
	double* a = &data_a[0];
	std::vector<double> data_b(DATASIZE, 1.0);
	double* b = &data_a[0];
	std::vector<double> data_c(DATASIZE, 0.0);
	double* c = &data_a[0];
	
	MPI_Barrier(MPI_COMM_WORLD);
	const double start_time = MPI_Wtime();
	for(int i = 0; i < ITERATIONS; ++i) {
		for(int n = 0; n < DATASIZE; ++n) {
			a[n] = (a[n] * b[n]) + c[n];
		}
	}
	const double end_time = MPI_Wtime();
	double duration = end_time - start_time;

	std::vector<double> timing(size, 0.0);
	MPI_Gather(&duration, 1, MPI_DOUBLE, &timing[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	if (rank == 0) {
		for (int r = 0; r < size; ++r) {
			std::cout << timing[r] << ' ';
		}
		std::cout << std::endl;
	}

	MPI_Finalize();
}


