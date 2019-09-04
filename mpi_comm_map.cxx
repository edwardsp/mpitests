#include <mpi.h>
#include <iostream>
#include <vector>

const int REPEATS = 1;
const int ITERATIONS = 1000;
const int DATASIZE = 1024;

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::vector<double> map(size*size, 0.0);
	std::vector<double> map_final(size*size, 0.0);
	std::vector<unsigned char> data_a(DATASIZE);
	std::vector<unsigned char> data_b(DATASIZE);
	
	for(int rep = 0; rep < REPEATS; ++rep) {
		for(int rank_src = 0; rank_src < size; ++rank_src) {
			for(int rank_dst = 0; rank_dst < size; ++rank_dst) {
				if (rank_src != rank_dst) {
					MPI_Barrier(MPI_COMM_WORLD);
					if(rank == rank_src) {
						double start_time = MPI_Wtime();
						for(int i = 0; i < ITERATIONS; ++i) {
							MPI_Send(&data_a[0], DATASIZE, MPI_BYTE, rank_dst, 101, MPI_COMM_WORLD);
							MPI_Recv(&data_b[0], DATASIZE, MPI_BYTE, rank_dst, 102, MPI_COMM_WORLD, 0);
						}
						double end_time = MPI_Wtime();
						map[rank_src*size + rank_dst] += end_time - start_time;
					} else if (rank == rank_dst) {
						for(int i = 0; i < ITERATIONS; ++i) {
							MPI_Recv(&data_a[0], DATASIZE, MPI_BYTE, rank_src, 101, MPI_COMM_WORLD, 0);
							MPI_Send(&data_b[0], DATASIZE, MPI_BYTE, rank_src, 102, MPI_COMM_WORLD);
						}
					}
				}

			}
		}
	}

	MPI_Reduce(&map[0], &map_final[0], size*size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		for(int rank_src = 0; rank_src < size; ++rank_src) {
			for(int rank_dst = 0; rank_dst < size; ++rank_dst) {
				std::cout << map_final[rank_src*size + rank_dst] << ' ';
			}
			std::cout << std::endl;
		}
	}

	MPI_Finalize();
}


