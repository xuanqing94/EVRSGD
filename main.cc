#include <mpi.h>
#include <stdio.h>


typedef struct Elements_ {
	int colIdx;
	double value;
} Elements;

typedef struct DataRow_ {
	int nLength;
	Elements *elements;
} DataRow;

typedef struct Data_ {
	int nCols;
	int nRows;
	DataRow *dataRows;
};

// load data from file, make sure data = NULL
void loadData(const char* file, Data* data) {

}

// clean memory
void rmData(Data* data) {

}

// server side
void server() {

}


// client side
void client(int clientId) {

}


int main(int argc, char** argv) {
  int rank, size;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	printf("I am %d of %d\n", rank + 1, size);

	MPI_Finalize();
	
	return 0;
}
