#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]){
	MPI_Comm intercomm;
	int msg, tag, dest;
	int size, rank;
	char port_name[MPI_MAX_PORT_NAME];
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	printf("Trying to connect to %s\n", argv[1]);
	gets(port_name);
        //MPI_Lookup_name("port", MPI_INFO_NULL, port_name);
	MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm);
        msg=42;
	tag=0;
	dest=0;	
	MPI_Send(&msg, 1, MPI_INT, dest, tag, intercomm);
	MPI_Comm_disconnect(&intercomm);
	MPI_Finalize();



}

