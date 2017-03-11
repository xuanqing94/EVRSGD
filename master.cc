#include <mpi.h>
#include <stdio.h>
int main(int argc, char **argv)
{
	MPI_Comm intercomm;
	MPI_Status status;
	char port_name[MPI_MAX_PORT_NAME];
	int size,rank, msg;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Open_port(MPI_INFO_NULL, port_name);
	//MPI_Publish_name("port", MPI_INFO_NULL, port_name);
	printf("server ready at %s\n", port_name);
	MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm);
	printf("client connected\n");
	MPI_Unpublish_name("port", MPI_INFO_NULL, port_name);
	MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, intercomm, &status);
	printf("msg: %d\n", msg);
	MPI_Comm_free(&intercomm);
	MPI_Close_port(port_name);
	MPI_Finalize();

}

