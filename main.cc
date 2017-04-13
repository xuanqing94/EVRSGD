#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include <atomic>
#include <string.h>

#include "command_line.h"
#include "loader.h"

#define MAX_TAG_DIGITS 16
#define MAX_IDX_DIGITS 16
#define MAX_VAL_DIGITS 16

#define ROW_TAG 1
#define ELE_TAG 2
#define NCOL_TAG 3
#define NROW_TAG 4
#define W_TAG 5
#define Z_TAG  6
#define GRADW_TAG 7

#define max(a,b) ((a) > (b)? (a):(b))

// from 0 to 1, 2, ..., nClients
//XXX should we use random assignment or scatter?
void distributeData(Data* data, int nClients) {
	int kthRow = 0;
	int nextReceiver = 0;
	MPI_Bcast(&data->nCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	while (kthRow != data->nRows) {
		nextReceiver = rand() % nClients + 1;
		DataRow *rowToSend = data->dataRows + kthRow;
		// send dataRow basic infomation, excluding elements
		if (MPI_Send(rowToSend, sizeof(DataRow), MPI_BYTE, 
					nextReceiver, ROW_TAG, MPI_COMM_WORLD) != MPI_SUCCESS) {
			fprintf(stderr, "Fail to send data row\n");
			exit(-1);
		}
		// send elements in that datarow
		Element *elements = rowToSend->elements;
		if (MPI_Send(elements, sizeof(Element) * rowToSend->nLength, MPI_BYTE, 
					nextReceiver, ELE_TAG, MPI_COMM_WORLD) != MPI_SUCCESS) {
			fprintf(stderr, "Fail to send elements\n");
			exit(-1);
		}
		kthRow++;
	}
	// broadcast end of data signal
	DataRow endOfData;
	endOfData.nLength = -1;
	for (int i = 1; i <= nClients; ++i) {
		if (MPI_Send(&endOfData, sizeof(DataRow), MPI_BYTE, i, ROW_TAG, MPI_COMM_WORLD) != MPI_SUCCESS) {
			fprintf(stderr, "Fail to send end-of-data signal\n");
			exit(-1);
		}
	}
}

void collectData(Data* data) {
	DataRow *dataRow = (DataRow*)malloc(sizeof(DataRow) * 1);
	Element *elementBuff = (Element*)malloc(sizeof(Element) * 1);
	long rowCapacity = 1;
	long rowCount = 0;
	long colCount = 0;
	long elemCapacity = 1;
	long elemCount = 0;
	MPI_Bcast(&colCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	while (true) {
		if (rowCount == rowCapacity) {
			rowCapacity <<= 1;
			dataRow = (DataRow*)realloc(dataRow, sizeof(DataRow) * rowCapacity);
		}
		DataRow *rowToRecv = dataRow + rowCount;
		MPI_Recv(rowToRecv, sizeof(DataRow), MPI_BYTE, 0, ROW_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (rowToRecv->nLength == -1) {
			break;
		}
		if (elemCount + rowToRecv->nLength >= elemCapacity) {
			elemCapacity += rowToRecv->nLength;
			elemCapacity <<= 1;
			elementBuff = (Element*)realloc(elementBuff, sizeof(Element) * elemCapacity);
		}
		Element *elemToRecv = elementBuff + elemCount;
		MPI_Recv(elemToRecv, sizeof(Element) * rowToRecv->nLength, MPI_BYTE, 0, ELE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		rowCount++;
		elemCount += rowToRecv->nLength;
	}
	long elemIdx = 0;
	for (int i = 0; i < rowCount; ++i) {
		dataRow[i].elements = elementBuff + elemIdx;
		elemIdx += dataRow[i].nLength;
	}
	data->nRows = rowCount;
	data->nCols = colCount;
	data->dataRows = dataRow;
}


// struct for thread arguments
typedef struct Arg_ {
	Data *data;
	double *z;
	pthread_mutex_t *lock_z;
	int interval;          // seconds every call
	double eps;            // minimal error
	std::atomic<bool> *stop;            // should stop or not
} Arg;

// Calculate the loss function given z
double functionVal_impl(Data* data, double* z) {
	double loss = 0.0;
	for (int i = 0; i < data->nRows; ++i) {
		DataRow *thisRow = data->dataRows + i;
		double innerProd = 0.0;
		// calculate w^Tx_i
		for (int j = 0; j < thisRow->nLength; ++j) {
			long colIdx = thisRow->elements[j].colIdx;
			innerProd += z[colIdx] * thisRow->elements[j].value;
		}
		if (thisRow->output > 0) {
			loss += log(1.0 + exp(-innerProd));
		}
		else {
			loss += log(1.0 + exp(innerProd));
		}
	}
	return loss;
}

void* functionVal(void* args) {
	Arg *arg = (Arg*)args;
	double *z_loc = (double*)calloc(arg->data->nCols, sizeof(double));
	double last_error = functionVal_impl(arg->data, z_loc);
	printf("Loss function: %lf\n", last_error);
	int counter = 0;
	while (true) {
		sleep(arg->interval);
		// acquire lock
		pthread_mutex_lock(arg->lock_z);
		// copy z
		memcpy(z_loc, arg->z, sizeof(double) * arg->data->nCols);
		// release lock
		pthread_mutex_unlock(arg->lock_z);
		// calculate the new value
		double this_error = functionVal_impl(arg->data, z_loc);
		printf("Time: %d, Loss function: %lf\n", ++counter, this_error);
		if (this_error - last_error > -arg->eps && this_error < last_error) {
			printf("Desired error attained, send stop signal\n");
			arg->stop->store(true);
			break;
		}
		else {
			this_error = last_error;
		}
	}
}

// calculate the gradient, suppose the original problem is logistic regression
// make sure len(w) == len(gradOut) == data->nCols
void grad(double *w, double* gradOut, Data* data) {
	memset(gradOut, 0, sizeof(double) * data->nCols);
	for (int i = 0; i < data->nRows; ++i) {
		DataRow *thisRow = data->dataRows + i;
		// calculate w^Tx
		double innerProd = 0.0;
		for (int j = 0; j < thisRow->nLength; ++j) {
			long idx = thisRow->elements[j].colIdx;
			innerProd += w[idx] * thisRow->elements[j].value;
		}
		if (thisRow->output == -1) {
			double coef = 1.0 - 1.0 / (exp(innerProd) + 1.0);
			for (int j = 0; j < thisRow->nLength; ++j) {
				long idx = thisRow->elements[j].colIdx;
				gradOut[idx] += thisRow->elements[j].value * coef;
			}
		}
		else { // thisRow->output == 1
			double coef = -1.0 / (exp(innerProd) + 1.0);
			for (int j = 0; j < thisRow->nLength; ++j) {
				long idx = thisRow->elements[j].colIdx;
				gradOut[idx] += thisRow->elements[j].value * coef;
			}
		}
	}
}

// server side for EVRSGD
void server_evrsgd(Data* data, int nClients, double eta, double rho, Arg* arg) {
	distributeData(data, nClients);
	printf("Data sent by server: %ld rows, %ld columns\n", data->nRows, data->nCols);
	int d = data -> nCols;
	double *bufferW = (double*)calloc(d * nClients, sizeof(double));
	double *sum_bufferW = (double*)calloc(d, sizeof(double));
	double *wk = (double*)calloc(d + 1, sizeof(double));
	while (!arg->stop->load()){
		MPI_Recv(wk, d + 1, MPI_DOUBLE, MPI_ANY_SOURCE, W_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int cur_k = wk[d];
		pthread_mutex_lock(arg->lock_z);
		for (int i = 0; i < d; i++){
			// update z
			arg->z[i] = (1 - eta * rho * nClients) * arg->z[i] + eta * rho * nClients * 
				(wk[i] - bufferW[(cur_k - 1) * d + i] + 1.0 / nClients * sum_bufferW[i]);
			// update summation of w
			sum_bufferW[i] += wk[i] - bufferW[(cur_k-1) * d + i];
			// update buffer W
			bufferW[(cur_k - 1) * d + i] = wk[i];
		}
		pthread_mutex_unlock(arg->lock_z);
		MPI_Send(arg->z, d, MPI_DOUBLE, cur_k, Z_TAG, MPI_COMM_WORLD);
	}
	printf("Server stopped.");
}

// client evrsgd side
void client_evrsgd(int clientId, double rho, double eta, double delta) {
	Data data;
	collectData(&data);
	printf("Data received by client: %ld rows, %ld columns\n", data.nRows, data.nCols);
	int d = data.nCols;
	double *wk = (double*)calloc(d + 1, sizeof(double));
	double *z = (double*)calloc(d, sizeof(double));
	double *gradwk = (double*)malloc(sizeof(double) * d);
	double *v = (double*)calloc(d, sizeof(double));
	time_t now;
	srand(clientId + (unsigned int)time(&now)); // set random seed
	while (1){
		// choose number of iterations randomly
		int n_iter = 10; //rand() % (max_iter - min_iter) + min_iter;
		for (int inner = 0; inner < n_iter; ++inner) {
			grad(wk, gradwk, &data);
			for (int i=0; i<d; i++){
				v[i] = delta * v[i] - eta *gradwk[i];
				wk[i] = wk[i] + v[i] - eta * rho * (wk[i] - z[i]);
			}
		}
		// send client ID to server
		wk[d]= (double) clientId;
		MPI_Send(wk, d + 1, MPI_DOUBLE, 0, W_TAG, MPI_COMM_WORLD);
		MPI_Recv(z, d, MPI_DOUBLE, 0, Z_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	rmData(&data);
}

// server side for EASGD
void server_easgd(Data* data, int nClients, double eta, double rho, Arg* arg) {
	distributeData(data, nClients);
	printf("Data sent by server: %ld rows, %ld columns\n", data->nRows, data->nCols);
	int d = data -> nCols;
	double *bufferW = (double*)calloc(d * nClients, sizeof(double));
	double *sum_bufferW = (double*)calloc(d, sizeof(double));
	double *wk = (double*)calloc(d + 1, sizeof(double));
	while (!arg->stop->load()){
		MPI_Recv(wk, d + 1, MPI_DOUBLE, MPI_ANY_SOURCE, W_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int cur_k = wk[d];
		pthread_mutex_lock(arg->lock_z);
		for (int i = 0; i < d; i++){
			// update summation of w
			sum_bufferW[i] += wk[i] - bufferW[(cur_k-1) * d + i];
			// update z
			arg->z[i] = (1 - eta * rho * nClients) * arg->z[i] + eta * rho * sum_bufferW[i];
			// update buffer W
			bufferW[(cur_k - 1) * d + i] = wk[i];
		}
		pthread_mutex_unlock(arg->lock_z);
		MPI_Send(arg->z, d, MPI_DOUBLE, cur_k, Z_TAG, MPI_COMM_WORLD);
	}
	printf("Server stopped.");
}

// server side for Downpour SGD
/* void server_dpsgd(Data* data, int nClients, double eta, double rho) {
	distributeData(data, nClients);
	printf("Data sent by server: %ld rows, %ld columns\n", data->nRows, data->nCols);
	int d = data -> nCols;
	double *bufferW = (double*)calloc(d, sizeof(double));
	double *buffergradW = (double*)calloc(d, sizeof(double));
	int block_size= max(d/nClients + 1, d-(nClients-1)*(d/nClients)+1);
	double *gradwk = (double*)calloc(block_size, sizeof(double));
	while (1){
		MPI_Recv(gradwk, block_size, MPI_DOUBLE, MPI_ANY_SOURCE, GRADW_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // have some problem with machine index
		int cur_k = gradwk[d/nClients];
		if (cur_k != nClients){
			for (int i = 0; i < d/nClients; i++){
				// update buffer W
				buffergradW[(cur_k - 1) * (d/nClients) + i] = wk[i];
			}
		}
		else{
			for (int i=0; i < (d-(nClients-1)*(d/nClients)); i++)
				buffergradW[(cur_k - 1) * (d/nClients) + i] = wk[i];
			//MPI_Send(bufferW, d, MPI_DOUBLE, cur_k, Z_TAG, MPI_COMM_WORLD);
		}
		for (int i=0; i <d; i++)
			bufferW[i] -= eta * buffergradW[i];

		MPI_Send(bufferW, d, MPI_DOUBLE, cur_k, Z_TAG, MPI_COMM_WORLD);
			if (cur_k != nClients){
			MPI_Send(bufferW+(cur_k-1)*d/nClient, d/nClient, MPI_DOUBLE, cur_k, Z_TAG, MPI_COMM_WORLD);
			}
			else {
			MPI_Send(bufferW+(cur_k-1)*d/nClient, d-(nClients-1)*(d/nClient), MPI_DOUBLE, cur_k, Z_TAG, MPI_COMM_WORLD);
			}
			}
		printf("Server stopped.");
	}
}
*/

// server side for downpour
void server_dpsgd(Data* data, int nClients, double eta, double rho, double delta, Arg* arg) {
	distributeData(data, nClients);
	printf("Data sent by server: %ld rows, %ld columns\n", data->nRows, data->nCols);
	int d = data -> nCols;
	double *v = (double*)calloc(d, sizeof(double));
	double *gradwk = (double*)calloc(d + 1, sizeof(double));
	while (!arg->stop->load()){
		MPI_Recv(gradwk, d + 1, MPI_DOUBLE, MPI_ANY_SOURCE, GRADW_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int cur_k = gradwk[d];
		pthread_mutex_lock(arg->lock_z);
		for (int i = 0; i < d; i++){
			// update z
			v[i] = delta * v[i] - eta * gradwk[i];
			arg->z[i] += v[i];
		}
		pthread_mutex_unlock(arg->lock_z);
		MPI_Send(arg->z, d, MPI_DOUBLE, cur_k, Z_TAG, MPI_COMM_WORLD);
	}
	printf("Server stopped.");
}


// client for eamsgd and easgd
void client_eamsgd(int clientId, double rho, double eta, double delta) {
	Data data;
	collectData(&data);
	printf("Data received by client: %ld rows, %ld columns\n", data.nRows, data.nCols);
	int d = data.nCols;
	double *wk = (double*)calloc(d + 1, sizeof(double));
	double *z = (double*)calloc(d, sizeof(double));
	double *gradwk = (double*)malloc(sizeof(double) * d);
	double *v = (double*)calloc(d, sizeof(double));
	time_t now;
	srand(clientId + (unsigned int)time(&now)); // set random seed
	while (1){
		// choose number of iterations randomly
		int n_iter = 1; //rand() % (max_iter - min_iter) + min_iter;
		for (int inner = 0; inner < n_iter; ++inner) {
			grad(wk, gradwk, &data);
			for (int i=0; i<d; i++){
				v[i] = delta * v[i] - eta * gradwk[i];
				wk[i] = wk[i] + v[i] - eta * rho * (wk[i] - z[i]);

			}
		}
		// send client ID to server
		wk[d]= (double) clientId;
		MPI_Send(wk, d + 1, MPI_DOUBLE, 0, W_TAG, MPI_COMM_WORLD);
		MPI_Recv(z, d, MPI_DOUBLE, 0, Z_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	rmData(&data);
}

/*
// client downpour sgd
void client_dpsgd(int clientId,  int size, double eta, double delta) {
	Data data;
	collectData(&data);
	printf("Data received by client: %ld rows, %ld columns\n", data.nRows, data.nCols);
	int d = data.nCols;
	double *wk = (double*)calloc(d + 1, sizeof(double));
	double *z = (double*)calloc(d, sizeof(double));
	double *gradwk = (double*)malloc(sizeof(double) * d);
	double *v = (double*)calloc(d, sizeof(double));
//	time_t now;
//	srand(clientId + (unsigned int)time(&now)); // set random seed
	while (1){
		// choose number of iterations randomly
		int n_iter = 1; //rand() % (max_iter - min_iter) + min_iter;
		for (int inner = 0; inner < n_iter; ++inner) {
			grad(wk, gradwk, &data);
		}
		// send client ID to server
		gradwk[d]= (double) clientId;
		if (clientId != size){
			MPI_Send(gradwk+(clientId-1)*(d/size), d/size, MPI_DOUBLE, 0, GRADW_TAG, MPI_COMM_WORLD);
		}
		else {
			MPI_Send(gradwk+(clientId-1)*(d/size), d- (clientId-1)*(d/size), MPI_DOUBLE, 0, GRADW_TAG, MPI_COMM_WORLD);
		}
		MPI_Recv(wk, d, MPI_DOUBLE, 0, Z_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	rmData(&data);
}
*/

void client_dpsgd(int clientId, double rho, double eta) {
	Data data;
	collectData(&data);
	printf("Data received by client: %ld rows, %ld columns\n", data.nRows, data.nCols);
	int d = data.nCols;
	double *z = (double*)calloc(d, sizeof(double));
	double *gradwk = (double*)malloc(sizeof(double) * (d + 1));
	while (1) {
		// choose number of iterations randomly
		int n_iter = 1; //rand() % (max_iter - min_iter) + min_iter;
		for (int inner = 0; inner < n_iter; ++inner) {
			grad(z, gradwk, &data);
		}
		// send client ID to server
		gradwk[d]= (double) clientId;
		MPI_Send(gradwk, d + 1, MPI_DOUBLE, 0, GRADW_TAG, MPI_COMM_WORLD);
		MPI_Recv(z, d, MPI_DOUBLE, 0, Z_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	rmData(&data);
}

int main(int argc, char** argv) {
	int rank, size;
	char input_file_name[1024] = {0};
	char test_file_name[1024] = {0};
	double eta;
	double rho;
	double eps = 1.0e-5;
	double delta;
	int method_flag;
	parse_command_line(argc, argv, input_file_name, test_file_name, &eta, &rho, &delta, &method_flag, &eps, rank); 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) { // for server
		// load data
		Data data;
		loadData(input_file_name, &data);
		printf("nRows: %ld, nFeature: %ld\n", data.nRows, data.nCols);

		std::atomic<bool> stop(false);
		int interval = 1;
		pthread_mutex_t lock_z;
		pthread_mutex_init(&lock_z, NULL);
		double *z = (double*)calloc(data.nCols, sizeof(double));
		// function arguments
		Arg arg;
		arg.data = &data;
		arg.z = z;
		arg.lock_z = &lock_z;
		arg.interval = interval;
		arg.eps = eps;
		arg.stop = &stop;

		pthread_t monitor_thread;
		pthread_create(&monitor_thread, NULL, &functionVal, &arg);
		if (method_flag == 0) {
			printf("Running EVRSGD method");
			server_evrsgd(&data, size - 1, eta, rho, &arg);
		}
		else if (method_flag == 1) {
			printf("Running EASGD method");
			server_easgd(&data, size - 1, eta, rho, &arg);
		}
		else if (method_flag == 2) {
			printf("Running EAMSGD method");
			server_easgd(&data, size -1,eta, rho, &arg);
		}
		else if (method_flag == 3) {
			printf("Running DOWNPOUR method");
			server_dpsgd(&data, size -1, eta, rho, 0, &arg);
		}
		else {
			fprintf(stderr, "Invalid method flag\n");
			exit(-1);
		}
		pthread_join(monitor_thread, NULL);
		rmData(&data);
	}
	else { // for clients
		if (method_flag == 0) {
			client_evrsgd(rank, rho, eta, 0);
		}
		else if (method_flag == 1) {
			client_eamsgd(rank, rho, eta, 0);
		}
		else if (method_flag == 2) {
			client_eamsgd(rank, rho, eta, delta);
		}
		else if (method_flag == 3) {
			client_dpsgd(rank, rho, eta);
		}

	}
	MPI_Finalize();
	return 0;
}
