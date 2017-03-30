#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include <stdatomic.h>
#include <string.h>
#include "command_line.h"

#define MAX_TAG_DIGITS 16
#define MAX_IDX_DIGITS 16
#define MAX_VAL_DIGITS 16

#define ROW_TAG 1
#define ELE_TAG 2
#define NCOL_TAG 3
#define NROW_TAG 4

typedef struct Elements_ {
	long colIdx;       // column index of that value
	double value;      // value in that element
} Element;

typedef struct DataRow_ {
	int output;        // y[i]
	long nLength;      // nnz in that row
	Element *elements; // data
} DataRow;

typedef struct Data_ {
	long nCols;        // number of columns
	long nRows;        // number of rows
	DataRow *dataRows; // data
} Data;

// load data from file, make sure data = NULL
void loadData(const char* file, Data* data) {
  FILE *fHandle = fopen(file, "r");
  if (fHandle == NULL) {
    fprintf(stderr, "Fail to open file: %s\n", file);
    data = NULL;
    return;
  }
  // read whole file, then parse it
  fseek(fHandle, 0, SEEK_END);
  unsigned long fSize = ftell(fHandle);
  fseek(fHandle, 0, SEEK_SET);
  char *content = (char*)malloc(fSize * sizeof(char));
  size_t count = fread(content, sizeof(char), fSize, fHandle);
  assert((count == fSize));

  long nSample = 0;
  long nFeature = 0;
  unsigned long nElem = 0;
  DataRow *dataRow = (DataRow*)malloc(sizeof(DataRow) * 1);
  Element *elemBuffer = (Element*)malloc(sizeof(Element) * 1); // we estimate #elem by fSize / 16.
  long rowBufferSize = 1; //1024
  long elemBufferSize = 1; // fSize >> 4
  int firstOutput = 0;
  for (unsigned long i = 0; i < fSize; ++i) {
    // allocate one row
    if (nSample == rowBufferSize) {
      // enlarge by 2 times
      rowBufferSize <<= 1;
      dataRow = (DataRow*)realloc(dataRow, sizeof(DataRow) * rowBufferSize);
    }
    DataRow *thisRow = dataRow + nSample;
    nSample++;
    // read class 
    char tag[MAX_TAG_DIGITS + 1] = {0};
    for(int tag_idx = 0; content[i] != ' ' && tag_idx < MAX_TAG_DIGITS; ++i, ++tag_idx) {
      tag[tag_idx] = content[i];
    }
    if (firstOutput == 0) {
      thisRow->output = 1;
      firstOutput = atoi(tag);
    }
    else {
      if (atoi(tag) == firstOutput) thisRow->output = 1;
      else thisRow->output = -1;
    }
    while (content[i] == ' ') i++;
    thisRow->nLength = 0;
    // read elements
    for(; content[i] != '\n';) {
        if (nElem == elemBufferSize) {
          elemBufferSize <<= 1;
          elemBuffer = (Element*)realloc(elemBuffer, sizeof(Element) * elemBufferSize);
        }
        Element *thisElem = elemBuffer + nElem;
        char indexStr[MAX_IDX_DIGITS + 1] = {0};
        char valStr[MAX_VAL_DIGITS + 1] = {0};
        for (int idx = 0; content[i] != ':'; ++i, ++idx) {
          indexStr[idx] = content[i];
        }
        i += 1;
        for (int idx = 0; content[i] != ' '; ++i, ++idx) {
          valStr[idx] = content[i];
        }
        while (content[i] == ' ') i++;
        thisElem->colIdx = atoi(indexStr) - 1;
        if (thisElem->colIdx + 1 > nFeature) 
          nFeature = thisElem->colIdx + 1;
        thisElem->value = atof(valStr);
        nElem += 1;
        thisRow->nLength += 1;
    }
  }
  // link elemBuffer with dataRow
  for (int i = 0; i < nSample; ++i) {
    dataRow[i].elements = elemBuffer;
    elemBuffer += dataRow[i].nLength;
  }
  data->dataRows = dataRow;
  data->nRows = nSample;
  data->nCols = nFeature;
  free(content);
  fclose(fHandle);
}

// clean memory
void rmData(Data* data) {
  free(data->dataRows->elements);
  free(data->dataRows);
}

// print
void printData(Data* data) {
  for (int i = 0; i < data->nRows; ++i) {
	  DataRow *dataRow = data->dataRows + i;
		printf("%d, length: %ld", dataRow->output, dataRow->nLength);
	  for (int j = 0; j < dataRow->nLength; ++j) {
		  printf("%ld:%lf ", (dataRow->elements[j]).colIdx, (dataRow->elements[j]).value);
		}
		printf("\n");
	}
}

// from 0 to 1, 2, ..., nClients
//XXX should we use random assignment or scatter?
void distributeData(Data* data, int nClients) {
  int kthRow = 0;
	int nextReceiver = 0;
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
	printf("Finished sending\n");
	// broadcast end of data signal
  DataRow endOfData;
	endOfData.nLength = 0;
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
	long elemCapacity = 1;
	long elemCount = 0;
	while (true) {
	  if (rowCount == rowCapacity) {
		  rowCapacity <<= 1;
		  dataRow = (DataRow*)realloc(dataRow, sizeof(DataRow) * rowCapacity);
		}
		DataRow *rowToRecv = dataRow + rowCount;
	  MPI_Recv(rowToRecv, sizeof(DataRow), MPI_BYTE, 0, ROW_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (rowToRecv->nLength == 0) {
		  rowCount--;
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
	data->dataRows = dataRow;
}


// struct for thread arguments
typedef struct Arg_ {
  Data *data;
	double *z;
	pthread_mutex_t *lock_z;
	int interval;          // seconds every call
  double eps;            // minimal error
	bool *stop;            // should stop or not
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

void functionVal(void* args) {
  Arg *arg = args;
	double *z_loc = (double*)calloc(sizeof(double) * arg->data->nCols);
	
	double last_error = functionVal_impl(arg->data, z_loc);
	printf("Loss function: %lf\n", last_error);
	while (true) {
	  sleep(args->interval);
		// acquire lock
		pthread_mutex_lock(args->lock_z);
    // copy z
		memcpy(z_loc, z, sizeof(double) * arg->data->nCols);
		// release lock
		pthread_mutex_unlock(args->lock_z);
		// calculate the new value
    double this_error = functionVal_impl(arg->data, z_loc);
    printf("Loss function: %lf\n", this_error);
		if (this_error - last_error > -eps && this_error < last_error) {
		  arg->stop = true;
			break;
		}
		else {
      this_error = last_error;
		}
	}
}

// server side
void server(Data* data, int nClients, double eta, double rho) {
	distributeData(data, nClients);
	int d = data -> nCols;
	double *bufferW = (double*)malloc(sizeof(double) * d * nClients);
	double *z = (double*)malloc(sizeof(double) * d);
	// initialize z,W
	for (int i=0; i < d; i++ )
		z[i]=0.0;
	for (int i=0; i < d * nClients; i++)
		bufferW[i]=0.0;
	double *sum_bufferW = (double*)malloc(sizeof(double)*d);
	for (int i=0; i<d; i++)
		sum_bufferW[i] = 0.0;
	while (1){
		double *wk = (double*)malloc(sizeof(double) * (d+1));
		MPI_Recv(wk, d+1, MPI_DOUBLE, MPI_ANY_SOURCE , MPI_ANY_TAG , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int cur_k = (int) wk[d];	
		for (int i=0; i <d; i++){
			z[i]= (1-eta*rho*nClients)*z[i] + eta*rho*nClients*(wk[i]-bufferW[(cur_k-1)*d+i]+1/nClients*sum_bufferW[i]);
			sum_bufferW[i] += wk[i]-bufferW[(cur_k-1)*d+i];
			bufferW[(cur_k-1)*d+i]=wk[i];
		}
		MPI_Send(z, d, MPI_DOUBLE, cur_k, MPI_ANY_TAG, MPI_COMM_WORLD);

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
		  double coef = -1.0 / (exp(innerProd) + 1.0);
      for (int j = 0; j < thisRow->nLength; ++j) {
			  long idx = thisRow->elements[j].colIdx;
				gradOut[idx] += thisRow->elements[j].value * coef;
			}
		}
		else { // thisRow->output == 1
		  double coef = 1.0 - 1.0 / (exp(innerProd) + 1.0);
			for (int j = 0; j < thisRow->nLength; ++j) {
			  long idx = thisRow->elements[j].colIdx;
				gradOut[idx] += thisRow->elements[j].value * coef;
			}
		}
	}
}

// client side
void client(int clientId, double rho, double eta) {
	Data data;
	//printData(&data);
  	collectData(&data);
	int d = (&data) -> nCols;
	double *wk = (double*)malloc(sizeof(double) * (d+1));
	double *z = (double*)malloc(sizeof(double) * d);
	double *gradwk = (double*)malloc(sizeof(double) * d);
	for (int i=0; i<d; i++){
		wk[i]=0.0;
		z[i]=0.0;
	}
	while (1){
	MPI_Recv(z, d, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	grad(wk, gradwk, &data);
	for (int i=0; i<d; i++)
		wk[i] = wk[i] - eta * (gradwk[i]+rho*(wk[i]-z[i]));
	//gradient
	wk[d]= (double) clientId;
	MPI_Send(wk, d+1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
	}
	rmData(&data);
}


int main(int argc, char** argv) {
  int rank, size;
  char input_file_name[1024];
  char test_file_name[1024];
  double eta;
	double rho;
	int *method_flag;
	parse_command_line(argc, argv, input_file_name, test_file_name, eta, rho, method_flag, rank); 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0) { // for server
	  // load data
		Data data;
		loadData("../test_realsim", &data);
    printf("nRows: %ld, nFeature: %ld\n", data.nRows, data.nCols);
    bool stop = false;
		double eps = 1.0e-6;
		interval = 1;
    pthread_mutex_t lock_z;
		pthread_mutex_init(&lock_z, NULL);
    double *z = (double*)calloc(sizeof(double) * data->nCols);
		// function arguments
		Arg arg;
		arg.data = data;
		arg.z = z;
		arg.lock_z = &lock_z;
		arg.interval = interval;
		arg.eps = eps;
		arg.stop = &stop;

		pthread_t monitor_thread;
		pthread_create(&monitor, NULL, (void*)(functionVal*)(void*), &arg);
		server(&data, size - 1, *rho, *eta);

		pthread_join(monitor_thread, NULL);
		rmData(&data);
	}
	else { // for clients
		client(rank, *rho, *eta);
	}
	MPI_Finalize();
	return 0;
}
