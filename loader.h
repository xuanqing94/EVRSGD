#ifndef _LOADER_H_
#define _LOADER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

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


void loadData(const char* file, Data* data);

// clean memory
void rmData(Data* data);

// print
void printData(Data* data, int headN=-1);

#endif // _LOADER_H_
