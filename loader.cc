#include "loader.h"

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
void printData(Data* data, int headN) {
  if (headN == -1) headN = data->nRows;
	if (headN < 0 || headN > data->nRows) {
	  fprintf(stderr, "Invalid headN parameter: %d\n", headN);
		exit(-1);
	}
  for (int i = 0; i < headN; ++i) {
	  DataRow *dataRow = data->dataRows + i;
		printf("%d, length: %ld\n", dataRow->output, dataRow->nLength);
	  for (int j = 0; j < dataRow->nLength; ++j) {
		  printf("%ld:%lf ", (dataRow->elements[j]).colIdx, (dataRow->elements[j]).value);
		}
		printf("\n");
	}
}

