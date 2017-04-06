#ifndef _COMMAND_LINE_H_
#define _COMMAND_LINE_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <mpi.h>

void exit_with_help();

void exit_input_error(int line_num);

void parse_command_line(int argc, char **argv, char *input_file_name, char *test_file_name, 
  double *eta, double *rho, double *delta, int *method_flag, double* eps, int rank);


#endif // _COMMAND_LINE_H_
