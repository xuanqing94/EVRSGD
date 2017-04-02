#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

void exit_with_help(){
 printf(
 "Usage: \n"   
 "options:\n"
 "-a algorithm: set type of solver (default 0)\n"
 " 0 -- EVRSGD\n"
 " 1 -- EASGD\n"
 " 2 -- EAMSGD\n"
 " 3 -- DOWNPOUR"
 "-r set parameter rho\n"
 "-s set stepsize eta\n"
 "-e minimum error: epsilon\n"
 "-f input data file name\n"
 "-d momentum parameter delta\n"
 );
 exit(1);
}

void exit_input_error(int line_num){
 fprintf(stderr, "Wrong input format at line %d\n", line_num);
 exit(1);
}

void parse_command_line(int argc, char **argv, char *input_file_name, char *test_file_name, double *eta, double *rho, double *delta, int *method_flag, double *eps, int rank){
 int i;
 for (i=1;i<argc;i++){
  if (argv[i][0] != '-') break;
  if (++i>=argc)
   exit_with_help();
  switch(argv[i-1][1]){
   case 'a':
    *method_flag = atoi(argv[i]);
    break;
   case 'r':
    *rho = atof(argv[i]);
    break;
   case 's':
    *eta = atof(argv[i]);
    break;
   case 'e':
     *eps = atof(argv[i]);
    break;
   case 'd':
     *delta = atof(argv[i]);
   case 'f':
     strcpy(input_file_name, argv[i]);
    break;
   default :
    if (rank == 0)
     fprintf(stderr, "unknown option: -%c\n", argv[i-1][1]);
    exit_with_help();
    break;
  }
 }
}
