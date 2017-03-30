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
	"-s type: set type of solver (default 0)\n"
	"	0 -- EVRSGD\n"
	"	1 -- EASGD\n"
	"	2 -- ESGD\n"
	"-r set parameter rho\n"
	"-e set stepsize eta\n"
	);
	exit(1);
}

void exit_input_error(int line_num){
	fprintf(stderr, "Wrong input format at line %d\n", line_num);
	exit(1);
}

void parse_command_line(int argc, char **argv, char *input_file_name, char *test_file_name, double *eta, double *rho, int *method_flag, int rank){
	int i;
	for (i=1;i<argc;i++){
		if (argv[i][0] != '-') break;
		if (++i>=argc)
			exit_with_help();
		switch(argv[i-1][1]){
			case 's':
				*method_flag = atoi(argv[i]);
				break;
			case 'r':
				*rho = atof(argv[i]);
				break;
			case 'e':
				*eta = atof(argv[i]);
				break;
			default :
				if (rank == 0)
					fprintf(stderr, "unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
				break;
		}
	}
	if (i>=argc -1)
		exit_with_help();
	strcpy(input_file_name, argv[i]);
	strcpy(test_file_name, argv[i+1]);
}
