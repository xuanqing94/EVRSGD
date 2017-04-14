USER=$(shell whoami)
MPI_PATH=/home/${USER}/mpich
CXX=${MPI_PATH}/bin/mpic++
FLAGS=-O2 -std=c++11 -pthread -L${MPI_PATH}/lib -I${MPI_PATH}/include 

.PHONY: clean hosts test

evrsgd.out: main.o command_line.o loader.o
	${CXX} main.o command_line.o loader.o -o evrsgd.out ${FLAGS}

test: hosts evrsgd.out	
	~/hydra/bin/mpiexec -n 6 --machinefile host_file ./evrsgd.out -a 1 -r 3 -s 0.1 -e 1.0e-5 -f ../real-sim 

hosts:
	echo thrall.cs.ucdavis.edu:6 > host_file
	echo jaina.cs.ucdavis.edu:1 > host_file
	echo illidan.cs.ucdavis.edu:6 >> host_file

main.o: main.cc
	${CXX} -c main.cc -o main.o ${FLAGS}

command_line.o: command_line.cc command_line.h 
	${CXX} -c command_line.cc -o command_line.o

loader.o: loader.cc loader.h
	${CXX} -c loader.cc -o loader.o

clean:
	rm -f *.o *.out
