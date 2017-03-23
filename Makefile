CXX=mpic++
CC=mpicc
GCXX=g++
FLAGS=-O2

evrsgd.out: main.o
	${CXX} main.o -o evrsgd.out ${FLAGS}

main.o: main.cc
	${CXX} -c main.cc -o main.o ${FLAGS}

clean:
	rm -f *.o *.out
