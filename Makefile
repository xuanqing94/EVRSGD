CXX=mpic++
FLAGS=-O2 -std=c++11

evrsgd.out: main.o command_line.o
	${CXX} main.o command_line.o -o evrsgd.out ${FLAGS}

main.o: main.cc
	${CXX} -c main.cc -o main.o ${FLAGS}

command_line.o: command_line.cc command_line.h 
	${CXX} -c command_line.cc -o command_line.o
clean:
	rm -f *.o *.out
