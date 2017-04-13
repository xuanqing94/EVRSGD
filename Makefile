CXX=/usr/local/openmpi/bin/mpic++
FLAGS=-lpthread -O2 -std=c++11

evrsgd.out: main.o command_line.o loader.o
	${CXX} main.o command_line.o loader.o -o evrsgd.out ${FLAGS}

main.o: main.cc
	${CXX} -c main.cc -o main.o ${FLAGS}

command_line.o: command_line.cc command_line.h 
	${CXX} -c command_line.cc -o command_line.o

loader.o: loader.cc loader.h
	${CXX} -c loader.cc -o loader.o

.PHONY: clean

clean:
	rm -f *.o *.out
