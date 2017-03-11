CXX = mpic++
CC = mpicc
GCXX = g++
EXECS=server worker

all: ${EXECS}
server: master.cc
	${CXX} -o master master.cc

worker: worker.cc
	${CXX} -o worker worker.cc

clean:
	rm -f ${EXECS}
