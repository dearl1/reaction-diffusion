CXX=mpicxx
CC=mpicxx
CXXFLAGS=-std=c++11 -O3
# LDLIBS = -lboost_program_options -llapack -lblas
LDLIBS = -llapack -lblas -lboost_program_options

default: myprog
# default: myprog run

main.o: main.cpp

ReactionDiffusion.o: ReactionDiffusion.cpp

myprog: main.o ReactionDiffusion.o
	mpicxx -O3 -o $@ $^ $(LDLIBS)

test1: myprog
	mpiexec --oversubscribe -np 1 myprog --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.75 --b 0.06 --eps 50 --mu1 5 --mu2 0

test2: myprog
	mpiexec --oversubscribe -np 1 myprog --dt 0.001 --T 100 --Nx 251 --Ny 251 --a 0.75 --b 0.06 --eps 13 --mu1 5 --mu2 0

test3: myprog
	mpiexec --oversubscribe -np 1 myprog --dt 0.001 --T 100 --Nx 101 --Ny 101 --a 0.5 --b 0.1 --eps 50 --mu1 5 --mu2 0

test4: myprog
	mpiexec --oversubscribe -np 1 myprog --dt 0.001 --T 100 --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --eps 12.5 --mu1 1 --mu2 0.01

.PHONY: clean run test1 test2 test3 test4
# .PHONY: clean run
clean:
	rm -f *.o myprog



