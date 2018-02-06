CC = mpic++
OPT = -Wall -g -O3 --std=c++11 -fopenmp -IParallel -DUSE_MPI
OBJ =  Parallel.o

all: $(OBJ) test.o
	$(CC) $(OPT) $(OBJ) test.o -o test.exe
%.o: %.cpp
	$(CC)  $(OPT) -c  $< -o $@
clean:
	rm -f *.o
	rm -f *.exe
