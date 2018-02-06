CC = h5c++
OPT = -Wall -g -O0 -lhdf5_cpp -lhdf5 -I$(HOME)/local_libraries/include -L$(HOME)/local_libraries/lib -I../LightParmParse
OBJ =  hdf5_generic.o ../LightParmParse/LParmParse.o utils.o

all: peeling nonpeeling
peeling: $(OBJ) text2hdf5.o
	$(CC) $(OPT) $(OBJ)  text2hdf5.o -o text2hdf5.exe

nonpeeling: $(OBJ) 2hdf5.o
	$(CC) $(OPT) $(OBJ)  2hdf5.o -o 2hdf5.exe
%.o: %.cpp
	$(CC)  $(OPT) -c  $< -o $@
clean:
	rm -f *.o
	rm -f *.exe
