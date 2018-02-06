CC = g++
OPT = -Wall -g -O0 --std=c++11 -fopenmp
OBJ = ../../BaseCell.o ../../BaseParticle.o ../../LymanAlphaLine.o ../../Utilities.o ../../rtsafe.o ../../RejectionMethod.o ../../RandomNumbers.o ../../BaseEmissionLine.o ../../LightParmParse/LParmParse.o main.o 

all: $(OBJ)
	$(CC) $(OPT) $(OBJ)  -o redistribution.exe
%.o: %.cpp
		$(CC)  $(OPT) -c  $< -o $@
clean:
		rm -f *.o
			rm -f *.exe

