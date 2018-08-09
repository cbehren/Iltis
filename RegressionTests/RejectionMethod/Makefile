CC = g++
OPT = -Wall -g -O3 --std=c++11 -I../../Parallel
OBJ =  ../../BaseCell.o ../../BaseDataset.o ../../BaseEmissionLine.o ../../BaseOutput.o ../../BaseParticle.o ../../BaseParticleVector.o ../../BaseSimulation.o ../../LymanAlphaLine.o ../../RandomNumbers.o ../../RejectionMethod.o ../../SphericalShellData.o ../../Utilities.o ../../LightParmParse/LParmParse.o ../../rtsafe.o ../../BaseEmissionModel.o ../../Parallel/Parallel.o ../../DustModule.o ../../ListEmissionModel.o ../../EmissionList.o ../../TraversalLength.o ../../Interpolate.o ../../NeutralFractionModule.o ../../InfiniteSlabData.o ../../UnigridDataset.o ../../Selectors.o
rejection: $(OBJ) main_rejection.o
	$(CC) $(OPT) $(OBJ)  main_rejection.o -o rejection.exe
%.o: %.cpp
		$(CC)  $(OPT) -c  $< -o $@
clean:
		rm -f *.o
			rm -f *.exe

