CC = g++
OPT = -Wall -g -O3 --std=c++11 -IParallel -fopenmp
# for MPI, use:
#CC = mpic++
#OPT = -Wall -g -O3 --std=c++11 -IParallel -fopenmp -DUSE_MPI
PLOTOBJ = BaseCell.o BaseDataset.o BaseEmissionLine.o BaseOutput.o BaseParticle.o BaseParticleVector.o PlottingInterface.o LymanAlphaLine.o RandomNumbers.o RejectionMethod.o SphericalShellData.o Utilities.o LightParmParse/LParmParse.o rtsafe.o BaseEmissionModel.o Parallel/Parallel.o DustModule.o ListEmissionModel.o EmissionList.o TraversalLength.o Interpolate.o NeutralFractionModule.o InfiniteSlabData.o UnigridDataset.o  PlottingOperators.o Selectors.o
OBJ =  BaseCell.o BaseDataset.o BaseEmissionLine.o BaseOutput.o BaseParticle.o BaseParticleVector.o BaseSimulation.o LymanAlphaLine.o RandomNumbers.o RejectionMethod.o SphericalShellData.o Utilities.o LightParmParse/LParmParse.o rtsafe.o BaseEmissionModel.o Parallel/Parallel.o DustModule.o ListEmissionModel.o EmissionList.o TraversalLength.o Interpolate.o NeutralFractionModule.o InfiniteSlabData.o UnigridDataset.o Selectors.o 
all: $(OBJ) main.o
	$(CC) $(OPT) $(OBJ) main.o -o Iltis.exe
plot: $(PLOTOBJ) plotter.o
	$(CC) $(OPT) $(PLOTOBJ) plotter.o -o plotter.exe
comm: $(OBJ) main_communication.o
	$(CC) $(OPT) $(OBJ) main_communication.o -o main_communication.exe
%.o: %.cpp
	$(CC)  $(OPT) -c  $< -o $@
clean:
	rm -f *.o
	rm -f *.exe
