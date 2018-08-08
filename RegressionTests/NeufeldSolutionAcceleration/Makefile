CC = mpic++
OPT = -Wall -g -O3 --std=c++11 -fopenmp -I../../Parallel -I ../../Octet -I../../RamsesMPI -I../../OctetMPI -I../../../BaseOctet -DUSE_MPI -I../../Ramses
OCTOBJ =  ../../Octet/RamsesOctet.o ../../Octet/RamsesCell.o  ../../OctetDataset.o
OCTMPIOBJ =  ../../OctetMPI/RamsesOctetMPI.o ../../OctetDatasetMPI.o
TPGOBJ = ../../TPGDataset.o ../../TPG/TPGObject.o
RAMSESMPIOBJ =  ../../RamsesMPI/RamsesReaderMPI.o ../../RamsesDatasetMPI.o
RAMSESOBJ = ../../Ramses/RamsesReader.o ../../Ramses/RAMSES_amr_data.o   ../../Ramses/RAMSES_hydro_data.o  ../../Ramses/RAMSES_info.o ../../Ramses/FortranUnformatted_IO.o ../../RamsesDataset.o
OBJ =  ../../BaseCell.o ../../BaseDataset.o ../../BaseEmissionLine.o ../../BaseOutput.o ../../BaseParticle.o ../../BaseParticleVector.o ../../BaseSimulation.o ../../LymanAlphaLine.o ../../RandomNumbers.o ../../RejectionMethod.o ../../SphericalShellData.o ../../Utilities.o ../../LightParmParse/LParmParse.o ../../rtsafe.o ../../BaseEmissionModel.o ../../Parallel/Parallel.o ../../DustModule.o ../../ListEmissionModel.o ../../EmissionList.o ../../TraversalLength.o ../../Interpolate.o ../../NeutralFractionModule.o ../../InfiniteSlabData.o ../../UnigridDataset.o ../../Selectors.o $(RAMSESOBJ) $(OCTOBJ) $(RAMSESMPIOBJ) $(OCTMPIOBJ) $(TPGOBJ)
all: $(OBJ) main.o
	$(CC) $(OPT) $(OBJ)  main.o -o neufeld.exe
%.o: %.cpp
		$(CC)  $(OPT) -c  $< -o $@
clean:
		rm -f *.o
			rm -f *.exe

