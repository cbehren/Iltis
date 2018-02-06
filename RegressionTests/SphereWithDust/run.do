export OMP_NUM_THREADS=6
mpirun -np 8 ../../LLTC.exe inputs_1e5
wc Output.txt > escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_5e5
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_1e6
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_5e6
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_1e7
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_5e7
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_1e5_t001
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_5e5_t001
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_1e6_t001
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_5e6_t001
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_1e7_t001
wc Output.txt >> escape_fractions.txt
mpirun -np 8 ../../LLTC.exe inputs_5e7_t001
wc Output.txt >> escape_fractions.txt
