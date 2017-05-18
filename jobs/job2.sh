#!/bin/bash -l

# job name
#SBATCH -J myjob
# account
#SBATCH -A edu17.SF2568
# email notification
#SBATCH --mail-type=BEGIN,END
# 10 minutes wall-clock time will be given to this job
#SBATCH -t 00:30:00
# Number of nodes
#SBATCH --nodes=1
# set tasks per node to 24 in order to disablr hyperthreading
#SBATCH --ntasks-per-node=24

module add i-compilers intelmpi


mpirun -l -np 1 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 2 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 3 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 4 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 5 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 6 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 7 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 8 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 9 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 10 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 11 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 12 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 13 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 14 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 15 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 16 ./bin/do_par ./data/lesmis.graph

mpirun -l -np 17 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 18 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 19 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 20 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 21 ./bin/do_par ./data/lesmis.graph
mpirun -l -np 22 ./bin/do_par ./data/lesmis.graph
