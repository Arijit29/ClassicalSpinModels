#! /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -q default
#PBS -j oe
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > pbsnodes

# Parameters for the 2d-Ising model

L=$l            # No. of sites
sweep=50000     # No. of MC sweeps
Nconf=500

time ./a.out $L $sweep $Nconf
