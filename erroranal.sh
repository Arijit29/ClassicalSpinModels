# Script for error analysis of Monte Carlo

L=48
sweep=100000
T=5

python ./visualise/MonteCarlo.py $L $sweep $T
python ./visualise/Statistics.py $L $sweep $T
