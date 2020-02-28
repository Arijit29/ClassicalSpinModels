#! /bin/bash

# Parameters for the 2d-Ising model

L=20            # No. of sites
J=1             # interaction
h=0             # External field
sweep=20000
Nconf=100

time ./a.out $L $J $h $sweep $Nconf
